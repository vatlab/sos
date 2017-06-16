#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import os
import sys
import re
import subprocess
import tempfile
import json
import platform
import shutil
import shlex
import docker

from io import BytesIO
from sos.utils import env
from sos.sos_eval import interpolate

#
# docker support
#
class SoS_DockerClient:
    '''A singleton class to ensure there is only one client'''
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(SoS_DockerClient, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        try:
            self.client = docker.from_env()
            self.client.info()
            # mount the /Volumes folder under mac, please refer to
            #    http://vatlab.github.io/SOS/doc/tutorials/SoS_Docker_Guide.html
            # for details.
            self.has_volumes = False
            if platform.system() == 'Darwin':
                try:
                    # this command log in to the docker machine, check if /Volumes has been mounted,
                    # and try to mount it if possible. This requires users to configure
                    subprocess.call("""docker-machine ssh "{}" 'mount | grep /Volumes || {{ echo "mounting /Volumes"; sudo mount  -t vboxsf Volumes /Volumes; }}' """.format(os.environ['DOCKER_MACHINE_NAME']),
                        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    env.logger.trace('Sucessfully mount /Volumes to virtual machine')
                    self.has_volumes = True
                except Exception as e:
                    env.logger.trace('Failed to mount /Volumes to virtual machine: {}'.format(e))
        except Exception as e:
            env.logger.debug('Docker client init fail: {}'.format(e))
            self.client = None

    def total_memory(self, image='ubuntu'):
        '''Get the available ram fo the docker machine in Kb'''
        try:
            ret = subprocess.check_output(
                '''docker run -t {} cat /proc/meminfo  | grep MemTotal'''.format(image),
                shell=True, stdin=subprocess.DEVNULL)
            # ret: MemTotal:       30208916 kB
            self.tot_mem = int(ret.split()[1])
        except:
            # some system does not have cat or grep
            self.tot_mem = None
        return self.tot_mem

    def _is_image_avail(self, image):
        images = sum([x.tags for x in self.client.images.list()], [])
        # some earlier version of docker-py returns docker.io/ for global repositories
        images = [x[10:] if x.startswith('docker.io/') else x for x in images]
        return (':' in image and image in images) or \
            (':' not in image and '{}:latest'.format(image) in images)

    def stream(self, line):
        # properly output streamed output
        try:
            sys.stdout.write(json.loads(line).get('stream', ''))
        except ValueError:
            # sometimes all the data is sent on a single line ????
            #
            # ValueError: Extra data: line 1 column 87 - line 1 column
            # 33268 (char 86 - 33267)
            # This ONLY works because every line is formatted as
            # {"stream": STRING}
            for obj in re.findall('{\s*"stream"\s*:\s*"[^"]*"\s*}', line):
                sys.stdout.write(json.loads(obj).get('stream', ''))

    def build(self, script, **kwargs):
        if not self.client:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        if script is not None:
            f = BytesIO(script.encode('utf-8'))
            image = self.client.images.build(fileobj=f, **kwargs)
            #self.stream(line.decode())
        else:
            image = self.client.images.build(**kwargs)
            #self.stream(line.decode())
        # if a tag is given, check if the image is built
        if 'tag' in kwargs and not self._is_image_avail(kwargs['tag']):
            raise RuntimeError('Image with tag {} is not created.'.format(kwargs['tag']))

    def load_image(self, image, **kwargs):
        if not self.client:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        env.logger.info('docker load {}'.format(image))
        self.client.images.load(image, **kwargs)

    def pull(self, image):
        if not self.client:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        # if image is specified, check if it is available locally. If not, pull it
        ret = 0
        if not self._is_image_avail(image):
            env.logger.info('docker pull {}'.format(image))
            # using subprocess instead of docker-py's pull function because this would have
            # much better progress bar display
            ret = subprocess.call('docker pull {}'.format(image), shell=True)
            #for line in self.client.pull(image, stream=True):
            #    self.stream(line)
        if not self._is_image_avail(image):
            raise RuntimeError('Failed to pull image {}'.format(image))
        return ret

    def commit(self, **kwargs):
        if not self.client:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        for line in self.client.commit(**kwargs):
            self.stream(line.decode())
        return 0

    def run(self, image, script='', interpreter='', args='', suffix='.sh', **kwargs):
        if self.client is None:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        #
        env.logger.debug('docker_run with keyword args {}'.format(kwargs))
        #
        # now, write a temporary file to a tempoary directory under the current directory, this is because
        # we need to share the directory to ...
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            # keep the temporary script for debugging purposes
            # tempdir = tempfile.mkdtemp(dir=os.getcwd())
            if script:
                tempscript = 'docker_run_{}{}'.format(os.getpid(), suffix)
                with open(os.path.join(tempdir, tempscript), 'w') as script_file:
                    script_file.write(script)
            #
            # if there is an interpreter and with args
            if not args:
                args = '${filename!q}'
            if not interpreter:
                interpreter = '/bin/bash'
                # if there is a shebang line, we ...
                if script.startswith('#!'):
                    # make the script executable
                    env.logger.warning('Shebang line in a docker-run script is ignored')
            elif not isinstance(interpreter, str):
                interpreter = interpreter[0]
            #
            binds = []
            if 'volumes' in kwargs:
                volumes = [kwargs['volumes']] if isinstance(kwargs['volumes'], str) else kwargs['volumes']
                for vol in volumes:
                    if not vol:
                        continue
                    if vol.count(':') != 1:
                        raise RuntimeError('Please specify columes in the format of host_dir:mnt_dir')
                    host_dir, mnt_dir = vol.split(':')
                    if platform.system() == 'Darwin':
                        # under Darwin, host_dir must be under /Users
                        if not os.path.abspath(host_dir).startswith('/Users') and not (self.has_volumes and os.path.abspath(host_dir).startswith('/Volumes')):
                            raise RuntimeError('hostdir ({}) under MacOSX must be under /Users or /Volumes (if properly configured, see https://github.com/vatlab/SOS/wiki/SoS-Docker-guide for details) to be usable in docker container'.format(host_dir))
                    binds.append('{}:{}'.format(os.path.abspath(host_dir), mnt_dir))
            #
            volumes_opt = ' '.join('-v {}'.format(x) for x in binds)
            # under mac, we by default share /Users within docker
            if platform.system() == 'Darwin':
                if not any(x.startswith('/Users:') for x in binds):
                    volumes_opt += ' -v /Users:/Users'
                if self.has_volumes:
                    volumes_opt += ' -v /Volumes:/Volumes'
            elif platform.system() == 'Linux':
                if not any(x.startswith('/home:') for x in binds):
                    volumes_opt += ' -v /home:/home'
            if not any(x.startswith('/tmp:') for x in binds):
                volumes_opt += ' -v /tmp:/tmp'
            #
            mem_limit_opt = ''
            if 'mem_limit' in kwargs:
                mem_limit_opt = '--memory={}'.format(kwargs['mem_limit'])
            #
            volumes_from_opt = ''
            if 'volumes_from' in kwargs:
                if isinstance(kwargs['volumes_from'], str):
                    volumes_from_opt = '--volumes_from={}'.format(kwargs['volumes_from'])
                elif isinstance(kwargs['volumes_from'], list):
                    volumes_from_opt = ' '.join('--volumes_from={}'.format(x) for x in kwargs['volumes_from'])
                else:
                    raise RuntimeError('Option volumes_from only accept a string or list of string'.format(kwargs['volumes_from']))
            # we also need to mount the script
            cmd_opt = ''
            if script and interpreter:
                volumes_opt += ' -v {}:{}'.format(shlex.quote(os.path.join(tempdir, tempscript)), '/var/lib/sos/{}'.format(tempscript))
                cmd_opt = interpolate('{} {}'.format(interpreter, args), '${ }',
                            {'filename': '/var/lib/sos/{}'.format(tempscript)})
            #
            working_dir_opt = '-w={}'.format(shlex.quote(os.path.abspath(os.getcwd())))
            if 'working_dir' in kwargs:
                if not os.path.isabs(kwargs['working_dir']):
                    env.logger.warning('An absolute path is needed for -w option of docker run command. "{}" provided, "{}" used.'
                        .format(kwargs['working_dir'], os.path.abspath(os.path.expanduser(kwargs['working_dir']))))
                    working_dir_opt = '-w={}'.format(os.path.abspath(os.path.expanduser(kwargs['working_dir'])))
                else:
                    working_dir_opt = '-w={}'.format(kwargs['working_dir'])
            #
            env_opt = ''
            if 'environment' in kwargs:
                if isinstance(kwargs['environment'], dict):
                    env_opt = ' '.join('-e {}={}'.format(x,y) for x,y in kwargs['environment'].items())
                elif isinstance(kwargs['environment'], list):
                    env_opt = ' '.join('-e {}'.format(x) for x in kwargs['environment'])
                elif isinstance(kwargs['environment'], str):
                    env_opt = '-e {}'.format(kwargs['environment'])
                else:
                    raise RuntimeError('Invalid value for option environment (str, list, or dict is allowd, {} provided)'.format(kwargs['environment']))
            #
            port_opt = '-P'
            if 'port' in kwargs:
                if isinstance(kwargs['port'], (str, int)):
                    port_opt = '-p {}'.format(kwargs['port'])
                elif isinstance(kwargs['port'], list):
                    port_opt = ' '.join('-p {}'.format(x) for x in kwargs['port'])
                else:
                    raise RuntimeError('Invalid value for option port (a list of intergers), {} provided'.format(kwargs['port']))
            #
            name_opt = ''
            if 'name' in kwargs:
                name_opt = '--name={}'.format(kwargs['name'])
            #
            stdin_opt = ''
            if 'stdin_open' in kwargs and kwargs['stdin_optn']:
                stdin_opt = '-i'
            #
            tty_opt = '-t'
            if 'tty' in kwargs and not kwargs['tty']:
                tty_opt = ''
            #
            user_opt = ''
            if 'user' in kwargs:
                user_opt = '-u {}'.format(kwargs['user'])
            #
            extra_opt = ''
            if 'extra_args' in kwargs:
                extra_opt = kwargs['extra_args']
            #
            security_opt = ''
            if platform.system() == 'Linux':
                # this is for a selinux problem when /var/sos/script cannot be executed
                security_opt = '--security-opt label:disable'
            command = 'docker run --rm {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
                security_opt,       # security option
                volumes_opt,        # volumes
                volumes_from_opt,   # volumes_from
                name_opt,           # name
                stdin_opt,          # stdin_optn
                tty_opt,            # tty
                port_opt,           # port
                working_dir_opt,    # working dir
                user_opt,           # user
                env_opt,            # environment
                mem_limit_opt,      # memory limit
                extra_opt,          # any extra parameters
                image,              # image
                cmd_opt
                )
            env.logger.info(command)
            ret = subprocess.call(command, shell=True)
            if ret != 0:
                msg = 'The script has been saved to .sos/{} so that you can execute it using the following command:\n{}'.format(
                    tempscript, command.replace(tempdir, os.path.abspath('./.sos')))
                shutil.copy(os.path.join(tempdir, tempscript), '.sos')
                if ret == 125:
                    raise RuntimeError('Docker daemon failed (exitcode=125). ' + msg)
                elif ret == 126:
                    raise RuntimeError('Failed to invoke specified command (exitcode=126). ' + msg)
                elif ret == 127:
                    raise RuntimeError('Failed to locate specified command (exitcode=127). ' + msg)
                elif ret == 137:
                    if not hasattr(self, 'tot_mem'):
                        self.tot_mem = self.total_memory(image)
                    if self.tot_mem is None:
                        raise RuntimeError('Script killed by docker. ' + msg)
                    else:
                        raise RuntimeError('Script killed by docker, probably because of lack of RAM (available RAM={:.1f}GB, exitcode=137). '.format(self.tot_mem/1024/1024) + msg)
                else:
                    raise RuntimeError('Executing script in docker returns an error (exitcode={}). '.format(ret) + msg)
        return 0
