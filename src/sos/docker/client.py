#!/usr/bin/env python3
#
# Copyright (C) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
from io import BytesIO

import docker
from sos.eval import interpolate
from sos.targets import sos_targets
from sos.utils import env

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
        except Exception:
            # some system does not have cat or grep
            self.tot_mem = None
        return self.tot_mem

    def _is_image_avail(self, image):
        try:
            images = sum([x.tags for x in self.client.images.list()], [])
        except AttributeError:
            raise RuntimeError(
                'Incompatible version of docker module detected. If you are using "docker-py", please uninstall it and install module "docker".')
        # some earlier version of docker-py returns docker.io/ for global repositories
        images = [x[10:] if x.startswith('docker.io/') else x for x in images]
        return (':' in image and image in images) or \
            (':' not in image and '{}:latest'.format(image) in images)

    def build(self, script, **kwargs):
        if not self.client:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        if script is not None:
            f = BytesIO(script.encode('utf-8'))
            self.client.images.build(fileobj=f, **kwargs)
            # self.stream(line.decode())
        else:
            self.client.images.build(**kwargs)
            # self.stream(line.decode())
        # if a tag is given, check if the image is built
        if 'tag' in kwargs and not self._is_image_avail(kwargs['tag']):
            raise RuntimeError('Image with tag {} is not created.'.format(kwargs['tag']))

    def load_image(self, image, **kwargs):
        if not self.client:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        env.logger.info('docker load {}'.format(image))
        self.client.images.load(image, **kwargs)

    def pull(self, image):
        if not self.client:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        # if image is specified, check if it is available locally. If not, pull it
        ret = 0
        if not self._is_image_avail(image):
            env.logger.info('docker pull {}'.format(image))
            # using subprocess instead of docker-py's pull function because this would have
            # much better progress bar display
            ret = subprocess.call('docker pull {}'.format(image), shell=True)
            # for line in self.client.pull(image, stream=True):
            #    self.stream(line)
        if not self._is_image_avail(image):
            raise RuntimeError('Failed to pull image {}'.format(image))
        return ret

    # def commit(self, **kwargs):
    #    if not self.client:
    #        raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
    #    for line in self.client.commit(**kwargs):
    #        self.stream(line.decode())
    #    return 0

    def run(self, image, script='', interpreter='', args='', suffix='.sh', **kwargs):
        if self.client is None:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
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
                args = '{filename:q}'
            #
            binds = []
            if 'volumes' in kwargs:
                volumes = [kwargs['volumes']] if isinstance(
                    kwargs['volumes'], str) else kwargs['volumes']
                for vol in volumes:
                    if not vol:
                        continue
                    if vol.count(':') != 1:
                        raise RuntimeError(
                            'Please specify columes in the format of host_dir:mnt_dir')
                    host_dir, mnt_dir = vol.split(':')
                    if platform.system() == 'Darwin':
                        # under Darwin, host_dir must be under /Users
                        if not os.path.abspath(host_dir).startswith('/Users') and not (self.has_volumes and os.path.abspath(host_dir).startswith('/Volumes')):
                            raise RuntimeError(
                                'hostdir ({}) under MacOSX must be under /Users or /Volumes (if properly configured, see https://github.com/vatlab/SOS/wiki/SoS-Docker-guide for details) to be usable in docker container'.format(host_dir))
                    binds.append('{}:{}'.format(os.path.abspath(host_dir), mnt_dir))
            #
            volumes_opt = ' '.join('-v {}'.format(x) for x in binds)
            # under mac, we by default share /Users within docker
            wdir = os.path.abspath(kwargs['working_dir']
                                   if 'working_dir' in kwargs else os.getcwd())
            if platform.system() == 'Darwin':
                if not any(x.startswith('/Users:') for x in binds):
                    volumes_opt += ' -v /Users:/Users'
                if self.has_volumes:
                    volumes_opt += ' -v /Volumes:/Volumes'
                if not wdir.startswith('/Users'):
                    volumes_opt += f' -v /{wdir}:/{wdir}'
            elif platform.system() == 'Linux':
                if not any(x.startswith('/home:') for x in binds):
                    volumes_opt += ' -v /home:/home'
                if not wdir.startswith('/home/'):
                    volumes_opt += f' -v /{wdir}:/{wdir}'
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
                    volumes_from_opt = f'--volumes_from={kwargs["volumes_from"]}'
                elif isinstance(kwargs['volumes_from'], list):
                    volumes_from_opt = ' '.join(
                        f'--volumes_from={x}' for x in kwargs['volumes_from'])
                else:
                    raise RuntimeError('Option volumes_from only accept a string or list of string: {} specified'.format(
                        kwargs['volumes_from']))

            # we also need to mount the script
            if script:
                volumes_opt += ' -v {}:{}'.format(shlex.quote(os.path.join(
                    tempdir, tempscript)), '/var/lib/sos/{}'.format(tempscript))
            cmd_opt = interpolate(f'{interpreter if isinstance(interpreter, str) else interpreter[0]} {args}', {
                'filename': sos_targets(f'/var/lib/sos/{tempscript}'),
                'script': script})
            #
            working_dir_opt = '-w={}'.format(shlex.quote(os.path.abspath(os.getcwd())))
            if 'working_dir' in kwargs:
                if not os.path.isabs(kwargs['working_dir']):
                    env.logger.warning('An absolute path is needed for -w option of docker run command. "{}" provided, "{}" used.'
                                       .format(kwargs['working_dir'], os.path.abspath(os.path.expanduser(kwargs['working_dir']))))
                    working_dir_opt = '-w={}'.format(os.path.abspath(
                        os.path.expanduser(kwargs['working_dir'])))
                else:
                    working_dir_opt = '-w={}'.format(kwargs['working_dir'])

            env_opt = ''
            if 'environment' in kwargs:
                if isinstance(kwargs['environment'], dict):
                    env_opt = ' '.join(f'-e {x}={y}' for x, y in kwargs['environment'].items())
                elif isinstance(kwargs['environment'], list):
                    env_opt = ' '.join(f'-e {x}' for x in kwargs['environment'])
                elif isinstance(kwargs['environment'], str):
                    env_opt = f'-e {kwargs["environment"]}'
                else:
                    raise RuntimeError('Invalid value for option environment (str, list, or dict is allowd, {} provided)'.format(
                        kwargs['environment']))
            #
            port_opt = '-P'
            if 'port' in kwargs:
                if isinstance(kwargs['port'], (str, int)):
                    port_opt = '-p {}'.format(kwargs['port'])
                elif isinstance(kwargs['port'], list):
                    port_opt = ' '.join('-p {}'.format(x) for x in kwargs['port'])
                else:
                    raise RuntimeError(
                        'Invalid value for option port (a list of intergers), {} provided'.format(kwargs['port']))
            #
            name_opt = ''
            if 'name' in kwargs:
                name_opt = f'--name={kwargs["name"]}'
            #
            stdin_opt = ''
            if 'stdin_open' in kwargs and kwargs['stdin_optn']:
                stdin_opt = '-i'
            #
            tty_opt = '-t'
            if 'tty' in kwargs and not kwargs['tty']:
                tty_opt = ''
            #
            if 'user' in kwargs:
                user_opt = f'-u {kwargs["user"]}'
            else:
                # Tocket #922
                user_opt = f'-u {os.getuid()}:{os.getgid()}'
            #
            extra_opt = ''
            if 'extra_args' in kwargs:
                extra_opt = kwargs['extra_args']
            #
            security_opt = ''
            if platform.system() == 'Linux':
                # this is for a selinux problem when /var/sos/script cannot be executed
                security_opt = '--security-opt label:disable'
            cmd = 'docker run --rm {} {} {} {} {} {} {} {} {} {} {} {} {} {}'.format(
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
            env.logger.debug(cmd)

            if env.config['run_mode'] == 'interactive':
                if 'stdout' in kwargs or 'stderr' in kwargs:
                    child = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE, bufsize=0)
                    out, err = child.communicate()
                    if 'stdout' in kwargs:
                        if kwargs['stdout'] is not False:
                            with open(kwargs['stdout'], 'ab') as so:
                                so.write(out)
                    else:
                        sys.stdout.write(out.decode())

                    if 'stderr' in kwargs:
                        if kwargs['stderr'] is not False:
                            with open(kwargs['stderr'], 'ab') as se:
                                se.write(err)
                    else:
                        sys.stderr.write(err.decode())
                    ret = child.returncode
                else:
                    # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
                    from .utils import pexpect_run
                    ret = pexpect_run(cmd)
            elif '__std_out__' in env.sos_dict and '__std_err__' in env.sos_dict:
                if 'stdout' in kwargs or 'stderr' in kwargs:
                    if 'stdout' in kwargs:
                        if kwargs['stdout'] is False:
                            so = subprocess.DEVNULL
                        else:
                            so = open(kwargs['stdout'], 'ab')
                    elif env.verbosity > 0:
                        so = open(env.sos_dict['__std_out__'], 'ab')
                    else:
                        so = subprocess.DEVNULL

                    if 'stderr' in kwargs:
                        if kwargs['stderr'] is False:
                            se = subprocess.DEVNULL
                        else:
                            se = open(kwargs['stderr'], 'ab')
                    elif env.verbosity > 1:
                        se = open(env.sos_dict['__std_err__'], 'ab')
                    else:
                        se = subprocess.DEVNULL

                    p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)
                    ret = p.wait()

                    if so != subprocess.DEVNULL:
                        so.close()
                    if se != subprocess.DEVNULL:
                        se.close()

                elif env.verbosity >= 1:
                    with open(env.sos_dict['__std_out__'], 'ab') as so, open(env.sos_dict['__std_err__'], 'ab') as se:
                        p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)
                        ret = p.wait()
                else:
                    p = subprocess.Popen(
                        cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                    ret = p.wait()
            else:
                if 'stdout' in kwargs:
                    if kwargs['stdout'] is False:
                        so = subprocess.DEVNULL
                    else:
                        so = open(kwargs['stdout'], 'ab')
                elif env.verbosity > 0:
                    so = None
                else:
                    so = subprocess.DEVNULL

                if 'stderr' in kwargs:
                    if kwargs['stderr'] is False:
                        se = subprocess.DEVNULL
                    else:
                        se = open(kwargs['stderr'], 'ab')
                elif env.verbosity > 1:
                    se = None
                else:
                    se = subprocess.DEVNULL

                p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)

                ret = p.wait()
                if so is not None and so != subprocess.DEVNULL:
                    so.close()
                if se is not None and se != subprocess.DEVNULL:
                    se.close()

            if ret != 0:
                debug_script_dir = os.path.join(env.exec_dir, '.sos')
                msg = 'The script has been saved to {}/{}. To reproduce the error please run:\n``{}``'.format(
                    debug_script_dir, tempscript, cmd.replace(tempdir, debug_script_dir))
                shutil.copy(os.path.join(tempdir, tempscript), debug_script_dir)
                if ret == 125:
                    msg = 'Docker daemon failed (exitcode=125). ' + msg
                elif ret == 126:
                    msg = 'Failed to invoke specified command (exitcode=126). ' + msg
                elif ret == 127:
                    msg = 'Failed to locate specified command (exitcode=127). ' + msg
                elif ret == 137:
                    if not hasattr(self, 'tot_mem'):
                        self.tot_mem = self.total_memory(image)
                    if self.tot_mem is None:
                        msg = 'Script killed by docker. ' + msg
                    else:
                        msg = 'Script killed by docker, probably because of lack of RAM (available RAM={:.1f}GB, exitcode=137). '.format(
                            self.tot_mem / 1024 / 1024) + msg
                else:
                    msg =  f"Executing script in docker returns an error (exitcode={ret}{', err=``%s``' % kwargs['stderr'] if 'stderr' in kwargs and os.path.isfile(kwargs['stderr']) else ''}).\n{msg}"
                raise subprocess.CalledProcessError(
                            returncode = ret,
                            cmd = cmd.replace(tempdir, debug_script_dir),
                            stderr = msg)
        return 0
