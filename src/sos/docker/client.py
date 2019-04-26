#!/usr/bin/env python3
#
# Copyright (C) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import platform
import shutil
import subprocess
import sys
import tempfile

from sos.eval import interpolate
from sos.targets import sos_targets, path
from sos.utils import env, pexpect_run

#
# docker support
#


class SoS_DockerClient:
    '''A singleton class to ensure there is only one client'''
    _instance = None

    client = shutil.which('docker')
    pulled_images = set()

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(SoS_DockerClient, cls).__new__(cls)
        return cls._instance

    def total_memory(self, image='ubuntu'):
        '''Get the available ram fo the docker machine in Kb'''
        try:
            ret = subprocess.check_output(
                f'''docker run -t {image} cat /proc/meminfo  | grep MemTotal''',
                shell=True,
                stdin=subprocess.DEVNULL)
            # ret: MemTotal:       30208916 kB
            self.tot_mem = int(ret.split()[1])
        except Exception:
            # some system does not have cat or grep
            self.tot_mem = None
        return self.tot_mem

    def _is_image_avail(self, image):
        # the command will return ID of the image if it exists
        try:
            return bool(
                subprocess.check_output(
                    f'''docker images {image} --no-trunc --format "{{{{.ID}}}}"''',
                    shell=True))
        except Exception as e:
            env.logger.warning(f'Failed to check image {image}: {e}')
            return False

    def _run_cmd(self, cmd, **kwargs):
        if env.config['run_mode'] == 'interactive':
            if 'stdout' in kwargs or 'stderr' in kwargs:
                child = subprocess.Popen(
                    cmd,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    bufsize=0)
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
                with open(env.sos_dict['__std_out__'],
                          'ab') as so, open(env.sos_dict['__std_err__'],
                                            'ab') as se:
                    p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)
                    ret = p.wait()
            else:
                p = subprocess.Popen(
                    cmd,
                    shell=True,
                    stderr=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL)
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
        return ret

    def build(self, script, **kwargs):
        if not self.client:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?'
            )
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            if script:
                docker_file = os.path.join(tempdir, 'Dockerfile')
                with open(docker_file, 'w') as df:
                    df.write(script)
                file_opt = ['-f', docker_file, '.']
            else:
                if 'file' not in kwargs:
                    raise RuntimeError(
                        'Docker file must be specified with option file if not directly included.'
                    )
                file_opt = ['--file', kwargs['file']]

            other_opts = []
            for arg, value in kwargs.items():
                # boolean args
                if arg in ('compress', 'disable_content_trust', 'force_rm',
                           'memory_swap', 'no_cache', 'pull', 'quiet', 'rm',
                           'squash', 'stream'):
                    if value is True:
                        other_opts.append(f'--{arg.replace("_", "-")}')
                    else:
                        env.logger.warning(
                            f'Boolean {arg} is ignored (True should be provided)'
                        )
                elif arg in ('add_host', 'build_arg', 'cache_from',
                             'cgroup_parent', 'cpu_period', 'cpu_quota',
                             'cpu-shares', 'cpuset_cpus', 'cpuset_mems',
                             'label', 'memory', 'network', 'platform',
                             'security_opt', 'shm_size', 'tag', 'target',
                             'ulimit'):
                    other_opts.extend([f'--{arg.replace("_", "-")}', value])

            cmd = subprocess.list2cmdline(['docker', 'build'] + file_opt +
                                          other_opts)

            env.logger.debug(cmd)
            if env.config['run_mode'] == 'dryrun':
                print(f'HINT: {cmd}')
                print(script)
                return 0

            ret = self._run_cmd(cmd, **kwargs)

            if ret != 0:
                if script:
                    debug_script_dir = os.path.join(env.exec_dir, '.sos')
                    msg = 'The Dockerfile has been saved to {}/Dockerfile. To reproduce the error please run:\n``{}``'.format(
                        debug_script_dir, cmd.replace(tempdir,
                                                      debug_script_dir))
                    shutil.copy(
                        os.path.join(tempdir, 'Dockerfile'), debug_script_dir)
                else:
                    msg = f'To reproduce this error please run {cmd}'
                raise subprocess.CalledProcessError(
                    returncode=ret, cmd=cmd, stderr=msg)
        # if a tag is given, check if the image is built
        if 'tag' in kwargs and not self._is_image_avail(kwargs['tag']):
            raise RuntimeError('Image with tag {} is not created.'.format(
                kwargs['tag']))

    def load_image(self, image, **kwargs):
        if not self.client:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?'
            )
        env.logger.info('docker load {}'.format(image))
        try:
            subprocess.call(f'''docker load -i {image} --quiet''', shell=True)
        except Exception as e:
            raise RuntimeError(f'Failed to load image {image}: {e}')

    def pull(self, image):
        if not self.client:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?'
            )
        if image in self.pulled_images:
            return
        # if image is specified, check if it is available locally. If not, pull it
        err_msg = ''
        try:
            print(f'HINT: Pulling docker image {image}')
            subprocess.check_output(
                'docker pull {}'.format(image),
                stderr=subprocess.STDOUT,
                shell=True,
                universal_newlines=True)
        except subprocess.CalledProcessError as exc:
            err_msg = exc.output
        if not self._is_image_avail(image):
            raise RuntimeError(
                f'Failed to pull docker image {image}:\n {err_msg}')
        self.pulled_images.add(image)

    def run(self,
            image,
            script='',
            interpreter='',
            args='',
            suffix='.sh',
            **kwargs):
        if self.client is None:
            raise RuntimeError(
                'Cannot connect to the Docker daemon. Is the docker daemon running on this host?'
            )
        #
        env.logger.debug('docker_run with keyword args {}'.format(kwargs))
        #
        # now, write a temporary file to a tempoary directory under the current directory, this is because
        # we need to share the directory to ...
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            # keep the temporary script for debugging purposes
            # tempdir = tempfile.mkdtemp(dir=os.getcwd())
            tempscript = 'docker_run_{}{}'.format(os.getpid(), suffix)
            if script:
                with open(os.path.join(tempdir, tempscript),
                          'w') as script_file:
                    # the input script might have windows new line but the container
                    # will need linux new line for proper execution #1023
                    script_file.write('\n'.join(script.splitlines()))
            #
            # if there is an interpreter and with args
            if not args:
                args = '{filename:pq}'
            #
            # under mac, we by default share /Users within docker
            wdir = os.path.abspath(os.getcwd())
            binds = []
            if 'volumes' in kwargs:
                volumes = [kwargs['volumes']] if isinstance(
                    kwargs['volumes'], str) else kwargs['volumes']
                has_wdir = False
                for vol in volumes:
                    if not vol:
                        continue
                    if isinstance(vol, (str, path)):
                        vol = str(vol)
                    else:
                        raise ValueError(
                            f'Unacceptable value {vol} for parameter volumes')
                    if vol.count(':') == 0:
                        host_dir, mnt_dir = vol, vol
                    elif vol.count(':') in (1, 2):
                        host_dir, mnt_dir = vol.split(':', 1)
                    else:
                        raise ValueError(
                            f'Invalid format for volume specification: {vol}')
                    binds.append(
                        f'{path(host_dir).resolve():p}:{path(mnt_dir):p}')
                    if wdir.startswith(
                            os.path.abspath(os.path.expanduser(host_dir))):
                        has_wdir = True
                volumes_opt = ' '.join('-v {}'.format(x) for x in binds)
                if not has_wdir:
                    volumes_opt += f' -v {path(wdir):p}:{path(wdir):p}'
            else:
                volumes_opt = f' -v {path(wdir):p}:{path(wdir):p}'

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
                    raise RuntimeError(
                        'Option volumes_from only accept a string or list of string: {} specified'
                        .format(kwargs['volumes_from']))

            # we also need to mount the script
            if script:
                volumes_opt += f' -v {path(tempdir)/tempscript:p}:/var/lib/sos/{tempscript}'
            cmd_opt = interpolate(
                f'{interpreter if isinstance(interpreter, str) else interpreter[0]} {args}',
                {
                    'filename': sos_targets(f'/var/lib/sos/{tempscript}'),
                    'script': script
                })
            #
            workdir_opt = ''
            if 'docker_workdir' in kwargs and kwargs[
                    'docker_workdir'] is not None:
                if not os.path.isabs(kwargs['docker_workdir']):
                    env.logger.warning(
                        'An absolute path is needed for -w option of docker run command. "{}" provided, "{}" used.'
                        .format(
                            kwargs['docker_workdir'],
                            os.path.abspath(
                                os.path.expanduser(kwargs['docker_workdir']))))
                    workdir_opt = f'-w={path(kwargs["docker_workdir"]).resolve():p}'
                else:
                    workdir_opt = f'-w={path(kwargs["docker_workdir"]):p}'
            elif 'docker_workdir' not in kwargs:
                # by default, map current working directoryself.
                workdir_opt = f'-w={path(wdir):p}'

            env_opt = ''
            if 'environment' in kwargs:
                if isinstance(kwargs['environment'], dict):
                    env_opt = ' '.join(
                        f'-e {x}={y}' for x, y in kwargs['environment'].items())
                elif isinstance(kwargs['environment'], list):
                    env_opt = ' '.join(f'-e {x}' for x in kwargs['environment'])
                elif isinstance(kwargs['environment'], str):
                    env_opt = f'-e {kwargs["environment"]}'
                else:
                    raise RuntimeError(
                        'Invalid value for option environment (str, list, or dict is allowd, {} provided)'
                        .format(kwargs['environment']))
            #
            port_opt = ''
            if 'port' in kwargs:
                if kwargs['port'] is True:
                    port_opt = '-P'
                elif isinstance(kwargs['port'], (str, int)):
                    port_opt = '-p {}'.format(kwargs['port'])
                elif isinstance(kwargs['port'], list):
                    port_opt = ' '.join(
                        '-p {}'.format(x) for x in kwargs['port'])
                else:
                    raise RuntimeError(
                        'Invalid value for option port (a list of intergers or True), {} provided'
                        .format(kwargs['port']))
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
            user_opt = ''
            if 'user' in kwargs:
                if kwargs['user'] is not None:
                    user_opt = f'-u {kwargs["user"]}'
            elif platform.system() != 'Windows':
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
                security_opt,  # security option
                volumes_opt,  # volumes
                volumes_from_opt,  # volumes_from
                name_opt,  # name
                stdin_opt,  # stdin_optn
                tty_opt,  # tty
                port_opt,  # port
                workdir_opt,  # working dir
                user_opt,  # user
                env_opt,  # environment
                mem_limit_opt,  # memory limit
                extra_opt,  # any extra parameters
                image,  # image
                cmd_opt)
            env.logger.debug(cmd)
            if env.config['run_mode'] == 'dryrun':
                print(f'HINT: {cmd}')
                print(script)
                return 0

            ret = self._run_cmd(cmd, **kwargs)

            if ret != 0:
                debug_script_dir = os.path.join(env.exec_dir, '.sos')
                msg = 'The script has been saved to {}/{}. To reproduce the error please run:\n``{}``'.format(
                    debug_script_dir, tempscript,
                    cmd.replace(f'{path(tempdir):p}',
                                f'{path(debug_script_dir):p}'))
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
                    out = f", stdout={kwargs['stdout']}" if 'stdout' in kwargs and os.path.isfile(
                        kwargs['stdout']) and os.path.getsize(
                            kwargs['stdout']) > 0 else ''
                    err = f", stderr={kwargs['stderr']}" if 'stderr' in kwargs and os.path.isfile(
                        kwargs['stderr']) and os.path.getsize(
                            kwargs['stderr']) > 0 else ''
                    msg = f"Executing script in docker returns an error (exitcode={ret}{err}{out}).\n{msg}"
                raise subprocess.CalledProcessError(
                    returncode=ret,
                    cmd=cmd.replace(tempdir, debug_script_dir),
                    stderr=msg)
        return 0
