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
# Singularity support
#


class SoS_SingularityClient:
    '''A singleton class to ensure there is only one client'''
    _instance = None

    pulled_images = set()

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(SoS_SingularityClient, cls).__new__(cls)
        return cls._instance

    def _ensure_singularity(self):
        if not shutil.which('singularity'):
            raise RuntimeError(f'Command singularity is not found')

    def _is_image_avail(self, image):
        # the command will return ID of the image if it exists
        try:
            return bool(subprocess.check_output(
                f'''Singularity images {image} --no-trunc --format "{{{{.ID}}}}"''',
                shell=True))
        except Exception as e:
            env.logger.warning(f'Failed to check image {image}: {e}')
            return False

    def _run_cmd(self, cmd, **kwargs):
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
                    p = subprocess.Popen(
                        cmd, shell=True, stderr=se, stdout=so)
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
        return ret

    def build(self, script=None, src=None, dest=None, **kwargs):
        self._ensure_singularity()
        if not dest:
            raise ValueError(f'Please specify result of sigularity build with option dest')
        if not script and not src:
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            if script:
                with open(os.path.join(tempdir, 'singularity.def'), 'w') as df:
                    df.write(df)
                file_opt = [dest, os.path.join(tempdir, 'singularity.def')]
            else:
                if not src:
                    raise ValueError(f'Please specify either a script file as script or a source url with option --src')
                file_opt = [dest, src]

            other_opts = []
            for arg, value in kwargs.items():
                # boolean args
                if arg in ('sandbox', 'writable', 'notest', 'checks', 'low', 'med', 'high'):
                    if value is True:
                        other_opts.append(f'--{arg.replace("_", "-")}')
                    else:
                        env.logger.warning(
                            f'Boolean {arg} is ignored (True should be provided)')
                elif arg in ('section', 'tag'):
                    other_opts.extend([f'--{arg.replace("_", "-")}', value])

            cmd = subprocess.list2cmdline(
                ['singularity', 'build'] + other_opts + file_opt)

            env.logger.debug(cmd)
            if env.config['run_mode'] == 'dryrun':
                print(f'HINT: {cmd}')
                print(script)
                return 0

            ret = self._run_cmd(cmd, **kwargs)

            if ret != 0:
                if script:
                    debug_script_dir = os.path.join(env.exec_dir, '.sos)
                    msg = 'The definition has been saved to {}/singularity.def. To reproduce the error please run:\n``{}``'.format(
                        debug_script_dir, cmd.replace(tempdir, path(debug_script_dir)))
                    shutil.copy(os.path.join(
                        tempdir, 'Singularityfile'), debug_script_dir)
                else:
                    msg = f'To reproduce this error please run {cmd}'
                raise subprocess.CalledProcessError(
                    returncode=ret, cmd=cmd, stderr=msg)

    def pull(self, image):
        self._ensure_singularity()

        if image in self.pulled_images:
            return
        # if image is specified, check if it is available locally. If not, pull it
        try:
            print(f'HINT: Pulling singularity image {image}')
            output = subprocess.check_output(
                'singularity pull {}'.format(image), stderr=subprocess.STDOUT, shell=True,
                universal_newlines=True)
            self.pulled_images.add(image)
        except subprocess.CalledProcessError as exc:
            env.logger.warning('Failed to pull {image}: {e}')

    def run(self, image, script='', interpreter='', args='', suffix='.sh', **kwargs):
        self._ensure_singularity()
        #
        env.logger.debug('singularity_run with keyword args {}'.format(kwargs))
        #
        # now, write a temporary file to a tempoary directory under the current directory, this is because
        # we need to share the directory to ...
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            # keep the temporary script for debugging purposes
            # tempdir = tempfile.mkdtemp(dir=os.getcwd())
            if script:
                tempscript = 'singularity_run_{}{}'.format(os.getpid(), suffix)
                with open(os.path.join(tempdir, tempscript), 'w') as script_file:
                    # the input script might have windows new line but the container
                    # will need linux new line for proper execution #1023
                    script_file.write('\n'.join(script.splitlines()))
            #
            # if there is an interpreter and with args
            if not args:
                args = '{filename:pq}'
            #
            # under mac, we by default share /Users within Singularity
            wdir = os.path.abspath(os.getcwd())
            binds = []
            if 'bind' in kwargs:
                volumes = [kwargs['bind']] if isinstance(
                    kwargs['bind'], str) else kwargs['bind']
                has_wdir = False
                for vol in volumes:
                    if not vol:
                        continue
                    if isinstance(vol, (str, path)):
                        vol = str(vol)
                    else:
                        raise ValueError(f'Unacceptable value {vol} for parameter bind')
                    if vol.count(':') == 0:
                        host_dir, mnt_dir = vol, vol
                    elif vol.count(':') in (1, 2):
                        host_dir, mnt_dir = vol.split(':', 1)
                    else:
                        raise ValueError(f'Invalid format for bind specification: {vol}')
                    binds.append(
                        f'{path(host_dir).resolve():p}:{path(mnt_dir):p}')
                    if wdir.startswith(os.path.abspath(os.path.expanduser(host_dir))):
                        has_wdir = True
                bind_opt = ' '.join('-v {}'.format(x) for x in binds)
                if not has_wdir:
                    bind_opt += f' -v {path(wdir):p}:{path(wdir):p}'
            else:
                bind_opt = f' -v {path(wdir):p}:{path(wdir):p}'


            # we also need to mount the script
            if script:
                bind_opt += f' -v {path(tempdir)/tempscript:p}:/var/lib/sos/{tempscript}'
            cmd_opt = interpolate(f'{interpreter if isinstance(interpreter, str) else interpreter[0]} {args}', {
                'filename': sos_targets(f'/var/lib/sos/{tempscript}'),
                'script': script})

            cmd = 'singularity exec {} {} {}'.format(
                bind_opt,        # volumes
                image,              # image
                cmd_opt
            )
            env.logger.debug(cmd)
            if env.config['run_mode'] == 'dryrun':
                print(f'HINT: {cmd}')
                print(script)
                return 0

            ret = self._run_cmd(cmd, **kwargs)

            if ret != 0:
                debug_script_dir = os.path.join(env.exec_dir, '.sos')
                msg = 'The script has been saved to {}/{}. To reproduce the error please run:\n``{}``'.format(
                    debug_script_dir, tempscript, cmd.replace(f'{path(tempdir):p}', f'{path(debug_script_dir):p}'))
                shutil.copy(os.path.join(tempdir, tempscript),
                            debug_script_dir)
                out = f", stdout={kwargs['stdout']}" if 'stdout' in kwargs and os.path.isfile(kwargs['stdout']) and os.path.getsize(kwargs['stdout']) > 0 else ''
                err = f", stderr={kwargs['stderr']}" if 'stderr' in kwargs and os.path.isfile(kwargs['stderr']) and os.path.getsize(kwargs['stderr']) > 0 else ''
                msg = f"Executing script in Singularity returns an error (exitcode={ret}{err}{out}).\n{msg}"
                raise subprocess.CalledProcessError(
                    returncode=ret,
                    cmd=cmd.replace(tempdir, debug_script_dir),
                    stderr=msg)
        return 0
