#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
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
import shlex
import json
import platform
import urllib
import blessings
from io import BytesIO
from docker import Client
from docker.utils import kwargs_from_env
from shutil import which
import multiprocessing as mp
from .utils import env, interpolate, glob_wildcards, downloadURL, fileMD5

__all__ = ['SoS_Action', 'execute_script', 'sos_run',
    'check_command', 'fail_if', 'warn_if', 'download',
    'run', 'bash', 'csh', 'tcsh', 'zsh', 'sh',
    'python', 'python3',
    'perl', 'ruby', 'node', 'JavaScript',
    'R', 'check_R_library',
    'docker_build', 'docker_commit'
    ]

from .sos_syntax import SOS_RUNTIME_OPTIONS
from .sos_script import SoS_Script

#
# docker support
#
class DockerClient:
    '''A singleton class to ensure there is only one client'''
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(DockerClient, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        self.client = Client(**kwargs_from_env(assert_hostname=False))
        try:
            self.client.info()
        except Exception as e:
            self.client = None

    def _is_image_avail(self, image):
        images = sum([x['RepoTags'] for x in self.client.images()], [])
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
            for line in self.client.build(fileobj=f, **kwargs):
                self.stream(line.decode())
        else:
            for line in self.client.build(**kwargs):
                self.stream(line.decode())
        # if a tag is given, check if the image is built
        if 'tag' in kwargs and not self._is_image_avail(kwargs['tag']):
            raise RuntimeError('Image with tag {} is not created.'.format(kwargs['tag']))

    def import_image(self, image, **kwargs):
        if not self.client:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        env.logger.info('docker import {}'.format(image))
        self.client.import_image(image, **kwargs)
       
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

    def run(self, image, script='', interpreter='', suffix='.sh', **kwargs):
        if self.client is None:
            raise RuntimeError('Cannot connect to the Docker daemon. Is the docker daemon running on this host?')
        env.logger.debug('docker_run with keyword args {}'.format(kwargs))
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
                        if not os.path.abspath(host_dir).startswith('/Users'):
                            raise RuntimeError('hostdir ({}) under MacOSX must be under /Users to be usable in docker container'.format(host_dir))
                    binds.append('{}:{}'.format(os.path.abspath(host_dir), mnt_dir))
            # we also need to mount the script
            if script and interpreter:
                binds.append('{}:{}'.format(os.path.join(tempdir, tempscript), '/var/lib/sos/{}'.format(tempscript)))
                cmd = interpreter.replace('{}', '/var/lib/sos/{}'.format(tempscript))
            else:
                cmd = ''
            command = 'docker run -t --rm {} {} {}'.format(' '.join('-v ' +x for x in binds), image, cmd)
            env.logger.info(command)
            ret = subprocess.call(command, shell=True)
            if ret != 0:
                raise RuntimeError('Executing script in docker returns an error')
        return 0

#
# A decoration function that allows SoS to replace all SoS actions
# with a null action.
#
def SoS_Action(run_mode='run'):
    run_mode = [run_mode] if isinstance(run_mode, str) else run_mode
    def runtime_decorator(func):
        def action_wrapper(*args, **kwargs):
            runtime_options = env.sos_dict.get('_runtime', {})
            #
            non_runtime_options = {k:v for k,v in kwargs.items() if k not in SOS_RUNTIME_OPTIONS}
            # docker files will be downloaded in run or prepare mode
            if 'docker_file' in runtime_options and env.run_mode in ('run', 'prepare'):
                docker = DockerClient()
                docker.import_image(runtime_options['docker_file'])
            # handle image
            if 'docker_image' in runtime_options and env.run_mode in ('run', 'prepare'):
                docker = DockerClient()
                docker.pull(runtime_options['docker_image'])
            if env.run_mode not in run_mode:
                return 0
            return func(*args, **non_runtime_options)
        action_wrapper.run_mode = run_mode
        return action_wrapper
    return runtime_decorator


class SoS_ExecuteScript:
    def __init__(self, script, interpreter, suffix):
        self.script = script
        self.interpreter = interpreter
        self.suffix = suffix

    def run(self, **kwargs):
        if '{}' not in self.interpreter:
            self.interpreter += ' {}'
        runtime_options = env.sos_dict.get('_runtime', {})
        if 'docker_image' in runtime_options:
            docker = DockerClient()
            docker.run(runtime_options['docker_image'], self.script, self.interpreter, self.suffix, 
                volumes=runtime_options.get('docker_volumes', []), **kwargs)
        else:
            self.script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=self.suffix, delete=False).name
            with open(self.script_file, 'w') as script_file:
                script_file.write(self.script)
            cmd = self.interpreter.replace('{}', shlex.quote(self.script_file))
            try:
                p = subprocess.Popen(cmd, shell=True)
                env.register_process(p.pid, 'Runing {}'.format(self.script_file))
                ret = p.wait()
            finally:
                env.deregister_process(p.pid)
            if ret != 0:
                raise RuntimeError('Failed to execute script')


@SoS_Action(run_mode=['dryrun', 'prepare', 'run'])
def sos_run(workflow, source=None):
    '''Execute a workflow from specified source, input, and output
    By default the workflow is defined in the existing SoS script, but
    extra sos files can be specified from paramter source. The workflow
    will be execute in the current step namespace with _input as workflow
    input. '''
    env.logger.info('Executing nested workflow {}'.format(workflow))
    script = SoS_Script(env.sos_dict['__step_context__'].content, env.sos_dict['__step_context__'].filename)
    wf = script.workflow(workflow, source=source)
    # if wf contains the current step or one of the previous one, this constitute
    # recusive nested workflow and should not be allowed
    if env.sos_dict['step_name'] in ['{}_{}'.format(x.name, x.index) for x in wf.sections if not x.is_parameters]:
        raise RuntimeError('Nested workflow {} contains the current step {}'.format(workflow, env.sos_dict['step_name']))
    return wf.run(args=env.sos_dict['__args__'], nested=True)

@SoS_Action(run_mode=['run'])
def execute_script(script, interpreter, suffix, **kwargs):
    return SoS_ExecuteScript(script, interpreter, suffix).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'prepare', 'run'])
def check_command(cmd, pattern = None):
    '''Raise an exception if output of `cmd` does not match specified `pattern`.
    Multiple patterns can be specified as a list of patterns.
    When pattern is None, check the existence of command `cmd`
    and raise an error if command does not exist.'''
    ret_val = 0
    if pattern is None and len(shlex.split(cmd)) == 1:
        name = which(cmd)
        if not name:
            raise RuntimeError('Command ``{}`` not found!'.format(cmd))
        env.logger.info('Command ``{}`` is located as ``{}``.'.format(cmd, name))
    else:
        try:
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, timeout=2).decode()
        except subprocess.TimeoutExpired as e:
            output = e.output
            ret_val = 1
            env.logger.warning(e)
            env.logger.warning(e.output.decode())
        except subprocess.CalledProcessError as e:
            ret_val = e.returncode
            output = e.output
            env.logger.warning(e)
            env.logger.warning(e.output.decode())
        #
        env.logger.trace('Output of command ``{}`` is ``{}``'.format(cmd, output))
        #
        if pattern:
            pattern = [pattern] if isinstance(pattern, str) else pattern
            if all([re.search(x, output, re.MULTILINE) is None for x in pattern]):
                raise RuntimeError('Output of command ``{}`` does not match specified regular expression ``{}``.'
                    .format(cmd, ' or '.join(pattern)))
    return ret_val

@SoS_Action(run_mode=['dryrun', 'prepare', 'run'])
def fail_if(expr, msg=''):
    '''Raise an exception with `msg` if condition `expr` is False'''
    if expr:
        raise RuntimeError(msg)
    return 0

@SoS_Action(run_mode=['dryrun', 'prepare', 'run'])
def warn_if(expr, msg=''):
    '''Yield an warning message `msg` if `expr` is False '''
    if expr:
        env.logger.warning(msg)
    return 0

@SoS_Action(run_mode=['prepare', 'run'])
def download(url_list, dest_dir='.', dest_file=None, decompress=False):
    '''Download files from specified URL, which should be space, tab or
    newline separated URLs. The files will be downloaded to specified
    destination. If `filename.md5` files are downloaded, they are used to 
    validate downloaded `filename`. Unless otherwise specified, compressed
    files are decompressed.
    '''
    urls = [x.strip() for x in url_list.split() if x.strip()]
    if dest_file is not None and len(urls) != 1:
        raise RuntimeError('Only one URL is allowed if a destination file is specified.')
    # 
    if dest_file is None:
        filenames = []
        for idx, url in enumerate(urls):
            token = urllib.parse.urlparse(url)
            # if no scheme or netloc, the URL is not acceptable
            if not all([getattr(token, qualifying_attr) for qualifying_attr in  ('scheme', 'netloc')]):
                filenames.append(None)
                continue
            filename = os.path.split(token.path)[-1]
            if not filename:
                filenames.append(None)
                continue
            filenames.append(os.path.join(dest_dir, filename))
    else:
        filenames = [dest_file]
    #
    succ = [False for x in urls]
    if len(succ) > 1:
        # first scroll several lines to reserve place for progress bar
        for url in urls:
            sys.stderr.write('\n')
        with mp.Pool(processes = env.sos_dict['CONFIG'].get('sos_download_processes', 5)) as pool:
            for idx, (url, filename) in enumerate(zip(urls, filenames)):
                if not filename:
                    continue
                succ[idx] = pool.apply_async(downloadURL, (url, filename,
                    decompress, len(urls) - idx))
            succ = [x.get() if isinstance(x, mp.pool.AsyncResult) else x for x in succ]
        #
        t = blessings.Terminal(stream=sys.stderr)
        sys.stderr.write(t.move( t.height, 0)) # + '\n')
    else:
        if dest_file is not None:
            succ[0] = downloadURL(urls[0], dest_file, decompress=decompress)
        else:
           if filenames[0]:
                succ[0] = downloadURL(urls[0], filenames[0], decompress=decompress)
    #
    for su, url in zip(succ, urls):
        if not su:
            env.logger.warning('Failed to download {}'.format(url))
    if not all(succ):
        raise RuntimeError('Not all files have been downloaded')
    # if downloaded files contains .md5 signature, use them to validate
    # downloaded files.
    for filename in [x for x in filenames if x.endswith('.md5')]:
        if filename[:-4] in filenames and os.path.isfile(filename[:-4]):
            with open(filename) as md5:
                rec_md5 = md5.readline().split()[0].strip()
                obs_md5 = fileMD5(filename[:-4], partial=False)
                if rec_md5 != obs_md5:
                    env.logger.warning('md5 signature mismatch for downloaded file {} (recorded {}, observed {})'
                        .format(filename[:-4], rec_md5, obs_md5))
                    return 1
    return 0

@SoS_Action(run_mode='run')
def run(script, **kwargs):
    return SoS_ExecuteScript(script, '/bin/bash', '.sh').run(**kwargs)

@SoS_Action(run_mode='run')
def bash(script, **kwargs):
    return SoS_ExecuteScript(script, '/bin/bash', '.sh').run(**kwargs)

@SoS_Action(run_mode='run')
def csh(script, **kwargs):
    return SoS_ExecuteScript(script, '/bin/csh', '.csh').run(**kwargs)

@SoS_Action(run_mode='run')
def tcsh(script, **kwargs):
    return SoS_ExecuteScript(script, '/bin/tcsh', '.sh').run(**kwargs)

@SoS_Action(run_mode='run')
def zsh(script, **kwargs):
    return SoS_ExecuteScript(script, '/bin/zsh', '.zsh').run(**kwargs)

@SoS_Action(run_mode='run')
def sh(script, **kwargs):
    return SoS_ExecuteScript(script, '/bin/sh', '.sh').run(**kwargs)

@SoS_Action(run_mode='run')
def python(script, **kwargs):
    return SoS_ExecuteScript(script, 'python', '.py').run(**kwargs)

@SoS_Action(run_mode='run')
def python3(script, **kwargs):
    return SoS_ExecuteScript(script, 'python3', '.py').run(**kwargs)

@SoS_Action(run_mode='run')
def perl(script, **kwargs):
    return SoS_ExecuteScript(script, 'perl', '.pl').run(**kwargs)

@SoS_Action(run_mode='run')
def ruby(script, **kwargs):
    return SoS_ExecuteScript(script, 'ruby', '.rb').run(**kwargs)

@SoS_Action(run_mode='run')
def node(script, **kwargs):
    return SoS_ExecuteScript(script, 'node', '.js').run(**kwargs)

@SoS_Action(run_mode='run')
def JavaScript(script, **kwargs):
    return SoS_ExecuteScript(script, 'node', '.js').run(**kwargs)

@SoS_Action(run_mode='run')
def R(script, **kwargs):
    return SoS_ExecuteScript(script, 'Rscript --default-packages=methods,utils,stats', '.R').run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'prepare', 'run'])
def check_R_library(name, version = None):
    '''Check existence and version match of R library.
    cran and bioc packages are unique yet might overlap with github.
    Therefore if the input name is {repo}/{pkg} the package will be
    installed from github if not available, else from cran or bioc
    '''
    output_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False).name
    if len(glob_wildcards('{repo}/{pkg}', [name])['repo']):
        # package is from github
        check_R_library('devtools')
        install_script = interpolate('''
        options(warn=-1)
        package_repo <- ${name!r}
        package <- basename(package_repo)
        if (require(package, character.only=TRUE, quietly=TRUE)) {
            write(paste(package, packageVersion(package), "AVAILABLE"), file="${output_file}")
        } else {
            devtools::install_github(package_repo)
            # if it still does not exist, write the package name to output
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "INSTALLED"), file="${output_file}")
            } else {
                write(paste(package, "NA", "MISSING"), file="${output_file}")
                quit("no")
            }
        }
        cur_version <- packageVersion(package)
        ''', '${ }', locals())
    else:
        # package is from cran or bioc
        install_script = interpolate('''
        options(warn=-1)
        package <- ${name!r}
        if (require(package, character.only=TRUE, quietly=TRUE)) {
            write(paste(package, packageVersion(package), "AVAILABLE"), file="${output_file}")
        } else {
            install.packages(package, repos="http://cran.us.r-project.org",
                quiet=FALSE)
            # if the package still does not exist
            if (!require(package, character.only=TRUE, quietly=TRUE)) {
                source("http://bioconductor.org/biocLite.R")
                biocLite(package, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
            }
            # if it still does not exist, write the package name to output
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "INSTALLED"), file="${output_file}")
            } else {
                write(paste(package, "NA", "MISSING"), file="${output_file}")
                quit("no")
            }
        }
        cur_version <- packageVersion(package)
        ''', '${ }', locals())
    version_script = ''
    if version is not None:
        version = [version] if isinstance(version, str) else version
        operators = []
        for idx, value in enumerate(version):
            value = str(value)
            if value.endswith('+'):
                operators.append('>=')
                version[idx] = value[:-1]
            elif value.endswith('-'):
                operators.append('<')
                version[idx] = value[:-1]
            else:
                operators.append('==')
        # check version and mark version mismatch
        # if current version satisfies any of the
        # requirement the check program quits
        for x, y in zip(version, operators):
            version_script += '''
            if (cur_version {1} {0}) {{
              quit("no")
            }}
            '''.format(repr(x), y)
        version_script += 'write(paste(package, cur_version, "VERSION_MISMATCH"), file = {})'.\
          format(repr(output_file))
    SoS_ExecuteScript(install_script + version_script, 'Rscript --default-packages=methods,utils,stats', '.R').run()
    ret_val = 0
    with open(output_file) as tmp:
        for line in tmp:
            lib, version, status = line.split()
            if status.strip() == "MISSING":
                raise RuntimeError('R Library {} is not available and cannot be installed.'.format(lib))
            elif status.strip() == 'AVAILABLE':
                env.logger.info('R library {} ({}) is available'.format(lib, version))
            elif status.strip() == 'INSTALLED':
                env.logger.info('R library {} ({}) has been installed'.format(lib, version))
            elif status.strip() == 'VERSION_MISMATCH':
                env.logger.warning('R library {} ({}) does not satisfy version requirement!'.format(lib, version))
                ret_val = 1
            else:
                raise RuntimeError('This should not happen: {}'.format(line))
    try:
        os.remove(self.output_file)
    except:
        pass
    return ret_val


@SoS_Action(run_mode='run')
def docker_build(dockerfile=None, **kwargs):
    '''docker build command. By default a script is sent to the docker build command but
    you can also specify different parameters defined inu//docker-py.readthedocs.org/en/stable/api/#build
    '''
    docker = DockerClient()
    docker.build(dockerfile, **kwargs)
    return 0


@SoS_Action(run_mode='run')
def docker_commit(**kwargs):
    docker = DockerClient()
    docker.commit(**kwargs)
    return 0

