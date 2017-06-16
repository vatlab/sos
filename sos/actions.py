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
import subprocess
import tempfile
import copy
import urllib
import urllib.request
import urllib.error
import urllib.parse
import shutil
import zipfile
import shlex
import gzip
import tarfile
import fasteners
import uuid

from functools import wraps

from collections.abc import Sequence
import multiprocessing as mp
from tqdm import tqdm as ProgressBar
from .utils import env, transcribe, AbortExecution, short_repr, get_traceback
from .sos_eval import Undetermined, interpolate
from .target import FileTarget, fileMD5, executable, UnknownTarget, BaseTarget

__all__ = ['SoS_Action', 'execute_script', 'sos_run',
    'fail_if', 'warn_if', 'stop_if',
    'download',
    'run', 'perl', 'ruby', 'node', 'JavaScript',
    'report', 'pandoc'
    ]

from .sos_syntax import SOS_RUNTIME_OPTIONS, SOS_ACTION_OPTIONS
from .sos_script import SoS_Script

def get_actions():
    # get the name of all actions, which are identified by an attribute
    # run_mode of the function
    return [k for k, v in globals().items() if hasattr(v, 'run_mode')]

#
# A decoration function that allows SoS to replace all SoS actions
# with a null action.
#
def SoS_Action(run_mode=['run', 'interactive'], acceptable_args=['*']):
    run_mode = [run_mode] if isinstance(run_mode, str) else run_mode
    def runtime_decorator(func):
        @wraps(func)
        def action_wrapper(*args, **kwargs):
            # if docker_image in args, a large number of docker-specific
            # args would be allowed.
            if '*' not in acceptable_args and 'docker_image' not in kwargs:
                for key in kwargs.keys():
                    if key not in acceptable_args and key not in SOS_ACTION_OPTIONS:
                        raise ValueError('Unrecognized option "{}" for action {}'.format(key, func))
            # docker files will be downloaded in run or prepare mode
            if 'docker_file' in kwargs and env.config['run_mode'] in ['run', 'interactive']:
                from .docker.client import SoS_DockerClient
                docker = SoS_DockerClient()
                docker.load_image(kwargs['docker_file'])
            # handle image
            if 'docker_image' in kwargs:
                from .docker.client import SoS_DockerClient
                docker = SoS_DockerClient()
                docker.pull(kwargs['docker_image'])
            if env.config['run_mode'] not in run_mode:
                # return dynamic expression when not in run mode, that is to say
                # the script logic cannot rely on the result of the action
                return Undetermined(func.__name__)
            if env.config['run_mode'] == 'interactive':
                for k,v in kwargs.items():
                    if k in SOS_RUNTIME_OPTIONS and k not in SOS_ACTION_OPTIONS:
                        env.logger.warning('Passing runtime option "{0}" to action is deprecated. Please use "task: {0}={1}" before action instead.'.format(k, v))
            if 'input' in kwargs and '__local_input__' in env.sos_dict:
                if isinstance(kwargs['input'], (str, BaseTarget)):
                    env.sos_dict['__local_input__'].append(kwargs['input'])
                elif isinstance(kwargs['input'], Sequence):
                    env.sos_dict['__local_input__'].extend(kwargs['input'])
                else:
                    env.logger.warning('Parameter input of action should be string or list of strings')
                for item in env.sos_dict['__local_input__']:
                    if isinstance(item, str):
                        if not FileTarget(item).exists('target'):
                            FileTarget(item).remove_sig()
                            raise UnknownTarget(item)
                    elif not item.exists('target'):
                        item.remove_sig()
                        raise UnknownTarget(item)
            if 'output' in kwargs and '__local_output__' in env.sos_dict:
                if isinstance(kwargs['output'], (str, BaseTarget)):
                    env.sos_dict['__local_output__'].append(kwargs['output'])
                elif isinstance(kwargs['output'], Sequence):
                    env.sos_dict['__local_output__'].extend(kwargs['output'])
                else:
                    env.logger.warnoutg('Parameter output of action should be stroutg or list of strings')
            if 'active' in kwargs:
                if isinstance(kwargs['active'], int):
                    if kwargs['active'] >= 0 and env.sos_dict['_index'] != kwargs['active']:
                        return None
                    if kwargs['active'] < 0 and env.sos_dict['_index'] != kwargs['active'] + env.sos_dict['__num_groups__']:
                        return None
                elif isinstance(kwargs['active'], Sequence):
                    allowed_index = list([x if x >= 0 else env.sos_dict['__num_groups__'] + x for x in kwargs['active']])
                    if env.sos_dict['_index'] not in allowed_index:
                        return None
                elif isinstance(kwargs['active'], slice):
                    allowed_index = list(range(env.sos_dict['__num_groups__']))[kwargs['active']]
                    if env.sos_dict['_index'] not in allowed_index:
                        return None
                else:
                    raise RuntimeError('Unacceptable value for option active: {}'.format(kwargs['active']))
            if 'workdir' in kwargs:
                if not kwargs['workdir'] or not isinstance(kwargs['workdir'], str):
                    raise RuntimeError('workdir option should be a path, {} provided'.format(kwargs['workdir']))
                if not os.path.isdir(os.path.expanduser(kwargs['workdir'])):
                    os.makedirs(os.path.expanduser(kwargs['workdir']))
                try:
                    olddir = os.getcwd()
                    os.chdir(os.path.expanduser(kwargs['workdir']))
                    try:
                        res = func(*args, **kwargs)
                    except Exception as e:
                        if 'allow_error' in kwargs and kwargs['allow_error']:
                            env.logger.warning(e)
                            res = None
                        else:
                            raise
                finally:
                    os.chdir(olddir)
            else:
                try:
                    res = func(*args, **kwargs)
                except Exception as e:
                    if 'allow_error' in kwargs and kwargs['allow_error']:
                        env.logger.warning(e)
                        res = None
                    else:
                        raise
            if '__local_output__' in env.sos_dict:
                for item in env.sos_dict['__local_output__']:
                    if isinstance(item, str):
                        if not FileTarget(item).exists():
                            raise RuntimeError("Target {} does not exist after execution of action".format(item))
                    elif not item.exists():
                            raise RuntimeError("Target {} does not exist after execution of action".format(item))

            return res
        action_wrapper.run_mode = run_mode
        return action_wrapper
    return runtime_decorator


class SoS_ExecuteScript:
    def __init__(self, script, interpreter, suffix, args=''):
        self.script = script
        self.interpreter = interpreter
        self.args = args
        if suffix:
            self.suffix = suffix
        elif sys.platform == 'win32':
            self.suffix = '.bat'
        else:
            self.suffix = '.sh'

    def run(self, **kwargs):
        if 'docker_image' in kwargs:
            if env.config['run_mode'] == 'dryrun':
                print('In docker image {}\n{}:\n{}\n'.format(kwargs['docker_image'], self.interpreter, self.script))
                return None
            from .docker.client import SoS_DockerClient
            docker = SoS_DockerClient()
            docker.run(kwargs['docker_image'], self.script, self.interpreter, self.args, self.suffix,
                **kwargs)
        else:
            if isinstance(self.interpreter, str):
                if self.interpreter and not shutil.which(shlex.split(self.interpreter)[0]):
                    raise RuntimeError('Failed to locate interpreter {}'.format(self.interpreter))
            elif isinstance(self.interpreter, Sequence):
                found = False
                for ip in self.interpreter:
                    if shutil.which(shlex.split(ip)[0]):
                        self.interpreter = ip
                        found = True
                        break
                if not found:
                    raise RuntimeError('Failed to locate any of the interpreters {}'
                        .format(', '.join(self.interpreter)))
            else:
                raise RuntimeError('Unacceptable interpreter {}'.format(self.interpreter))

            transcribe(self.script, action=self.interpreter)
            debug_script_file = os.path.join(env.exec_dir, '.sos', '{}_{}_{}{}'.format(env.sos_dict['step_name'],
                env.sos_dict['_index'], str(uuid.uuid4())[:8], self.suffix))
            env.logger.debug('Script for step {} is saved to {}'.format(env.sos_dict['step_name'], debug_script_file))
            with open(debug_script_file, 'w') as sfile:
                sfile.write(self.script)
            env.logger.trace(self.script)

            try:
                p = None
                script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=self.suffix, delete=False).name
                with open(script_file, 'w') as sfile:
                    sfile.write(self.script)
                if self.interpreter:
                    # if there is an interpreter and with args
                    if not self.args:
                        self.args = '${filename!q}'
                else:
                    if sys.platform == 'win32':
                        # in the case there is no interpreter, we put the script
                        # at first (this is the case for windows)
                        #
                        # and we donot add default args.
                        self.interpreter = '${filename!q}'
                    else:
                        # if there is a shebang line, we ...
                        if self.script.startswith('#!'):
                            # make the script executable
                            os.chmod(script_file, 0o775)
                            self.interpreter = '${filename!q}'
                        else:
                            self.interpreter = '/bin/bash'
                            if not self.args:
                                self.args = '-ev ${filename!q}'
                #
                if env.config['run_mode'] == 'dryrun':
                    print('{}:\n{}\n'.format(self.interpreter, self.script))
                    return None
                cmd = interpolate('{} {}'.format(self.interpreter, self.args), '${ }', {'filename': script_file, 'script': self.script})
                #
                if env.config['run_mode'] == 'interactive':
                    # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
                    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                    out, err = p.communicate()
                    if out and env.verbosity > 1:
                        sys.stdout.write(out.decode())
                    if err and env.verbosity > 0:
                        sys.stderr.write(err.decode())
                    ret = p.returncode
                    sys.stdout.flush()
                    sys.stderr.flush()
                elif '__std_out__' in env.sos_dict and '__std_err__' in env.sos_dict:
                    if env.verbosity > 1:
                        with open(env.sos_dict['__std_out__'], 'ab') as so, open(env.sos_dict['__std_err__'], 'ab') as se:
                            p = subprocess.Popen(cmd, shell=True, stderr=se, stdout=so)
                    else:
                        p = subprocess.Popen(cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                    ret = p.wait()
                else:
                    p = subprocess.Popen(cmd, shell=True,
                                         stderr=None if env.verbosity > 0 else subprocess.DEVNULL,
                                         stdout=None if env.verbosity > 1 else subprocess.DEVNULL)
                    ret = p.wait()
                if ret != 0:
                    with open(debug_script_file, 'w') as sfile:
                        sfile.write(self.script)
                    if self.interpreter == '/bin/bash' and self.args == '-ev ${filename!q}':
                        debug_args = '${filename!q}'
                    else:
                        debug_args = self.args
                    cmd = interpolate('{} {}'.format(self.interpreter, debug_args), '${ }', {'filename': debug_script_file, 'script': self.script})
                    raise RuntimeError('Failed to execute script (ret={}).\nPlease use command\n\t``{}``\nunder "{}" to test it.'
                        .format(ret, ' \\\n\t  '.join(cmd.split()), os.getcwd()))
            except RuntimeError:
                raise
            except Exception as e:
                env.logger.error('Failed to execute script: {}'.format(e))
                raise
            finally:
                os.remove(script_file)


@SoS_Action(run_mode=['run', 'interactive'])
def sos_run(workflow=None, targets=None, shared=[], args={}, **kwargs):
    '''Execute a workflow from specified source, input, and output
    By default the workflow is defined in the existing SoS script, but
    extra sos files can be specified from paramter source. The workflow
    will be execute in the current step namespace with _input as workflow
    input. '''
    script = SoS_Script(env.sos_dict['__step_context__'].content, env.sos_dict['__step_context__'].filename)
    wf = script.workflow(workflow, use_default=not targets)
    # if wf contains the current step or one of the previous one, this constitute
    # recusive nested workflow and should not be allowed
    if env.sos_dict['step_name'] in ['{}_{}'.format(x.name, x.index) for x in wf.sections]:
        raise RuntimeError('Nested workflow {} contains the current step {}'.format(workflow, env.sos_dict['step_name']))
    # args can be specified both as a dictionary or keyword arguments
    args.update(kwargs)
    args['__args__'] = env.sos_dict['__args__']
    if isinstance(shared, str):
        shared = [shared]

    # for nested workflow, _input would becomes the input of workflow.
    env.sos_dict.set('__step_output__', copy.deepcopy(env.sos_dict.get('_input', None)))
    shared.append('__step_output__')
    try:
        my_name = env.sos_dict['step_name']
        args_output = ', '.join('{}={}'.format(x, short_repr(y)) for x,y in args.items() if not x.startswith('__'))
        env.logger.info('Executing workflow ``{}`` with input ``{}`` and {}'
            .format(workflow, short_repr(env.sos_dict.get('_input', None), True),
            'no args' if not args_output else args_output))

        if not hasattr(env, '__pipe__'):
            # if env has no __pipe__, this means the nested workflow is executed from 
            # within a task, and we can just use an executor to execute it.
            #         # tell the master process to receive a workflow
            shared = {x: (env.sos_dict[x] if x in env.sos_dict else None) for x in shared}
            if env.config['run_mode'] == 'run':
                from sos.sos_executor import Base_Executor
                executor = Base_Executor(wf, args=args, shared=shared, config=env.config)
                if shared:
                    q = mp.Pipe()
                else:
                    q = None
                p = mp.Process(target=executor.run, kwargs={'targets': targets, 'parent_pipe': q[1], 'my_workflow_id': None})
                p.start()
                if shared:
                    res = q[0].recv()
                    if res is None:
                        # worker is killed
                        sys.exit(0)
                    elif isinstance(res, Exception):
                        raise res
                    env.sos_dict.quick_update(res['shared'])
                else:
                    res = None
                p.join()
            else:
                from sos.jupyter.sos_executor import Interactive_Executor
                executor = Interactive_Executor(wf, args=args, shared=shared, config=env.config)
                res = executor.run(targets=targets)
            return res
        
        else:
            # tell the master process to receive a workflow
            env.__pipe__.send('workflow {}'.format(uuid.uuid4()))
            # really send the workflow
            shared = {x: (env.sos_dict[x] if x in env.sos_dict else None) for x in shared}
            env.__pipe__.send((wf, targets, args, shared, env.config))
            res = env.__pipe__.recv()
            if res is None:
                sys.exit(0)
            elif isinstance(res, Exception):
                raise res
            else:
                env.sos_dict.quick_update(res['shared'])
                return res
    finally:
        # restore step_name in case the subworkflow re-defines it
        env.sos_dict.set('step_name', my_name)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script', 'interpreter', 'suffix', 'args'])
def execute_script(script, interpreter, suffix, args='', **kwargs):
    '''Execute specified script using specified interpreter.'''
    if env.config['run_mode'] == 'dryrun':
        print('{}:\n{}\n'.format(interpreter, script))
        return None
    return SoS_ExecuteScript(script, interpreter, suffix, args).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['expr', 'msg'])
def fail_if(expr, msg=''):
    '''Raise an exception with `msg` if condition `expr` is False'''
    if expr:
        raise RuntimeError(msg)
    return 0

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['expr', 'msg'])
def warn_if(expr, msg=''):
    '''Yield an warning message `msg` if `expr` is False '''
    if expr:
        env.logger.warning(msg)
    return 0

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['expr', 'msg'])
def stop_if(expr, msg=''):
    '''Abort the execution of the current step or loop and yield
    an warning message `msg` if `expr` is False '''
    if expr:
        raise AbortExecution(msg)
    return 0

#
# download file with progress bar
#
def downloadURL(URL, dest, decompress=False, index=None):
    dest = os.path.abspath(os.path.expanduser(dest))
    dest_dir, filename = os.path.split(dest)
    #
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    if not os.path.isdir(dest_dir):
        raise RuntimeError('Failed to create destination directory to download {}'.format(URL))
    #
    message = filename
    if len(message) > 30:
        message = message[:10] + '...' + message[-16:]
    #
    dest_tmp = dest + '.tmp_{}'.format(os.getpid())
    term_width = shutil.get_terminal_size((80, 20)).columns
    try:
        env.logger.debug('Download {} to {}'.format(URL, dest))
        sig = FileTarget(dest)
        if os.path.isfile(dest):
            prog = ProgressBar(desc=message + ': \033[32m validating\033[0m', disable=env.verbosity <= 1,
                position=index, leave=True, bar_format='{desc}', total=10000000)
            if env.config['sig_mode'] == 'build':
                if decompress:
                    prog.set_description(message + ': \033[32m scanning decompressed files\033[0m')
                    prog.update()
                    if zipfile.is_zipfile(dest):
                        zip = zipfile.ZipFile(dest)
                        names = zip.namelist()
                        for name in names:
                            dest_file = os.path.join(dest_dir, name)
                            if not os.path.isfile(dest_file):
                                env.logger.warning('Missing decompressed file {}'.format(dest_file))
                            else:
                                sig.add(dest_file)
                    elif tarfile.is_tarfile(dest):
                        with tarfile.open(dest, 'r:*') as tar:
                            # only extract files
                            file_count = 0
                            for tarinfo in tar:
                                if tarinfo.isfile():
                                    dest_file = os.path.join(dest_dir, tarinfo.name)
                                    if not os.path.isfile(dest_file):
                                        env.logger.warning('Missing decompressed file {}'.format(dest_file))
                                    else:
                                        sig.add(dest_file)
                                        file_count += 1
                                # sometimes the file is very large but we do not need to decompress all
                                # and track all files.
                                if file_count > 10:
                                    break
                    elif dest.endswith('.gz'):
                        decomp = dest[:-3]
                        if not os.path.isfile(decomp):
                            env.logger.warning('Missing decompressed file {}'.format(decomp))
                        sig.add(decomp)
                prog.set_description(message + ': \033[32m writing signature\033[0m')
                prog.update()
                sig.write_sig()
                prog.set_description(message + ': \033[32m signature calculated\033[0m')
                prog.update()
                prog.close()
                return True
            elif env.config['sig_mode'] == 'ignore':
                prog.set_description(message + ': \033[32m use existing\033[0m')
                prog.update()
                prog.close()
                return True
            else:
                prog.update()
                if sig.validate():
                    prog.set_description(message + ': \033[32m Validated\033[0m')
                    prog.update()
                    prog.close()
                    return True
                else:
                    prog.set_description(message + ':\033[91m Signature mismatch\033[0m')
                    prog.update()
        else:
            prog = ProgressBar(desc=message, disable=env.verbosity <= 1, position=index,
                leave=True, bar_format='{desc}', total=10000000)

        #
        # Stop using pycurl because of libcurl version compatibility problems
        # that happen so often and difficult to fix. Error message looks like
        #
        # Reason: Incompatible library version: pycurl.cpython-35m-darwin.so
        # requires version 9.0.0 or later, but libcurl.4.dylib provides version 7.0.0
        #
        #with open(dest_tmp, 'wb') as f:
        #    c = pycurl.Curl()
        #    c.setopt(pycurl.URL, str(URL))
        #    c.setopt(pycurl.WRITEFUNCTION, f.write)
        #    c.setopt(pycurl.SSL_VERIFYPEER, False)
        #    c.setopt(pycurl.NOPROGRESS, False)
        #    c.setopt(pycurl.PROGRESSFUNCTION, prog.curlUpdate)
        #    c.perform()
        #if c.getinfo(pycurl.HTTP_CODE) == 404:
        #    prog.set_description(message + ':\033[91m 404 Error {}\033[0m'.format(' '*(term_width - len(message) - 12)))
        #    try:
        #        os.remove(dest_tmp)
        #    except OSError:
        #        pass
        #    return False
        with open(dest_tmp, 'wb') as f:
            try:
                u = urllib.request.urlopen(str(URL))
                try:
                    file_size = int(u.getheader("Content-Length"))
                    prog = ProgressBar(total=file_size, desc=message, position=index, leave=False)
                except:
                    file_size = None
                file_size_dl = 0
                block_sz = 8192
                while True:
                    buffer = u.read(block_sz)
                    if not buffer:
                        break
                    file_size_dl += len(buffer)
                    f.write(buffer)
                    prog.update(len(buffer))
            except urllib.error.HTTPError as e:
                prog.set_description(message + ':\033[91m {} Error\033[0m'.format(e.code))
                prog.update()
                prog.close()
                try:
                    os.remove(dest_tmp)
                except OSError:
                    pass
                return False
            except Exception as e:
                prog.set_description(message + ':\033[91m {}\033[0m'.format(e))
                prog.update()
                prog.close()
                try:
                    os.remove(dest_tmp)
                except OSError:
                    pass
                return False
        #
        os.rename(dest_tmp, dest)
        decompressed = 0
        if decompress:
            if zipfile.is_zipfile(dest):
                prog.set_description(message + ':\033[91m Decompressing\033[0m')
                prog.update()
                prog.close()
                zip = zipfile.ZipFile(dest)
                zip.extractall(dest_dir)
                names = zip.namelist()
                for name in names:
                    if not os.path.isfile(os.path.join(dest_dir, name)):
                        return False
                    else:
                        sig.add(os.path.join(dest_dir, name))
                        decompressed += 1
            elif tarfile.is_tarfile(dest):
                prog.set_description(message + ':\033[91m Decompressing\033[0m')
                prog.update()
                prog.close()
                with tarfile.open(dest, 'r:*') as tar:
                    tar.extractall(dest_dir)
                    # only extract files
                    files = [x.name for x in tar.getmembers() if x.isfile()]
                    for name in files:
                        if not os.path.isfile(os.path.join(dest_dir, name)):
                            return False
                        else:
                            sig.add(os.path.join(dest_dir, name))
                            decompressed += 1
            elif dest.endswith('.gz'):
                prog.set_description(message + ':\033[91m Decompressing\033[0m')
                prog.update()
                prog.close()
                decomp = dest[:-3]
                with gzip.open(dest, 'rb') as fin, open(decomp, 'wb') as fout:
                    buffer = fin.read(100000)
                    while buffer:
                        fout.write(buffer)
                        buffer = fin.read(100000)
                sig.add(decomp)
                decompressed += 1
        decompress_msg = '' if not decompressed else ' ({} file{} decompressed)'.format(
            decompressed, '' if decompressed <= 1 else 's')
        prog.set_description(message + ':\033[32m downloaded{} {}\033[0m'.format(decompress_msg,
            ' '*(term_width - len(message) - 13 - len(decompress_msg))))
        prog.update()
        prog.close()
        # if a md5 file exists
        # if downloaded files contains .md5 signature, use them to validate
        # downloaded files.
        if os.path.isfile(dest + '.md5'):
            prog.set_description(message + ':\033[91m Verifying md5 signature\033[0m')
            prog.update()
            prog.close()
            with open(dest + '.md5') as md5:
                rec_md5 = md5.readline().split()[0].strip()
                obs_md5 = fileMD5(dest, partial=False)
                if rec_md5 != obs_md5:
                    prog.set_description(message + ':\033[91m MD5 signature mismatch\033[0m')
                    prog.update()
                    prog.close()
                    env.logger.warning('md5 signature mismatch for downloaded file {} (recorded {}, observed {})'
                        .format(filename[:-4], rec_md5, obs_md5))
            prog.set_description(message + ':\033[91m MD5 signature verified\033[0m')
            prog.update()
            prog.close()
    except Exception as e:
        if env.verbosity > 2:
             sys.stderr.write(get_traceback())
        env.logger.error('Failed to download: {}'.format(e))
        return False
    finally:
        # if there is something wrong still remove temporary file
        if os.path.isfile(dest_tmp):
            os.remove(dest_tmp)
    sig.write_sig()
    return os.path.isfile(dest)


@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['URLs', 'dest_dir', 'dest_file', 'decompress'])
def download(URLs, dest_dir='.', dest_file=None, decompress=False):
    '''Download files from specified URL, which should be space, tab or
    newline separated URLs. The files will be downloaded to specified
    destination. If `filename.md5` files are downloaded, they are used to
    validate downloaded `filename`. Unless otherwise specified, compressed
    files are decompressed.
    '''
    if env.config['run_mode'] == 'dryrun':
        print('download\n{}\n'.format(URLs))
        return None
    if isinstance(URLs, str):
        urls = [x.strip() for x in URLs.split() if x.strip()]
    else:
        urls = list(URLs)
    #
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
        with mp.Pool(processes = env.sos_dict['CONFIG'].get('sos_download_processes', 5)) as pool:
            for idx, (url, filename) in enumerate(zip(urls, filenames)):
                if not filename:
                    continue
                succ[idx] = pool.apply_async(downloadURL, (url, filename,
                    decompress, idx))
            succ = [x.get() if isinstance(x, mp.pool.AsyncResult) else x for x in succ]
    else:
        if dest_file is not None:
            succ[0] = downloadURL(urls[0], dest_file, decompress=decompress)
        else:
           if filenames[0]:
                succ[0] = downloadURL(urls[0], filenames[0], decompress=decompress)
    #
    #for su, url in zip(succ, urls):
    #    if not su:
    #        env.logger.warning('Failed to download {}'.format(url))
    if not all(succ):
        raise RuntimeError('Not all files have been downloaded')
    return 0

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script', 'args'])
def run(script, args='', **kwargs):
    '''Execute specified script using bash.'''
    return SoS_ExecuteScript(script, '', '', args).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script', 'args'])
def perl(script, args='', **kwargs):
    '''Execute specified script using perl.'''
    return SoS_ExecuteScript(script, 'perl', '.pl', args).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script', 'args'])
def ruby(script, args='', **kwargs):
    '''Execute specified script using ruby.'''
    return SoS_ExecuteScript(script, 'ruby', '.rb', args).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script', 'args'])
def node(script, args='', **kwargs):
    '''Execute specified script using node.'''
    return SoS_ExecuteScript(script, 'node', '.js', args).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script', 'args'])
def JavaScript(script, args='', **kwargs):
    '''Execute specified script using node.'''
    return SoS_ExecuteScript(script, 'node', '.js', args).run(**kwargs)


def collect_input(script, input):
    # determine file extension
    if input is not None:
        if isinstance(input, str):
            ext = os.path.splitext(input)[-1]
        elif isinstance(input, Sequence) and len(input) > 0:
            ext = os.path.splitext(input[0])[-1]
        else:
            raise ValueError('Unknown input file for action pandoc')
    else:
        ext = '.md'

    input_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=ext, delete=False).name
    with open(input_file, 'w') as tmp:
        if script is not None and script.strip():
            tmp.write(script.rstrip() + '\n\n')
        if isinstance(input, str):
            try:
                with open(input) as ifile:
                    tmp.write(ifile.read() + '\n\n')
            except  Exception as e:
                    raise ValueError('Failed to read input file {}: {}'.format(input, e))
        elif isinstance(input, Sequence):
            for ifile in input:
                try:
                    with open(ifile) as itmp:
                        tmp.write(itmp.read().rstrip() + '\n\n')
                except  Exception as e:
                        raise ValueError('Failed to read input file {}: {}'.format(ifile, e))
    return input_file


@SoS_Action(run_mode=['dryrun', 'run', 'interactive'], acceptable_args=['script'])
def report(script=None, input=None, output=None, **kwargs):
    '''Write script to an output file specified by `output`, which can be
    a filename to which the content of the script will be written,
    any object with a "write" attribute (e.g. a file handle) for which the "write"
    function will be called with the content. If output is unspecified, the content
    will be written to standard output or appended to a file specified with command
    line option `-r`. '''
    if env.config['run_mode'] == 'dryrun':
        print('report:\n{}'.format('' if script is None else script))
        if input is not None:
            for ifile in input:
                print('  from file: {}'.format(ifile))
        return
    file_handle = None
    if isinstance(output, str):
        if not output or output == '-':
            writer = sys.stdout.write
        elif output.startswith('>>'):
            file_handle = open(os.path.expanduser(output[2:]), 'a')
            writer = file_handle.write
        else:
            file_handle = open(os.path.expanduser(output), 'w')
            writer = file_handle.write
    elif hasattr(output, 'write'):
        writer = output.write
    elif output is None or output == '':
        writer = sys.stdout.write
    else:
        raise ValueError('Invalid output {}.'.format(output))

    # file lock to prevent race condition
    with fasteners.InterProcessLock('/tmp/report_lock'):
        if isinstance(script, str) and script.strip():
            writer(script.rstrip() + '\n\n')
        if input is not None:
            if isinstance(input, str):
                env.logger.debug('Loading report from {}'.format(input))
                with open(input) as ifile:
                    writer(ifile.read().rstrip() + '\n\n')
            elif isinstance(input, Sequence):
                for ifile in input:
                    try:
                        env.logger.debug('Loading report from {}'.format(ifile))
                        with open(ifile) as itmp:
                            writer(itmp.read().rstrip() + '\n\n')
                    except Exception as e:
                        raise ValueError('Failed to read input file {}: {}'.format(ifile, e))
            else:
                raise ValueError('Unknown input file for action report')
    #
    if file_handle:
        file_handle.close()


@SoS_Action(run_mode=['run', 'interactive'], acceptable_args=['script', 'args'])
def pandoc(script=None, input=None, output=None, args='${input!q} --output ${output!q}', **kwargs):
    '''Convert input file to output using pandoc

    The input can be specified in three ways:

    1. instant script, which is assumed to be in md format

    pandoc:   output='report.html'
      script

    2. one or more input files. The format is determined by extension of input file

    pandoc(input, output='report.html')

    3. input file specified by command line option `-r` .
    pandoc(output='report.html')

    If no output is specified, it is assumed to be in html format
    and is written to standard output.

    You can specify more options such as "from" and "to" by customizing
    the args parameter of the action. The default value of args is
    `${input!q} --output ${output!q}'
    '''
#
#     # this is output format
#     pandoc [OPTIONS] [FILES]
# Input formats:  commonmark, docbook, docx, epub, haddock, html, json*, latex,
#                 markdown, markdown_github, markdown_mmd, markdown_phpextra,
#                 markdown_strict, mediawiki, native, odt, opml, org, rst, t2t,
#                 textile, twiki
#                 [ *only Pandoc's JSON version of native AST]
# Output formats: asciidoc, beamer, commonmark, context, docbook, docx, dokuwiki,
#                 dzslides, epub, epub3, fb2, haddock, html, html5, icml, json*,
#                 latex, man, markdown, markdown_github, markdown_mmd,
#                 markdown_phpextra, markdown_strict, mediawiki, native, odt,
#                 opendocument, opml, org, pdf**, plain, revealjs, rst, rtf, s5,
#                 slideous, slidy, tei, texinfo, textile
#                 [**for pdf output, use latex or beamer and -o FILENAME.pdf]
# Options:
#   -f FORMAT, -r FORMAT  --from=FORMAT, --read=FORMAT
#   -t FORMAT, -w FORMAT  --to=FORMAT, --write=FORMAT
#   -o FILENAME           --output=FILENAME
#                         --data-dir=DIRECTORY
#   -R                    --parse-raw
#   -S                    --smart
#
# IGNORED
#
    if not executable('pandoc').exists():
        raise UnknownTarget(executable('pandoc'))

    input_file = collect_input(script, input)

    write_to_stdout = False
    if output is None:
        write_to_stdout = True
        output_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.html', delete=False).name
    elif isinstance(output, str):
        output_file = os.path.expanduser(output)
    else:
        raise RuntimeError('A filename is expected, {} provided'.format(output))
    #
    ret = 1
    try:
        p = None
        cmd = interpolate('pandoc {}'.format(args), '${ }', {'input': input_file, 'output': output_file})
        env.logger.trace('Running command "{}"'.format(cmd))
        if env.config['run_mode'] == 'interactive':
            # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
            p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            out, err = p.communicate()
            sys.stdout.write(out.decode())
            sys.stderr.write(err.decode())
            ret = p.returncode
        else:
            p = subprocess.Popen(cmd, shell=True)
            ret = p.wait()
    except Exception as e:
        env.logger.error(e)
    if ret != 0:
        temp_file = os.path.join('.sos', '{}_{}.md'.format('pandoc', os.getpid()))
        shutil.copyfile(input_file, temp_file)
        cmd = interpolate('pandoc {}'.format(args), '${ }', {'input': temp_file, 'output': output_file})
        raise RuntimeError('Failed to execute script. Please use command \n{}\nunder {} to test it.'
            .format(cmd, os.getcwd()))
    if write_to_stdout:
        with open(output_file) as out:
            sys.stdout.write(out.read())
    else:
        env.logger.info('Report saved to {}'.format(output))
    try:
        os.remove(input_file)
    except:
        pass

