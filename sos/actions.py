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
import subprocess
import tempfile
import shlex
import glob
import copy
import urllib
import urllib.request
import urllib.error
import urllib.parse
import shutil
import zipfile
import gzip
import tarfile
from functools import wraps

if sys.platform != 'win32':
    import blessings

from collections.abc import Sequence
import multiprocessing as mp
from .utils import env, ProgressBar, natural_keys, transcribe, AbortExecution, short_repr, get_traceback
from .sos_eval import Undetermined
from .target import FileTarget, fileMD5
from .monitor import ProcessMonitor, summarizeExecution

__all__ = ['SoS_Action', 'execute_script', 'sos_run',
    'fail_if', 'warn_if', 'stop_if',
    'download',
    'run', 'bash', 'csh', 'tcsh', 'zsh', 'sh',
    'python', 'python3',
    'perl', 'ruby', 'node', 'JavaScript',
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
def SoS_Action(run_mode=['run', 'interactive']):
    run_mode = [run_mode] if isinstance(run_mode, str) else run_mode
    def runtime_decorator(func):
        @wraps(func)
        def action_wrapper(*args, **kwargs):
            # docker files will be downloaded in run or prepare mode
            if 'docker_file' in kwargs and env.run_mode in ['run', 'interactive']:
                from .docker.client import DockerClient
                docker = DockerClient()
                docker.import_image(kwargs['docker_file'])
            # handle image
            if 'docker_image' in kwargs:
                from .docker.client import DockerClient
                docker = DockerClient()
                docker.pull(kwargs['docker_image'])
            if env.run_mode not in run_mode:
                # return dynamic expression when not in run mode, that is to say
                # the script logic cannot rely on the result of the action
                return Undetermined(func.__name__)
            if env.run_mode == 'interactive':
                for k,v in kwargs.items():
                    if k in SOS_RUNTIME_OPTIONS and k not in SOS_ACTION_OPTIONS:
                        env.logger.warning('Passing runtime option "{0}" to action is deprecated. Please use "task: {0}={1}" before action instead.'.format(k, v))
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
                    res = func(*args, **kwargs)
                finally:
                    os.chdir(olddir)
            else:
                res = func(*args, **kwargs)
            return res
        action_wrapper.run_mode = run_mode
        return action_wrapper
    return runtime_decorator


class SoS_ExecuteScript:
    def __init__(self, script, interpreter, suffix, args=''):
        self.script = script
        self.interpreter = interpreter
        if args:
            self.interpreter += ' ' + args
        self.suffix = suffix

    def run(self, **kwargs):
        transcribe(self.script, action=self.interpreter)
        if '{}' not in self.interpreter:
            self.interpreter += ' {}'
        debug_script_file = os.path.join(env.exec_dir, '.sos', '{}_{}{}'.format(env.sos_dict['step_name'],
            env.sos_dict['_index'], self.suffix))
        env.logger.debug('Script for step {} is saved to {}'.format(env.sos_dict['step_name'], debug_script_file))
        with open(debug_script_file, 'w') as sfile:
            sfile.write(self.script)
        if 'docker_image' in kwargs:
            from .docker.client import DockerClient
            docker = DockerClient()
            docker.run(kwargs['docker_image'], self.script, self.interpreter, self.suffix,
                **kwargs)
        else:
            try:
                script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=self.suffix, delete=False).name
                with open(script_file, 'w') as sfile:
                    sfile.write(self.script)
                cmd = self.interpreter.replace('{}', shlex.quote(script_file))
                #
                if env.run_mode == 'interactive':
                    # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
                    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                    pid = p.pid
                    env.register_process(p.pid, 'Runing {}'.format(script_file))
                    out, err = p.communicate()
                    sys.stdout.write(out.decode())
                    sys.stderr.write(err.decode())
                    ret = p.returncode
                    sys.stdout.flush()
                    sys.stderr.flush()
                else:
                    p = subprocess.Popen(cmd, shell=True)
                    pid = p.pid
                    if '__step_sig__' in env.sos_dict and env.sos_dict['__step_sig__'] is not None:
                        m = ProcessMonitor(pid, msg=self.script, sig=env.sos_dict['__step_sig__'])
                        m.start()
                    env.register_process(pid, 'Runing {}'.format(script_file))
                    ret = p.wait()
                    if '__step_sig__' in env.sos_dict and env.sos_dict['__step_sig__'] is not None:
                        summarizeExecution(pid)

                if ret != 0:
                    with open(debug_script_file, 'w') as sfile:
                        sfile.write(self.script)
                    cmd = self.interpreter.replace('{}', shlex.quote(debug_script_file))
                    raise RuntimeError('Failed to execute script. The script is saved to {}. Please use command "{}" under {} to test it.'
                        .format(debug_script_file, cmd, os.getcwd()))
            except Exception as e:
                env.logger.error('Failed to execute script: {}'.format(e))
                raise
            finally:
                env.deregister_process(p.pid)
                os.remove(script_file)


@SoS_Action(run_mode=['run', 'interactive'])
def sos_run(workflow, **kwargs):
    '''Execute a workflow from specified source, input, and output
    By default the workflow is defined in the existing SoS script, but
    extra sos files can be specified from paramter source. The workflow
    will be execute in the current step namespace with _input as workflow
    input. '''
    from .sos_executor import Base_Executor, MP_Executor
    script = SoS_Script(env.sos_dict['__step_context__'].content, env.sos_dict['__step_context__'].filename)
    wf = script.workflow(workflow)
    # if wf contains the current step or one of the previous one, this constitute
    # recusive nested workflow and should not be allowed
    if env.sos_dict['step_name'] in ['{}_{}'.format(x.name, x.index) for x in wf.sections]:
        raise RuntimeError('Nested workflow {} contains the current step {}'.format(workflow, env.sos_dict['step_name']))
    # for nested workflow, _input would becomes the input of workflow.
    for k,v in kwargs.items():
        env.sos_dict.set(k, v)
    env.sos_dict.set('__step_output__', copy.deepcopy(env.sos_dict['_input']))
    try:
        my_name = env.sos_dict['step_name']
        if env.run_mode == 'dryrun':
            env.logger.info('Checking nested workflow {}'.format(workflow))
            return Base_Executor(wf, args=env.sos_dict['__args__'], nested=True).dryrun()
        elif env.run_mode in ('run', 'interactive'):
            env.logger.info('Executing workflow ``{}`` with input ``{}``'
                .format(workflow, short_repr(env.sos_dict['_input'], True)))

            if env.__task_engine__:
                import pkg_resources
                # import all executors
                executor_class = None
                for entrypoint in pkg_resources.iter_entry_points(group='sos_executors'):
                    # Grab the function that is the actual plugin.
                    name = entrypoint.name
                    if name == env.__task_engine__:
                        try:
                            executor_class = entrypoint.load()
                        except Exception as e:
                            print('Failed to load queue executor {}: {}'.format(entrypoint.name, e))

                if not executor_class:
                    sys.exit('Could not locate specified queue executor {}'.format(env.__task_engine__))
            else:
                if env.max_jobs == 1:
                    executor_class = Base_Executor
                else:
                    executor_class = MP_Executor

            return executor_class(wf, args=env.sos_dict['__args__'], nested=True).run()
    finally:
        # restore step_name in case the subworkflow re-defines it
        env.sos_dict.set('step_name', my_name)

@SoS_Action(run_mode=['run', 'interactive'])
def execute_script(script, interpreter, suffix, args='', **kwargs):
    '''Execute specified script using specified interpreter.'''
    return SoS_ExecuteScript(script, interpreter, suffix, args).run(**kwargs)

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'])
def fail_if(expr, msg=''):
    '''Raise an exception with `msg` if condition `expr` is False'''
    if expr:
        raise RuntimeError(msg)
    return 0

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'])
def warn_if(expr, msg=''):
    '''Yield an warning message `msg` if `expr` is False '''
    if expr:
        env.logger.warning(msg)
    return 0

@SoS_Action(run_mode=['dryrun', 'run', 'interactive'])
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
        prog = ProgressBar(message, disp=env.verbosity > 1)
        sig = FileTarget(dest)
        if os.path.isfile(dest):
            if env.sig_mode == 'construct':
                prog.done(message + ': \033[32m writing signature\033[0m')
                sig.write_sig()
                prog.done(message + ': \033[32m signature calculated\033[0m')
                return True
            elif env.sig_mode == 'ignore':
                prog.done(message + ': \033[32m use existing\033[0m')
                return True
            else:
                prog.done(message + ': \033[32m Validating signature\033[0m')
                if sig.validate():
                    prog.done(message + ': \033[32m Validated\033[0m')
                    return True
                else:
                    prog.done(message + ':\033[91m Signature mismatch\033[0m')
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
        #    prog.done(message + ':\033[91m 404 Error {}\033[0m'.format(' '*(term_width - len(message) - 12)))
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
                    prog.urllibUpdate(file_size, file_size_dl)
            except urllib.error.HTTPError as e:
                prog.done(message + ':\033[91m {} Error\033[0m'.format(e.code))
                try:
                    os.remove(dest_tmp)
                except OSError:
                    pass
                return False
            except Exception as e:
                prog.done(message + ':\033[91m {}\033[0m'.format(e))
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
                prog.done(message + ':\033[91m Decompressing\033[0m')
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
                prog.done(message + ':\033[91m Decompressing\033[0m')
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
                prog.done(message + ':\033[91m Decompressing\033[0m')
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
        prog.done(message + ':\033[32m downloaded{} {}\033[0m'.format(decompress_msg,
            ' '*(term_width - len(message) - 13 - len(decompress_msg))))
        # if a md5 file exists
        # if downloaded files contains .md5 signature, use them to validate
        # downloaded files.
        if os.path.isfile(dest + '.md5'):
            prog.done(message + ':\033[91m Verifying md5 signature\033[0m')
            with open(dest + '.md5') as md5:
                rec_md5 = md5.readline().split()[0].strip()
                obs_md5 = fileMD5(dest, partial=False)
                if rec_md5 != obs_md5:
                    prog.done(message + ':\033[91m MD5 signature mismatch\033[0m')
                    env.logger.warning('md5 signature mismatch for downloaded file {} (recorded {}, observed {})'
                        .format(filename[:-4], rec_md5, obs_md5))
            prog.done(message + ':\033[91m MD5 signature verified\033[0m')
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


@SoS_Action(run_mode=['run', 'interactive'])
def download(URLs, dest_dir='.', dest_file=None, decompress=False):
    '''Download files from specified URL, which should be space, tab or
    newline separated URLs. The files will be downloaded to specified
    destination. If `filename.md5` files are downloaded, they are used to
    validate downloaded `filename`. Unless otherwise specified, compressed
    files are decompressed.
    '''
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
        if sys.platform != 'win32':
            t = blessings.Terminal(stream=sys.stderr)
            sys.stderr.write(t.move( t.height, 0)) # + '\n')
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

@SoS_Action(run_mode=['run', 'interactive'])
def run(script, args='', **kwargs):
    '''Execute specified script using bash.'''
    return SoS_ExecuteScript(script, '/bin/bash', '.sh', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def bash(script, args='', **kwargs):
    '''Execute specified script using bash.'''
    return SoS_ExecuteScript(script, '/bin/bash', '.sh', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def csh(script, args='', **kwargs):
    '''Execute specified script using csh.'''
    return SoS_ExecuteScript(script, '/bin/csh', '.csh', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def tcsh(script, args='', **kwargs):
    '''Execute specified script using tcsh.'''
    return SoS_ExecuteScript(script, '/bin/tcsh', '.sh', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def zsh(script, args='', **kwargs):
    '''Execute specified script using zsh.'''
    return SoS_ExecuteScript(script, '/bin/zsh', '.zsh', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def sh(script, args='', **kwargs):
    '''Execute specified script using sh.'''
    return SoS_ExecuteScript(script, '/bin/sh', '.sh', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def python(script, args='', **kwargs):
    '''Execute specified script using python (which can be python 2 or 3 depending on system configuration.'''
    return SoS_ExecuteScript(script, 'python', '.py', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def python3(script, args='', **kwargs):
    '''Execute specified script using python 3.'''
    return SoS_ExecuteScript(script, 'python3', '.py', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def perl(script, args='', **kwargs):
    '''Execute specified script using perl.'''
    return SoS_ExecuteScript(script, 'perl', '.pl', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def ruby(script, args='', **kwargs):
    '''Execute specified script using ruby.'''
    return SoS_ExecuteScript(script, 'ruby', '.rb', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def node(script, args='', **kwargs):
    '''Execute specified script using node.'''
    return SoS_ExecuteScript(script, 'node', '.js', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def JavaScript(script, args='', **kwargs):
    '''Execute specified script using node.'''
    return SoS_ExecuteScript(script, 'node', '.js', args).run(**kwargs)

@SoS_Action(run_mode=['run', 'interactive'])
def report(script=None, from_file=None, to_file=None, mode='a', **kwargs):
    if to_file is not None:
        if not to_file:
            raise RuntimeError('Invalid parameter to_file "{}"'.format(to_file))
        report_file = to_file
    else:
        report_file = '.sos/report.md'
    #
    content = ''
    if script is not None:
        content = script
    if from_file is not None:
        try:
            with open(from_file) as rep:
                content += rep.read().decode()
        except Exception as e:
            raise RuntimeError('Failed to import report from {}: {}'.format(from_file, e))
    #
    # write report file (the ${} expressions must have been interpolated.
    if report_file == '__STDERR__':
        sys.stderr.write(content)
    else:
        with open(report_file, mode)as md:
            md.write(content)


@SoS_Action(run_mode=['run', 'interactive'])
def pandoc(script=None, output=None, **kwargs):
    '''This action can be used in three ways

    pandoc:   outputfile='report.html'
      script

    pandoc(filename='report.sos', outputfile='report.html')

    pandoc(outputfile='report.html')

    '''
    # in run mode, collect report and call pandoc
    sos_script = env.sos_dict['__step_context__'].filename
    # this is the case for stirng input (test only)
    if sos_script is None:
        sos_script = 'string_input'
    if script is not None:
        # get a temporary file with the content
        script_file = '{}.md'.format(os.path.basename(sos_script))
        with open(script_file, 'w') as report:
            report.write(script)
    elif 'filename' in kwargs:
        script_file = kwargs['filename']
    elif env.run_mode == 'interactive' and '__summary_report__' in env.sos_dict:
        script_file = env.sos_dict['__summary_report__']
    else:
        step_reports = glob.glob(os.path.join('.sos', 'report', '*'))
        step_reports.sort(key=natural_keys)
        # merge the files
        script_file = '{}.md'.format(os.path.basename(sos_script))
        env.logger.trace('Gathering reports {} to {}'.format(', '.join(step_reports), script_file))
        with open(script_file, 'w') as combined:
            for step_report in step_reports:
                with open(step_report, 'r') as md:
                    combined.write(md.read())
# $ pandoc -h
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
    arg_from = '--from={}'.format(kwargs['from']) if 'from' in kwargs else ''
    arg_to = '--to={}'.format(kwargs['to']) if 'to' in kwargs else ''
    if output is None:
        raise RuntimeError('Parameter output is required for action pandoc')
    elif not isinstance(output, str):
        raise RuntimeError('A filename is expected, {} provided'.format(output))
    arg_output = '--output={}'.format(shlex.quote(output))
    extra_args = ''
    if 'extra_args' in kwargs:
        ea = kwargs['extra_args']
        if isinstance(ea, str):
            extra_args = ea
        elif isinstance(ea, Sequence):
            extra_args = ' '.join(list(ea))
        elif isinstance(ea, dict):
            extra_args = ' '.join('--{}={}'.format(k,v) for k,v in ea.items())
    #
    command = 'pandoc {} {} {} {} {{}}'.format(arg_from, arg_to, arg_output, extra_args)
    try:
        cmd = command.replace('{}', shlex.quote(script_file))
        env.logger.trace('Running command "{}"'.format(cmd))
        if env.run_mode == 'interactive':
            # need to catch output and send to python output, which will in trun be hijacked by SoS notebook
            p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            pid = p.pid
            env.register_process(p.pid, 'Runing {}'.format(script_file))
            out, err = p.communicate()
            sys.stdout.write(out.decode())
            sys.stderr.write(err.decode())
            ret = p.returncode
        else:
            p = subprocess.Popen(cmd, shell=True)
            pid = p.pid
            env.register_process(pid, 'Runing {}'.format(script_file))
            ret = p.wait()
    except Exception as e:
        env.logger.error(e)
    finally:
        env.deregister_process(p.pid)
    if ret != 0:
        temp_file = os.path.join('.sos', '{}_{}.md'.format('pandoc', os.getpid()))
        shutil.copyfile(script_file, temp_file)
        cmd = command.replace('{}', shlex.quote(temp_file))
        raise RuntimeError('Failed to execute script. The script is saved to {}. Please use command "{}" to test it.'
            .format(temp_file, cmd))
    env.logger.info('Report saved to {}'.format(output))

