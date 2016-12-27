#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
#
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
import fnmatch
import contextlib
import subprocess
import argparse
import pkg_resources
import pydoc

from ipykernel.ipkernel import IPythonKernel


import pickle
from types import ModuleType
from sos.utils import env, WorkflowDict, short_repr, pretty_size
from sos._version import __sos_version__, __version__
from sos.sos_eval import SoS_exec, SoS_eval, interpolate, get_default_global_sigil
from sos.sos_syntax import SOS_SECTION_HEADER

from IPython.lib.clipboard import ClipboardEmpty, osx_clipboard_get, tkinter_clipboard_get
from IPython.core.error import UsageError
from IPython.core.display import HTML
from IPython.utils.tokenutil import line_at_cursor, token_at_cursor
from jupyter_client import manager, find_connection_file

from textwrap import dedent
from io import StringIO

from .completer import SoS_Completer
from .inspector import SoS_Inspector

from .sos_executor import runfile

class FlushableStringIO(StringIO):
    '''This is a string buffer for output, but it will only
    keep the first 200 lines and the last 10 lines.
    '''
    def __init__(self, kernel, name, *args, **kwargs):
        StringIO.__init__(self, *args, **kwargs)
        self.kernel = kernel
        self.name = name
        self.nlines = 0

    def write(self, content):
        if self.nlines > 1000:
            return
        nlines = content.count('\n')
        if self.nlines + nlines > 1000:
            StringIO.write(self, '\n'.join(content.split('\n')[:200]))
        else:
            StringIO.write(self, content)
        self.nlines += nlines

    def flush(self):
        content = self.getvalue()
        if self.nlines > 1000:
            lines = content.split('\n')
            content = '\n'.join(lines[:180]) + \
                '\n-- {} more lines --\n'.format(self.nlines - 180)
        elif self.nlines > 200:
            lines = content.split('\n')
            content = '\n'.join(lines[:180]) + \
                '\n-- {} lines --\n'.format(self.nlines - 190) + \
                '\n'.join(lines[-10:])
        if content.strip():
            self.kernel.send_response(self.kernel.iopub_socket, 'stream',
                {'name': self.name, 'text': content})
        self.truncate(0)
        self.seek(0)
        self.nlines = 0
        return len(content.strip())

__all__ = ['SoS_Kernel']

def clipboard_get():
    """ Get text from the clipboard.
    """
    if sys.platform == 'darwin':
        try:
            return osx_clipboard_get()
        except:
            return tkinter_clipboard_get()
    else:
        return tkinter_clipboard_get()


def get_previewers():
    # Note: data is zest.releaser specific: we want to pass
    # something to the plugin
    group = 'sos_previewers'
    result = []
    for entrypoint in pkg_resources.iter_entry_points(group=group):
        # if ':' in entry point name, it should be a function
        try:
            name, priority = entrypoint.name.split(',', 1)
            priority = int(priority)
        except Exception as e:
            env.logger.warning('Ignore incorrect previewer entry point {}: {}'.format(entrypoint, e))
            continue
        # If name points to a function in a module. Let us try to import the module
        if ':' in name:
            import importlib
            try:
                mod, func = name.split(':')
                imported = importlib.import_module(mod)
                result.append((getattr(imported, func), entrypoint, priority))
            except ImportError:
                env.logger.warning('Failed to load function {}:{}'.format(mod, func))
        else:
            result.append((name, entrypoint, priority))
    #
    result.sort(key=lambda x: -x[2])
    return result


class SoS_Kernel(IPythonKernel):
    implementation = 'SOS'
    implementation_version = __version__
    language = 'sos'
    language_version = __sos_version__
    language_info = {
        'mimetype': 'text/x-sos',
        'name': 'sos',
        'file_extension': '.sos',
        'pygments_lexer': 'sos',
        'codemirror_mode': 'sos',
        'nbconvert_exporter': 'sos.jupyter.converter.SoS_Exporter',
    }
    banner = "SoS kernel - script of scripts"

    ALL_MAGICS = {
        'dict',
        'matplotlib',
        'cd',
        'set',
        'restart',
        'with',
        'use',
        'get',
        'put',
        'paste',
        'run',
        'preview',
        'sandbox',
    }
    MAGIC_DICT = re.compile('^%dict(\s|$)')
    MAGIC_CONNECT_INFO = re.compile('^%connect_info(\s|$)')
    MAGIC_MATPLOTLIB = re.compile('^%matplotlib(\s|$)')
    MAGIC_CD = re.compile('^%cd(\s|$)')
    MAGIC_SET = re.compile('^%set(\s|$)')
    MAGIC_RESTART = re.compile('^%restart(\s|$)')
    MAGIC_WITH = re.compile('^%with(\s|$)')
    MAGIC_USE = re.compile('^%use(\s|$)')
    MAGIC_GET = re.compile('^%get(\s|$)')
    MAGIC_PUT = re.compile('^%put(\s|$)')
    MAGIC_PASTE = re.compile('^%paste(\s|$)')
    MAGIC_RUN = re.compile('^%run(\s|$)')
    MAGIC_RERUN = re.compile('^%rerun(\s|$)')
    MAGIC_PREVIEW = re.compile('^%preview(\s|$)')
    MAGIC_SANDBOX = re.compile('^%sandbox(\s|$)')

    def get_use_parser(self):
        parser = argparse.ArgumentParser(prog='%use',
            description='''Switch to a specified subkernel.''')
        parser.add_argument('kernel', nargs='?', default='',
            help='Kernel to switch to.')
        parser.add_argument('-i', '--in', nargs='*', dest='in_vars',
            help='Input variables (variables to get from SoS kernel)')
        parser.add_argument('-o', '--out', nargs='*', dest='out_vars',
            help='''Output variables (variables to put back to SoS kernel
            before switching back to the SoS kernel''')
        parser.error = self._parse_error
        return parser

    def get_with_parser(self):
        parser = argparse.ArgumentParser(prog='%with',
            description='''Use specified the subkernel to evaluate current
            cell''')
        parser.add_argument('kernel', nargs='?', default='',
            help='Kernel to switch to.')
        parser.add_argument('-i', '--in', nargs='*', dest='in_vars',
            help='Input variables (variables to get from SoS kernel)')
        parser.add_argument('-o', '--out', nargs='*', dest='out_vars',
            help='''Output variables (variables to put back to SoS kernel
            before switching back to the SoS kernel''')
        parser.error = self._parse_error
        return parser

    def get_preview_parser(self):
        parser = argparse.ArgumentParser(prog='%preview',
            description='''Preview files, sos variables, or expressions''')
        parser.add_argument('items', nargs='*', 
            help='''filename, variable name, or expression''')
        parser.add_argument('--off', action='store_true',
            help='''Turn off file preview''')
        parser.error = self._parse_error
        return parser

    def get_set_parser(self):
        parser = argparse.ArgumentParser(prog='%set',
            description='''Set persistent command line options for SoS runs.''')
        parser.error = self._parse_error
        return parser

    def get_run_parser(self):
        parser = argparse.ArgumentParser(prog='%run',
            description='''Execute the current cell with specified command line
            arguments. Arguments set by magic %set will be appended at the
            end of command line''')
        parser.error = self._parse_error
        return parser

    def get_rerun_parser(self):
        parser = argparse.ArgumentParser(prog='%rerun',
            description='''Re-execute the last executed code, most likely with
            different command line options''')
        parser.error = self._parse_error
        return parser

    def get_get_parser(self):
        parser = argparse.ArgumentParser(prog='%get',
            description='''Get specified variables from the SoS kernel
            to the existing subkernel.''')
        parser.add_argument('vars', nargs='?',
            help='''Names of SoS variables''')
        parser.error = self._parse_error
        return parser

    def get_put_parser(self):
        parser = argparse.ArgumentParser(prog='%put',
            description='''Put specified variables in the subkernel to the
            SoS kernel.''')
        parser.add_argument('vars', nargs='?',
            help='''Names of SoS variables''')
        parser.error = self._parse_error
        return parser

    def get_supported_languages(self):
        if self._supported_languages is not None:
            return self._supported_languages
        group = 'sos_languages'
        self._supported_languages = {}
        for entrypoint in pkg_resources.iter_entry_points(group=group):
            # Grab the function that is the actual plugin.
            name = entrypoint.name
            try:
                plugin = entrypoint.load()()
                self._supported_languages[name] = plugin
                # for convenience, we create two entries for, e.g. R and ir
                if name != plugin.kernel_name:
                    self._supported_languages[plugin.kernel_name] = plugin
            except Exception as e:
                raise RuntimeError('Failed to load language {}: {}'.format(entrypoint.name, e))
        return self._supported_languages

    supported_languages = property(lambda self:self.get_supported_languages())

    def get_completer(self):
        if self._completer is None:
            self._completer = SoS_Completer(self)
        return self._completer

    completer = property(lambda self:self.get_completer())

    def get_inspector(self):
        if self._inspector is None:
            self._inspector = SoS_Inspector(self)
        return self._inspector

    inspector = property(lambda self:self.get_inspector())

    def __init__(self, **kwargs):
        super(SoS_Kernel, self).__init__(**kwargs)
        self.options = ''
        self.kernel = 'sos'
        self.banner = self.banner + '\nConnection file {}'.format(os.path.basename(find_connection_file()))
        # FIXME: this should in theory be a MultiKernelManager...
        self.kernels = {}
        #self.shell = InteractiveShell.instance()
        self.format_obj = self.shell.display_formatter.format

        self.previewers = None
        self.original_keys = None
        self._supported_languages = None
        self._completer = None
        self._inspector = None
        self._execution_count = 1

    def _reset_dict(self):
        env.sos_dict = WorkflowDict()
        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        SoS_exec("run_mode = 'interactive'", None)
        self.original_keys = set(env.sos_dict._dict.keys()) | {'SOS_VERSION', 'CONFIG', \
            'step_name', '__builtins__', 'input', 'output', 'depends'}

    @contextlib.contextmanager
    def redirect_sos_io(self):
        save_stdout = sys.stdout
        save_stderr = sys.stderr
        sys.stdout = FlushableStringIO(self, 'stdout')
        sys.stderr = FlushableStringIO(self, 'stderr')
        yield
        sys.stdout = save_stdout
        sys.stderr = save_stderr

    def do_is_complete(self, code):
        '''check if new line is in order'''
        code = code.strip()
        if not code:
            return {'status': 'complete', 'indent': ''}
        if any(code.startswith(x) for x in ['%dict', '%paste', '%edit', '%cd', '!']):
            return {'status': 'complete', 'indent': ''}
        if code.endswith(':') or code.endswith(','):
            return {'status': 'incomplete', 'indent': '  '}
        lines = code.split('\n')
        if lines[-1].startswith(' ') or lines[-1].startswith('\t'):
            # if it is a new line, complte
            empty = [idx for idx,x in enumerate(lines[-1]) if x not in (' ', '\t')][0]
            return {'status': 'incomplete', 'indent': lines[-1][:empty]}
        #
        if SOS_SECTION_HEADER.match(lines[-1]):
            return {'status': 'incomplete', 'indent': ''}
        #
        return {'status': 'incomplete', 'indent': ''}


    def do_inspect(self, code, cursor_pos, detail_level=0):
        line, offset = line_at_cursor(code, cursor_pos)
        name = token_at_cursor(code, cursor_pos)
        data = self.inspector.inspect(name, line, cursor_pos - offset)
        
        reply_content = {'status' : 'ok'}
        reply_content['metadata'] = {}
        reply_content['found'] = True if data else False
        reply_content['data'] = data
        return reply_content

    def do_complete(self, code, cursor_pos):
        text, matches = self.completer.complete_text(code, cursor_pos)
        return {'matches' : matches,
                'cursor_end' : cursor_pos,
                'cursor_start' : cursor_pos - len(text),
                'metadata' : {},
                'status' : 'ok'}

    def warn(self, message):
        if message.strip():
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stderr', 'text': message})

    def get_magic_and_code(self, code, warn_remaining=False):
        lines = code.split('\n')
        pieces = lines[0].strip().split(None, 1)
        if len(pieces) == 2:
            command_line = pieces[1]
        else:
            command_line = ''
        remaining_code = '\n'.join(lines[1:])
        if warn_remaining and remaining_code.strip():
            self.warn('Statement {} ignored'.format(short_repr(remaining_code)))
        return command_line, remaining_code

    def run_cell(self, code, store_history):
        #
        if not self.KM.is_alive():
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': 'Restarting kernel "{}"\n'.format(self.kernel)})
            self.KM.restart_kernel(now=False)
            self.KC = self.KM.client()
        # flush stale replies, which could have been ignored, due to missed heartbeats
        while self.KC.shell_channel.msg_ready():
            self.KC.shell_channel.get_msg()
        # executing code in another kernel
        self.KC.execute(code, silent=False, store_history=not store_history)

        # first thing is wait for any side effects (output, stdin, etc.)
        _execution_state = "busy"
        while _execution_state != 'idle':
            # display intermediate print statements, etc.
            while self.KC.iopub_channel.msg_ready():
                sub_msg = self.KC.iopub_channel.get_msg()
                msg_type = sub_msg['header']['msg_type']
                if msg_type == 'status':
                    _execution_state = sub_msg["content"]["execution_state"]
                else:
                    if msg_type in ('execute_input', 'execute_result'):
                        # override execution count with the master count,
                        # not sure if it is needed
                        sub_msg['content']['execution_count'] = self._execution_count
                    #
                    # NOTE: we do not send status of sub kernel alone because
                    # these are generated automatically during the execution of
                    # "this cell" in SoS kernel
                    #
                    self.send_response(self.iopub_socket, msg_type, sub_msg['content'])
        # now get the real result
        reply = self.KC.get_shell_msg(timeout=10)
        reply['content']['execution_count'] = self._execution_count
        return reply['content']

    def switch_kernel(self, kernel, ret_vars=[]):
        if kernel and kernel != 'sos':
            if kernel != self.kernel:
                if kernel in self.kernels:
                    # switch to an active kernel
                    self.KM, self.KC = self.kernels[kernel]
                    self.RET_VARS = ret_vars
                    self.kernel = kernel
                else:
                    try:
                        # start a new kernel
                        self.kernels[kernel] = manager.start_new_kernel(startup_timeout=60, 
                            kernel_name=self.supported_languages[kernel].kernel_name if kernel in self.supported_languages else kernel)
                        self.KM, self.KC = self.kernels[kernel]
                        self.RET_VARS = ret_vars
                        self.kernel = kernel
                        if self.kernel in self.supported_languages:
                            init_stmts = self.supported_languages[self.kernel].init_statements
                            if init_stmts:
                                self.run_cell(init_stmts, False)
                    except Exception as e:
                        self.warn('Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}\n'.format(kernel, e))
        else:
            # kernel == '' or kernel == 'sos'
            if not kernel:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Kernel "{}" is used.\n'.format(self.kernel)})
            elif kernel == 'sos':
                # if we are switching back to sos, the return variables need to be
                # returned
                self.handle_magic_put(self.RET_VARS)
                self.RET_VARS = []
                self.kernel = 'sos'

    def restart_kernel(self, kernel):
        if kernel == 'sos':
            # cannot restart myself ...
            self.warn('Cannot restart sos kernel from within sos.')
        elif kernel:
            if kernel in self.kernels:
                try:
                    self.kernels[kernel][0].shutdown_kernel(restart=False)
                except Exception as e:
                    self.warn('Failed to shutdown kernel {}: {}\n'.format(kernel, e))
            #
            try:
                self.kernels[kernel] = manager.start_new_kernel(startup_timeout=60, kernel_name=kernel)
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Kernel {} {}started\n'.format(kernel, 're' if kernel in self.kernels else '')})
                #self.send_response(self.iopub_socket, 'stream',
                #    {'name': 'stdout', 'text': 'Kernel "{}" started\n'.format(kernel)})
                if kernel == self.kernel:
                    self.KM, self.KC = self.kernels[kernel]
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}\n'.format(kernel, e)})
        else:
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': 'Specify one of the kernels to restart: sos{}\n'
                    .format(''.join(', {}'.format(x) for x in self.kernels))})

    def _parse_error(self, msg):
        self.warn(msg)

    def get_dict_parser(self):
        parser = argparse.ArgumentParser(prog='%dict',
            description='Inspect or reset SoS dictionary')
        parser.add_argument('vars', nargs='*')
        parser.add_argument('-k', '--keys', action='store_true',
            help='Return only keys')
        parser.add_argument('-r', '--reset', action='store_true',
            help='Rest SoS dictionary (clear all user variables)')
        parser.add_argument('-a', '--all', action='store_true',
            help='Return all variales, including system functions and variables')
        parser.add_argument('-d', '--del', nargs='+', metavar='VAR', dest='__del__',
            help='Remove specified variables from SoS dictionary')
        parser.error = self._parse_error
        return parser

    def get_sandbox_parser(self):
        parser = argparse.ArgumentParser(prog='%sandbox',
            description='''Execute content of a cell in a temporary directory
                with fresh dictionary (by default).''')
        parser.add_argument('-d', '--dir',
            help='''Execute workflow in specified directory. The directory
                will be created if does not exist, and will not be removed
                after the completion. ''')
        parser.add_argument('-k', '--keep-dict', action='store_true',
            help='''Keep current sos dictionary.''')
        parser.add_argument('-e', '--expect-error', action='store_true',
            help='''If set, expect error from the excution and report
                success if an error occurs.''')
        parser.error = self._parse_error
        return parser

    def handle_magic_dict(self, line):
        'Magic that displays content of the dictionary'
        # do not return __builtins__ beacuse it is too long...
        import shlex
        parser = self.get_dict_parser()
        args = parser.parse_args(shlex.split(line))

        for x in args.vars:
            if not x in env.sos_dict:
                self.warn('Unrecognized sosdict option or variable name {}'.format(x))
                return

        if args.reset:
            self._reset_dict()
            return

        if args.__del__:
            for x in args.__del__:
                if x in env.sos_dict:
                    env.sos_dict.pop(x)
            return

        if args.keys:
            if args.all:
                self.send_result(env.sos_dict._dict.keys())
            elif args.vars:
                self.send_result(set(args.vars))
            else:
                self.send_result({x for x in env.sos_dict._dict.keys() if not x.startswith('__')} - self.original_keys)
        else:
            if args.all:
                self.send_result(env.sos_dict._dict)
            elif args.vars:
                self.send_result({x:y for x,y in env.sos_dict._dict.items() if x in args.vars})
            else:
                self.send_result({x:y for x,y in env.sos_dict._dict.items() if x not in self.original_keys and not x.startswith('__')})

    def handle_magic_set(self, options):
        if options.strip():
            #self.send_response(self.iopub_socket, 'stream',
            #    {'name': 'stdout', 'text': 'sos options set to "{}"\n'.format(options)})
            if not options.strip().startswith('-'):
                self.warn('Magic %set cannot set positional argument, {} provided.'.format(options))
            else:
                self.options = options.strip()
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'sos options is set to "{}"\n'.format(self.options)})
        else:
            if self.options:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'sos options "{}" reset to ""\n'.format(self.options)})
                self.options = ''
            else:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Usage: set persistent sos options such as workflow name and -v 3 (debug output)\n'})

    def handle_magic_get(self, items):
        for item in items:
            if item not in env.sos_dict:
                self.warn('Variable {} does not exist'.format(item))
                return
        if self.kernel == 'python3':
            # if it is a python kernel, passing specified SoS variables to it
            sos_data = pickle.dumps({x:env.sos_dict[x] for x in items})
            # this can fail if the underlying python kernel is python 2
            self.KC.execute("import pickle\nglobals().update(pickle.loads({!r}))".format(sos_data),
                silent=True, store_history=False)
        elif self.kernel in self.supported_languages:
            lan = self.supported_languages[self.kernel]
            try:
                statements = []
                for item in items:
                    new_name, py_repr = lan.py_to_lan(item, env.sos_dict[item])
                    if new_name != item:
                        self.warn('Variable {} is passed from SoS to kernel {} as {}'
                            .format(item, self.kernel, new_name))
                    statements.append(py_repr)
            except Exception as e:
                self.warn('Failed to get variable: {}\n'.format(e))
                return
            self.KC.execute('\n'.join(statements), silent=True, store_history=False)
        elif self.kernel == 'sos':
            self.warn('Magic %get can only be executed by subkernels')
            return
        else:
            self.warn('Language {} does not support magic %get.'.format(self.kernel))
            return
        # first thing is wait for any side effects (output, stdin, etc.)
        _execution_state = "busy"
        while _execution_state != 'idle':
            # display intermediate print statements, etc.
            while self.KC.iopub_channel.msg_ready():
                sub_msg = self.KC.iopub_channel.get_msg()
                msg_type = sub_msg['header']['msg_type']
                if msg_type == 'status':
                    _execution_state = sub_msg["content"]["execution_state"]
                else:
                    self.send_response(self.iopub_socket, msg_type,
                        sub_msg['content'])

    def handle_magic_put(self, items):
        if not items:
            return
        if self.kernel == 'python':
            # if it is a python kernel, passing specified SoS variables to it
            self.KC.execute('import pickle\npickle.dumps({{ {} }})'.format(','.join('"{0}":{0}'.format(x) for x in items)),
                silent=False, store_history=False)
            # first thing is wait for any side effects (output, stdin, etc.)
            _execution_state = "busy"
            while _execution_state != 'idle':
                # display intermediate print statements, etc.
                while self.KC.iopub_channel.msg_ready():
                    sub_msg = self.KC.iopub_channel.get_msg()
                    msg_type = sub_msg['header']['msg_type']
                    if msg_type == 'status':
                        _execution_state = sub_msg["content"]["execution_state"]
                    else:
                        if msg_type == 'execute_result':
                            #self.warn(repr(sub_msg['content']['data']))
                            try:
                                env.sos_dict.update(
                                    pickle.loads(eval(sub_msg['content']['data']['text/plain']))
                                    )
                            except Exception as e:
                                self.warn('Failed to put variable {}: {}\n'.format(', '.join(items), e))
                            break
                        else:
                            self.send_response(self.iopub_socket, msg_type,
                                sub_msg['content'])
            # verify
            for item in items:
                if not item in env.sos_dict:
                    self.warn('Failed to put variable {} to SoS namespace\n'.format(item))
        elif self.kernel in self.supported_languages:
            lan = self.supported_languages[self.kernel]
            new_names, py_repr = lan.lan_to_py(items)
            for n, nn in zip(items, new_names):
                if n != nn:
                    self.warn('Variable {} is put to SoS as {}'.format(n, nn))
            self.KC.execute(py_repr, silent=False, store_history=False)
            # first thing is wait for any side effects (output, stdin, etc.)
            _execution_state = "busy"
            while _execution_state != 'idle':
                # display intermediate print statements, etc.
                while self.KC.iopub_channel.msg_ready():
                    sub_msg = self.KC.iopub_channel.get_msg()
                    msg_type = sub_msg['header']['msg_type']
                    if msg_type == 'status':
                        _execution_state = sub_msg["content"]["execution_state"]
                    else:
                        # irkernel (since the new version) does not produce execute_result, only
                        # display_data
                        if msg_type in ('display_data', 'execute_result'):
                            try:
                                env.sos_dict.update(
                                    lan.py_to_dict(sub_msg['content']['data']['text/plain'])
                                    )
                            except Exception as e:
                                 self.warn(str(e))
                            break
                        else:
                            self.send_response(self.iopub_socket, msg_type,
                                sub_msg['content'])
            # verify
            for n, nn in zip(items, new_names):
                if not nn in env.sos_dict:
                    self.warn('Failed to put variable {} to SoS namespace\n'.format(n))
        elif self.kernel == 'sos':
            self.warn('Magic %put can only be executed by subkernels')
        else:
            self.warn('Language {} does not support magic %put.'.format(self.kernel))

    def _interpolate_option(self, option, quiet=False):
        # interpolate command
        try:
            new_option = interpolate(option, sigil='${ }', local_dict=env.sos_dict._dict)
            if new_option != option and not quiet:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text':
                    new_option.strip() + '\n## -- End interpolated command --\n'})
            return new_option
        except Exception as e:
            self.warn('Failed to interpolate {}: {}\n'.format(short_repr(option), e))
            return None

    def handle_magic_preview(self, options):
        # we do string interpolation here because the execution of 
        # statements before it can change the meanings of them.
        options = self._interpolate_option(options, quiet=True)
        if options is None:
            return
        # find filenames and quoted expressions
        import shlex
        parser = self.get_preview_parser()
        args = parser.parse_args(shlex.split(options, posix=False))
        if not args.items or args.off:
            return
        self.send_response(self.iopub_socket, 'display_data',
            {
              'source': 'SoS',
              'metadata': {},
              'data': { 'text/html': HTML('<pre><font color="green">## %preview {}</font></pre>'.format(options)).data}
            })
        # expand items
        for item in args.items:
            try:
                if os.path.isfile(item):
                    self.preview_file(item)
                    continue
            except Exception as e:
                self.warn('\n> Failed to preview file {}: {}'.format(item, e))
                continue
            try:
                self.send_response(self.iopub_socket, 'display_data',
                    {'metadata': {},
                    'data': {'text/plain': '>>> ' + item + ':\n',
                        'text/html': HTML('<pre><font color="green">> {}:</font></pre>'.format(item)).data
                        }
                    })
                format_dict, md_dict = self.preview_var(item)
                self.send_response(self.iopub_socket, 'display_data',
                    {'execution_count': self._execution_count, 'data': format_dict,
                    'metadata': md_dict})
            except Exception as e:
                self.warn('\n> Failed to preview file or expression {}'.format(item))

    def handle_magic_cd(self, option):
        # interpolate command
        option = self._interpolate_option(option, quiet=True)
        if option is None:
            return
        to_dir = option.strip()
        try:
            os.chdir(os.path.expanduser(to_dir))
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': os.getcwd()})
        except Exception as e:
            self.warn('Failed to change dir to {}: {}'.format(os.path.expanduser(to_dir), e))

    def handle_shell_command(self, cmd):
        # interpolate command
        cmd = self._interpolate_option(cmd, quiet=False)
        if cmd is None:
            return
        with self.redirect_sos_io():
            try:
                p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                out, err = p.communicate()
                sys.stdout.write(out.decode())
                sys.stderr.write(err.decode())
                # ret = p.returncode
                sys.stderr.flush()
                sys.stdout.flush()
            except Exception as e:
                sys.stderr.flush()
                sys.stdout.flush()
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': str(e)})

    def run_sos_code(self, code, silent):
        code = dedent(code)
        if os.path.isfile('.sos/report.md'):
            os.remove('.sos/report.md')
        with self.redirect_sos_io():
            try:
                # record input and output
                res = runfile(code=code, args=self.options)
                self.send_result(res, silent)
            except Exception:
                sys.stderr.flush()
                sys.stdout.flush()
                #self.send_response(self.iopub_socket, 'display_data',
                #    {
                #        'source': 'SoS',
                #        'metadata': {},
                #        'data': { 'text/html': HTML('<hr color="black" width="60%">').data}
                #    })
                raise
            except KeyboardInterrupt:
                self.warn('Keyboard Interrupt\n')
                return {'status': 'abort', 'execution_count': self._execution_count}
            finally:
                sys.stderr.flush()
                sys.stdout.flush()
        #
        if not silent and (not hasattr(self, 'preview_output') or self.preview_output):
            # Send standard output
            #if os.path.isfile('.sos/report.md'):
            #    with open('.sos/report.md') as sr:
            #        sos_report = sr.read()
                #with open(self.report_file, 'a') as summary_report:
                #    summary_report.write(sos_report + '\n\n')
            #    if sos_report.strip():
            #        self.send_response(self.iopub_socket, 'display_data',
            #            {
            #                'source': 'SoS',
            #                'metadata': {},
            #                'data': {'text/markdown': sos_report}
            #            })
            #
            if 'input' in env.sos_dict:
                input_files = env.sos_dict['input']
                if input_files is None:
                    input_files = []
            else:
                input_files = []
            if 'output' in env.sos_dict:
                output_files = env.sos_dict['output']
                if output_files is None:
                    output_files = []
            else:
                output_files = []
            # use a table to list input and/or output file if exist
            if output_files:
                self.send_response(self.iopub_socket, 'display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html': HTML('<pre><font color="green">## -- Preview output --</font></pre>').data}
                        })
                if hasattr(self, 'in_sandbox') and self.in_sandbox:
                    # if in sand box, do not link output to their files because these
                    # files will be removed soon.
                    self.send_response(self.iopub_socket, 'display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html':
                                HTML('''<pre> input: {}\noutput: {}\n</pre>'''.format(
                                ', '.join(x for x in input_files),
                                ', '.join(x for x in output_files))).data
                            }
                        })
                else:
                    self.send_response(self.iopub_socket, 'display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html':
                                HTML('''<pre> input: {}\noutput: {}\n</pre>'''.format(
                                ', '.join('<a target="_blank" href="{0}">{0}</a>'.format(x) for x in input_files),
                                ', '.join('<a target="_blank" href="{0}">{0}</a>'.format(x) for x in output_files))).data
                            }
                        })
                for filename in output_files:
                    self.preview_file(filename)

    def preview_var(self, item):
        if item in env.sos_dict:
            obj = env.sos_dict[item]
        else:
            obj = SoS_eval(item, sigil=get_default_global_sigil())
        if callable(obj) or isinstance(obj, ModuleType):
            return {'text/plain': pydoc.render_doc(obj, title='SoS Documentation: %s')}, {}
        else:
            return self.format_obj(obj)

    def preview_file(self, filename):
        if not os.path.isfile(filename):
            self.warn('\n> ' + filename + ' does not exist')
            return
        self.send_response(self.iopub_socket, 'display_data',
             {'metadata': {},
             'data': {
                 'text/plain': '\n> {} ({}):'.format(filename, pretty_size(os.path.getsize(filename))),
                 'text/html': HTML('<pre><font color="green">> {} ({}):</font></pre>'.format(filename, pretty_size(os.path.getsize(filename)))).data,
                }
             })
        previewer_func = None
        # lazy import of previewers
        if self.previewers is None:
            self.previewers = get_previewers()
        for x, y, _ in self.previewers:
            if isinstance(x, str):
                if fnmatch.fnmatch(os.path.basename(filename), x):
                    # we load entrypoint only before it is used. This is to avoid
                    # loading previewers that require additional external modules
                    # we can cache the loaded function but there does not seem to be
                    # a strong performance need for this.
                    previewer_func = y.load()
                    break
            else:
                # it should be a function
                try:
                    if x(filename):
                        try:
                            previewer_func = y.load()
                        except Exception as e:
                            self.warn('Failed to load previewer {}: {}'.format(y, e))
                            continue
                        break
                except Exception as e:
                    self.warn(e)
                    continue
        #
        # if no previewer can be found
        if previewer_func is None:
            return
        try:
            result = previewer_func(filename)
            if not result:
                return
            if isinstance(result, str):
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': result})
            elif isinstance(result, dict):
                self.send_response(self.iopub_socket, 'display_data',
                    {'source': filename, 'data': result, 'metadata': {}})
            else:
                self.warn('Unrecognized preview content: {}'.format(result))
                return
        except Exception as e:
            self.warn('Failed to preview {}: {}'.format(filename, e))

    def send_result(self, res, silent=False):
        # this is Ok, send result back
        if not silent and res is not None:
            format_dict, md_dict = self.format_obj(res)
            self.send_response(self.iopub_socket, 'execute_result',
                {'execution_count': self._execution_count, 'data': format_dict,
                'metadata': md_dict})

    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        # evaluate user expression
        ret = self._do_execute(code=code, silent=silent, store_history=store_history,
            user_expressions=user_expressions, allow_stdin=allow_stdin)

        out = {}
        for key, expr in (user_expressions or {}).items():
            try:
                #value = self.shell._format_user_obj(SoS_eval(expr, sigil=get_default_global_sigil()))
                value = SoS_eval(expr, sigil=get_default_global_sigil())
                value = self.shell._format_user_obj(value)
            except Exception as e:
                self.warn('Failed to evaluate user expression {}: {}'.format(expr, e))
                value = self.shell._user_obj_error()
            out[key] = value
        ret['user_expressions'] = out
        #
        self._execution_count += 1
        # make sure post_executed is triggered after the completion of all cell content
        self.shell.user_ns.update(env.sos_dict._dict)
        # trigger post processing of object and display matplotlib figures
        self.shell.events.trigger('post_execute')
        return ret

    def remove_leading_comments(self, code):
        lines = code.splitlines()
        try:
            idx = [x.startswith('#') or not x.strip() for x in lines].index(False)
            return os.linesep.join(lines[idx:])
        except Exception:
            # if all line is empty
            return ''

    def _do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        # if the kernel is SoS, remove comments and newlines
        code = self.remove_leading_comments(code)

        if self.original_keys is None:
            self._reset_dict()
        if code == 'import os\n_pid = os.getpid()':
            # this is a special probing command from vim-ipython. Let us handle it specially
            # so that vim-python can get the pid.
            return
        if self.MAGIC_DICT.match(code):
            # %dict should be the last magic
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_dict(options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_CONNECT_INFO.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            cfile = find_connection_file()
            with open(cfile) as conn:
                conn_info = conn.read()
            self.send_response(self.iopub_socket, 'stream',
                  {'name': 'stdout', 'text': 'Connection file: {}\n{}'.format(cfile, conn_info)})
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_MATPLOTLIB.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                self.shell.enable_gui = lambda gui: None
                gui, backend = self.shell.enable_matplotlib(options)
            except Exception as e:
                self.warn('Failed to set matplotlib backnd {}: {}'.format(options, e))
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_SET.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_set(options)
            # self.options will be set to inflence the execution of remaing_code
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_RESTART.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.restart_kernel(options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_WITH.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_with_parser()
                args = parser.parse_args(options.split())
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            original_kernel = self.kernel
            self.switch_kernel(args.kernel, args.out_vars)
            if args.in_vars:
                self.handle_magic_get(args.in_vars)
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.switch_kernel(original_kernel)
        elif self.MAGIC_USE.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_use_parser()
                args = parser.parse_args(options.split())
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'abort',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            self.switch_kernel(args.kernel, args.out_vars)
            if args.in_vars:
                self.handle_magic_get(args.in_vars)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_GET.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_get(options.split())
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_PUT.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_put(options.split())
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_PASTE.match(code):
            options, remaining_code = self.get_magic_and_code(code, True)
            try:
                old_options = self.options
                self.options = options + ' ' + self.options
                try:
                    code = clipboard_get()
                except ClipboardEmpty:
                    raise UsageError("The clipboard appears to be empty")
                except Exception as e:
                    env.logger.error('Could not get text from the clipboard: {}'.format(e))
                    return
                #
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': code.strip() + '\n## -- End pasted text --\n'})
                return self._do_execute(code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.options = old_options
        elif self.MAGIC_RUN.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            old_options = self.options
            self.options = options + ' ' + self.options
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.options = old_options
        elif self.MAGIC_RERUN.match(code):
            options, remaining_code = self.get_magic_and_code(code, True)
            old_options = self.options
            self.options = options + ' ' + self.options
            try:
                if not hasattr(self, 'last_executed_code'):
                    self.warn('No saved script')
                    self.last_executed_code = ''
                return self._do_execute(self.last_executed_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.options = old_options
        elif self.MAGIC_SANDBOX.match(code):
            import tempfile
            import shutil
            import shlex
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_sandbox_parser()
            args = parser.parse_args(shlex.split(options))
            self.in_sandbox = True
            try:
                old_dir = os.getcwd()
                if args.dir:
                    args.dir = os.path.expanduser(args.dir)
                    if not os.path.isdir(args.dir):
                        os.makedirs(args.dir)
                    env.exec_dir = os.path.abspath(args.dir)
                    os.chdir(args.dir)
                else:
                    new_dir = tempfile.mkdtemp()
                    env.exec_dir = os.path.abspath(new_dir)
                    os.chdir(new_dir)
                if not args.keep_dict:
                    old_dict = env.sos_dict
                    self._reset_dict()
                ret = self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
                if args.expect_error and ret['status'] == 'error':
                    #self.warn('\nSandbox execution failed.')
                    return {'status': 'ok', 
                        'payload': [], 'user_expressions': {},
                        'execution_count': self._execution_count}
                else:
                    return ret
            finally:
                if not args.keep_dict:
                    env.sos_dict = old_dict
                os.chdir(old_dir)
                if not args.dir:
                    shutil.rmtree(new_dir)
                self.in_sandbox = False
                #env.exec_dir = old_dir
        elif self.MAGIC_PREVIEW.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            import shlex
            parser = self.get_preview_parser()
            args = parser.parse_args(shlex.split(options, posix=False))
            if args.off:
                self.preview_output = False
            else:
                self.preview_output = True
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                if not args.off:
                    self.handle_magic_preview(options)
        elif self.MAGIC_CD.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_cd(options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('!'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_shell_command(code.split(' ')[0][1:] + ' ' + options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.kernel != 'sos':
            # handle string interpolation before sending to the underlying kernel
            if code:
                self.last_executed_code = code
            code = self._interpolate_option(code, quiet=False)
            if code is None:
                return
            try:
                return self.run_cell(code, store_history)
            except KeyboardInterrupt:
                self.warn('Keyboard Interrupt\n')
                return {'status': 'abort', 'execution_count': self._execution_count}
        else:
            if code:
                self.last_executed_code = code
            # run sos
            try:
                self.run_sos_code(code, silent)
                return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self._execution_count}
            except Exception as e:
                self.warn(str(e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            finally:
                # even if something goes wrong, we clear output so that the "preview"
                # will not be viewed by a later step.
                env.sos_dict.pop('input', None)
                env.sos_dict.pop('output', None)

    def do_shutdown(self, restart):
        #
        for name, (km, kv) in self.kernels.items():
            try:
                km.shutdown_kernel(restart=restart)
            except Exception as e:
                self.warn('Failed to shutdown kernel {}: {}'.format(name, e))


if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
