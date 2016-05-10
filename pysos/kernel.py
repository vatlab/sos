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

import sys
from .utils import env, WorkflowDict
from ._version import __sos_version__, __version__
from .sos_eval import SoS_exec, SoS_eval
from .sos_executor import Interactive_Executor
from .sos_syntax import SOS_SECTION_HEADER

from IPython.lib.clipboard import ClipboardEmpty, osx_clipboard_get, tkinter_clipboard_get
from IPython.core.error import UsageError
from ipykernel.kernelbase import Kernel

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


class SoS_Kernel(Kernel):
    implementation = 'SoS'
    implementation_version = __version__
    language = 'sos'
    language_version = __sos_version__
    language_info = {
        'mimetype': 'text/x-sos',
        'name': 'sos',
        'file_extension': '.sos',
        'pygments_lexer': 'sos',
        'codemirror_mode': 'sos',
    }
    banner = "SoS kernel - script of scripts"

    def __init__(self, **kwargs):
        super(SoS_Kernel, self).__init__(**kwargs)
        self._start_sos()

    def _start_sos(self):
        env.sos_dict = WorkflowDict()
        SoS_exec('from pysos import *')
        env.sos_dict.set('__interactive__', True)
        self.executor = Interactive_Executor()
        self.original_keys = set(env.sos_dict._dict.keys())
        self.original_keys.add('__builtins__')
        self.options = ''

    def eval_code(self, code, line):
        code = code.strip()
        try:
            # is it an expression?
            compile(code, '<string>', 'eval')
            #if line.strip():
            #    env.logger.warning('{} ignored for expression evaluation'.format(line))
            return SoS_eval(code)
        except:
            # is it a list of statement?
            try:
                compile(code, '<string>', 'exec')
                #if line.strip():
                #    env.logger.warning('{} ignored for statement execution'.format(line))
                return SoS_exec(code)
            except:
                return self.executor.run_interactive(code, command_line=line)

    def do_inspect(self, code, cursor_pos, detail_level=0):
        'Inspect code'
        return {
            'status': 'ok',
            'found': 'true',
            'data': {x:y for x,y in env.sos_dict._dict.items() if x not in self.original_keys and not x.startswith('__')},
            'metadata':''}

    def sosdict(self, line):
        'Magic that displays content of the dictionary'
        # do not return __builtins__ beacuse it is too long...
        actions = line.strip().split()
        for action in actions:
            if action not in ['reset', 'all', 'keys']:
                raise RuntimeError('Unrecognized sosdict option {}'.format(action))
        if 'reset' in actions:
            return self._reset()
        if 'keys' in actions:
            if 'all' in actions:
                return env.sos_dict._dict.keys()
            else:
                return {x for x in env.sos_dict._dict.keys() if not x.startswith('__')} - self.original_keys
        else:
            if 'all' in actions:
                return env.sos_dict._dict
            else:
                return {x:y for x,y in env.sos_dict._dict.items() if x not in self.original_keys and not x.startswith('__')}

    def do_is_complete(self, code):
        '''check if new line is in order'''
        code = code.strip()
        if not code:
            return {'status': 'complete', 'indent': ''}
        if any(code.startswith(x) for x in ['#sosdict', '#sospaste']):
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
        # check syntax??
        try:
            compile(code, '<string>', 'exec')
            return {'status': 'complete', 'indent': ''}
        except:
            try:
                self.executor.parse_script(code)
                return {'status': 'complete', 'indent': ''}
            except:
                return {'status': 'unknown', 'indent': ''}

    def get_magic_option(self, code):
        lines = code.split('\n')
        pieces = lines[0].strip().split(None, 1)
        if len(pieces) == 2:
            command_line = pieces[1]
        else:
            command_line = ''
        return command_line

    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        mode = 'code'
        if code.startswith('#sosdict'):
            mode = 'dict'
            command_line = self.get_magic_option(code)
        elif code.startswith('#sosset'):
            options = self.get_magic_option(code)
            if options.strip():
                print('sos options set to "{}"'.format(options))
                self.options = options.strip()
            else:
                if self.options:
                    print('sos options "{}" reset to ""'.format(self.options))
                    self.options = ''
                else:
                    print('Usage: set persistent sos options such as -v 3 (debug output) -p (prepare) and -t (transcribe)')
        elif code.startswith('#sospaste'):
            command_line = self.options + ' ' + self.get_magic_option(code)
            try:
                code = clipboard_get()
            except ClipboardEmpty:
                raise UsageError("The clipboard appears to be empty")
            except Exception as e:
                env.logger.error('Could not get text from the clipboard: {}'.format(e))
                return
            #
            print(code.strip())
            print('## -- End pasted text --')
        elif code.startswith('#sosrun'):
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options + ' ' + self.get_magic_option(code)
        else:
            command_line = self.options
        #
        try:
            if mode == 'dict':
                res = self.sosdict(command_line)
            else:
                res = self.eval_code(code, command_line)
        except Exception as e:
            stream_content = {'name': 'stderr', 'text': repr(e)}
            self.send_response(self.iopub_socket, 'stream', stream_content)
            return  {'status': 'error',
                'ename': e.__class__.__name__,
                'evalue': repr(e),
                'traceback': [],
                'execution_count': self.execution_count,
               }

        # this is Ok, send result back
        if not silent and res is not None:
            stream_content = {'name': 'stdout', 'text': repr(res)}
            self.send_response(self.iopub_socket, 'stream', stream_content)
        #
        return {'status': 'ok',
                # The base class increments the execution count
                'execution_count': self.execution_count,
                'payload': [],
                'user_expressions': {}, #res,
               }

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
