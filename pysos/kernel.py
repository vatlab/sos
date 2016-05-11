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
import time
from .utils import env, WorkflowDict
from ._version import __sos_version__, __version__
from .sos_eval import SoS_exec, SoS_eval
from .sos_executor import Interactive_Executor
from .sos_syntax import SOS_SECTION_HEADER

from IPython.lib.clipboard import ClipboardEmpty, osx_clipboard_get, tkinter_clipboard_get
from IPython.core.error import UsageError
from ipykernel.kernelbase import Kernel
from jupyter_client import manager
from IPython.utils import io

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
        self.kernel = None
        self.kernels = {}

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
        if any(code.startswith(x) for x in ['#dict', '#paste']):
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

    def handle_iopub(self, msg_id=''):
        """Process messages on the IOPub channel

           This method consumes and processes messages on the IOPub channel,
           such as stdout, stderr, execute_result and status.

           It only displays output that is caused by this session.
        """
        while self.KC.iopub_channel.msg_ready():
            sub_msg = self.KC.iopub_channel.get_msg()
            msg_type = sub_msg['header']['msg_type']
            parent = sub_msg["parent_header"]
            if msg_type == 'status':
                self._execution_state = sub_msg["content"]["execution_state"]
            elif msg_type == 'stream':
                if sub_msg["content"]["name"] == "stdout":
                    print(sub_msg["content"]["text"], file=io.stdout, end="")
                    io.stdout.flush()
                elif sub_msg["content"]["name"] == "stderr":
                    print(sub_msg["content"]["text"], file=io.stderr, end="")
                    io.stderr.flush()
            elif msg_type == 'execute_result':
                self.execution_count = int(sub_msg["content"]["execution_count"])
                data = sub_msg["content"]["data"]
                if 'text/plain' in data:
                    print(data['text/plain'])
                else:
                    print('{} not handled'.format(', '.join(data.keys())))
                #self.handle_rich_data(format_dict)
                # taken from DisplayHook.__call__:
                # This is currently not handled.
                #hook = self.displayhook
                #hook.start_displayhook()
                #hook.write_output_prompt()
                #hook.write_format_data(format_dict)
                #hook.log_output(format_dict)
                #hook.finish_displayhook()
            elif msg_type == 'display_data':
                data = sub_msg["content"]["data"]
                #handled = self.handle_rich_data(data)
                #if not handled:
                if 'text/plain' in data:
                    print(data['text/plain'])
            elif msg_type == 'execute_input':
                pass
            elif msg_type == 'clear_output':
                print("\r", file=io.stdout, end="")
            elif msg_type == 'error':
                for frame in sub_msg["content"]["traceback"]:
                    print(frame, file=io.stderr)

    def run_cell(self, code, store_history):
        #
        # flush stale replies, which could have been ignored, due to missed heartbeats
        while self.KC.shell_channel.msg_ready():
            self.KC.shell_channel.get_msg()
        # executing code in another kernel
        msg_id = self.KC.execute(code, silent=False, store_history=not store_history)

        # first thing is wait for any side effects (output, stdin, etc.)
        self._executing = True
        self._execution_state = "busy"
        while self._execution_state != 'idle':
            # display intermediate print statements, etc.
            self.handle_iopub(msg_id)
            time.sleep(0.01)

        reply = self.KC.get_shell_msg(timeout=10)
        return reply['content']


    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        if code == 'import os\n_pid = os.getpid()':
            # this is a special probing command from vim-ipython. Let us handle it specially
            # so that vim-python can get the pid.
            return {'status': 'ok',
                # The base class increments the execution count
                'execution_count': None,
                'payload': [],
                'user_expressions': {'_pid': {'data': {'text/plain': os.getpid()}}}
               }
        mode = 'code'
        if code.startswith('#dict'):
            mode = 'dict'
            command_line = self.get_magic_option(code)
        elif code.startswith('#set'):
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
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options
        elif code.startswith('#kernel'):
            options = self.get_magic_option(code)
            if options == 'sos':
                options = ''
            #
            if options.strip():
                if options != self.kernel:
                    if options in self.kernels:
                        print('Using kernel "{}"'.format(options))
                        self.KM, self.KC = self.kernels[options]
                        self.kernel = options
                        self.handle_iopub()
                    else:
                        try:
                            self.kernels[options] = manager.start_new_kernel(startup_timeout=60, kernel_name=options)
                            print('Kernel "{}" started'.format(options))
                            self.KM, self.KC = self.kernels[options]
                            self.kernel = options
                            self.handle_iopub()
                        except:
                            print('Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed.'.format(options))
            else:
                if self.kernel:
                    print('switching back to sos kernel')
                    self.kernel = ''
                else:
                    print('Usage: switch current kernel to another Jupyter kernel (e.g. ir for R)')
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options
        elif code.startswith('#paste'):
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
        elif code.startswith('#run'):
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options + ' ' + self.get_magic_option(code)
        else:
            command_line = self.options
        #
        try:
            if mode == 'dict':
                res = self.sosdict(command_line)
            elif self.kernel:
                return self.run_cell(code, store_history)
            else:
                res = self.executor.run_interactive(code, command_line)
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
