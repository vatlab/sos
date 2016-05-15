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
import base64
import imghdr
import contextlib

from .utils import env, WorkflowDict, shortRepr
from .signature import textMD5, FileInfo
from ._version import __sos_version__, __version__
from .sos_eval import SoS_exec, interpolate
from .sos_executor import Interactive_Executor
from .sos_syntax import SOS_SECTION_HEADER

from IPython.lib.clipboard import ClipboardEmpty, osx_clipboard_get, tkinter_clipboard_get
from IPython.core.error import UsageError
from IPython.core.display import HTML
from ipykernel.kernelbase import Kernel
from jupyter_client import manager

from io import StringIO

from nbconvert.exporters import Exporter

@contextlib.contextmanager
def redirect_sos_io():
    save_stdout = sys.stdout
    save_stderr = sys.stderr
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    yield
    sys.stdout = save_stdout
    sys.stderr = save_stderr

__all__ = ['SoS_Exporter', 'SoS_Kernel']

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

def is_image(filename):
    with open(filename, 'rb') as f:
        image = f.read()
    return imghdr.what(None, image) is not None


def display_data_for_image(filename):
    with open(filename, 'rb') as f:
        image = f.read()

    image_type = imghdr.what(None, image)
    if image_type is None:
        raise ValueError("Not a valid image: %s" % image)

    image_data = base64.b64encode(image).decode('ascii')
    content = {
        'data': {
            'image/' + image_type: image_data
        },
        'metadata': {}
    }
    return content

class SoS_Exporter(Exporter):
    def __init__(self, config=None, **kwargs):
        self.output_extension = '.sos'
        self.output_mimetype = 'text/x-sos'
        Exporter.__init__(self, config, **kwargs)

    def from_notebook_cell(self, cell, fh):
        if not hasattr(cell, 'execution_count') or cell.execution_count is None:
            fh.write('\n#cell {}\n'.format(cell.cell_type))
        else:
            fh.write('\n#cell {} {}\n'.format(cell.cell_type, cell.execution_count))
        if cell.cell_type == 'code':
            if any(cell.source.startswith(x) for x in ('%run', '%restart', '%dict', '%use', '%with', '%set', '%paste')):
                env.logger.warning('SoS magic "{}" has to remove them before executing the script with sos command.'.format(cell.source.split('\n')[0]))
            fh.write(cell.source + '\n')
        elif cell.cell_type == "markdown":
            fh.write('\n'.join('! ' + x for x in cell.source.split('\n')) + '\n')

    def from_notebook_node(self, nb, resources, **kwargs):
        #
        with StringIO() as fh:
            fh.write('#!/usr/bin/env sos-runner\n')
            fh.write('#fileformat=SOSNB1.0\n')
            for cell in nb.cells:
                self.from_notebook_cell(cell, fh)
            content = fh.getvalue()
        resources['output_extension'] = '.sos'
        return content, resources

class SoS_Kernel(Kernel):
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
        'nbconvert_exporter': 'pysos.kernel.SoS_Exporter',
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
        self.kernel = 'sos'
        # FIXME: this should in theory be a MultiKernelManager...
        self.kernels = {}
        self.original_kernel = None

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
        if any(code.startswith(x) for x in ['%dict', '%paste']):
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
        msg_id = self.KC.execute(code, silent=False, store_history=not store_history)

        # first thing is wait for any side effects (output, stdin, etc.)
        _execution_state = "busy"
        while _execution_state != 'idle':
            # display intermediate print statements, etc.
            while self.KC.iopub_channel.msg_ready():
                sub_msg = self.KC.iopub_channel.get_msg()
                msg_type = sub_msg['header']['msg_type']
                if msg_type == 'execute_input':
                    continue
                elif msg_type == 'status':
                    _execution_state = sub_msg["content"]["execution_state"]
                # pass along all messages except for execute_input
                #self.send_response(self.iopub_socket, 'stream',
                #            {'name': 'stdout', 'text': sub_msg + '\n'})
                self.send_response(self.iopub_socket, msg_type, sub_msg['content'])
        # now get the real result
        reply = self.KC.get_shell_msg(timeout=10)
        # FIXME: not sure if other part of the reply is useful...
        return reply['content']

    def switch_kernel(self, kernel):
        if kernel and kernel != 'sos':
            if kernel != self.kernel:
                if kernel in self.kernels:
                    #self.send_response(self.iopub_socket, 'stream',
                    #    {'name': 'stdout', 'text': 'Using kernel "{}"\n'.format(kernel)})
                    self.KM, self.KC = self.kernels[kernel]
                    self.kernel = kernel
                else:
                    try:
                        self.kernels[kernel] = manager.start_new_kernel(startup_timeout=60, kernel_name=kernel)
                        #self.send_response(self.iopub_socket, 'stream',
                        #    {'name': 'stdout', 'text': 'Kernel "{}" started\n'.format(kernel)})
                        self.KM, self.KC = self.kernels[kernel]
                        self.kernel = kernel
                    except Exception as e:
                        self.send_response(self.iopub_socket, 'stream',
                            {'name': 'stdout', 'text': 'Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}\n'.format(kernel, e)})
        else:
            # kerl is None ('') or kernel == 'sos'
            if self.kernel != 'sos':
                #self.send_response(self.iopub_socket, 'stream',
                #    {'name': 'stdout', 'text': 'Switching back to sos kernel\n'})
                self.kernel = 'sos'
            elif kernel == '':
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Usage: switch current kernel to another Jupyter kernel (e.g. R or ir for R)\n'})

    def restart_kernel(self, kernel):
        if kernel == 'sos':
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': 'Please use Kernel -> Restart to restart SoS kernel'})
        elif kernel:
            if kernel in self.kernels:
                try:
                    self.kernels[kernel][0].shutdown_kernel(restart=False)
                except Exception as e:
                    env.logger.warning('Failed to shutdown kernel {}: {}'.format(kernel, e))
            #
            try:
                self.kernels[kernel] = manager.start_new_kernel(startup_timeout=60, kernel_name=kernel)
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Kernel {} {}started'.format(kernel, 're' if kernel in self.kernels else '')})
                #self.send_response(self.iopub_socket, 'stream',
                #    {'name': 'stdout', 'text': 'Kernel "{}" started\n'.format(kernel)})
                if kernel == self.kernel:
                    self.KM, self.KC = self.kernels[kernel]
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}\n'.format(kernel, e)})


    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        self.execution_count += 1
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
        if code.startswith('%dict'):
            mode = 'dict'
            command_line = self.get_magic_option(code)
        elif code.startswith('%set'):
            options = self.get_magic_option(code)
            if options.strip():
                #self.send_response(self.iopub_socket, 'stream',
                #    {'name': 'stdout', 'text': 'sos options set to "{}"\n'.format(options)})
                self.options = options.strip()
            else:
                if self.options:
                    #self.send_response(self.iopub_socket, 'stream',
                    #    {'name': 'stdout', 'text': 'sos options "{}" reset to ""\n'.format(self.options)})
                    self.options = ''
                else:
                    self.send_response(self.iopub_socket, 'stream',
                        {'name': 'stdout', 'text': 'Usage: set persistent sos options such as -v 3 (debug output) -p (prepare) and -t (transcribe)\n'})
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options
        elif code.startswith('%restart'):
            options = self.get_magic_option(code)
            if options == 'R':
                options = 'ir'
            self.restart_kernel(options)
            return {'status': 'ok',
                    # The base class increments the execution count
                    'execution_count': self.execution_count,
                    'payload': [],
                    'user_expressions': {},
                   }
        elif code.startswith('%with') or code.startswith('%use'):
            options = self.get_magic_option(code)
            if options == 'R':
                options = 'ir'
            #
            if code.startswith('%with'):
                self.original_kernel = self.kernel
            self.switch_kernel(options)
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options
        elif code.startswith('%paste'):
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
        elif code.startswith('%run'):
            lines = code.split('\n')
            code = '\n'.join(lines[1:])
            command_line = self.options + ' ' + self.get_magic_option(code)
        else:
            command_line = self.options
        #
        try:
            try:
                if mode == 'dict':
                    res = self.sosdict(command_line)
                elif self.kernel != 'sos':
                    # handle string interpolation before sending to the underlying kernel
                    try:
                        new_code = interpolate(code, sigil='${ }', local_dict=env.sos_dict._dict)
                        if new_code != code:
                            print(new_code.strip())
                            code = new_code
                            self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stdout', 'text':
                                new_code.strip() + '\n## -- End interpolated text --\n'})
                    except Exception as e:
                        self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stdout', 'text': 'Failed to interpolate {}: {}'.format(shortRepr(code), e)})
                    return self.run_cell(code, store_history)
                else:
                    with redirect_sos_io():
                        try:
                            res = self.executor.run_interactive(code, command_line)
                            sos_out = sys.stdout.getvalue()
                            sos_err = sys.stderr.getvalue()
                        except KeyboardInterrupt:
                            return {'status': 'abort', 'execution_count': self.execution_count}
                    #
                    if not silent:
                        # Send standard output
                        if sos_out:
                            self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stdout', 'text': sos_out})
                        if sos_err:
                            self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stderr', 'text': sos_err})
                        #
                        if '__step_input__' in env.sos_dict:
                            input_files = env.sos_dict['__step_input__']
                        else:
                            input_files = []
                        if '__step_output__' in env.sos_dict:
                            output_files = env.sos_dict['__step_output__']
                            if output_files is None:
                                output_files = []
                        else:
                            output_files = []
                        # use a table to list input and/or output file if exist
                        start_out = True
                        if input_files or output_files:
                            self.send_response(self.iopub_socket, 'display_data',
                                    {'data': { 'text/html':
                                        HTML('''<pre> input: {}\noutput: {}\n</pre>'''.format(
                                            ', '.join('<a href="{0}">{0}</a>'.format(x) for x in input_files),
                                            ', '.join('<a href="{0}">{0}</a>'.format(x) for x in output_files))).data
                                        }
                                    })
                        # Send images, if any
                        for filename in output_files:
                            if is_image(filename):
                                self.send_response(self.iopub_socket, 'stream',
                                     {'name': 'stdout', 'text': '> ' + filename + ' ({:.1f} KB)'.format(os.path.getsize(filename) / 1024)})
                                data = display_data_for_image(filename)
                                self.send_response(self.iopub_socket, 'display_data', data)
                                start_out = False
                            elif filename.lower().endswith('.pdf'):
                                self.send_response(self.iopub_socket, 'stream',
                                         {'name': 'stdout', 'text': '> ' + filename + ' ({:.1f} KB)'.format(os.path.getsize(filename) / 1024)})
                                self.send_response(self.iopub_socket, 'display_data',
                                    {'source': filename,
                                     'data': { 'text/html': HTML('<iframe src={0} width="100%"></iframe>'.format(filename)).data}})
                                start_out = False
                            elif filename.lower().endswith('.csv') or filename.lower().endswith('.tsv'):
                                try:
                                    import pandas
                                    data = pandas.read_csv(filename, nrows=10)
                                    html = data.to_html()
                                    self.send_response(self.iopub_socket, 'stream',
                                         {'name': 'stdout', 'text': '> ' +  filename + '(first 10 rows)\n'})
                                    self.send_response(self.iopub_socket, 'display_data',
                                        {'source': filename,
                                         'data': { 'text/html':
                                         HTML(html).data}})
                                    start_out = False
                                except Exception as e:
                                    #env.logger.warning('Cannot display {}: {}'.format(filename, e))
                                    self.send_response(self.iopub_socket, 'stream',
                                         {'name': 'stdout', 'text': repr(e)})
                            else:
                                fi = FileInfo(filename)
                                desc = fi.describe()
                                if desc:
                                    self.send_response(self.iopub_socket, 'stream',
                                         {'name': 'stdout', 'text': '> ' + filename + ' ({:.1f} KB)'.format(os.path.getsize(filename) / 1024)})
                                    self.send_response(self.iopub_socket, 'display_data',
                                        {'source': filename,
                                         'data': {'text/html':
                                         HTML('<pre>{}</pre>'.format(desc)).data}})
                                    start_out = False

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
                    'user_expressions': {},
                   }
        finally:
            if self.original_kernel is not None:
                self.switch_kernel(self.original_kernel)
                self.original_kernel = None

    def do_shutdown(self, restart):
        #
        for name, (km, kv) in self.kernels.items():
            try:
                km.shutdown_kernel(restart=restart)
            except Exception as e:
                env.logger.warning('Failed to shutdown kernel {}: {}'.format(name, e))

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
