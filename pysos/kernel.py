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
import copy
import fnmatch
import contextlib
import subprocess

import zipfile
import tarfile
import gzip
import pickle
from collections.abc import Sequence

from .utils import env, WorkflowDict, short_repr, pretty_size, dehtml
from ._version import __sos_version__, __version__
from .sos_eval import SoS_exec, interpolate
from .sos_executor import Interactive_Executor
from .sos_syntax import SOS_SECTION_HEADER
from .converter import SoS_Exporter

from IPython.core.interactiveshell import InteractiveShell
from IPython.lib.clipboard import ClipboardEmpty, osx_clipboard_get, tkinter_clipboard_get
from IPython.core.error import UsageError
from IPython.core.display import HTML
from ipykernel.kernelbase import Kernel
from jupyter_client import manager, find_connection_file
from ipykernel.zmqshell import ZMQDisplayPublisher

from textwrap import dedent
from io import StringIO


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
        self.kernel.send_response(self.kernel.iopub_socket, 'stream',
            {'name': self.name, 'text': content})
        self.truncate(0)
        self.seek(0)
        self.nlines = 0
        return len(content.strip())

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

class SoS_FilePreviewer():
    def __init__(self):
        pass

    def display_data_for_image(self, filename):
        with open(filename, 'rb') as f:
            image = f.read()

        image_type = imghdr.what(None, image)
        image_data = base64.b64encode(image).decode('ascii')
        if image_type != 'png':
            try:
                from wand.image import Image
                img = Image(filename=filename)
                return { 'image/' + image_type: image_data,
                    'image/png': base64.b64encode(img._repr_png_()).decode('ascii') }
            except Exception:
                return { 'image/' + image_type: image_data }
        else:
            return { 'image/' + image_type: image_data }

    def preview(self, filename):
        if imghdr.what(filename) is not None:
            # image
            return 'display_data', self.display_data_for_image(filename)
        elif filename.lower().endswith('.pdf'):
            try:
                # this import will fail even if wand is installed
                # if imagemagick is not installed properly.
                from wand.image import Image
                img = Image(filename=filename)
                return 'display_data', {
                    'text/html': HTML('<iframe src={0} width="100%"></iframe>'.format(filename)).data,
                    'image/png': base64.b64encode(img._repr_png_()).decode('ascii') }
            except Exception as e:
                env.logger.error(e)
                return 'display_data', { 'text/html':
                    HTML('<iframe src={0} width="100%"></iframe>'.format(filename)).data}
        elif filename.lower().endswith('.html'):
            with open(filename) as html:
                content = html.read()
            return 'display_data', { 'text/html': content,
                'text/plain': dehtml(content) }
        elif filename.lower().endswith('.csv') or filename.lower().endswith('.tsv'):
            try:
                import pandas
                data = pandas.read_csv(filename, nrows=10)
                html = data.to_html()
                return 'display_data', { 'text/html': HTML(html).data}
            except Exception:
                pass
        # is it a compressed file?
        if zipfile.is_zipfile(filename):
            zip = zipfile.ZipFile(filename)
            names = zip.namelist()
            return '{} files\n'.format(len(names)) + '\n'.join(names[:5]) + ('\n...' if len(names) > 5 else '')
        elif tarfile.is_tarfile(filename):
            with tarfile.open(filename, 'r:*') as tar:
                # only extract files
                names = [x.name for x in tar.getmembers() if x.isfile()]
            return '{} files\n'.format(len(names)) + '\n'.join(names[:5]) + ('\n...' if len(names) > 5 else '')
        elif filename.endswith('.gz'):
            content = b''
            with gzip.open(filename, 'rb') as fin:
                for line in range(5):
                    content += fin.readline()
            try:
                return content.decode()
            except:
                return 'binary data'
        else:
            content = b''
            with open(filename, 'rb') as fin:
                for line in range(5):
                    content += fin.readline()
            try:
                return content.decode()
            except:
                pass
        return 'binary data'

class BioPreviewer(SoS_FilePreviewer):
    def  __init__(self):
        SoS_FilePreviewer.__init__(self)

    def previewBam(self, filename):
        try:
            import pysam
        except ImportError:
            return 'pysam is needed to preview bam format'
        try:
            res = ''
            with pysam.AlignmentFile(filename, 'rb') as bam:
                headers = bam.header
                for record_type in ('RG', 'PG', 'SQ'):
                    if record_type not in headers:
                        continue
                    else:
                        records = headers[record_type]
                    res += record_type + ':\n'
                    for i, record in enumerate(records):
                        if type(record) == str:
                            res += '  ' + short_repr(record) + '\n'
                        elif type(record) == dict:
                            res += '  '
                            for idx, (k, v) in enumerate(record.items()):
                                if idx < 4:
                                    res += '{}: {}    '.format(k, short_repr(v))
                                elif idx == 4:
                                    res += '...'
                                    break
                        if i > 4:
                            res += '\n  ...\n'
                            break
                        else:
                            res += '\n'
            return res
        except Exception as e:
            return 'failed to preview {}'.format(e)

    def preview(self, filename):
        if filename.lower().endswith('.bam'):
            return self.previewBam(filename)
        else:
            return SoS_FilePreviewer.preview(filename)

def R_repr(obj):
    if isinstance(obj, (int, float, str)):
        return repr(obj)
    elif isinstance(obj, Sequence):
        return 'c(' + ','.join(R_repr(x) for x in obj) + ')'
    else:
        import pandas
        if isinstance(obj, pandas.DataFrame):
            try:
                import feather
            except ImportError:
                raise UsageError('The feather-format module is required to pass pandas DataFrame as R data.frame'
                    'See https://github.com/wesm/feather/tree/master/python for details.')
            feather_tmp_ = os.path.join(env.exec_dir, '.sos', '__data.feather')
            if os.path.isfile(feather_tmp_):
                os.remove(feather_tmp_)
            feather.write_dataframe(obj, feather_tmp_)
            return 'read_feather("{}")'.format(feather_tmp_)
        else:
            raise UsageError('Passing {} from Python to R is not yet supported.'.format(short_repr(obj)))

# an R function that tries to convert an R object to Python repr
dputToString_func = '''
dputToString <- function (obj) {
  con <- textConnection(NULL,open="w")
  tryCatch({dput(obj,con);
           textConnectionValue(con)},
           finally=close(con))
}
'''

def from_R_repr(expr):
    '''
    Convert text in the following format to a Python object

    '[1] "structure(list(a = 1, b = \\"a\\"), .Names = c(\\"a\\", \\"b\\"))"'
    '''
    try:
        text = expr[20:-1].split('), .Names')[0].replace('\\', '')
        convert_funcs = {
            'c': lambda *args: list(args)
        }
        return eval('dict({})'.format(text), convert_funcs)
    except Exception as e:
        raise UsageError('Failed to convert {} to Python object: {}'.format(expr, e))

def python2R(name):
    #
    # return equivalent R representation for python object.
    # This is limited to select python objects but we will obtain
    # native R object using this method.
    if name not in env.sos_dict:
        raise UsageError('{} not exist'.format(name))
    try:
        r_repr = R_repr(env.sos_dict[name])
        # there can be some more general solution but right now let us just
        # do this.
        if 'read_feather' in r_repr:
            return 'library(feather)\n{} <- {}'.format(name, r_repr)
        else:
            return '{} <- {}'.format(name, r_repr)
    except Exception as e:
        raise UsageError('Failed to convert variable {} to R: {}'.format(name, e))


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
    shell = InteractiveShell.instance()

    def __init__(self, **kwargs):
        super(SoS_Kernel, self).__init__(**kwargs)
        self.options = ''
        self.kernel = 'sos'
        self.banner = self.banner + '\nConnection file {}'.format(os.path.basename(find_connection_file()))
        # FIXME: this should in theory be a MultiKernelManager...
        self.kernels = {}
        self.original_kernel = None
        self.format_obj = self.shell.display_formatter.format

        # InteractiveShell uses a default publisher that only displays text/plain
        # using the ZMQDisplayPublisher will display matplotlib inline figures
        self.shell.display_pub = ZMQDisplayPublisher()
        self.shell.display_pub.session = self.session
        self.shell.display_pub.pub_socket = self.iopub_socket

        self.shell.enable_gui = lambda x: False
        self.previewer = {'*': SoS_FilePreviewer().preview, '*.bam': BioPreviewer().preview }

        self.report_file = os.path.join(env.exec_dir, 'summary_report.md')
        if os.path.isfile(self.report_file):
            os.remove(self.report_file)
        # touch the file
        with open(self.report_file, 'w'):
            pass
        self.original_keys = None

    def _reset_dict(self):
        env.sos_dict = WorkflowDict()
        SoS_exec('import os, sys, glob')
        SoS_exec('from pysos.runtime import *')
        SoS_exec("run_mode = 'interactive'")
        self.executor = Interactive_Executor()
        self.original_keys = set(env.sos_dict._dict.keys())
        self.original_keys.add('__builtins__')
        env.sos_dict.set('__summary_report__', self.report_file)

    @contextlib.contextmanager
    def redirect_sos_io(self):
        save_stdout = sys.stdout
        save_stderr = sys.stderr
        sys.stdout = FlushableStringIO(self, 'stdout')
        sys.stderr = FlushableStringIO(self, 'stderr')
        yield
        sys.stdout = save_stdout
        sys.stderr = save_stderr

    #
    # Right now we are not sure what to return for do_inspect
    # http://jupyter-client.readthedocs.io/en/latest/messaging.html
    # returnning sos_dict is wasteful and is prone to error (unpickleable objects
    # causing error messages)
    #
    #def do_inspect(self, code, cursor_pos, detail_level=0):
    #    'Inspect code'
    #    # x:y for x,y in env.sos_dict._dict.items() if x not in self.original_keys and not x.startswith('__')},
    #    return {
    #        'status': 'ok',
    #        'found': 'true',
    #        'data': {}, #
    #        'metadata': {}}

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

    def get_magic_and_code(self, code, warn_remaining=False):
        lines = code.split('\n')
        pieces = lines[0].strip().split(None, 1)
        if len(pieces) == 2:
            command_line = pieces[1]
        else:
            command_line = ''
        remaining_code = '\n'.join(lines[1:])
        if warn_remaining and remaining_code.strip():
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stderr', 'text': 
                'Statement {} ignored'.format(short_repr(remaining_code))})
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
                elif msg_type in ('execute_input', 'execute_result'):
                    # override execution count with the master count,
                    # not sure if it is needed
                    sub_msg['content']['execution_count'] = self.execution_count
                self.send_response(self.iopub_socket, msg_type, sub_msg['content'])
        # now get the real result
        reply = self.KC.get_shell_msg(timeout=10)
        reply['content']['execution_count'] = self.execution_count
        return reply['content']

    def switch_kernel(self, kernel):
        if kernel and kernel != 'sos':
            if kernel != self.kernel:
                if kernel in self.kernels:
                    self.KM, self.KC = self.kernels[kernel]
                    self.kernel = kernel
                else:
                    try:
                        self.kernels[kernel] = manager.start_new_kernel(startup_timeout=60, kernel_name=kernel)
                        self.KM, self.KC = self.kernels[kernel]
                        self.kernel = kernel
                    except Exception as e:
                        self.send_response(self.iopub_socket, 'stream',
                            {'name': 'stderr', 'text': 'Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}\n'.format(kernel, e)})
        else:
            # kernel == '' or kernel == 'sos'
            if kernel == '':
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Kernel "{}" is used.\n'.format(self.kernel)})
            elif kernel == 'sos':
                #self.send_response(self.iopub_socket, 'stream',
                #    {'name': 'stdout', 'text': 'Switching back to sos kernel\n'})
                self.kernel = 'sos'

    def restart_kernel(self, kernel):
        if kernel == 'sos':
            # cannot restart myself ...
            self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Cannot restart sos kernel from within sos.'})
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

    def handle_magic_dict(self, line):
        'Magic that displays content of the dictionary'
        # do not return __builtins__ beacuse it is too long...
        actions = line.strip().split()
        for action in actions:
            if action not in ['reset', 'all', 'keys']:
                self.send_response(self.iopub_socket, 'stream',
                      {'name': 'stderr', 'text': 'Unrecognized sosdict option {}'.format(action)})
                return
        if 'reset' in actions:
            self.send_result(self._reset_dict())
        if 'keys' in actions:
            if 'all' in actions:
                self.send_result(env.sos_dict._dict.keys())
            else:
                self.send_result({x for x in env.sos_dict._dict.keys() if not x.startswith('__')} - self.original_keys)
        else:
            if 'all' in actions:
                self.send_result(env.sos_dict._dict)
            else:
                self.send_result({x:y for x,y in env.sos_dict._dict.items() if x not in self.original_keys and not x.startswith('__')})

    def handle_magic_set(self, options):
        if options.strip():
            #self.send_response(self.iopub_socket, 'stream',
            #    {'name': 'stdout', 'text': 'sos options set to "{}"\n'.format(options)})
            self.options = options.strip()
        else:
            if self.options:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'sos options "{}" reset to ""\n'.format(self.options)})
                self.options = ''
            else:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Usage: set persistent sos options such as -v 3 (debug output) -p (prepare) and -t (transcribe)\n'})
  
            self.switch_kernel(options)

    def handle_magic_get(self, options):
        items = options.split()
        for item in items:
            if item not in env.sos_dict:
                self.send_response(self.iopub_socket, 'stream',
                     {'name': 'stderr', 'text': 'Variable {} does not exist'.format(item)})
                return
        if self.kernel == 'python':
            # if it is a python kernel, passing specified SoS variables to it
            sos_data = pickle.dumps({x:env.sos_dict[x] for x in items})
            # this can fail if the underlying python kernel is python 2
            self.KC.execute("import pickle\nglobals().update(pickle.loads({!r}))".format(sos_data),
                silent=True, store_history=False)
        elif self.kernel == 'ir':
            try:
                sos_data = '\n'.join(python2R(x) for x in items)
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Failed to get variable: {}'.format(e)})
                return
            self.KC.execute(sos_data, silent=True, store_history=False)
        else:
            self.send_response(self.iopub_socket, 'stream',
                 {'name': 'stderr', 'text': 'Can not pass variables to kernel {}'.format(self.kernel)})
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

    def handle_magic_put(self, options):
        items = options.split()
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
                    elif msg_type == 'execute_result':
                        #self.send_response(self.iopub_socket, 'stream',
                        #    {'name': 'stderr', 'text': repr(sub_msg['content']['data'])})
                        try:
                            env.sos_dict.update(
                                pickle.loads(eval(sub_msg['content']['data']['text/plain']))
                                )
                        except Exception as e:
                            self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stderr', 'text': 'Failed to put variable {}: {}'.format(', '.join(items), e)})
                        break
        elif self.kernel == 'ir':
            # if it is a python kernel, passing specified SoS variables to it
            self.KC.execute('{}\ndputToString(list({}))'.format(dputToString_func, 
                ','.join('{0}={0}'.format(x) for x in items)),
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
                    elif msg_type == 'execute_result':
                        #self.send_response(self.iopub_socket, 'stream',
                        #    {'name': 'stderr', 'text': repr(sub_msg['content']['data'])})
                        try:
                            env.sos_dict.update(
                                from_R_repr(sub_msg['content']['data']['text/plain'])
                                )
                        except Exception as e:
                             self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stderr', 'text': str(e)})
                        break
        else:
            raise UsageError('Can only pass variables to python kernel')

    def handle_shell_command(self, cmd):
        # interpolate command
        try:
            new_cmd = interpolate(cmd, sigil='${ }', local_dict=env.sos_dict._dict)
            if new_cmd != cmd:
                cmd = new_cmd
                if not cmd.startswith('cd ') and not cmd.startswith('cd\t'):
                    self.send_response(self.iopub_socket, 'stream',
                        {'name': 'stdout', 'text':
                        new_cmd.strip() + '\n## -- End interpolated command --\n'})
        except Exception as e:
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': 'Failed to interpolate {}: {}'.format(short_repr(cmd), e)})
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': str(e)})
        # command cd is handled differently because it is the only one that
        # has effect on sos.
        if cmd.startswith('cd ') or cmd.startswith('cd\t'):
            to_dir = cmd[3:].strip()
            try:
                os.chdir(os.path.expanduser(os.path.expandvars(to_dir)))
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': repr(e)})
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': os.getcwd()})
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
        last_input = [] if '__step_input__' not in env.sos_dict else copy.deepcopy(env.sos_dict['__step_input__'])
        last_output = [] if '__step_output__' not in env.sos_dict else copy.deepcopy(env.sos_dict['__step_output__'])
        if os.path.isfile('.sos/report.md'):
            os.remove('.sos/report.md')
        with self.redirect_sos_io():
            try:
                # record input and output
                res = self.executor.run(code, self.options)
                self.send_result(res, silent)
            except Exception:
                sys.stderr.flush()
                sys.stdout.flush()
                self.send_response(self.iopub_socket, 'display_data',
                    {
                        'source': 'SoS',
                        'metadata': {},
                        'data': { 'text/html': HTML('<hr color="black" width="60%">').data}
                    })
                raise
            except KeyboardInterrupt:
                return {'status': 'abort', 'execution_count': self.execution_count}
            finally:
                sys.stderr.flush()
                sys.stdout.flush()
        #
        if not silent:
            start_output = True
            # Send standard output
            if os.path.isfile('.sos/report.md'):
                with open('.sos/report.md') as sr:
                    sos_report = sr.read()
                with open(self.report_file, 'a') as summary_report:
                    summary_report.write(sos_report + '\n\n')
                if sos_report.strip():
                    self.send_response(self.iopub_socket, 'display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': {'text/markdown': sos_report}
                        })
                    start_output = False
            #
            if '__step_input__' in env.sos_dict and env.sos_dict['__step_input__'] != last_input:
                input_files = env.sos_dict['__step_input__']
            else:
                input_files = []
            if '__step_output__' in env.sos_dict and env.sos_dict['__step_output__'] != last_output:
                output_files = env.sos_dict['__step_output__']
                if output_files is None:
                    output_files = []
            else:
                output_files = []
            # use a table to list input and/or output file if exist
            if input_files or output_files:
                if not start_output:
                    self.send_response(self.iopub_socket, 'display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html': HTML('<hr color="black" width="60%">').data}
                        })
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
            # Send images, if any
            for filename in output_files:
                self.send_response(self.iopub_socket, 'stream',
                     {'name': 'stdout', 'text': '\n> ' + filename + ' ({})'.format(pretty_size(os.path.getsize(filename)))})
                previewer = [x for x in self.previewer.keys() if fnmatch.fnmatch(os.path.basename(filename), x)]
                if not previewer:
                    continue
                else:
                    # choose the longest matching pattern (e.g. '*' and '*.pdf', choose '*.pdf')
                    previewer_name = max(previewer, key=len)
                    previewer_func = self.previewer[previewer_name]
                    if not previewer_func:
                        continue
                    if not callable(previewer_func):
                        raise RuntimeError('Previewer {} is not callable'.format(previewer_name))
                    try:
                        result = previewer_func(filename)
                        if not result:
                            continue
                        if isinstance(result, str):
                            self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stdout', 'text': '\n'+result})
                        else:
                            msg_type, msg_data = result
                            self.send_response(self.iopub_socket, msg_type,
                                {'source': filename, 'data': msg_data, 'metadata': {}})
                    except Exception as e:
                        self.send_response(self.iopub_socket, 'stream',
                            {'name': 'stderr', 'text': 'Failed to preview {}: {}'.format(filename, e) })

    def send_result(self, res, silent=False):
        # this is Ok, send result back
        if not silent and res is not None:
            format_dict, md_dict = self.format_obj(res)
            self.send_response(self.iopub_socket, 'execute_result',
                {'execution_count': self.execution_count, 'data': format_dict,
                'metadata': md_dict})

    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        if self.original_keys is None:
            self._reset_dict()
        if code == 'import os\n_pid = os.getpid()':
            # this is a special probing command from vim-ipython. Let us handle it specially
            # so that vim-python can get the pid.
            return 
        if code.startswith('%dict'):
            # %dict should be the last magic
            options, remaining_code = self.get_magic_and_code(code, True)
            self.handle_magic_dict(options)
            return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
        elif code.startswith('%connect_info'):
            options, remaining_code = self.get_magic_and_code(code, False)
            cfile = find_connection_file()
            with open(cfile) as conn:
                conn_info = conn.read()
            self.send_response(self.iopub_socket, 'stream',
                  {'name': 'stdout', 'text': 'Connection file: {}\n{}'.format(cfile, conn_info)})
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%matplotlib'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.shell.enable_gui = lambda gui: None
            gui, backend = self.shell.enable_matplotlib(options)
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%set'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_set(options)
            # self.options will be set to inflence the execution of remaing_code
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%restart'):
            options, remaining_code = self.get_magic_and_code(code, True)
            if options == 'R':
                options = 'ir'
            self.restart_kernel(options)
            return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
        elif code.startswith('%with'):
            options, remaining_code = self.get_magic_and_code(code, False)
            if options == 'R':
                options = 'ir'
            original_kernel = self.kernel
            self.switch_kernel(options)
            try:
                return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.switch_kernel(original_kernel)
        elif code.startswith('%use'):
            options, remaining_code = self.get_magic_and_code(code, False)
            if options == 'R':
                options = 'ir'
            self.switch_kernel(options)
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%get'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_get(options)
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%put'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_put(options)
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%paste'):
            options, remaining_code = self.get_magic_and_code(code, True)
            self.options += ' ' + options.strip()
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
            return self.do_execute(code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('%run'):
            options, remaining_code = self.get_magic_and_code(code, True)
            old_options = self.options
            self.options += ' ' + options
            try:
                return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.options = old_options
        elif code.startswith('!'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_shell_command(code.split(' ')[0][1:] + ' ' + options)
            return self.do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.kernel != 'sos':
            # handle string interpolation before sending to the underlying kernel
            try:
                new_code = interpolate(code, sigil='${ }', local_dict=env.sos_dict._dict)
                if new_code != code:
                    code = new_code
                    self.send_response(self.iopub_socket, 'stream',
                        {'name': 'stdout', 'text':
                        new_code.strip() + '\n## -- End interpolated text --\n'})
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Failed to interpolate {}: {}'.format(short_repr(code), e)})
            return self.run_cell(code, store_history)
        else:
            # run sos
            try:
                self.run_sos_code(code, silent)
                self.shell.user_ns.update(env.sos_dict._dict)
                # trigger post processing of object and display matplotlib figures
                self.shell.events.trigger('post_execute')
                return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
            except Exception as e:
                stream_content = {'name': 'stderr', 'text': str(e)}
                self.send_response(self.iopub_socket, 'stream', stream_content)
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self.execution_count,
                   }

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
