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
import tempfile
import argparse
import pkg_resources

import pickle
from collections.abc import Sequence

from sos.utils import env, WorkflowDict, short_repr, pretty_size, _parse_error
from sos._version import __sos_version__, __version__
from sos.sos_eval import SoS_exec, SoS_eval, interpolate
from sos.sos_syntax import SOS_SECTION_HEADER
from sos.converter import SoS_Exporter

from IPython.core.interactiveshell import InteractiveShell
from IPython.lib.clipboard import ClipboardEmpty, osx_clipboard_get, tkinter_clipboard_get
from IPython.core.error import UsageError
from IPython.core.display import HTML
from ipykernel.kernelbase import Kernel
from jupyter_client import manager, find_connection_file
from ipykernel.zmqshell import ZMQDisplayPublisher

from textwrap import dedent
from io import StringIO

from .sos_executor import Interactive_Executor

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

def get_supported_languages():
    group = 'sos_languages'
    result = {}
    for entrypoint in pkg_resources.iter_entry_points(group=group):
        # Grab the function that is the actual plugin.
        name = entrypoint.name
        try:
            plugin = entrypoint.load()()
            result[name] = plugin
            # for convenience, we create two entries for, e.g. R and ir
            if name != plugin.kernel_name:
                result[plugin.kernel_name] = plugin
        except Exception as e:
            raise RuntimeError('Failed to load language {}: {}'.format(entrypoint.name, e))
    return result

def get_previewers():
    # Note: data is zest.releaser specific: we want to pass
    # something to the plugin
    group = 'sos_previewers'
    result = []
    for entrypoint in pkg_resources.iter_entry_points(group=group):
        # Grab the function that is the actual plugin.
        plugin = entrypoint.load()
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
                result.append((getattr(imported, func), plugin, priority))
            except ImportError:
                env.logger.warning('Failed to load function {}:{}'.format(mod, func))
        else:
            result.append((name, plugin, priority))
    #
    result.sort(key=lambda x: -x[2])
    with open('a.txt', 'w') as ss:
        ss.write(repr(result))
    return result

def homogeneous_type(seq):
    iseq = iter(seq)
    first_type = type(next(iseq))
    if first_type in (int, float):
        return True if all( (type(x) in (int, float)) for x in iseq ) else False
    else:
        return True if all( (type(x) is first_type) for x in iseq ) else False

def R_repr(obj):
    if isinstance(obj, bool):
        return 'TRUE' if obj else 'FALSE'
    elif isinstance(obj, (int, float, str)):
        return repr(obj)
    elif isinstance(obj, Sequence):
        if len(obj) == 0:
            return 'c()'

        # if the data is of homogeneous type, let us use c()
        # otherwise use list()
        # this can be confusion but list can be difficult to handle
        if homogeneous_type(obj):
            return 'c(' + ','.join(R_repr(x) for x in obj) + ')'
        else:
            return 'list(' + ','.join(R_repr(x) for x in obj) + ')'
    elif obj is None:
        return 'NULL'
    elif isinstance(obj, dict):
        return 'list(' + ','.join('{}={}'.format(x, R_repr(y)) for x,y in obj.items()) + ')'
    elif isinstance(obj, set):
        return 'list(' + ','.join(R_repr(x) for x in obj) + ')'
    else:
        import numpy
        import pandas
        if isinstance(obj, (numpy.intc, numpy.intp, numpy.int8, numpy.int16, numpy.int32, numpy.int64,\
                numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.float16, numpy.float32, \
                numpy.float64)):
            return repr(obj)
        elif isinstance(obj, numpy.matrixlib.defmatrix.matrix):
            try:
                import feather
            except ImportError:
                raise UsageError('The feather-format module is required to pass numpy matrix as R matrix'
                    'See https://github.com/wesm/feather/tree/master/python for details.')
            feather_tmp_ = tempfile.NamedTemporaryFile(suffix='.feather', delete=False).name
            feather.write_dataframe(pandas.DataFrame(obj).copy(), feather_tmp_)
            return 'data.matrix(read_feather("{}"))'.format(feather_tmp_)
        elif isinstance(obj, numpy.ndarray):
            return 'c(' + ','.join(R_repr(x) for x in obj) + ')'
        elif isinstance(obj, pandas.DataFrame):
            try:
                import feather
            except ImportError:
                raise UsageError('The feather-format module is required to pass pandas DataFrame as R data.frame'
                    'See https://github.com/wesm/feather/tree/master/python for details.')
            feather_tmp_ = tempfile.NamedTemporaryFile(suffix='.feather', delete=False).name
            try:
                data = obj.copy()
                feather.write_dataframe(data, feather_tmp_)
            except:
                # if data cannot be written, we try to manipulate data
                # frame to have consistent types and try again
                for c in data.columns:
                    if not homogeneous_type(data[c]):
                        data[c] = [str(x) for x in data[c]]
                feather.write_dataframe(data, feather_tmp_)
            return 'read_feather("{}")'.format(feather_tmp_)
        else:
            return repr('Unsupported datatype {}'.format(short_repr(obj)))


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
        'nbconvert_exporter': 'sos.kernel.SoS_Exporter',
    }
    banner = "SoS kernel - script of scripts"
    shell = InteractiveShell.instance()

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
    MAGIC_PREVIEW = re.compile('^%preview(\s|$)')

    def __init__(self, **kwargs):
        super(SoS_Kernel, self).__init__(**kwargs)
        self.options = ''
        self.kernel = 'sos'
        self.banner = self.banner + '\nConnection file {}'.format(os.path.basename(find_connection_file()))
        # FIXME: this should in theory be a MultiKernelManager...
        self.kernels = {}
        self.supported_languages = None
        self.format_obj = self.shell.display_formatter.format

        # InteractiveShell uses a default publisher that only displays text/plain
        # using the ZMQDisplayPublisher will display matplotlib inline figures
        self.shell.display_pub = ZMQDisplayPublisher()
        self.shell.display_pub.session = self.session
        self.shell.display_pub.pub_socket = self.iopub_socket

        self.shell.enable_gui = lambda x: False
        self.previewers = None

        self.report_file = os.path.join(env.exec_dir, 'summary_report.md')
        if os.path.isfile(self.report_file):
            os.remove(self.report_file)
        # touch the file
        with open(self.report_file, 'w'):
            pass
        self.original_keys = None

    def _reset_dict(self):
        env.sos_dict = WorkflowDict()
        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        SoS_exec("run_mode = 'interactive'", None)
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

    def warn(self, message):
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
                        sub_msg['content']['execution_count'] = self.execution_count
                    #
                    # NOTE: we do not send status of sub kernel alone because
                    # these are generated automatically during the execution of
                    # "this cell" in SoS kernel
                    #
                    self.send_response(self.iopub_socket, msg_type, sub_msg['content'])
        # now get the real result
        reply = self.KC.get_shell_msg(timeout=10)
        reply['content']['execution_count'] = self.execution_count
        return reply['content']

    def switch_kernel(self, kernel, ret_vars=[]):
        if self.supported_languages is None:
            self.supported_languages = get_supported_languages()
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

    def handle_magic_dict(self, line):
        'Magic that displays content of the dictionary'
        # do not return __builtins__ beacuse it is too long...
        actions = line.strip().split()
        keys = [x for x in actions if not x.startswith('-')]
        for x in keys:
            if not x in env.sos_dict:
                self.warn('Unrecognized sosdict option or variable name {}'.format(x))
                return
        for x in [x for x in actions if x.startswith('-')]:
            if not x in ['-r', '--reset', '-k', '--keys', '-a', '--all']:
                raise RuntimeError('Unrecognized option {} for magic %dict'.format(x))
        if '--reset' in actions or '-r' in actions:
            self.send_result(self._reset_dict())
        if '--keys' in actions or '-k' in actions:
            if '--all' in actions or '-a' in actions:
                self.send_result(env.sos_dict._dict.keys())
            elif keys:
                self.send_result(set(keys))
            else:
                self.send_result({x for x in env.sos_dict._dict.keys() if not x.startswith('__')} - self.original_keys)
        else:
            if '--all' in actions or '-a' in actions:
                self.send_result(env.sos_dict._dict)
            elif keys:
                self.send_result({x:y for x,y in env.sos_dict._dict.items() if x in keys})
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
                    new_name, py_repr = lan.repr_of_py_obj(item, env.sos_dict[item])
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
            new_names, py_repr = lan.py_repr_of_obj(items)
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
                                    lan.py_from_repr_of_obj(sub_msg['content']['data']['text/plain'])
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
        options = self._interpolate_option(options, quiet=True)
        if options is None:
            return
        # find filenames and quoted expressions
        import shlex
        items = shlex.split(options)
        if not items:
            return
        self.send_response(self.iopub_socket, 'display_data',
            {
              'source': 'SoS',
              'metadata': {},
              'data': { 'text/html': HTML('<pre><font color="green">## %preview {}</font></pre>'.format(options)).data}
            })
        # expand items
        for item in items:
            try:
                if os.path.isfile(item):
                    self.preview(item)
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
                if item in env.sos_dict:
                    obj = env.sos_dict[item]
                else:
                    # we disallow other sigil
                    obj = SoS_eval(item, sigil='${ }')
                format_dict, md_dict = self.format_obj(obj)
                self.send_response(self.iopub_socket, 'display_data',
                    {'execution_count': self.execution_count, 'data': format_dict,
                    'metadata': md_dict})
            except Exception as e:
                self.warn('\n> Failed to preview expression {}: {}'.format(item, e))

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
                self.warn('Keyboard Interrupt\n')
                return {'status': 'abort', 'execution_count': self.execution_count}
            finally:
                sys.stderr.flush()
                sys.stdout.flush()
        #
        if not silent:
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
            #
            if 'input' in env.sos_dict:
                input_files = env.sos_dict['input']
            else:
                input_files = []
            if 'output' in env.sos_dict:
                output_files = env.sos_dict['output']
                if output_files is None:
                    output_files = []
            else:
                output_files = []
            # use a table to list input and/or output file if exist
            if input_files or output_files:
                self.send_response(self.iopub_socket, 'display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html': HTML('<pre><font color="green">## -- Preview output --</font></pre>').data}
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
                self.preview(filename)

    def preview(self, filename):
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
        for x,y,_ in self.previewers:
            if isinstance(x, str):
                if fnmatch.fnmatch(os.path.basename(filename), x):
                    previewer_func = y
                    break
            else:
                # it should be a function
                try:
                    if x(filename):
                        previewer_func = y
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
                {'execution_count': self.execution_count, 'data': format_dict,
                'metadata': md_dict})

    def parse_in_out_vars(self, args):
        parser = argparse.ArgumentParser()
        parser.add_argument('kernel', nargs='?', default='')
        parser.add_argument('-i', '--in', nargs='*', dest='in_vars')
        parser.add_argument('-o', '--out', nargs='*', dest='out_vars')
        parser.error = _parse_error
        return parser.parse_args(args)

    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        ret = self._do_execute(code=code, silent=silent, store_history=store_history,
            user_expressions=user_expressions, allow_stdin=allow_stdin)
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
        except Exception as e:
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
            options, remaining_code = self.get_magic_and_code(code, True)
            self.handle_magic_dict(options)
            return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
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
            self.shell.enable_gui = lambda gui: None
            gui, backend = self.shell.enable_matplotlib(options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_SET.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_set(options)
            # self.options will be set to inflence the execution of remaing_code
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_RESTART.match(code):
            options, remaining_code = self.get_magic_and_code(code, True)
            self.restart_kernel(options)
            return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
        elif self.MAGIC_WITH.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                args = self.parse_in_out_vars(options.split())
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self.execution_count,
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
                args = self.parse_in_out_vars(options.split())
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'abort',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self.execution_count,
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
            return self._do_execute(code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_RUN.match(code):
            options, remaining_code = self.get_magic_and_code(code, True)
            old_options = self.options
            self.options += ' ' + options
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.options = old_options
        elif self.MAGIC_PREVIEW.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
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
            code = self._interpolate_option(code, quiet=False)
            if code is None:
                return
            try:
                return self.run_cell(code, store_history)
            except KeyboardInterrupt:
                self.warn('Keyboard Interrupt\n')
                return {'status': 'abort', 'execution_count': self.execution_count}
        else:
            # run sos
            try:
                self.run_sos_code(code, silent)
                return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
            except Exception as e:
                self.warn(str(e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self.execution_count,
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



#
# The following is a hack to let SoS kernel work with spyder version > 3.0
# A complete rewritten will be needed after spyder officially supports
# third-party kernel.
#
# Note that this kernel is only used by Spyder, not by jupyter notebook
# and qtconsole.
#
from ipykernel.comm import CommManager

class SoS_SpyderKernel(SoS_Kernel):
    """Spyder kernel for Jupyter"""

    MAGIC_EDIT = re.compile('^%edit(\s|$)')

    def __init__(self, *args, **kwargs):
        super(SoS_SpyderKernel, self).__init__(*args, **kwargs)

        self.namespace_view_settings = {}
        self._pdb_obj = None
        self._pdb_step = None
        # ???
        self.shell.kernel = self

        self.comm_manager = CommManager(parent=self, kernel=self)

        self.shell.configurables.append(self.comm_manager)
        comm_msg_types = [ 'comm_open', 'comm_msg', 'comm_close' ]
        for msg_type in comm_msg_types:
            self.shell_handlers[msg_type] = getattr(self.comm_manager, msg_type)

    def parse_edit_magic(self, args):
        import shlex
        parser = argparse.ArgumentParser()
        parser.add_argument('filenames', nargs='+')
        parser.add_argument('--cd', action='store_true', dest='__switch_dir__')
        parser.error = _parse_error
        return parser.parse_args(shlex.split(args))

    def handle_magic_edit(self, options):
        import subprocess
        options = self._interpolate_option(options)
        if options is None:
            return
        args = self.parse_edit_magic(options)
        import1 = "import sys"
        import2 = "from spyder.app.start import send_args_to_spyder"
        code = "send_args_to_spyder([{}])".format(','.join('"{}"'.format(x) for x in args.filenames))
        cmd = "{0} -c '{1}; {2}; {3}'".format(sys.executable,
            import1, import2, code)
        subprocess.call(cmd, shell=True)
        if args.__switch_dir__:
            script_dir = os.path.dirname(os.path.abspath(args.filenames[-1]))
            os.chdir(script_dir)
            self.send_response(self.iopub_socket, 'stream',
                  {'name': 'stdout', 'text': 'Current working directory is set to {}\n'.format(script_dir)})

    # add an additional magic that only useful for spyder
    def _do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        code = self.remove_leading_comments(code)

        if self.MAGIC_EDIT.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_edit(options)
            # self.options will be set to inflence the execution of remaing_code
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        else:
            return super(SoS_SpyderKernel, self)._do_execute(code, silent, store_history, user_expressions, allow_stdin)

    def _reset_dict(self):
        super(SoS_SpyderKernel, self)._reset_dict()
        # spyder 3 executes commands such as
        #
        # get_ipython().kernel....
        #
        # so we will have to expose the shell to the SoS dictionary
        env.sos_dict.set('__ipython__', self.shell)
        SoS_exec('''
def get_ipython():
    return __ipython__
''', None)
        self.original_keys.add('__ipython__')
        self.original_keys.add('get_ipython')

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
