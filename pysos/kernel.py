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
import re
import fnmatch
import contextlib
import subprocess
import tempfile
import argparse

import zipfile
import tarfile
import gzip
import pickle
from collections.abc import Sequence

from .utils import env, WorkflowDict, short_repr, pretty_size, dehtml, _parse_error
from ._version import __sos_version__, __version__
from .sos_eval import SoS_exec, SoS_eval, interpolate
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
                data = pandas.read_csv(filename)
                html = data._repr_html_()
                return 'display_data', { 'text/html': HTML(html).data}
            except Exception:
                pass
        elif filename.lower().endswith('.xlsx') or filename.lower().endswith('.xls'):
            try:
                import pandas
                data = pandas.read_excel(filename)
                html = data._repr_html_()
                return 'display_data', { 'text/html': HTML(html).data}
            except Exception:
                pass
        # is it a compressed file?
        elif zipfile.is_zipfile(filename):
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

# R    length (n)    Python
# NULL        None
# logical    1    boolean
# integer    1    integer
# numeric    1    double
# character    1    unicode
# logical    n > 1    array
# integer    n > 1    array
# numeric    n > 1    list
# character    n > 1    list
# list without names    n > 0    list
# list with names    n > 0    dict
# matrix    n > 0    array
# data.frame    n > 0    DataFrame

kernel_init_command = {
    'ir': r'''
..py.repr.logical.1 <- function(obj) {
    if(obj)
        'True'
    else
        'False'
}
..py.repr.integer.1 <- function(obj) {
    as.character(obj)
}
..py.repr.double.1 <- function(obj) {
    as.character(obj)
}
..py.repr.character.1 <- function(obj) {
    options(useFancyQuotes=FALSE)
    dQuote(obj)
}
..py.repr.dataframe <- function(obj) {
    if (!require("feather")) {
        install.packages('feather', repos='http://cran.stat.ucla.edu/')
        }
    library(feather)
    tf = tempfile('feather')
    write_feather(obj, tf)
    paste0("read_dataframe('", tf, "')")
}
..py.repr.matrix <- function(obj) {
    if (!require("feather")) {
        install.packages('feather', repos='http://cran.stat.ucla.edu/')
        }
    library(feather)
    tf = tempfile('feather')
    write_feather(as.data.frame(obj), tf)
    paste0("read_dataframe('", tf, "').as_matrix()")
}
..py.repr.n <- function(obj) {
    paste("[",
        paste(sapply(obj, ..py.repr), collapse=','),
        "]")
}
..py.repr <- function(obj) {
    if (is.matrix(obj)) {
        ..py.repr.matrix(obj)
    } else if (is.data.frame(obj)) {
        ..py.repr.dataframe(obj)
    } else if (is.null(obj)) {
        'None'
    } else if (is.integer(obj)) {
        if (length(obj) == 1)
            ..py.repr.integer.1(obj)
        else
            ..py.repr.n(obj)
    } else if (is.double(obj)){
        if (length(obj) == 1)
            ..py.repr.double.1(obj)
        else
            ..py.repr.n(obj)
    } else if (is.character(obj)) {
        if (length(obj) == 1)
            ..py.repr.character.1(obj)
        else
            ..py.repr.n(obj)
    } else if (is.logical(obj)) {
        if (length(obj) == 1)
            ..py.repr.logical.1(obj)
        else
            ..py.repr.n(obj)
    } else if (is.list(obj)) {
        # if the list has no name
        if (is.null(names(obj)))
            ..py.repr.n(obj)
        else {
            options(useFancyQuotes=FALSE)
            paste("{",
                  paste(sapply(names(obj), function (x)
                      paste(dQuote(gsub("\\.", "_", as.character(x))), ":", ..py.repr(obj[[x]]))),
                      collapse=','),
                  "}")
        }
    } else {
        "'Untransferrable variable'"
    }
}
'''}

def from_R_repr(expr):
    '''
    Convert expression returned from R to python
    '''
    try:
        if 'read_dataframe' in expr:
            from feather import read_dataframe
        # the result is something like
        # [1] "{'a': 1}"
        return eval(eval(expr.split(' ', 1)[-1]))
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
            return '''if (!require("feather")) {{
                install.packages('feather', repos='http://cran.stat.ucla.edu/')
                }}
                library(feather)
                {} <- {}'''.format('.' + name[1:] if name.startswith('_') else name, r_repr)
        else:
            return '{} <- {}'.format('.' + name[1:] if name.startswith('_') else name, r_repr)
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

    MAGIC_DICT = re.compile('^%dict(\s|$)')
    MAGIC_CONNECT_INFO = re.compile('^%connect_info(\s|$)')
    MAGIC_MATPLOTLIB = re.compile('^%matplotlib(\s|$)')
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
        SoS_exec('import os, sys, glob', None)
        SoS_exec('from pysos.runtime import *', None)
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
        #
        return {'status': 'incomplete', 'indent': ''}

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
                        self.kernels[kernel] = manager.start_new_kernel(startup_timeout=60, kernel_name=kernel)
                        self.KM, self.KC = self.kernels[kernel]
                        self.RET_VARS = ret_vars
                        self.kernel = kernel
                        if self.kernel in kernel_init_command:
                            self.run_cell(kernel_init_command[self.kernel], False)
                    except Exception as e:
                        self.send_response(self.iopub_socket, 'stream',
                            {'name': 'stderr', 'text': 'Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}\n'.format(kernel, e)})
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
            self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Cannot restart sos kernel from within sos.'})
        elif kernel:
            if kernel in self.kernels:
                try:
                    self.kernels[kernel][0].shutdown_kernel(restart=False)
                except Exception as e:
                    env.logger.warning('Failed to shutdown kernel {}: {}\n'.format(kernel, e))
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
                self.send_response(self.iopub_socket, 'stream',
                      {'name': 'stderr', 'text': 'Unrecognized sosdict option or variable name {}'.format(x)})
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
            for item in items:
                if item.startswith('_'):
                    self.send_response(self.iopub_socket, 'stream',
                        {'name': 'stderr', 'text': 'Variable {} is imported as {}\n'.format(item, '.' + item[1:])})
            try:
                sos_data = '\n'.join(python2R(x) for x in items)
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Failed to get variable: {}\n'.format(e)})
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
                            #self.send_response(self.iopub_socket, 'stream',
                            #    {'name': 'stderr', 'text': repr(sub_msg['content']['data'])})
                            try:
                                env.sos_dict.update(
                                    pickle.loads(eval(sub_msg['content']['data']['text/plain']))
                                    )
                            except Exception as e:
                                self.send_response(self.iopub_socket, 'stream',
                                    {'name': 'stderr', 'text': 'Failed to put variable {}: {}\n'.format(', '.join(items), e)})
                            break
                        else:
                            self.send_response(self.iopub_socket, msg_type,
                                sub_msg['content'])
            # verify
            for item in items:
                if not item in env.sos_dict:
                    self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stderr', 'text': 'Failed to put variable {} to SoS namespace\n'.format(item)})
        elif self.kernel == 'ir':
            for item in items:
                if '.' in item:
                    self.send_response(self.iopub_socket, 'stream',
                        {'name': 'stderr', 'text': 'Variable {} is exported as {}\n'.format(item, item.replace('.', '_'))})
            # if it is a python kernel, passing specified SoS variables to it
            self.KC.execute('..py.repr(list({}))'.format(
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
                    else:
                        # irkernel (since the new version) does not produce execute_result, only
                        # display_data
                        if msg_type in ('display_data', 'execute_result'):
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
                            self.send_response(self.iopub_socket, msg_type,
                                sub_msg['content'])
            # verify
            for item in items:
                if not item.replace('.', '_') in env.sos_dict:
                    self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stderr', 'text': 'Failed to put variable {} to SoS namespace\n'.format(item)})
        else:
            self.send_response(self.iopub_socket, 'stream',
                                {'name': 'stderr', 'text': 'Can only pass variables to python kernel'})

    def handle_magic_preview(self, options):
        try:
            options = interpolate(options, sigil='${ }', local_dict=env.sos_dict._dict)
        except Exception as e:
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': 'Failed to interpolate {}: {}\n'.format(short_repr(options), e)})
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': str(e)})
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
              'data': { 'text/html': HTML('<pre>## %preview {}</pre>'.format(options)).data}
            })
        # expand items
        for item in items:
            try:
                if os.path.isfile(item):
                    self.preview(item)
                    continue
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': item + ':\n'})
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
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': '\n> Failed to find file or evaluate expression {}: {}'.format(item, e)})

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
                {'name': 'stdout', 'text': 'Failed to interpolate {}: {}\n'.format(short_repr(cmd), e)})
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': str(e)})
            return
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
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Keyboard Interrupt\n'})
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
                            'data': { 'text/html': HTML('<pre>## -- Preview output --</pre>').data}
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
            self.send_response(self.iopub_socket, 'stream',
                 {'name': 'stderr', 'text': '\n> ' + filename + ' does not exist'})
            return
        self.send_response(self.iopub_socket, 'stream',
             {'name': 'stdout', 'text': '\n> ' + filename + ' ({})'.format(pretty_size(os.path.getsize(filename)))})
        previewer = [x for x in self.previewer.keys() if fnmatch.fnmatch(os.path.basename(filename), x)]
        if not previewer:
            return
        # choose the longest matching pattern (e.g. '*' and '*.pdf', choose '*.pdf')
        previewer_name = max(previewer, key=len)
        previewer_func = self.previewer[previewer_name]
        if not previewer_func:
            return
        if not callable(previewer_func):
            raise RuntimeError('Previewer {} is not callable'.format(previewer_name))
        try:
            result = previewer_func(filename)
            if not result:
                return
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

    def _do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        if code.startswith('\n') or code.startswith(' '):
            code = re.sub('^\s*\n', '', code, re.M)
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
            if options == 'R':
                options = 'ir'
            self.restart_kernel(options)
            return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self.execution_count}
        elif self.MAGIC_WITH.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                args = self.parse_in_out_vars(options.split())
            except Exception as e:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Invalid option "{}": {}\n'.format(options, e)})
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self.execution_count,
                   }
            if args.kernel == 'R':
                args.kernel = 'ir'
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
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Invalid option "{}": {}\n'.format(options, e)})
                return {'status': 'abort',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self.execution_count,
                   }
            if args.kernel == 'R':
                args.kernel = 'ir'
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
        elif code.startswith('!'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_shell_command(code.split(' ')[0][1:] + ' ' + options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
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
            try:
                return self.run_cell(code, store_history)
            except KeyboardInterrupt:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stderr', 'text': 'Keyboard Interrupt\n'})
                return {'status': 'abort', 'execution_count': self.execution_count}
        else:
            # run sos
            try:
                self.run_sos_code(code, silent)
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

    #
    # The following is copied from spyder_kernel.py, which are needed for
    # spyder to work. I cannot however derive this class also from
    # SpyderKernel because that kernel assumes a ipkernel kernel.
    #
    @property
    def _pdb_frame(self):
        """Return current Pdb frame if there is any"""
        if self._pdb_obj is not None and self._pdb_obj.curframe is not None:
            return self._pdb_obj.curframe

    @property
    def _pdb_locals(self):
        """
        Return current Pdb frame locals if available. Otherwise
        return an empty dictionary
        """
        if self._pdb_frame:
            return self._pdb_obj.curframe_locals
        else:
            return {}

    # -- Public API ---------------------------------------------------
    # For the Variable Explorer
    def get_namespace_view(self):
        """
        Return the namespace view

        This is a dictionary with the following structure

        {'a': {'color': '#800000', 'size': 1, 'type': 'str', 'view': '1'}}

        Here:
        * 'a' is the variable name
        * 'color' is the color used to show it
        * 'size' and 'type' are self-evident
        * and'view' is its value or the text shown in the last column
        """
        settings = self.namespace_view_settings
        if settings:
            ns = self._get_current_namespace()
            more_excluded_names = ['In', 'Out']
            view = make_remote_view(ns, settings, more_excluded_names)
            return view

    def get_var_properties(self):
        """
        Get some properties of the variables in the current
        namespace
        """
        settings = self.namespace_view_settings
        if settings:
            ns = self._get_current_namespace()
            data = get_remote_data(ns, settings, mode='editable',
                                   more_excluded_names=['In', 'Out'])

            properties = {}
            for name, value in list(data.items()):
                properties[name] = {
                    'is_list':  isinstance(value, (tuple, list)),
                    'is_dict':  isinstance(value, dict),
                    'len': self._get_len(value),
                    'is_array': self._is_array(value),
                    'is_image': self._is_image(value),
                    'is_data_frame': self._is_data_frame(value),
                    'is_series': self._is_series(value),
                    'array_shape': self._get_array_shape(value),
                    'array_ndim': self._get_array_ndim(value)
                }

            return properties
        else:
            return {}

    def get_value(self, name):
        """Get the value of a variable"""
        ns = self._get_current_namespace()
        value = ns[name]
        publish_data({'__spy_data__': value})

    def set_value(self, name, value):
        """Set the value of a variable"""
        ns = self._get_reference_namespace(name)
        value = deserialize_object(value)[0]
        if isinstance(value, CannedObject):
            value = value.get_object()
        ns[name] = value

    def remove_value(self, name):
        """Remove a variable"""
        ns = self._get_reference_namespace(name)
        ns.pop(name)

    def copy_value(self, orig_name, new_name):
        """Copy a variable"""
        ns = self._get_reference_namespace(orig_name)
        ns[new_name] = ns[orig_name]

    def load_data(self, filename, ext):
        """Load data from filename"""
        glbs = self._mglobals()

        load_func = iofunctions.load_funcs[ext]
        data, error_message = load_func(filename)

        if error_message:
            return error_message

        for key in list(data.keys()):
            new_key = fix_reference_name(key, blacklist=list(glbs.keys()))
            if new_key != key:
                data[new_key] = data.pop(key)

        try:
            glbs.update(data)
        except Exception as error:
            return str(error)

        return None

    def save_namespace(self, filename):
        """Save namespace into filename"""
        ns = self._get_current_namespace()
        settings = self.namespace_view_settings
        data = get_remote_data(ns, settings, mode='picklable',
                               more_excluded_names=['In', 'Out']).copy()
        return iofunctions.save(data, filename)

    # --- For Pdb
    def get_pdb_step(self):
        """Return info about pdb current frame"""
        return self._pdb_step

    # --- For the Help plugin
    def is_defined(self, obj, force_import=False):
        """Return True if object is defined in current namespace"""
        ns = self._get_current_namespace(with_magics=True)
        return isdefined(obj, force_import=force_import, namespace=ns)

    def get_doc(self, objtxt):
        """Get object documentation dictionary"""
        obj, valid = self._eval(objtxt)
        if valid:
            return getdoc(obj)

    def get_source(self, objtxt):
        """Get object source"""
        obj, valid = self._eval(objtxt)
        if valid:
            return getsource(obj)

    # -- Private API ---------------------------------------------------
    # --- For the Variable Explorer
    def _get_current_namespace(self, with_magics=False):
        """
        Return current namespace

        This is globals() if not debugging, or a dictionary containing
        both locals() and globals() for current frame when debugging
        """
        ns = {}
        glbs = self._mglobals()

        if self._pdb_frame is None:
            ns.update(glbs)
        else:
            ns.update(glbs)
            ns.update(self._pdb_locals)

        # Add magics to ns so we can show help about them on the Help
        # plugin
        if with_magics:
            line_magics = self.shell.magics_manager.magics['line']
            cell_magics = self.shell.magics_manager.magics['cell']
            ns.update(line_magics)
            ns.update(cell_magics)

        return ns

    def _get_reference_namespace(self, name):
        """
        Return namespace where reference name is defined

        It returns the globals() if reference has not yet been defined
        """
        glbs = self._mglobals()
        if self._pdb_frame is None:
            return glbs
        else:
            lcls = self._pdb_locals
            if name in lcls:
                return lcls
            else:
                return glbs

    def _mglobals(self):
        """Return current globals -- handles Pdb frames"""
        if self._pdb_frame is not None:
            return self._pdb_frame.f_globals
        else:
            return self.shell.user_ns

    def _get_len(self, var):
        """Return sequence length"""
        try:
            return len(var)
        except TypeError:
            return None

    def _is_array(self, var):
        """Return True if variable is a NumPy array"""
        try:
            import numpy
            return isinstance(var, numpy.ndarray)
        except ImportError:
            return False

    def _is_image(self, var):
        """Return True if variable is a PIL.Image image"""
        try:
            from PIL import Image
            return isinstance(var, Image.Image)
        except ImportError:
            return False

    def _is_data_frame(self, var):
        """Return True if variable is a DataFrame"""
        try:
            from pandas import DataFrame
            return isinstance(var, DataFrame)
        except:
            return False

    def _is_series(self, var):
        """Return True if variable is a Series"""
        try:
            from pandas import Series
            return isinstance(var, Series)
        except:
            return False

    def _get_array_shape(self, var):
        """Return array's shape"""
        try:
            if self._is_array(var):
                return var.shape
            else:
                return None
        except AttributeError:
            return None

    def _get_array_ndim(self, var):
        """Return array's ndim"""
        try:
            if self._is_array(var):
                return var.ndim
            else:
                return None
        except AttributeError:
            return None

    # --- For Pdb
    def _register_pdb_session(self, pdb_obj):
        """Register Pdb session to use it later"""
        self._pdb_obj = pdb_obj

    # --- For the Help plugin
    def _eval(self, text):
        """
        Evaluate text and return (obj, valid)
        where *obj* is the object represented by *text*
        and *valid* is True if object evaluation did not raise any exception
        """
        assert is_text_string(text)
        ns = self._get_current_namespace(with_magics=True)
        try:
            return eval(text, ns), True
        except:
            return None, False

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
