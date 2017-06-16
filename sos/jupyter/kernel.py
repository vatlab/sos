#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
import time
import base64
import shlex
import fnmatch
import contextlib
import subprocess
import argparse
import pkg_resources
import pydoc

from ipykernel.ipkernel import IPythonKernel
from collections import Sized, defaultdict, OrderedDict

from types import ModuleType
from sos.utils import env, WorkflowDict, short_repr, pretty_size, PrettyRelativeTime
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
from .sos_step import PendingTasks

# I have not figured out how to log message to console so I have to use a file
# to send debug message.
def log_to_file(msg):
    with open(os.path.join(os.path.expanduser('~'), 'jupyter_debug.txt'), 'a') as log:
        log.write('{}\n'.format(msg))

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
            if self.name == 'stdout' and self.kernel._render_result:
                format_dict, md_dict = self.kernel.format_obj(self.kernel.render_result(content))
                self.kernel.send_response(self.kernel.iopub_socket, 'display_data',
                    {'source': 'SoS', 'metadata': md_dict,
                     'data': format_dict
                    })
            else:
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
    # we put string infront of functions so that filename matching will be used before
    # content matching. For example, the xlsx file is actually a zip file so it could be
    # previewed as a zip file first.
    return [x for x in result if isinstance(x[0], str)] + [x for x in result if not isinstance(x[0], str)]

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
        'cd',
        'debug',
        'dict',
        'get',
        'matplotlib',
        'paste',
        'preview',
        'put',
        'render',
        'rerun',
        'run',
        'sandbox',
        'set',
        'sessioninfo',
        'sosrun',
        'sossave',
        'shutdown',
        'taskinfo',
        'tasks',
        'toc',
        'use',
        'with',
    }
    MAGIC_DICT = re.compile('^%dict(\s|$)')
    MAGIC_CONNECT_INFO = re.compile('^%connect_info(\s|$)')
    MAGIC_MATPLOTLIB = re.compile('^%matplotlib(\s|$)')
    MAGIC_CD = re.compile('^%cd(\s|$)')
    MAGIC_SET = re.compile('^%set(\s|$)')
    MAGIC_SHUTDOWN = re.compile('^%shutdown(\s|$)')
    MAGIC_WITH = re.compile('^%with(\s|$)')
    MAGIC_FRONTEND = re.compile('^%frontend(\s|$)')
    MAGIC_USE = re.compile('^%use(\s|$)')
    MAGIC_GET = re.compile('^%get(\s|$)')
    MAGIC_PUT = re.compile('^%put(\s|$)')
    MAGIC_PASTE = re.compile('^%paste(\s|$)')
    MAGIC_RUN = re.compile('^%run(\s|$)')
    MAGIC_SOSRUN = re.compile('^%sosrun(\s|$)')
    MAGIC_SOSSAVE = re.compile('^%sossave(\s|$)')
    MAGIC_SAVE = re.compile('^%save(\s|$)')
    MAGIC_RERUN = re.compile('^%rerun(\s|$)')
    MAGIC_PREVIEW = re.compile('^%preview(\s|$)')
    MAGIC_SANDBOX = re.compile('^%sandbox(\s|$)')
    MAGIC_DEBUG = re.compile('^%debug(\s|$)')
    MAGIC_TASKINFO = re.compile('^%taskinfo(\s|$)')
    MAGIC_TASKS = re.compile('^%tasks(\s|$)')
    MAGIC_SKIP = re.compile('^%skip(\s|$)')
    MAGIC_TOC = re.compile('^%toc(\s|$)')
    MAGIC_RENDER = re.compile('^%render(\s|$)')
    MAGIC_SESSIONINFO = re.compile('^%sessioninfo(\s|$)')

    def get_use_parser(self):
        parser = argparse.ArgumentParser(prog='%use',
            description='''Switch to a specified subkernel.''')
        parser.add_argument('name', nargs='?', default='',
            help='''Displayed name of kernel to start (if no kernel with name is
            specified) or switch to (if a kernel with this name is already started).
            The name is usually a kernel name (e.g. %%use ir) or a language name
            (e.g. %%use R) in which case the language name will be used. One or
            more parameters --language or --kernel will need to be specified
            if a new name is used to start a separate instance of a kernel.''')
        parser.add_argument('-k', '--kernel',
            help='''kernel name as displayed in the output of jupyter kernelspec
            list. Default to the default kernel of selected language (e.g. ir for
            language R.''')
        parser.add_argument('-l', '--language',
            help='''Language extension that enables magics such as %%get and %%put
            for the kernel, which should be in the name of a registered language
            (e.g. R), or a specific language module in the format of
            package.module:class. SoS maitains a list of languages and kernels
            so this option is only needed for starting a new instance of a kernel.
            ''')
        parser.add_argument('-c', '--color',
            help='''Background color of new or existing kernel, which overrides
            the default color of the language. A special value "default" can be
            used to reset color to default.''')
        parser.add_argument('-r', '--restart', action='store_true',
            help='''Restart the kernel if it is running.''')
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
        parser.add_argument('name', nargs='?', default='',
            help='''Displayed name of kernel to start (if no kernel with name is
            specified) or switch to (if a kernel with this name is already started).
            The name is usually a kernel name (e.g. %use ir) or a language name
            (e.g. %use R) in which case the language name will be used. One or
            more parameters --language or --kernel will need to be specified
            if a new name is used to start a separate instance of a kernel.''')
        parser.add_argument('-k', '--kernel',
            help='''kernel name as displayed in the output of jupyter kernelspec
            list. Default to the default kernel of selected language (e.g. ir for
            language R.''')
        parser.add_argument('-l', '--language',
            help='''Language extension that enables magics such as %get and %put
            for the kernel, which should be in the name of a registered language
            (e.g. R), or a specific language module in the format of
            package.module:class. SoS maitains a list of languages and kernels
            so this option is only needed for starting a new instance of a kernel.
            ''')
        parser.add_argument('-c', '--color',
            help='''Background color of existing or new kernel, which overrides
            the default color of the language. A special value "default" can be
            used to reset color to default.''')
        parser.add_argument('-r', '--restart', action='store_true',
            help='''Restart the kernel if it is running.''')
        parser.add_argument('-i', '--in', nargs='*', dest='in_vars',
            help='Input variables (variables to get from SoS kernel)')
        parser.add_argument('-o', '--out', nargs='*', dest='out_vars',
            help='''Output variables (variables to put back to SoS kernel
            before switching back to the SoS kernel''')
        parser.error = self._parse_error
        return parser

    def get_frontend_parser(self):
        parser = argparse.ArgumentParser(prog='%frontend',
            description='''Use specified the subkernel to evaluate current
            cell. soft with tells the start kernel of the cell. If no other
            switch happens the kernel will switch back. However, a %use inside
            the cell will still switch the global kernel. In contrast, a hard
            %with magic will absorb the effect of %use.''')
        parser.add_argument('--list-kernel', action='store_true',
            help='List kernels')
        parser.add_argument('--default-kernel',
            help='Default global kernel')
        parser.add_argument('--cell-kernel',
            help='Kernel to switch to.')
        parser.add_argument('--use-panel', action='store_true',
            help='If panel is open')
        # pass cell index from notebook so that we know which cell fired
        # the command. Use to set metadata of cell through frontend message
        parser.add_argument('--resume', action='store_true',
            help='''If the cell is automatically reresumed by frontend, in which
            case -s force should be handled differently.'''),
        parser.add_argument('--cell', dest='cell_idx', type=int,
            help='Index of cell')
        parser.add_argument('--workflow', const='', nargs='?',
            help='Workflow defined in the notebook')
        parser.add_argument('--filename',
            help='filename of the current notebook')
        parser.error = self._parse_error
        return parser

    def get_preview_parser(self):
        parser = argparse.ArgumentParser(prog='%preview',
            description='''Preview files, sos variables, or expressions in the
                side panel, or notebook if side panel is not opened, unless
                options --panel or --notebook is specified.''')
        parser.add_argument('items', nargs='*',
            help='''Filename, variable name, or expression. Wildcard characters
                such as '*' and '?' are allowed for filenames.''')
        parser.add_argument('-k', '--kernel',
            help='''kernel in which variables will be previewed. By default
            the variable will be previewed in the current kernel of the cell.''')
        parser.add_argument('-w', '--workflow', action='store_true',
            help='''Preview notebook workflow''')
        # this option is currently hidden
        parser.add_argument('-s', '--style', choices=['table', 'scatterplot'],
            help='''Option to visualize dataframe, which by default is "table".
            The %%preview magic also accepts arbitrary additional keyword
            arguments, which would be interpreted by individual style. Passing
            '-h' with '--style' would display the usage information of particular
            style.''')
        parser.add_argument('--off', action='store_true',
            help='''Turn off file preview''')
        loc = parser.add_mutually_exclusive_group()
        loc.add_argument('-p', '--panel', action='store_true',
            help='''Preview in side panel even if the panel is currently closed''')
        loc.add_argument('-n', '--notebook', action='store_true',
            help='''Preview in the main notebook.''')
        parser.error = self._parse_error
        return parser

    def get_set_parser(self):
        parser = argparse.ArgumentParser(prog='%set',
            description='''Set persistent command line options for SoS runs.''')
        parser.error = self._parse_error
        return parser

    def get_shutdown_parser(self):
        parser = argparse.ArgumentParser(prog='%shutdown',
            description='''Shutdown or restart specified subkernel''')
        parser.add_argument('kernel', nargs='?',
            help='''Name of the kernel to be restarted, default to the
            current running kernel.''')
        parser.add_argument('-r', '--restart', action='store_true',
            help='''Restart the kernel''')
        parser.error = self._parse_error
        return parser

    def get_run_parser(self):
        parser = argparse.ArgumentParser(prog='%run',
            description='''Execute the current cell with specified command line
            arguments. Arguments set by magic %set will be appended at the
            end of command line''')
        parser.error = self._parse_error
        return parser

    def get_sosrun_parser(self):
        parser = argparse.ArgumentParser(prog='%sosrun',
            description='''Execute the entire notebook with steps consisting of SoS
            cells (cells with SoS kernel) with section header, with specified command
            line arguments. Arguments set by magic %set will be appended at the
            end of command line''')
        parser.error = self._parse_error
        return parser

    def get_save_parser(self):
        parser = argparse.ArgumentParser(prog='%save',
            description='''Save the jupyter notebook as workflow (consisting of all sos
            steps defined in cells starting with section header) or a HTML report to
            specified file.''')
        parser.add_argument('filename', nargs='?',
            help='''Filename of saved report or script. Default to notebookname with file
            extension determined by option --to.''')
        parser.add_argument('-f', '--force', action='store_true',
            help='''If destination file already exists, overwrite it.''')
        parser.add_argument('-a', '--append', action='store_true',
            help='''If destination file already exists, append to it.''')
        parser.add_argument('-x', '--set-executable', dest = "setx", action='store_true',
            help='''Set `executable` permission to saved script.''')
        parser.error = self._parse_error
        return parser

    def get_sossave_parser(self):
        parser = argparse.ArgumentParser(prog='%sossave',
            description='''Save the jupyter notebook as workflow (consisting of all sos
            steps defined in cells starting with section header) or a HTML report to
            specified file.''')
        parser.add_argument('filename', nargs='?',
            help='''Filename of saved report or script. Default to notebookname with file
            extension determined by option --to.''')
        parser.add_argument('-t', '--to', dest='__to__', choices=['sos', 'html'],
            help='''Destination format, default to sos.''')
        parser.add_argument('-c', '--commit', action='store_true',
            help='''Commit the saved file to git directory using command
            git commit FILE''')
        parser.add_argument('-a', '--all', action='store_true',
            help='''The --all option for sos convert script.ipynb script.sos, which
            saves all cells and their metadata to a .sos file, that contains all input
            information of the notebook but might not be executable in batch mode.''')
        parser.add_argument('-m', '--message',
            help='''Message for git commit. Default to "save FILENAME"''')
        parser.add_argument('-p', '--push', action='store_true',
            help='''Push the commit with command "git push"''')
        parser.add_argument('-f', '--force', action='store_true',
            help='''If destination file already exists, overwrite it.''')
        parser.add_argument('-x', '--set-executable', dest = "setx", action='store_true',
            help='''Set `executable` permission to saved script.''')
        parser.add_argument('--template', default='sos-report',
            help='''Template to generate HTML output, default to sos-report,
            which uses a control panel to control the display of contents.''')
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
            description='''Get specified variables from another kernel, which is
                by default the SoS kernel.''')
        parser.add_argument('--from', dest='__from__',
            help='''Name of kernel from which the variables will be obtained.
                Default to the SoS kernel.''')
        parser.add_argument('vars', nargs='*',
            help='''Names of SoS variables''')
        parser.error = self._parse_error
        return parser

    def get_put_parser(self):
        parser = argparse.ArgumentParser(prog='%put',
            description='''Put specified variables in the subkernel to another
            kernel, which is by default the SoS kernel.''')
        parser.add_argument('--to', dest='__to__',
            help='''Name of kernel from which the variables will be obtained.
                Default to the SoS kernel.''')
        parser.add_argument('vars', nargs='*',
            help='''Names of SoS variables''')
        parser.error = self._parse_error
        return parser

    def get_debug_parser(self):
        parser = argparse.ArgumentParser(prog='%debug',
            description='''Turn on or off debug information''')
        parser.add_argument('status', choices=['on', 'off'],
            help='''Turn on or off debugging''')
        parser.error = self._parse_error
        return parser

    def get_taskinfo_parser(self):
        parser = argparse.ArgumentParser(prog='%taskinfo',
            description='''Get information on specified task. By default
            sos would query against all running task queues but it would
            start a task queue and query status if option -q is specified.
            ''')
        parser.add_argument('task', help='ID of task')
        parser.add_argument('-q', '--queue',
            help='''Task queue on which the task is executed.''')
        parser.error = self._parse_error
        return parser

    def get_tasks_parser(self):
        parser = argparse.ArgumentParser(prog='%tasks',
            description='''Show a list of tasks from specified queue''')
        parser.add_argument('tasks', nargs='*', help='ID of tasks')
        parser.add_argument('-s', '--status', nargs='*',
            help='''Display tasks of specified status. Default to all.'''),
        parser.add_argument('-q', '--queue',
            help='''Task queue on which the tasks are retrived.''')
        parser.add_argument('--age', help='''Limit to tasks that is created more than
            (default) or within specified age. Value of this parameter can be in units
            s (second), m (minute), h (hour), or d (day, default), with optional
            prefix + for older (default) and - for younder than specified age.''')
        parser.error = self._parse_error
        return parser

    def get_skip_parser(self):
        parser = argparse.ArgumentParser(prog='%skip',
            description='''Skip the rest of the cell content and return''')
        parser.error = self._parse_error
        return parser

    def get_toc_parser(self):
        parser = argparse.ArgumentParser(prog='%toc',
            description='''Display toc in the side panel''')
        parser.error = self._parse_error
        return parser

    def get_render_parser(self):
        parser = argparse.ArgumentParser(prog='%render',
            description='''Treat the output of a SoS cell as another format, default to markdown.''')
        parser.add_argument('format', default='Markdown', nargs='?',
            help='''Format to render output of cell, default to Markdown, but can be any
            format that is supported by the IPython.display module such as HTML, Math, JSON,
            JavaScript and SVG.''')
        parser.error = self._parse_error
        return parser

    def get_sessioninfo_parser(self):
        parser = argparse.ArgumentParser(prog='%sessioninfo',
            description='''List the session info of all subkernels, and information
            stored in variable sessioninfo''')
        parser.error = self._parse_error
        return parser

    def find_kernel(self, name, kernel=None, language=None, color=None, notify_frontend=True):
        # find from subkernel name
        def update_existing(idx):
            x = self._kernel_list[idx]
            if (kernel is not None and kernel != x[1]) or (language is not None and language != x[2]):
                raise ValueError('Cannot change kernel or language of predefined subkernel {}'.format(name))
            if color is not None:
                if color == 'default':
                    if self._kernel_list[idx][2]:
                        self._kernel_list[idx][3] = self._supported_languages[self._kernel_list[idx][2]].background_color
                    else:
                        self._kernel_list[idx][3] = ''
                else:
                    self._kernel_list[idx][3] = color
                if notify_frontend:
                    self.send_frontend_msg('kernel-list', self.get_kernel_list())

        # find from language name (subkernel name, which is usually language name)
        for idx,x in enumerate(self.get_kernel_list()):
            if x[0] == name:
                if x[0] == 'SoS' or x[2] or language is None:
                    update_existing(idx)
                    return x
                else:
                    if not kernel:
                        kernel = name
                    break
        # find from kernel name
        for idx,x in enumerate(self._kernel_list):
            if x[1] == name:
                # if exist language or no new language defined.
                if x[2] or language is None:
                    update_existing(idx)
                    return x
                else:
                    # otherwise, try to use the new language
                    kernel = name
                    break
        # now, no kernel is found, name has to be a new name and we need some definition
        # if kernel is defined
        def add_or_replace(kdef):
            for idx, x in enumerate(self._kernel_list):
                if x[0] == kdef[0]:
                    self._kernel_list[idx] = kdef
                    return self._kernel_list[idx]
            else:
                self._kernel_list.append(kdef)
                return self._kernel_list[-1]

        if kernel is not None:
            # in this case kernel should have been defined in kernel list
            if kernel not in [x[1] for x in self._kernel_list]:
                raise ValueError('Unrecognized Jupyter kernel name {}. Please make sure it is properly installed and appear in the output of command "jupyter kenelspec list"'.format(kernel))
            # now this a new instance for an existing kernel
            kdef = [x for x in self._kernel_list if x[1] == kernel][0]
            if not language:
                if color == 'default':
                    if kdef[2]:
                        color = self._supported_languages[kdef[2]].background_color
                    else:
                        color = kdef[3]
                new_def = add_or_replace([name, kdef[1], kdef[2], kdef[3] if color is None else color])
                if notify_frontend:
                    self.send_frontend_msg('kernel-list', self.get_kernel_list())
                return new_def
            else:
                # if language is defined,
                if ':' in language:
                    # if this is a new module, let us create an entry point and load
                    from pkg_resources import EntryPoint
                    mn, attr = language.split(':', 1)
                    ep = EntryPoint(name=kernel, module_name=mn, attrs=tuple(attr.split('.')))
                    try:
                        plugin = ep.resolve()(self)
                        self._supported_languages[name] = plugin
                        # for convenience, we create two entries for, e.g. R and ir
                        # but only if there is no existing definition
                        if name != plugin.kernel_name and plugin.kernel_name not in self._supported_languages:
                            self._supported_languages[plugin.kernel_name] = plugin
                    except Exception as e:
                        raise RuntimeError('Failed to load language {}: {}'.format(language, e))
                    #
                    if color == 'default':
                        color = plugin.background_color
                    new_def = add_or_replace([name, kdef[1], kernel, kdef[3] if color is None else color])
                else:
                    # if should be defined ...
                    if language not in self._supported_languages:
                        raise RuntimeError('Unrecognized language definition {}, which should be a known language name or a class in the format of package.module:class'.format(language))
                    #
                    self._supported_languages[name] = self._supported_languages[language]
                    if color == 'default':
                        color = self._supported_languages[name].background_color
                    new_def = add_or_replace([name, kdef[1], language, kdef[3] if color is None else color])
                if notify_frontend:
                    self.send_frontend_msg('kernel-list', self.get_kernel_list())
                return new_def
        elif language is not None:
            # kernel is not defined and we only have language
            if ':' in language:
                # if this is a new module, let us create an entry point and load
                from pkg_resources import EntryPoint
                mn, attr = language.split(':', 1)
                ep = EntryPoint(name='__unknown__', module_name=mn, attrs=tuple(attr.split('.')))
                try:
                    plugin = ep.resolve()(self)
                    self._supported_languages[name] = plugin
                except Exception as e:
                    raise RuntimeError('Failed to load language {}: {}'.format(language, e))
                #
                if plugin.kernel_name not in [x[1] for x in self._kernel_list]:
                    raise ValueError('Unrecognized Jupyter kernel name {} defined by language {}. Please make sure it is properly installed and appear in the output of command "jupyter kenelspec list"'.format(
                        plugin.kernel_name, language))

                if color == 'default':
                    color = plugin.background_color
                new_def = add_or_replace([name, plugin.kernel_name, plugin.kernel_name, plugin.background_color if color is None else color])
            else:
                # if should be defined ...
                if language not in self._supported_languages:
                    raise RuntimeError('Unrecognized language definition {}'.format(language))
                #
                if self._supported_languages[language].kernel_name not in [x[1] for x in self._kernel_list]:
                    raise ValueError('Unrecognized Jupyter kernel name {} defined by language {}. Please make sure it is properly installed and appear in the output of command "jupyter kenelspec list"'.format(
                        self._supported_languages[language].kernel_name, language))

                new_def = add_or_replace([
                    name, self._supported_languages[language].kernel_name, language,
                        self._supported_languages[language].background_color if color is None or color == 'default' else color])

            self.send_frontend_msg('kernel-list', self.get_kernel_list())
            return new_def
        else:
            raise ValueError('No pre-defined subkernel named {} is found. Please define it with one or both of parameters --kernel and --language'.format(name))

    def get_supported_languages(self):
        if self._supported_languages is not None:
            return self._supported_languages
        group = 'sos_languages'
        self._supported_languages = {}
        for entrypoint in pkg_resources.iter_entry_points(group=group):
            # Grab the function that is the actual plugin.
            name = entrypoint.name
            try:
                plugin = entrypoint.load()(self)
                self._supported_languages[name] = plugin
                # for convenience, we create two entries for, e.g. R and ir
                if name != plugin.kernel_name:
                    self._supported_languages[plugin.kernel_name] = plugin
            except Exception as e:
                self.log.error('Failed to load language {}: {}'.format(entrypoint.name, e))
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
        self.kernel = 'SoS'
        # a dictionary of started kernels, with the format of
        #
        # 'R': ['ir', 'sos.R.sos_R', '#FFEEAABB']
        #
        # Note that:
        #
        # 'R' is the displayed name of the kernel.
        # 'ir' is the kernel name.
        # 'sos.R.sos_R' is the language module.
        # '#FFEEAABB' is the background color
        #
        self.kernels = {}
        #self.shell = InteractiveShell.instance()
        self.format_obj = self.shell.display_formatter.format

        self.previewers = None
        self.original_keys = None
        self._supported_languages = None
        self._completer = None
        self._inspector = None
        self._real_execution_count = 1
        self._execution_count = 1
        self._debug_mode = False
        self._use_panel = None
        self._resume_execution = False

        self.frontend_comm = None
        self.comm_manager.register_target('sos_comm', self.sos_comm)
        self.cell_idx = None
        self.my_tasks = {}
        #
        self._workflow_mode = False
        self._render_result = False
        env.__task_notifier__ = self.notify_task_status

    def handle_taskinfo(self, task_id, task_queue, side_panel=None):
        # requesting information on task
        from sos.hosts import Host
        host = Host(task_queue)
        result = host._task_engine.query_tasks([task_id], verbosity=2, html=True)
        #log_to_file(result)
        if side_panel is True:
            self.send_frontend_msg('display_data',
                {'metadata': {},
                 'data': {'text/plain': result,
                 'text/html': HTML(result).data
                 }})
        else:
            self.send_response(self.iopub_socket, 'display_data',
                {'metadata': {},
                 'data': {'text/plain': result,
                 'text/html': HTML(result).data
                 }})
        # now, there is a possibility that the status of the task is different from what
        # task engine knows (e.g. a task is rerun outside of jupyter). In this case, since we
        # already get the status, we should update the task engine...
        #
        # <tr><th align="right"  width="30%">Status</th><td align="left">completed</td></tr>
        status = result.split('>Status<', 1)[-1].split('</td',1)[0].split('>')[-1]
        host._task_engine.update_task_status(task_id, status)


    def handle_tasks(self, tasks, queue='localhost', status=None, age=None):
        from sos.hosts import Host
        try:
            host = Host(queue)
        except Exception as e:
            self.warn('Invalid task queu {}: {}'.format(queue, e))
            return
        # get all tasks
        for tid, tst, tdt in host._task_engine.monitor_tasks(tasks, status=None, age=age):
            self.notify_task_status(['new-status', queue, tid, tst, tdt])
        self.send_frontend_msg('update-duration', {})

    def handle_sessioninfo(self):
        #
        from sos.utils import loaded_modules
        result = OrderedDict()
        #
        from sos._version import __version__
        result['SoS'] = [('SoS Version', __version__)]
        result['SoS'].extend(loaded_modules(env.sos_dict))
        #
        cur_kernel = self.kernel
        try:
            for kernel in self.kernels.keys():
                kinfo = self.find_kernel(kernel)
                self.switch_kernel(kernel)
                result[kernel] = [
                        ('Kernel', kinfo[1]),
                        ('Language', kinfo[2])
                ]
                if kernel not in self.supported_languages:
                    continue
                lan = self.supported_languages[kernel]
                if hasattr(lan, 'sessioninfo'):
                    try:
                        sinfo = lan.sessioninfo()
                        if isinstance(sinfo, str):
                            result[kernel].append([sinfo])
                        elif isinstance(sinfo, dict):
                            result[kernel].extend(list(sinfo.items()))
                        elif isinstance(sinfo, list):
                            result[kernel].extend(sinfo)
                        else:
                            self.warn('Unrecognized session info: {}'.format(sinfo))
                    except Exception as e:
                        self.warn('Failed to obtain sessioninfo of kernel {}: {}'.format(kernel, e))
        finally:
            self.switch_kernel(cur_kernel)
        #
        if 'sessioninfo' in env.sos_dict:
            result.update(env.sos_dict['sessioninfo'])
        #
        res = ''
        for key, items in result.items():
            res += '<p class="session_section">{}</p>\n'.format(key)
            res += '<table class="session_info">\n'
            for item in items:
                res += '<tr>\n'
                if isinstance(item, str):
                    res += '<td colspan="2"><pre>{}</pre></td>\n'.format(item)
                elif len(item) == 1:
                    res += '<td colspan="2"><pre>{}</pre></td>\n'.format(item[0])
                elif len(item) == 2:
                    res += '<th>{}</th><td><pre>{}</pre></td>\n'.format(item[0], item[1])
                else:
                    self.warn('Invalid session info item of type {}: {}'.format(item.__class__.__name__, short_repr(item)))
                res += '</tr>\n'
            res += '</table>\n'
        self.send_response(self.iopub_socket, 'display_data',
            { 'source': 'SoS', 'metadata': {},
              'data': {'text/html': HTML(res).data} })

    def sos_comm(self, comm, msg):
        # record frontend_comm to send messages
        self.frontend_comm = comm

        @comm.on_msg
        def handle_frontend_msg(msg):
            content = msg['content']['data']
            #log_to_file(msg)
            for k,v in content.items():
                if k == 'list-kernel':
                    self.send_frontend_msg('kernel-list', self.get_kernel_list(v))
                elif k == 'kill-task':
                    # kill specified task
                    from sos.hosts import Host
                    Host(v[1])._task_engine.kill_tasks([v[0]])
                    self.notify_task_status(['change-status', v[1], v[0], 'aborted'])
                elif k == 'resume-task':
                    # kill specified task
                    from sos.hosts import Host
                    Host(v[1])._task_engine.resume_task(v[0])
                    self.notify_task_status(['change-status', v[1], v[0], 'pending'])
                elif k == 'task-info':
                    self.handle_taskinfo(v[0], v[1], side_panel=True)
                elif k == 'update-task-status':
                    if not isinstance(v, list):
                        continue
                    # split by host ...
                    host_status = defaultdict(list)
                    for name in v:
                        if not name.startswith('status_'):
                            continue
                        try:
                            tqu, tid = name[7:].rsplit('_', 1)
                        except:
                            # incorrect ID...
                            continue
                        host_status[tqu].append(tid)
                    #log_to_file(host_status)
                    #
                    from sos.hosts import Host
                    for tqu, tids in host_status.items():
                        try:
                            h = Host(tqu)
                        except:
                            continue
                        for tid in tids:
                            tst = h._task_engine.check_task_status(tid, unknown='unknown')
                            self.notify_task_status(['change-status', tqu, tid, tst])
                    self.send_frontend_msg('update-duration', {})
                else:
                    # this somehow does not work
                    self.warn('Unknown message {}: {}'.format(k, v))

    def notify_task_status(self, task_status):
        status_class = {
            'pending': 'fa-square-o',
            'submitted': 'fa-spinner',
            'running': 'fa-spinner fa-pulse fa-spin',
            'result-ready': 'fa-files-o',
            'completed': 'fa-check-square-o',
            'failed': 'fa-times-circle-o',
            'aborted': 'fa-frown-o',
            'result-mismatch': 'fa-question',
            'unknown': 'fa-question',
            }

        action_class = {
            'pending': 'fa-stop',
            'submitted': 'fa-stop',
            'running': 'fa-stop',
            'result-ready': 'fa-play',
            'completed': 'fa-play',
            'failed':  'fa-play',
            'aborted':  'fa-play',
            'result-mismatch': 'fa-play',
            'unknown': 'fa-question',
        }

        action_func = {
            'pending': 'kill_task',
            'submitted': 'kill_task',
            'running': 'kill_task',
            'result-ready': 'resume_task',
            'completed': 'resume_task',
            'failed':  'resume_task',
            'aborted':  'resume_task',
            'result-mismatch': 'resume_task',
            'unknown': 'function(){}',
        }

        if task_status[0] == 'new-status':
            tqu, tid, tst, tdt = task_status[1:]
            self.send_response(self.iopub_socket, 'display_data',
                {
                    'source': 'SoS',
                    'metadata': {},
                    'data': { 'text/html':
                        HTML('''<table id="table_{0}_{1}" style="border: 0px"><tr style="border: 0px">
                        <td style="border: 0px">
                        <i id="status_{0}_{1}"
                            class="fa fa-2x fa-fw {2}"
                            onmouseover="$('#status_{0}_{1}').addClass('{3}').removeClass('{2}')"
                            onmouseleave="$('#status_{0}_{1}').addClass('{2}').removeClass('{3}')"
                            onclick="{4}('{1}', '{0}')"
                        ></i> </td>
                        <td style="border:0px"><a onclick="task_info('{1}', '{0}')"><pre>{1}</pre></a></td>
                        <td style="border:0px">&nbsp;</td>
                        <td style="border:0px;text-align=right;">
                        <pre><time id="duration_{0}_{1}" datetime="{5}">{6}</time></pre></td>
                        </tr>
                        </table>'''.format(tqu, tid, status_class[tst], action_class[tst], action_func[tst], tdt*1000,
                            PrettyRelativeTime(time.time() - tdt))).data
                        }
                })
            # keep tracks of my tasks to avoid updating status of
            # tasks that does not belong to the notebook
            self.my_tasks[(tqu, tid)] = time.time()
        elif task_status[0] == 'remove-task':
            tqu, tid = task_status[1:]
            if (tqu, tid) in self.my_tasks:
                self.send_frontend_msg('remove-task', [tqu, tid])
        elif task_status[0] == 'change-status':
            tqu, tid, tst = task_status[1:]
            self.send_frontend_msg('task-status', [tqu, tid, tst, status_class[tst], action_class[tst], action_func[tst]])
            self.my_tasks[(tqu, tid)] = time.time()
        elif task_status[0] == 'pulse-status':
            tqu, tid, tst = task_status[1:]
            if (tqu, tid) in self.my_tasks:
                if time.time() - self.my_tasks[(tqu, tid)] < 20:
                    # if it has been within the first 20 seconds of new or updated message
                    # can confirm to verify it has been successfully delivered. Otherwise
                    # ignore such message
                    self.send_frontend_msg('task-status', [tqu, tid, tst, status_class[tst], action_class[tst], action_func[tst]])
        else:
            raise RuntimeError('Unrecognized status change message {}'.format(task_status))

    def send_frontend_msg(self, msg_type, msg={}):
        # if comm is never created by frontend, the kernel is in test mode without frontend
        if not self.frontend_comm:
            return
        if self._use_panel is False and msg_type in ('display_data', 'stream', 'preview-input'):
            if msg_type in ('display_data', 'stream'):
                self.send_response(self.iopub_socket, msg_type, msg)
            elif msg_type == 'preview-input':
                self.send_response(self.iopub_socket, 'display_data',
                    {
                        'source': 'SoS',
                        'metadata': {},
                        'data': { 'text/html': HTML('<div class="sos_hint">{}</div>'.format(msg)).data}
                    })
        else:
            self.frontend_comm.send(msg, {'msg_type': msg_type})

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
        message = str(message).rstrip() + '\n'
        if message.strip():
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stderr', 'text': message})

    def get_magic_and_code(self, code, warn_remaining=False):
        if code.startswith('%') or code.startswith('!'):
            lines = re.split(r'(?<!\\)\n', code, 1)
            # remove lines joint by \
            lines[0] = lines[0].replace('\\\n', '')
        else:
            lines[0] = code.split('\n', 1)

        pieces = lines[0].strip().split(None, 1)
        if len(pieces) == 2:
            command_line = pieces[1]
        else:
            command_line = ''
        remaining_code = lines[1] if len(lines) > 1 else ''
        if warn_remaining and remaining_code.strip():
            self.warn('Statement {} ignored'.format(short_repr(remaining_code)))
        return command_line, remaining_code

    def run_cell(self, code, silent, store_history, on_error=None):
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
        self.KC.execute(code, silent=silent, store_history=store_history)

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
                    if silent and msg_type in ['display_data', 'stream', 'execute_result']:
                        continue
                    # NOTE: we do not send status of sub kernel alone because
                    # these are generated automatically during the execution of
                    # "this cell" in SoS kernel
                    if self._render_result and msg_type == 'stream' and sub_msg['content']['name'] == 'stdout':
                        format_dict, md_dict = self.format_obj(self.render_result(sub_msg['content']['text']))
                        self.send_response(self.iopub_socket, 'display_data',
                            {'source': 'SoS', 'metadata': md_dict,
                             'data': format_dict
                            })
                    elif msg_type == 'error':
                        if on_error and not self._debug_mode:
                            self.warn(on_error)
                        else:
                            self.send_response(self.iopub_socket, msg_type, sub_msg['content'])
                    elif msg_type == 'stream' and sub_msg['content']['name'] == 'stderr':
                        self.warn(sub_msg['content']['text'])
                    else:
                        self.send_response(self.iopub_socket, msg_type, sub_msg['content'])
        #
        # now get the real result
        reply = self.KC.get_shell_msg(timeout=10)
        reply['content']['execution_count'] = self._execution_count
        return reply['content']

    def switch_kernel(self, kernel, in_vars=[], ret_vars=[], kernel_name=None, language=None, color=None):
        # switching to a non-sos kernel
        if not kernel:
            kinfo = self.find_kernel(self.kernel)
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': '''\
Current subkernel: {} (kernel={}, language={}, color="{}")
Active subkernels: {}
Available subkernels:\n{}'''.format(
        kinfo[0], kinfo[1], kinfo[2] if kinfo[2] else "undefined", kinfo[3],
        ', '.join(self.kernels.keys()),
        '\n'.join(['    {} ({})'.format(x[0], x[1]) for x in self.get_kernel_list()]))})
            return
        kinfo = self.find_kernel(kernel, kernel_name, language, color)
        if kinfo[0] == self.kernel:
            # the same kernel, do nothing?
            # but the senario can be
            #
            # kernel in SoS
            # cell R
            # %use R -i n
            #
            # SoS get:
            #
            # %softwidth --default-kernel R --cell-kernel R
            # %use R -i n
            #
            # Now, SoS -> R without variable passing
            # R -> R should honor -i n

            # or, when we randomly jump cells, we should more aggreessively return
            # automatically shared variables to sos (done by the following) (#375)
            if kinfo[0] != 'SoS':
                self.switch_kernel('SoS')
                self.switch_kernel(kinfo[0], in_vars, ret_vars)
        elif kinfo[0] == 'SoS':
            # switch from non-sos to sos kernel
            self.handle_magic_put(self.RET_VARS)
            self.RET_VARS = []
            self.kernel = 'SoS'
        elif self.kernel != 'SoS':
            # not to 'sos' (kernel != 'sos'), see if they are the same kernel under
            self.switch_kernel('SoS', in_vars, ret_vars)
            self.switch_kernel(kinfo[0], in_vars, ret_vars)
        else:
            if self._debug_mode:
                self.warn('Switch from {} to {}'.format(self.kernel, kinfo[0]))
            # case when self.kernel == 'sos', kernel != 'sos'
            # to a subkernel
            new_kernel = False
            if kinfo[0] not in self.kernels:
                # start a new kernel
                try:
                    self.kernels[kinfo[0]] = manager.start_new_kernel(
                            startup_timeout=60, kernel_name=kinfo[1], cwd=os.getcwd())
                    new_kernel = True
                except Exception as e:
                    self.warn('Failed to start kernel "{}". Use "jupyter kernelspec list" to check if it is installed: {}'.format(kernel, e))
                    return
            self.KM, self.KC = self.kernels[kinfo[0]]
            self.RET_VARS = ret_vars
            self.kernel = kinfo[0]
            if new_kernel and self.kernel in self.supported_languages:
                init_stmts = self.supported_languages[self.kernel].init_statements
                if init_stmts:
                    self.run_cell(init_stmts, True, False)
            #
            self.handle_magic_get(in_vars)

    def shutdown_kernel(self, kernel, restart=False):
        kernel = self.find_kernel(kernel)[0]
        if kernel == 'SoS':
            # cannot restart myself ...
            self.warn('Cannot restart SoS kernel from within SoS.')
        elif kernel:
            if kernel not in self.kernels:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': '{} is not running'.format(kernel) })
            elif restart:
                orig_kernel = self.kernel
                try:
                    # shutdown
                    self.shutdown_kernel(kernel)
                    # switch back to kernel (start a new one)
                    self.switch_kernel(kernel)
                finally:
                    # finally switch to starting kernel
                    self.switch_kernel(orig_kernel)
            else:
                # shutdown
                if self.kernel == kernel:
                    self.switch_kernel('SoS')
                try:
                    self.kernels[kernel][0].shutdown_kernel(restart=False)
                except Exception as e:
                    self.warn('Failed to shutdown kernel {}: {}\n'.format(kernel, e))
                finally:
                    self.kernels.pop(kernel)
        else:
            self.send_response(self.iopub_socket, 'stream',
                {'name': 'stdout', 'text': 'Specify one of the kernels to shutdown: SoS{}\n'
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
        parser = self.get_dict_parser()
        try:
            args = parser.parse_args(shlex.split(line))
        except SystemExit:
            return

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
                self.warn('Magic %set cannot set positional argument, {} provided.\n'.format(options))
            else:
                self.options = options.strip()
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Set sos options to "{}"\n'.format(self.options)})
        else:
            if self.options:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Reset sos options from "{}" to ""\n'.format(self.options)})
                self.options = ''
            else:
                self.send_response(self.iopub_socket, 'stream',
                    {'name': 'stdout', 'text': 'Usage: set persistent sos command line options such as "-v 3" (debug output)\n'})

    def handle_magic_get(self, items, from_kernel=None, explicit=False):
        if from_kernel is None or from_kernel.lower() == 'sos':
            # autmatically get all variables with names start with 'sos'
            default_items = [x for x in env.sos_dict.keys() if x.startswith('sos') and x not in self.original_keys]
            items = default_items if not items else items + default_items
            for item in items:
                if item not in env.sos_dict:
                    self.warn('Variable {} does not exist'.format(item))
                    return
            if not items:
                return
            if self.kernel in self.supported_languages:
                lan = self.supported_languages[self.kernel]
                try:
                    lan.get_vars(items)
                except Exception as e:
                    self.warn('Failed to get variable: {}\n'.format(e))
                    return
            elif self.kernel == 'SoS':
                self.warn('Magic %get without option --kernel can only be executed by subkernels')
                return
            else:
                if explicit:
                    self.warn('Language {} does not support magic %get.'.format(self.kernel))
                return
        elif self.kernel.lower() == 'sos':
            # if another kernel is specified and the current kernel is sos
            # we get from subkernel
            try:
                self.switch_kernel(from_kernel)
                self.handle_magic_put(items)
            finally:
                self.switch_kernel('SoS')
        else:
            # if another kernel is specified, we should try to let that kernel pass
            # the variables to this one directly
            try:
                my_kernel = self.kernel
                self.switch_kernel(from_kernel)
                # put stuff to sos or my_kernel directly
                self.handle_magic_put(items, to_kernel=my_kernel, explicit=explicit)
            finally:
                # then switch back
                self.switch_kernel(my_kernel)

    def get_response(self, statement, msg_types, name=None):
        # get response of statement of specific msg types.
        responses = []
        self.KC.execute(statement, silent=False, store_history=False)
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
                    if msg_type in msg_types and (name is None or sub_msg['content'].get('name', None) in name):
                        responses.append([msg_type, sub_msg['content']])
                    else:
                        if self._debug_mode:
                            self.warn('{}: {}'.format(msg_type, sub_msg['content']))
                        self.send_response(self.iopub_socket, msg_type,
                            sub_msg['content'])
        if not responses and self._debug_mode:
            self.warn('Failed to get a response from message type {}'.format(msg_types))

        return responses

    def handle_magic_put(self, items, to_kernel=None, explicit=False):
        if self.kernel.lower() == 'sos':
            if to_kernel is None:
                self.warn('Magic %put without option --kernel can only be executed by subkernels')
                return
            # if another kernel is specified and the current kernel is sos
            try:
                # switch to kernel and bring in items
                self.switch_kernel(to_kernel, in_vars=items)
            finally:
                # switch back
                self.switch_kernel('SoS')
        else:
            # put to sos kernel or another kernel
            #
            # items can be None if unspecified
            if not items:
                # we do not simply return because we need to return default variables (with name startswith sos
                items = []
            if self.kernel not in self.supported_languages:
                if explicit:
                    self.warn('Subkernel {} does not support magic %put.'.format(self.kernel))
                return
            #
            lan = self.supported_languages[self.kernel]
            # pass language name to to_kernel
            if to_kernel:
                objects = lan.put_vars(items, to_kernel=self.find_kernel(to_kernel)[2])
            else:
                objects = lan.put_vars(items, to_kernel='SoS')
            if isinstance(objects, dict):
                # returns a SOS dictionary
                try:
                    env.sos_dict.update(objects)
                except Exception as e:
                    self.warn('Failed to execute %put: {}'.format(e))
                    return

                if to_kernel is None:
                    return
                # if another kernel is specified and the current kernel is not sos
                # we need to first put to sos then to another kernel
                try:
                    my_kernel = self.kernel
                    # switch to the destination kernel and bring in vars
                    self.switch_kernel(to_kernel, in_vars=items)
                finally:
                    # switch back to the original kernel
                    self.switch_kernel(my_kernel)
            elif isinstance(objects, str):
                # an statement that will be executed in the destination kernel
                if to_kernel is None or to_kernel == 'SoS':
                    # evaluate in SoS, this should not happen or rarely happen
                    # because the subkernel should return a dictionary for SoS kernel
                    try:
                        exec(objects, env.sos_dict._dict)
                    except Exception as e:
                        self.warn('Failed to put variables {} to SoS kernel: {}'.format(items, e))
                        return
                try:
                    my_kernel = self.kernel
                    # switch to the destination kernel
                    self.switch_kernel(to_kernel)
                    # execute the statement to pass variables directly to destination kernel
                    self.run_cell(objects, True, False)
                finally:
                    # switch back to the original kernel
                    self.switch_kernel(my_kernel)
            else:
                self.warn('Unrecognized return value of type {} for action %put'.format(object.__class__.__name__))
                return


    def _interpolate_option(self, option, quiet=False):
        # interpolate command
        try:
            new_option = interpolate(option, sigil='${ }', local_dict=env.sos_dict._dict)
            if new_option != option and not quiet:
                self.send_response(self.iopub_socket, 'display_data',
                    {
                        'source': 'SoS',
                        'metadata': {},
                        'data': {
                            'text/html': HTML('<div class="sos_hint">{}</div>'.format(
                            new_option.strip() + '<br>## -- End interpolated command --<br>')).data }
                        })
            return new_option
        except Exception as e:
            self.warn('Failed to interpolate {}: {}\n'.format(short_repr(option), e))
            return None

    def handle_magic_preview(self, items, kernel=None, style=None):
        # find filenames and quoted expressions
        self.send_frontend_msg('preview-input', '%preview {}'.format(' '.join(items)))
        # expand items
        handled = [False for x in items]
        for idx, item in enumerate(items):
            try:
                if os.path.isfile(item):
                    handled[idx] = True
                    self.preview_file(item, style)
                    continue
                else:
                    import glob
                    files = glob.glob(item)
                    if files:
                        for pfile in files:
                            self.preview_file(pfile, style)
                        handled[idx] = True
                        continue
            except Exception as e:
                self.warn('\n> Failed to preview file {}: {}'.format(item, e))
                continue

        # all are files
        if all(handled):
            return

        # non-sos kernel
        use_sos = kernel in ('sos', 'SoS') or (kernel is None and self.kernel == 'SoS')
        orig_kernel = self.kernel
        if kernel is not None and self.kernel != self.find_kernel(kernel)[0]:
            self.switch_kernel(kernel)
        if self._use_panel:
            self.send_frontend_msg('preview-kernel', self.kernel)
        try:
            for item in (x for x,y in zip(items, handled) if not y):
                try:
                    # quoted
                    if (item.startswith('"') and item.endswith('"')) or \
                        (item.startswith("'") and item.endswith("'")):
                        try:
                            item = eval(item)
                        except:
                            pass
                    if use_sos:
                        obj_desc, preview = self.preview_var(item, style)
                        if preview is None:
                            continue
                        else:
                            format_dict, md_dict = preview
                        self.send_frontend_msg('display_data',
                            {'metadata': {},
                            'data': {'text/plain': '>>> ' + item + ':\n',
                                'text/html': HTML('<div class="sos_hint">> {}: {}</div>'
                                    .format(item, obj_desc)).data
                                }
                            })
                        self.send_frontend_msg('display_data',
                            {'execution_count': self._execution_count, 'data': format_dict,
                            'metadata': md_dict})
                    else:
                        self.send_frontend_msg('display_data',
                            {'metadata': {},
                            'data': {'text/plain': '>>> ' + item + ':\n',
                                'text/html': HTML('<div class="sos_hint">> {}:</div>'.format(item)).data
                                }
                            })
                        # evaluate
                        responses = self.get_response(item, ['stream', 'display_data', 'execution_result'])
                        for response in responses:
                            self.send_frontend_msg(response[0], response[1])
                except Exception as e:
                    self.send_frontend_msg('stream',
                        {'name': 'stderr', 'text': '> Failed to preview file or expression {}: {}'.format(item, e)})
        finally:
            self.switch_kernel(orig_kernel)

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
                if isinstance(cmd, str):
                    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                else:
                    p = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
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
        with self.redirect_sos_io():
            try:
                # record input and output
                fopt = ''
                res = runfile(code=code, args=self.options + fopt, kernel=self)
                self.send_result(res, silent)
            except PendingTasks as e:
                # send cell index and task IDs to frontend
                self.send_frontend_msg('tasks-pending', [self.cell_idx, e.tasks])
                return
            except Exception as e:
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
                    input_files = [x for x in input_files if isinstance(x, str)]
            else:
                input_files = []
            if 'output' in env.sos_dict:
                output_files = env.sos_dict['output']
                if output_files is None:
                    output_files = []
                else:
                    output_files = [x for x in output_files if isinstance(x, str)]
            else:
                output_files = []
            # use a table to list input and/or output file if exist
            if output_files:
                self.send_frontend_msg('preview-input',
                       '%preview {}'.format(' '.join(output_files)))
                if hasattr(self, 'in_sandbox') and self.in_sandbox:
                    # if in sand box, do not link output to their files because these
                    # files will be removed soon.
                    self.send_frontend_msg('display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html':
                                HTML('''<div class="sos_hint"> input: {}<br>output: {}\n</div>'''.format(
                                ', '.join(x for x in input_files),
                                ', '.join(x for x in output_files))).data
                            }
                        })
                else:
                    self.send_frontend_msg('display_data',
                        {
                            'source': 'SoS',
                            'metadata': {},
                            'data': { 'text/html':
                                HTML('''<div class="sos_hint"> input: {}<br>output: {}\n</div>'''.format(
                                ', '.join('<a target="_blank" href="{0}">{0}</a>'.format(x) for x in input_files),
                                ', '.join('<a target="_blank" href="{0}">{0}</a>'.format(x) for x in output_files))).data
                            }
                        })
                for filename in output_files:
                    self.preview_file(filename, style=None)

    def preview_var(self, item, style=None):
        if item in env.sos_dict:
            obj = env.sos_dict[item]
        else:
            obj = SoS_eval(item, sigil=get_default_global_sigil())
        # get the basic information of object
        txt = type(obj).__name__
        # we could potentially check the shape of data frame and matrix
        # but then we will need to import the numpy and pandas libraries
        if hasattr(obj, 'shape'):
            txt += ' of shape {}'.format(getattr(obj, 'shape'))
        elif isinstance(obj, Sized):
            txt += ' of length {}'.format(obj.__len__())
        if callable(obj) or isinstance(obj, ModuleType):
            return txt, ({'text/plain': pydoc.render_doc(obj, title='SoS Documentation: %s')}, {})
        elif hasattr(obj, 'to_html'):
            try:
                from .visualize import Visualizer
                result = Visualizer(self, style).preview(obj)
                if isinstance(result, (list, tuple)) and len(result) == 2:
                    return txt, result
                elif isinstance(result, dict):
                    return txt, (result, {})
                elif result is None:
                    return txt, None
                else:
                    raise ValueError('Unrecognized return value from visualizer: {}.'.format(short_repr(result)))
            except Exception as e:
                self.warn(e)
                return txt, self.format_obj(obj)
        else:
            return txt, self.format_obj(obj)

    def preview_file(self, filename, style=None):
        if not os.path.isfile(filename):
            self.warn('\n> ' + filename + ' does not exist')
            return
        self.send_frontend_msg('display_data',
             {'metadata': {},
             'data': {
                 'text/plain': '\n> {} ({}):'.format(filename, pretty_size(os.path.getsize(filename))),
                 'text/html': HTML('<div class="sos_hint">> {} ({}):</div>'.format(filename, pretty_size(os.path.getsize(filename)))).data,
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
                            self.send_frontend_msg('stream', {
                                'name': 'stderr',
                                'text': 'Failed to load previewer {}: {}'.format(y, e) })
                            continue
                        break
                except Exception as e:
                    self.send_frontend_msg('stream', {
                                'name': 'stderr',
                                'text': str(e)})
                    continue
        #
        # if no previewer can be found
        if previewer_func is None:
            return
        try:
            result = previewer_func(filename, self, style)
            if not result:
                return
            if isinstance(result, str):
                self.send_frontend_msg('stream',
                    {'name': 'stdout', 'text': result})
            elif isinstance(result, dict):
                self.send_frontend_msg('display_data',
                    {'source': filename, 'data': result, 'metadata': {}})
            elif isinstance(result, [list, tuple]) and len(result) == 2:
                self.send_frontend_msg('display_data',
                    {'source': filename, 'data': result[0], 'metadata': result[1]})
            else:
                self.send_frontend_msg('stream', {
                    'name': 'stderr',
                    'text': 'Unrecognized preview content: {}'.format(result)})
                return
        except Exception as e:
            self.send_frontend_msg('stream', {
                'name': 'stderr',
                'text': 'Failed to preview {}: {}'.format(filename, e)})

    def render_result(self, res):
        if self._render_result is False:
            return res
        if not isinstance(res, str):
            self.warn('Cannot render result {} in type {} as {}.'.format(short_repr(res),
                res.__class__.__name__, self._render_result))
        else:
            # import the object from IPython.display
            mod = __import__('IPython.display')
            if not hasattr(mod.display, self._render_result):
                self.warn('Unrecognized render format {}'.format(self._render_result))
            else:
                func = getattr(mod.display, self._render_result)
                res = func(res)
        return res

    def send_result(self, res, silent=False):
        # this is Ok, send result back
        if not silent and res is not None:
            format_dict, md_dict = self.format_obj(self.render_result(res))
            self.send_response(self.iopub_socket, 'execute_result',
                {'execution_count': self._execution_count, 'data': format_dict,
                'metadata': md_dict})

    def get_kernel_list(self, notebook_kernel_list=None):
        if not hasattr(self, '_kernel_list'):
            from jupyter_client.kernelspec import KernelSpecManager
            km = KernelSpecManager()
            specs = km.find_kernel_specs()
            # get supported languages
            self._kernel_list = []
            lan_map = {self.supported_languages[x].kernel_name:(x, self.supported_languages[x].background_color) for x in self.supported_languages.keys()
                    if x != self.supported_languages[x].kernel_name}
            for spec in specs.keys():
                if spec == 'sos':
                    # the SoS kernel will be default theme color.
                    self._kernel_list.append(['SoS', 'sos', '', ''])
                elif spec in lan_map:
                    # e.g. ir ==> R
                    self._kernel_list.append([lan_map[spec][0], spec, lan_map[spec][0], lan_map[spec][1]])
                else:
                    # undefined language also use default theme color
                    self._kernel_list.append([spec, spec, '', ''])
        # now, using a list of kernels sent from the kernel, we might need to adjust
        # our list or create new kernels.
        if notebook_kernel_list is not None:
            for [name, kernel, lan, color] in notebook_kernel_list:
                try:
                    # if we can find the kernel, fine...
                    self.find_kernel(name, kernel, lan, color, notify_frontend=False)
                except Exception as e:
                    # otherwise do not worry about it.
                    env.logger.warning('Failed to locate subkernel {} with kernerl "{}" and language "{}": {}'.format(
                        name, kernel, lan, e))
        return self._kernel_list

    def do_execute(self, code, silent, store_history=True, user_expressions=None,
                   allow_stdin=False):
        if self._debug_mode:
            self.warn(code)
        # a flag for if the kernel is hard switched (by %use)
        self.hard_switch_kernel = False
        # evaluate user expression
        try:
            ret = self._do_execute(code=code, silent=silent, store_history=store_history,
                user_expressions=user_expressions, allow_stdin=allow_stdin)
        except Exception as e:
            return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }

        if ret is None:
            ret = {'status': 'ok',
                   'payload': [], 'user_expressions': {},
                   'execution_count': self._execution_count}

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
        if not silent and store_history:
            self._real_execution_count += 1
        self._execution_count = self._real_execution_count
        # make sure post_executed is triggered after the completion of all cell content
        self.shell.user_ns.update(env.sos_dict._dict)
        # trigger post processing of object and display matplotlib figures
        self.shell.events.trigger('post_execute')
        # tell the frontend the kernel for the "next" cell
        if store_history:
            self.send_frontend_msg('default-kernel', self.kernel)
        return ret

    def remove_leading_comments(self, code):
        lines = code.splitlines()
        try:
            idx = [x.startswith('#') or not x.strip() for x in lines].index(False)
            return '\n'.join(lines[idx:])
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
        if self.MAGIC_SKIP.match(code):
            self.warn('The %skip magic is deprecated and will be removed later.')
            return {'status': 'ok', 'payload': [], 'user_expressions': {}, 'execution_count': self._execution_count}
        elif self.MAGIC_RENDER.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_render_parser()
            try:
                args = parser.parse_args(shlex.split(options))
            except SystemExit:
                return
            try:
                self._render_result = args.format
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self._render_result = False
        elif self.MAGIC_SESSIONINFO.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_sessioninfo_parser()
            try:
                args = parser.parse_known_args(shlex.split(options))
            except SystemExit:
                return
            self.handle_sessioninfo()
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_TOC.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.send_frontend_msg('show_toc')
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_DICT.match(code):
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
        elif self.MAGIC_SHUTDOWN.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_shutdown_parser()
            try:
                args = parser.parse_args(shlex.split(options))
            except SystemExit:
                return
            self.shutdown_kernel(args.kernel if args.kernel else self.kernel, args.restart)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_FRONTEND.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_frontend_parser()
                try:
                    args = parser.parse_args(shlex.split(options))
                except SystemExit:
                    return
                self.cell_idx = args.cell_idx
                # for panel cell, we return a non-informative execution count
                if self.cell_idx is None or self.cell_idx < 0:
                    self._execution_count = '-'
                self._notebook_name = args.filename if args.filename else 'Untitled'
                if args.workflow is not None:
                    self._workflow = '#!/usr/bin/env sos-runner\n#fileformat=SOS1.0\n\n' + \
                        base64.b64decode(args.workflow).decode()
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            self._use_panel = args.use_panel is True
            if args.list_kernel:
                # https://github.com/jupyter/help/issues/153#issuecomment-289026056
                #
                # when the frontend is refreshed, cached comm would be lost and
                # communication would be discontinued. However, a kernel-list
                # request would be sent by the new-connection so we reset the
                # frontend_comm to re-connect to the frontend.
                self.comm_manager.register_target('sos_comm', self.sos_comm)

            # args.default_kernel should be valid
            if self.find_kernel(args.default_kernel)[0] != self.find_kernel(self.kernel)[0]:
                self.switch_kernel(args.default_kernel)
            #
            if args.cell_kernel == 'undefined':
                args.cell_kernel = args.default_kernel
            #
            original_kernel = self.kernel
            if self.find_kernel(args.cell_kernel)[0] != self.find_kernel(self.kernel)[0]:
                self.switch_kernel(args.cell_kernel)
            try:
                if args.resume:
                    self._resume_execution = True
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self._resume_execution = False
                if not self.hard_switch_kernel:
                    self.switch_kernel(original_kernel)
        elif self.MAGIC_WITH.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_with_parser()
                try:
                    args = parser.parse_args(shlex.split(options))
                except SystemExit:
                    return
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            original_kernel = self.kernel
            if args.restart and args.name in self.kernels:
                self.shutdown_kernel(args.name)
            try:
                self.switch_kernel(args.name, args.in_vars, args.out_vars,
                    args.kernel, args.language, args.color)
            except Exception as e:
                self.warn('Failed to switch to subkernel {} (kernel {}, language {}): {}'.format(args.name,
                    args.kernel, args.language, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.switch_kernel(original_kernel)
        elif self.MAGIC_USE.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_use_parser()
                try:
                    args = parser.parse_args(shlex.split(options))
                except SystemExit:
                    return
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'abort',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            if args.restart and args.name in self.kernels:
                self.shutdown_kernel(args.name)
                self.warn('{} is shutdown'.format(args.name))
            try:
                self.switch_kernel(args.name, args.in_vars, args.out_vars,
                    args.kernel, args.language, args.color)
                self.hard_switch_kernel = True
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            except Exception as e:
                self.warn('Failed to switch to subkernel {} (kernel {}, language {}): {}'.format(args.name,
                    args.kernel, args.language, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
        elif self.MAGIC_GET.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_get_parser()
                try:
                    args = parser.parse_args(options.split())
                except SystemExit:
                    return
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            self.handle_magic_get(args.vars, args.__from__, explicit=True)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_PUT.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_put_parser()
                try:
                    args = parser.parse_args(options.split())
                except SystemExit:
                    return
            except Exception as e:
                self.warn('Invalid option "{}": {}\n'.format(options, e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                   }
            self.handle_magic_put(args.vars, args.__to__, explicit=True)
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
                    env.logger.warn('Failed to get text from the clipboard: {}'.format(e))
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
                self._workflow_mode = True
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self._workflow_mode = False
                self.options = old_options
        elif self.MAGIC_SOSRUN.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            options = self._interpolate_option(options, quiet=False)
            old_options = self.options
            self.options = options + ' ' + self.options
            try:
                self._workflow_mode = True
                #self.send_frontend_msg('preview-workflow', self._workflow)
                self._do_execute(self._workflow, silent, store_history, user_expressions, allow_stdin)
            finally:
                self._workflow_mode = False
                self.options = old_options
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_SAVE.match(code):
            if self.kernel != 'SoS':
                # pass the %save magic to underlying kernel
                return self._do_execute(code, silent, store_history, user_expressions, allow_stdin)
            # if sos kernel ...
            options, remaining_code = self.get_magic_and_code(code, False)
            # parse options
            # save
            # execute the rest
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_SOSSAVE.match(code):
            # get the saved filename
            options, remaining_code = self.get_magic_and_code(code, False)
            try:
                parser = self.get_sossave_parser()
                try:
                    args = parser.parse_args(shlex.split(options))
                except SystemExit:
                    return
                if args.filename:
                    filename = self._interpolate_option(args.filename, quiet=False).strip()
                    if filename.lower().endswith('.html'):
                        if args.__to__ is None:
                            ftype = 'html'
                        else:
                            self.warn('%sossave to an .html file in {} format'.format(args.__to__))
                    else:
                        ftype = 'sos'
                else:
                    ftype = args.__to__ if args.__to__ else 'sos'
                    filename = self._notebook_name + '.' + ftype

                if os.path.isfile(filename) and not args.force:
                    raise ValueError('Cannot overwrite existing output file {}'.format(filename))
                #self.send_frontend_msg('preview-workflow', self._workflow)
                if ftype == 'sos':
                    if not args.all:
                        with open(filename, 'w') as script:
                            script.write(self._workflow)
                    else:
                        # convert to sos report
                        from sos.jupyter.converter import notebook_to_script
                        arg = argparse.Namespace()
                        arg.all = True
                        notebook_to_script(self._notebook_name + '.ipynb', filename, sargs=arg, unknown_args=[])
                    if args.setx:
                        import stat
                        os.chmod(filename, os.stat(filename).st_mode | stat.S_IEXEC)
                else:
                    # convert to sos report
                    from sos.jupyter.converter import notebook_to_html
                    arg = argparse.Namespace()
                    arg.template = args.template
                    notebook_to_html(self._notebook_name + '.ipynb', filename, sargs=arg, unknown_args=[])

                self.send_response(self.iopub_socket, 'display_data',
                    {'source': 'SoS', 'metadata': {},
                     'data': {
                         'text/plain': 'Workflow saved to {}\n'.format(filename),
                         'text/html': HTML('Workflow saved to <a href="{0}" target="_blank">{0}</a>'.format(filename)).data
                          }
                     })
                #
                if args.commit:
                    self.handle_shell_command(['git', 'commit', filename, '-m',
                        args.message if args.message else 'save {}'.format(filename)])
                if args.push:
                    self.handle_shell_command(['git', 'push'])
            except Exception as e:
                self.warn('Failed to save workflow: {}'.format(e))
                return {'status': 'error',
                    'ename': e.__class__.__name__,
                    'evalue': str(e),
                    'traceback': [],
                    'execution_count': self._execution_count,
                }
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_RERUN.match(code):
            options, remaining_code = self.get_magic_and_code(code, True)
            old_options = self.options
            self.options = options + ' ' + self.options
            try:
                self._workflow_mode = True
                if not hasattr(self, 'last_executed_code') or not self.last_executed_code:
                    self.warn('No saved script')
                    self.last_executed_code = ''
                return self._do_execute(self.last_executed_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                self.options = old_options
                self._workflow_mode = False
        elif self.MAGIC_SANDBOX.match(code):
            import tempfile
            import shutil
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_sandbox_parser()
            try:
                args = parser.parse_args(shlex.split(options))
            except SystemExit:
                return
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
            options = self._interpolate_option(options, quiet=True)
            parser = self.get_preview_parser()
            options = shlex.split(options, posix=False)
            help_option = []
            if ('-s' in options or '--style' in options) and '-h' in options:
                # defer -h to subparser
                options.remove('-h')
                help_option = ['-h']
            try:
                args, style_options = parser.parse_known_args(options)
            except SystemExit:
                return
            #
            style_options.extend(help_option)
            style = {'style': args.style, 'options': style_options }
            #
            if args.off:
                self.preview_output = False
            else:
                self.preview_output = True
            #
            if args.panel:
                self._use_panel = True
            elif args.notebook:
                self._use_panel = False
            # else, use default _use_panel
            try:
                return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
            finally:
                # preview workflow
                if args.workflow:
                    self.send_frontend_msg('stream',
                        {'name': 'stdout', 'text': self._workflow})
                if not args.off and args.items:
                    self.handle_magic_preview(args.items, args.kernel, style)
        elif self.MAGIC_CD.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_magic_cd(options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_DEBUG.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_debug_parser()
            try:
                args = parser.parse_args(options.split())
            except SystemExit:
                return
            self._debug_mode = args.status == 'on'
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_TASKINFO.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_taskinfo_parser()
            try:
                args = parser.parse_args(options.split())
            except SystemExit:
                return
            self.handle_taskinfo(args.task, args.queue)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.MAGIC_TASKS.match(code):
            options, remaining_code = self.get_magic_and_code(code, False)
            parser = self.get_tasks_parser()
            try:
                args = parser.parse_args(options.split())
            except SystemExit:
                return
            self.handle_tasks(args.tasks, args.queue if args.queue else 'localhost', args.status, args.age)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif code.startswith('!'):
            options, remaining_code = self.get_magic_and_code(code, False)
            self.handle_shell_command(code.split(' ')[0][1:] + ' ' + options)
            return self._do_execute(remaining_code, silent, store_history, user_expressions, allow_stdin)
        elif self.kernel != 'SoS':
            # handle string interpolation before sending to the underlying kernel
            if code:
                self.last_executed_code = code
            code = self._interpolate_option(code, quiet=False)
            if self.cell_idx is not None:
                self.send_frontend_msg('cell-kernel', [self.cell_idx, self.kernel])
                self.cell_idx = None
            if code is None:
                return
            try:
                return self.run_cell(code, silent, store_history)
            except KeyboardInterrupt:
                self.warn('Keyboard Interrupt\n')
                return {'status': 'abort', 'execution_count': self._execution_count}
        else:
            if code:
                self.last_executed_code = code
            # run sos
            try:
                self.run_sos_code(code, silent)
                if self.cell_idx is not None:
                    self.send_frontend_msg('cell-kernel', [self.cell_idx, 'SoS'])
                    self.cell_idx = None
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

    def __del__(self):
        # upon releasing of sos kernel, kill all subkernels. This I thought would be
        # called by the Jupyter cleanup code or the OS (because subkernels are subprocesses)
        # but they are not.
        self.do_shutdown(False)

if __name__ == '__main__':
    from ipykernel.kernelapp import IPKernelApp
    IPKernelApp.launch_instance(kernel_class=SoS_Kernel)
