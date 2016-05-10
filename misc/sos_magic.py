#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)

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

# This file must have been copied to ~/.ipython/extensions/ during
# the installation of SoS (python setup.py install), with a sos
# profile with pre-loaded sos magics. You can use
#
#    ipython --profile sos
#
# to use sos_magic in a ipython session with sos profile, or run
#
#    %load_ext sos_magic
#
# after you load ipython, or you can add 'sos_magic' to
#
#     c.InteractiveShellApp.extensions
#
# in ~/.ipython/profile_default/ipython_config.py, something like:
#
# c.InteractiveShellApp.extensions = [
#    'autoreload',
#    'sos_magic'
# ]
#

from pysos.utils import env, WorkflowDict
from pysos.sos_eval import SoS_exec, SoS_eval
from pysos.sos_executor import Interactive_Executor

from IPython.core.error import TryNext
from IPython.lib.clipboard import ClipboardEmpty
from IPython.core.magic import Magics, magics_class, line_magic, cell_magic, line_cell_magic

# The class MUST call this class decorator at creation time
@magics_class
class SoS_Magics(Magics):
    '''Magics that works with Script of Scripts'''
    def __init__(self, shell):
        super(SoS_Magics, self).__init__(shell)
        self._reset()

    def _reset(self):
        env.sos_dict = WorkflowDict()
        SoS_exec('from pysos import *')
        env.sos_dict.set('__interactive__', True)
        self.executor = Interactive_Executor()
        self.original_keys = set(env.sos_dict._dict.keys())
        self.original_keys.add('__builtins__')
        self.options = ''

    @line_cell_magic
    def sos(self, line, cell=None):
        'Magic execute sos expression and statements'
        # if in line mode, no command line
        if cell is None:
            try:
                # is it an expression?
                compile(line, '<string>', 'eval')
                return SoS_eval(line)
            except: # if it is satement
                return SoS_exec(line)
        else:
            try:
                # is it an expression?
                compile(cell, '<string>', 'eval')
                #if line.strip():
                #    env.logger.warning('{} ignored for expression evaluation'.format(line))
                return SoS_eval(cell)
            except:
                # is it a list of statement?
                try:
                    compile(cell, '<string>', 'exec')
                    #if line.strip():
                    #    env.logger.warning('{} ignored for statement execution'.format(line))
                    return SoS_exec(cell)
                except:
                    return self.executor.run_interactive(cell, command_line=self.options + line.strip())

    @line_magic
    def sospaste(self, line):
        'Magic that execute sos expression and statements from clipboard'
        # get and print clipboard content
        try:
            block = self.shell.hooks.clipboard_get()
        except TryNext as clipboard_exc:
            message = getattr(clipboard_exc, 'args')
            if message:
                error(message[0])
            else:
                error('Could not get text from the clipboard.')
            return
        except ClipboardEmpty:
            raise UsageError("The clipboard appears to be empty")
        #
        print(block.strip())
        print('## -- End pasted text --')
        try:
            # is it an expression?
            compile(block, '<string>', 'eval')
            return SoS_eval(block)
        except:
            # is it a list of statement?
            try:
                compile(block, '<string>', 'exec')
                return SoS_exec(block)
            except:
                return self.executor.run_interactive(block, command_line=self.options + line.strip())

    @line_magic
    def sosset(self, line):
        'Magic that set perminant options for sos and sospaste'
        # do not return __builtins__ beacuse it is too long...
        if line.strip():
            print('sos options set to "{}"'.format(line.strip()))
            self.options = line.strip() + ' '
        else:
            if self.options:
                print('sos options "{}" is reset to ""').format(self.options.strip())
                self.options = ''
            else:
                print('Usage: set persistent sos options such as -v3 -i (inspect) -p (prepare) -t (transcribe)')

    @line_magic
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

def load_ipython_extension(ipython):
    ipython.register_magics(SoS_Magics)
