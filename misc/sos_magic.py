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

# cp sos_magic.py to ~/.ipython/extensions/
#
# Run
#    %load_ext sos_magic
#
# after you load ipython, or you can add c.InteractiveShellApp.extensions
# in ~/.ipython/profile_default/ipython_config.py, something like:
#
# c.InteractiveShellApp.extensions = [
#    'autoreload',
#    'sos_magic'
# ]
#  


from pysos.utils import env, WorkflowDict
from pysos.sos_eval import SoS_exec, SoS_eval

env.sos_dict = WorkflowDict()
SoS_exec('from pysos import *')


from IPython.core.error import TryNext
from IPython.lib.clipboard import ClipboardEmpty
from IPython.core.magic import Magics, magics_class, line_magic, cell_magic, line_cell_magic

# The class MUST call this class decorator at creation time
@magics_class
class SoS_Magic(Magics):

    @line_cell_magic
    def sos(self, line, cell=None):
        'Magic execute sos expression and statements'
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
                compile(line, '<string>', 'eval')
                return SoS_eval(line)
            except: # if it is satement
                return SoS_exec(line)


    @line_cell_magic
    def sospaste(self, line, cell=None):
        'Magic that execute sos expression and statements from clipboard'
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
        # block
        try:
            # is it an expression?
            compile(block, '<string>', 'eval')
            return SoS_eval(block)
        except: # if it is satement
            return SoS_exec(block)

    @line_cell_magic
    def sosdict(self, line, cell=None):
        'Magic that displays content of the dictionary'
        return env.sos_dict._dict

def load_ipython_extension(ipython):
    ipython.register_magics(SoS_Magic)
