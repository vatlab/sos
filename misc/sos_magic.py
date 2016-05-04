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

import os
import shlex
import argparse
from pysos.utils import env, WorkflowDict
from pysos.sos_eval import SoS_exec, SoS_eval
from pysos.sos_script import SoS_Script
from pysos.sos_executor import Sequential_Executor

env.sos_dict = WorkflowDict()
SoS_exec('from pysos import *')
env.sos_dict.set('__interactive__', True)


from IPython.core.error import TryNext
from IPython.lib.clipboard import ClipboardEmpty
from IPython.core.magic import Magics, magics_class, line_magic, cell_magic, line_cell_magic

def parse_args(command_line):
    parser = argparse.ArgumentParser()
    # ignored
    parser.add_argument('-j', type=int, metavar='JOBS', default=1, dest='__max_jobs__')
    # ignored
    parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE')
    # ignored
    parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE')
    runmode = parser.add_argument_group(title='Run mode options')
    runmode.add_argument('-d', action='store_true', dest='__dryrun__')
    runmode.add_argument('-p', action='store_true', dest='__prepare__')
    runmode.add_argument('-f', action='store_true', dest='__rerun__')
    runmode.add_argument('-F', action='store_true', dest='__construct__')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2)
    #
    args, workflow_args = parser.parse_known_args(shlex.split(command_line))
    return args, workflow_args

def sos_run_block(block, command_line='', pasted=False):
    '''Execute a block of SoS script that is sent by iPython.'''
    # first, we try to parse it
    #
    if not block.strip():
        return
    if pasted:
        print(block)
        print('\n## -- End pasted text --')
    # will exit if there is parsing error
    script = SoS_Script(content=block)
    # if there is only a global section
    if not script.sections:
        # if there is nothing, return
        env.logger.warning('No SoS step is defined. Execute as statements.')
        return SoS_exec(block)
    #
    if command_line.strip() and not command_line.strip().startswith('-'):
        wf_and_args = command_line.strip().split(None, 1)
        wf_name = wf_and_args[0]
        if len(wf_and_args) > 1:
            command_line = wf_and_args[1]
    else:
        wf_name = None
    #
    args, workflow_args = parse_args(command_line)
    # special execution mode
    if args.__dryrun__:
        env.run_mode = 'dryrun'
    elif args.__prepare__:
        env.run_mode = 'prepare'
    else:
        env.run_mode = 'run'
    if args.__rerun__:
        env.sig_mode = 'ignore'
    elif args.__construct__:
        env.sig_mode = 'construct'
    #
    workflow = script.workflow(wf_name)
    if args.__report__:
        executor = Sequential_Executor(workflow, report=args.__report__)
    else:
        executor = Sequential_Executor(workflow, report='.sos/ipython.md')
    executor.run(workflow_args, cmd_name='<script> {}'.format(wf_name), config_file=args.__config__)

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
            if line.strip():
                env.logger.warning('line used as part of cell')
            try:
                block = line + '\n' + cell
                # is it an expression?
                compile(block, '<string>', 'eval')
                return SoS_eval(block)
            except: # if it is satement
                return SoS_exec(block)

    def get_clipboard(self):
        try:
            return self.shell.hooks.clipboard_get()
        except TryNext as clipboard_exc:
            message = getattr(clipboard_exc, 'args')
            if message:
                error(message[0])
            else:
                error('Could not get text from the clipboard.')
            return
        except ClipboardEmpty:
            raise UsageError("The clipboard appears to be empty")


    @line_cell_magic
    def sospaste(self, line, cell=None):
        'Magic that execute sos expression and statements from clipboard'
        block = self.get_clipboard()
        if line.strip():
            env.logger.warning('line {} ignored'.format(line))
        if cell:
            env.logger.warning('cell content {} ignored'.format(cell))
        try:
            # is it an expression?
            compile(block, '<string>', 'eval')
            return SoS_eval(block)
        except: # if it is satement
            return SoS_exec(block)

    @line_cell_magic
    def sosdict(self, line, cell=None):
        'Magic that displays content of the dictionary'
        # do not return __builtins__ beacuse it is too long...
        env.sos_dict._dict.pop('__builtins__', None)
        return env.sos_dict._dict

    @cell_magic
    def sosrun(self, line, cell):
        block = line + '\n' + cell
        env.run_mode = 'run'
        sos_run_block(block, command_line=line)

    @line_cell_magic
    def sosrunpaste(self, line, cell=None):
        block = self.get_clipboard()
        env.run_mode = 'run'
        sos_run_block(block, command_line=line, pasted=True)

    @line_cell_magic
    def sosreset(self, line, cell=None):
        if line.strip():
            env.logger.warning('line {} ignored'.format(line))
        if cell:
            env.logger.warning('cell content {} ignored'.format(cell))
        env.sos_dict = WorkflowDict()
        SoS_exec('from pysos import *')
        env.sos_dict.set('__interactive__', True)

def load_ipython_extension(ipython):
    ipython.register_magics(SoS_Magic)
