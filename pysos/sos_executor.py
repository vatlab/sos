#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
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
import yaml
import glob
import shlex
import shutil
import argparse

import multiprocessing as mp

from . import __version__
from .sos_step import Inspect_Step_Executor, Prepare_Step_Executor, Run_Step_Executor, Interactive_Step_Executor
from .utils import env, Error, WorkflowDict,  get_traceback, ProgressBar, \
    frozendict, natural_keys, dict_merge, ArgumentError
from .sos_eval import Undetermined, SoS_eval, SoS_exec
from .sos_script import SoS_Script, SoS_Step, SoS_ScriptContent
from .sos_syntax import SOS_SECTION_HEADER

__all__ = []


def __null_func__(*args, **kwargs):
    '''This function will be passed to SoS's namespace and be executed
    to evaluate functions of input, output, and depends directives.'''
    return args, kwargs

class ExecuteError(Error):
    """Raised when there are errors in inspect mode. Such errors are not raised
    immediately, but will be collected and raised at the end """

    def __init__(self, workflow):
        Error.__init__(self, 'SoS workflow contains errors: %s' % workflow)
        self.workflow = workflow
        self.errors = []
        self.traces = []
        self.args = (workflow, )

    def append(self, line, error):
        lines = [x for x in line.split('\n') if x.strip()]
        if not lines:
            short_line = '<empty>'
        else:
            short_line = lines[0][:40] if len(lines[0]) > 40 else lines[0]
        self.errors.append(short_line)
        self.traces.append(get_traceback())
        if isinstance(error, Exception):
            self.message += '\n[%s] %s:\n\t%s' % (short_line, error.__class__.__name__, error)
        else:
            self.message += '\n[%s]:\n\t%s' % (short_line, error)

class Base_Executor:
    '''This is the base class of all executor that provides common
    set up and tear functions for all executors.'''
    def __init__(self, workflow, args=[], config_file=None, new_dict=True):
        self.workflow = workflow

        # if creating a new dictionary, set it up with some basic varibles
        # and functions
        if new_dict:
            env.sos_dict = WorkflowDict()
            env.sos_dict['__execute_errors__'] = ExecuteError(self.workflow.name)
            # inject a few things
            env.sos_dict.set('__null_func__', __null_func__)
            env.sos_dict.set('__args__', args)
            env.sos_dict.set('__unknown_args__', args)
            # initial values
            env.sos_dict.set('SOS_VERSION', __version__)
            env.sos_dict.set('SOS_SCRIPT', self.workflow.sections[0].context.filename)
            env.sos_dict.set('__step_output__', [])
            SoS_exec('import os, sys, glob')
            SoS_exec('from pysos import *')

        # load configuration files
        cfg = {}
        sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yaml')
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        #
        # local config file
        sos_config_file = os.path.join('.sos', 'config.yaml')
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse local sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        #
        if config_file is not None:
            if not os.path.isfile(config_file):
                raise RuntimeError('Config file {} not found'.format(config_file))
            try:
                with open(config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format? ({})'.format(config_file, e))
        #
        env.sos_dict.set('CONFIG', frozendict(cfg))

    def inspect(self, nested=False):
        '''Run the script in inspect mode to check for errors.'''
        env.run_mode = 'inspect'
        # passing run_mode to SoS dict so that users can execute blocks of
        # python statements in different run modes.
        env.sos_dict.set('run_mode', env.run_mode)

        # process steps of the pipeline
        for idx, section in enumerate(self.workflow.sections):
            # handle skip, which might have to be evaluated till now.
            #
            # the global section has to be executed here because step options might need
            # information from it. Also, the variables in the global section should be
            # global. In addition, the global section has to be executed multiple times
            # because sections can come from different scripts (nested workflows).
            if section.global_def:
                try:
                    SoS_exec(section.global_def)
                except Exception as e:
                    if env.verbosity > 2:
                        sys.stderr.write(get_traceback())
                    raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                        section.global_def, e))
            #
            # Here we require that skip to be evaluatable at inspect and prepare mode
            # up until this step. There is to say, it cannot rely on the result of
            # an action that is only available till run time (action would return
            # Undetermined when executed not in the specified runmode)
            if 'skip' in section.options:
                if isinstance(section.options['skip'], Undetermined):
                    try:
                        val_skip = section.options['skip'].value(section.sigil)
                        if val_skip is None:
                            val_skip = False
                    except Exception as e:
                        raise RuntimeError('Failed to evaluate value of section option skip={}: {}'.format(section.options['skip'], e))
                else:
                    val_skip = section.options['skip']
                if val_skip is None or val_skip is True:
                    continue
                elif val_skip is not False:
                    raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
            #
            # execute section with specified input
            queue = mp.Queue()
            executor = Inspect_Step_Executor(section, queue)
            proc = mp.Process(target=executor.run)
            proc.start()
            res = queue.get()
            proc.join()
            # if the job is failed
            if isinstance(res, Exception):
                raise RuntimeError(res)
            #
            for k, v in res.items():
                env.sos_dict.set(k, v)
        # at the end
        if not nested:
            exception = env.sos_dict['__execute_errors__']
            if exception.errors:
                # if there is any error, raise it
                raise exception

    def prepare(self, nested=False):
        '''Run the script in prepare mode to prepare resources.'''
        env.run_mode = 'prepare'
        # passing run_mode to SoS dict so that users can execute blocks of
        # python statements in different run modes.
        env.sos_dict.set('run_mode', env.run_mode)
        
        # process step of the pipelinp
        #
        # the steps can be executed in the pool (Not implemented)
        # if nested = true, start a new progress bar
        DAG = {}
        for idx, section in enumerate(self.workflow.sections):
            # handle skip, which might have to be evaluated till now.
            #
            # the global section has to be executed here because step options might need
            # infomration from it. Also, the variables in the global section should be
            # global. In addition, the global section has to be executed multiple times
            # because sections can come from different scripts (nested workflows).
            if section.global_def:
                try:
                    SoS_exec(section.global_def)
                except Exception as e:
                    if env.verbosity > 2:
                        sys.stderr.write(get_traceback())
                    raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                        section.global_def, e))
            #
            # Important:
            #
            # Here we require that skip to be evaluatable at inspect and prepare mode
            # up until this step. There is to say, it cannot rely on the result of
            # an action that is only available till run time (action would return
            # Undetermined when executed not in the specified runmode)
            #
            if 'skip' in section.options:
                if isinstance(section.options['skip'], Undetermined):
                    try:
                        val_skip = section.options['skip'].value(section.sigil)
                        if val_skip is None:
                            val_skip = False
                    except Exception as e:
                        raise RuntimeError('Failed to evaluate value of section option skip={}: {}'.format(section.options['skip'], e))
                else:
                    val_skip = section.options['skip']
                if val_skip is None or val_skip is True:
                    continue
                elif val_skip is not False:
                    raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
            #
            # execute section with specified input
            queue = mp.Queue()
            executor = Prepare_Step_Executor(section, queue, DAG)
            proc = mp.Process(target=executor.run)
            proc.start()
            res = queue.get()
            proc.join()
            # if the job is failed
            if isinstance(res, Exception):
                raise RuntimeError(res)
            #
            for k, v in res.items():
                if k == '__dag__':
                    DAG.update(v)
                else:
                    env.sos_dict.set(k, v)
            prog.progress(1)
        prog.done()
        # at the end
        if not nested and env.run_mode != 'run':
            exception = env.sos_dict['__execute_errors__']
            if exception.errors:
                # if there is any error, raise it
                raise exception


class Sequential_Executor(Base_Executor):
    #
    # Execute a workflow sequentially in batch mode
    def __init__(self, workflow, args=[], config_file=None):
        Base_Executor.__init__(self, workflow, args, config_file, new_dict=True)

    def run(self, DAG=None, nested=False):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        # process step of the pipelinp
        #
        # the steps can be executed in the pool (Not implemented)
        # if nested = true, start a new progress bar
        prog = ProgressBar(self.workflow.name, len(self.workflow.sections),
            disp=len(self.workflow.sections) > 1 and env.verbosity == 1)
        for idx, section in enumerate(self.workflow.sections):
            # handle skip, which might have to be evaluated till now.
            #
            # the global section has to be executed here because step options might need
            # infomration from it. Also, the variables in the global section should be
            # global. In addition, the global section has to be executed multiple times
            # because sections can come from different scripts (nested workflows).
            if section.global_def:
                try:
                    SoS_exec(section.global_def)
                except Exception as e:
                    if env.verbosity > 2:
                        sys.stderr.write(get_traceback())
                    raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                        section.global_def, e))
            #
            if 'skip' in section.options:
                if isinstance(section.options['skip'], Undetermined):
                    try:
                        val_skip = section.options['skip'].value(section.sigil)
                        if val_skip is None:
                            val_skip = False
                    except Exception as e:
                        raise RuntimeError('Failed to evaluate value of section option skip={}: {}'.format(section.options['skip'], e))
                else:
                    val_skip = section.options['skip']
                if val_skip is None or val_skip is True:
                    continue
                elif val_skip is not False:
                    raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
            #
            # execute section with specified input
            # 1. for first step of workflow, _step.input=[]
            # 2. for subworkflow, _step.input = _input
            # 3. for second to later step, _step.input = _step.output
            # each section can use a separate process
            queue = mp.Queue()
            executor = Run_Step_Executor(section, queue, DAG)
            proc = mp.Process(target=executor.run)
            proc.start()
            res = queue.get()
            proc.join()
            # if the job is failed
            if isinstance(res, Exception):
                raise RuntimeError(res)
            #
            for k, v in res.items():
                if k == '__dag__':
                    DAG.update(v)
                else:
                    env.sos_dict.set(k, v)
            prog.progress(1)
        prog.done()
        # at the end
        if not nested and env.run_mode != 'run':
            exception = env.sos_dict['__execute_errors__']
            if exception.errors:
                # if there is any error, raise it
                raise exception


class Interactive_Executor(Base_Executor):
    '''Interactive executor called from by iPython Jupyter or Spyder'''
    def __init__(self):
        # we actually do not have our own workflow, everything is passed from ipython
        Base_Executor.__init__(self, None, None, None, False)

    def parse_command_line(self, command_line):
        parser = argparse.ArgumentParser()
        # no default workflow so it will execute any workflow if the code piece
        # defines only one workflow
        parser.add_argument('workflow', metavar='WORKFLOW', nargs='?')
        parser.add_argument('-j', type=int, metavar='JOBS', default=1, dest='__max_jobs__')
        parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE')
        parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE', 
            default=os.path.join('.sos', '__step_report.md'))
        parser.add_argument('-t', dest='__transcript__', nargs='?', const='__STDERR__',
            metavar='TRANSCRIPT')
        runmode = parser.add_argument_group(title='Run mode options')
        runmode.add_argument('-f', action='store_true', dest='__rerun__')
        runmode.add_argument('-F', action='store_true', dest='__construct__')
        # default to 1 to avoid output env.logger.info to notebook
        parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=1)
        #
        args, workflow_args = parser.parse_known_args(shlex.split(command_line))
        return args, workflow_args

    def parse_script(self, code):
        '''Used by the kernel to judge if the code is complete'''
        return SoS_Script(content=code)

    def run(self, block, command_line=''):
        '''Execute a block of SoS script that is sent by iPython/Jupyer/Spyer
        The code can be simple SoS/Python statements, one SoS step, or more
        or more SoS workflows with multiple steps. This executor,
        1. adds a section header to the script if there is no section head
        2. execute the workflow in interactive mode, which is different from
           batch mode in a number of ways, which most notably without support
           for nested workflow.
        3. Optionally execute the workflow in inspection or preparation mode
           for debugging purposes.
        '''
        # if there is no valid code do nothing
        if not block.strip():
            return
        # if there is no section header, add a header so that the block
        # appears to be a SoS script with one section
        if not any([SOS_SECTION_HEADER.match(line) for line in block.split()]):
            block = '[interactive_0]\n' + block

        script = SoS_Script(content=block)
        env.run_mode = 'interactive'
        try:
            args, workflow_args = self.parse_command_line(command_line)
            self.workflow = script.workflow(args.workflow)

            if args.__rerun__:
                sig_mode = 'ignore'
            elif args.__construct__:
                sig_mode = 'construct'
            else:
                sig_mode = 'default'

            if os.path.isfile(args.__report__):
                os.remove(args.__report__)

            return self.execute_workflow(args=[], config_file=args.__config__)
        finally:
            env.verbosity = 1
            env.sig_mode = 'default'

    def execute_workflow(self, args=[], config_file=None):
        '''Execute a workflow with specified command line args. '''
        #
        # process step of the pipelinp
        self.setup(args, False, config_file)
        # this is the result returned by the workflow, if the
        # last stement is an expression.
        last_res = None
        #
        for idx, section in enumerate(self.workflow.sections):
            if 'skip' in section.options:
                val_skip = section.options['skip']
                if val_skip is None or val_skip is True:
                    continue
                elif val_skip is not False:
                    raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
            #
            last_res = Step_Executor(section).run_interactive()
            # if the step is failed
            if isinstance(last_res, Exception):
                raise RuntimeError(last_res)
        return last_res



