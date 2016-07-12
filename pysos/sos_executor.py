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
from .sos_step import Step_Executor
from .utils import env, Error, WorkflowDict,  get_traceback, ProgressBar, \
    frozendict, natural_keys, dict_merge, ArgumentError
from .sos_eval import Undetermined, SoS_eval, SoS_exec
from .sos_script import SoS_Script, SoS_Step, SoS_ScriptContent
from .sos_syntax import SOS_SECTION_HEADER

__all__ = []

def __null_func__(*args, **kwargs):
    '''This is a utility function for the parser'''
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
    def __init__(self, workflow, report, transcript, debug):
        self.workflow = workflow
        self.report = report
        self.transcript = transcript
        self.debug = debug

    def load_config(self, config_file=None):
        '''load global, local and user-specified config files'''
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

    def setup(self, args=[], nested=False, cmd_name='', config_file=None):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        if nested:
            # if this is a subworkflow, we use _input as step input the workflow
            env.sos_dict.set('__step_output__', env.sos_dict['_input'])
        else:
            # Because this workflow might belong to a combined workflow, we do not clear
            # locals before the execution of workflow.
            # Need to choose what to inject to globals
            if env.run_mode != 'interactive':
                env.sos_dict = WorkflowDict()
            if env.run_mode != 'run':
                env.sos_dict['__execute_errors__'] = ExecuteError(self.workflow.name)
            # inject a few things
            env.sos_dict.set('__null_func__', __null_func__)
            env.sos_dict.set('__args__', args)
            env.sos_dict.set('__unknown_args__', args)
            # initial values
            env.sos_dict.set('SOS_VERSION', __version__)
            env.sos_dict.set('SOS_SCRIPT', self.workflow.sections[0].context.filename)
            # passing run_mode to SoS dict so that users can execute blocks of python statements
            # in different run modes.
            env.sos_dict.set('run_mode', env.run_mode)
            # load global, local and user-specified config file
            self.load_config(config_file)
            # there is no default input for the first step...
            env.sos_dict.set('__step_output__', [])
            SoS_exec('import os, sys, glob')
            SoS_exec('from pysos import *')
        #
        if os.path.isdir(os.path.join('.sos', 'report')):
            shutil.rmtree(os.path.join('.sos', 'report'))
        os.makedirs(os.path.join('.sos', 'report'))
        env.sos_dict.set('__transcript__', None)

    def finalize(self):
        # collect reports and write to a file
        if env.run_mode != 'run':
            return
        step_reports = glob.glob(os.path.join('.sos', 'report', '*'))
        step_reports.sort(key=natural_keys)
        # merge the files
        if step_reports and self.report:
            if self.report == '__STDOUT__':
                combined = sys.stdout
            else:
                combined = open(self.report, 'w')
            for step_report in step_reports:
                with open(step_report, 'r') as md:
                    combined.write(md.read())
            if self.report != '__STDOUT__':
                combined.close()
                env.logger.info('Report saved to {}'.format(self.report))

class Sequential_Executor(Base_Executor):
    #
    # A SoS workflow with multiple steps
    #
    def __init__(self, workflow, report=None, transcript=None, debug=False):
        Base_Executor.__init__(self, workflow, report, transcript, debug)


    def run(self, args=[], nested=False, cmd_name='', config_file=None,
        run_mode='run', sig_mode='default', verbosity=2):
        env.verbosity = verbosity
        if not self.workflow.sections:
            env.logger.trace('Skip because no section is defined')
            return
        DAG = {}
        if not nested or run_mode == 'inspect':
            env.run_mode = 'inspect'
            try:
                self.setup(args, nested, cmd_name, config_file)
                if run_mode == 'inspect':
                    env.sos_dict.set('__transcript__', self.transcript)
                self.execute(args, nested, cmd_name, config_file, DAG=DAG)
            except Exception:
                if verbosity and verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise
            if '__unknown_args__' in env.sos_dict and env.sos_dict['__unknown_args__']:
                raise ArgumentError('Unhandled command line argument {}'.format(' '.join(env.sos_dict['__unknown_args__'])))
        if run_mode in ['prepare', 'run'] and (not nested or run_mode == 'prepare'):
            env.run_mode = 'prepare'
            env.sig_mode = sig_mode
            try:
                self.setup(args, nested, cmd_name, config_file)
                if run_mode == 'prepare':
                    env.sos_dict.set('__transcript__', self.transcript)
                self.execute(args, nested, cmd_name, config_file, DAG=DAG)
            except Exception:
                if verbosity and verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise
        if run_mode == 'run' and (not nested or run_mode == 'run'):
            env.run_mode = 'run'
            env.sig_mode = sig_mode
            try:
                self.setup(args, nested, cmd_name, config_file)
                if run_mode == 'run':
                    env.sos_dict.set('__transcript__', self.transcript)
                self.execute(args, nested, cmd_name, config_file, DAG=DAG)
            except Exception:
                if verbosity and verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise
        self.finalize()

    def execute(self, args=[], nested=False, cmd_name='', config_file=None, DAG={}):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        #
        # process step of the pipelinp
        #
        # the steps can be executed in the pool (Not implemented)
        # if nested = true, start a new progress bar
        prog = ProgressBar(self.workflow.name, len(self.workflow.sections),
            disp=len(self.workflow.sections) > 1 and env.verbosity == 1 and env.run_mode == 'run')
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
            # 1. for first step of workflow, _step.input=[]
            # 2. for subworkflow, _step.input = _input
            # 3. for second to later step, _step.input = _step.output
            # each section can use a separate process
            if self.debug or env.run_mode == 'interactive':
                res = Step_Executor(section).run(DAG)
            else:
                queue = mp.Queue()
                proc = mp.Process(target=Step_Executor(section).run_with_queue,
                    args=(queue, DAG))
                proc.start()
                res = queue.get()
                proc.join()
            # if the job is failed
            if isinstance(res, Exception):
                # error must have been displayed.
                #if env.verbosity > 2 and hasattr(res, 'traces'):
                #    env.logger.error(res.traces)
                raise RuntimeError(res)
            #res = section.run()
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
        # ignored
        parser.add_argument('-j', type=int, metavar='JOBS', default=1, dest='__max_jobs__')
        # ignored
        parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE')
        # ignored
        parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE', 
            default=os.path.join('.sos', '__step_report.md'))
        parser.add_argument('-t', dest='__transcript__', nargs='?', const='__STDERR__',
            metavar='TRANSCRIPT')
        runmode = parser.add_argument_group(title='Run mode options')
        runmode.add_argument('-i', action='store_true', dest='__inspect__')
        runmode.add_argument('-p', action='store_true', dest='__prepare__')
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
        '''Execute a block of SoS script that is sent by iPython.'''
        # first, we try to parse it
        #
        if not block.strip():
            for key in ['__step_input__', '__step_output__']:
                env.sos_dict.pop(key, None)
            return
        # first check if block contains section header
        # if there is no section header, add a header so that the block
        # appears to be a SoS script with one section
        if not any([SOS_SECTION_HEADER.match(line) for line in block.split()]):
            block = '[0]\n' + block

        script = SoS_Script(content=block)
        #
        if command_line.strip() and not command_line.startswith('-'):
            wf_and_args = command_line.split(None, 1)
            wf_name = wf_and_args[0]
            command_line = command_line[len(wf_name)+1:]
        else:
            wf_name = None
        #
        try:
            args, workflow_args = self.parse_command_line(command_line)
            if os.path.isfile(args.__report__):
                os.remove(args.__report__)
            # if there is only a global section
            #
            sig_mode = 'default'
            run_mode = 'interactive'
            #if args.__rerun__:
            #    sig_mode = 'ignore'
            #if args.__prepare__:
            #    run_mode = 'prepare'
            #if args.__inspect__:
            #    run_mode = 'inspect'
            if args.__construct__:
                sig_mode = 'construct'
            #
            self.workflow = script.workflow(wf_name)
            return self.execute_workflow(args=[], nested=False, cmd_name='', config_file=None)
        finally:
            env.verbosity = 2
            env.run_mode = 'interactive'
            env.sig_mode = 'default'

    def execute_workflow(self, args=[], nested=False, cmd_name='', config_file=None):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        #
        # process step of the pipelinp
        self.setup(args, nested, cmd_name, config_file)
        last_res = None
        #
        # the steps can be executed in the pool (Not implemented)
        # if nested = true, start a new progress bar
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
                        val_sk =usection.options['skip'].value(section.sigil)
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
            DAG = {}
            last_res = Step_Executor(section).run_interactive()
            # if the job is failed
            if isinstance(last_res, Exception):
                # error must have been displayed.
                #if env.verbosity > 2 and hasattr(res, 'traces'):
                #    env.logger.error(res.traces)
                raise RuntimeError(last_res)
        # at the end
        if not nested and env.run_mode != 'run':
            exception = env.sos_dict['__execute_errors__']
            if exception.errors:
                # if there is any error, raise it
                raise exception
        self.finalize()
        return last_res



