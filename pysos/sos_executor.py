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

from collections.abc import Sequence

from . import __version__
from .sos_step import Step_Executor
from .utils import env, Error, WorkflowDict,  get_traceback, ProgressBar, frozendict, natural_keys
from .sos_eval import Undetermined, SoS_eval, SoS_exec
from .sos_script import SoS_Script, SoS_Step, SoS_ScriptContent

__all__ = []

def __null_func__(*args, **kwargs):
    '''This is a utility function for the parser'''
    return args, kwargs

class ArgumentError(Error):
    """Raised when an invalid argument is passed."""
    def __init__(self, msg):
        Error.__init__(self, msg)
        self.args = (msg, )

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
    def __init__(self, workflow, report, transcript):
        self.workflow = workflow
        self.report = report
        self.transcript = transcript

    def _parse_error(self, msg):
        '''This function will replace error() function in argparse module so that SoS
        can hijack errors raised from it.'''
        raise ArgumentError(msg)

    def parse_args(self, section, args, check_unused=False, cmd_name=''):
        '''Parse command line arguments and set values to parameters section'''
        env.logger.debug('Execute ``{}_parameters``'.format(section.name))
        env.sos_dict.set('step_name', '{}_parameters'.format(section.name))
        if section.global_def:
            try:
                SoS_exec(section.global_def)
            except Exception as e:
                if env.verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(section.global_def, e))

        def str2bool(v):
            if v.lower() in ('yes', 'true', 't', '1'):
                return True
            elif v.lower() in ('no', 'false', 'f', '0'):
                return False
            else:
                raise ArgumentError('Invalid value for bool argument "{}" (only yes,no,true,false,t,f,0,1 are allowed)'.format(v))
        #
        parser = argparse.ArgumentParser(prog='sos-runner {}'.format(cmd_name))
        parser.register('type', 'bool', str2bool)
        arguments = {}
        for key, defvalue, comment in section.parameters:
            try:
                defvalue = SoS_eval(defvalue, section.sigil)
                arguments[key] = defvalue
            except Exception as e:
                raise RuntimeError('Incorrect default value {} for parameter {}: {}'.format(defvalue, key, e))
            if isinstance(defvalue, type):
                if defvalue == bool:
                    parser.add_argument('--{}'.format(key), type='bool', help=comment, required=True, nargs='?')
                else:
                    # if only a type is specified, it is a required document of required type
                    parser.add_argument('--{}'.format(key), type=str if hasattr(defvalue, '__iter__') else defvalue,
                        help=comment, required=True, nargs='+' if hasattr(defvalue, '__iter__') else '?')
            else:
                if isinstance(defvalue, bool):
                    parser.add_argument('--{}'.format(key), type='bool', help=comment,
                        nargs='?', default=defvalue)
                else:
                    if isinstance(defvalue, str):
                        deftype = str
                    elif isinstance(defvalue, Sequence):
                        if len(defvalue) > 0:
                            deftype = type(defvalue[0])
                        else:
                            deftype = str
                    else:
                        deftype = type(defvalue)
                    parser.add_argument('--{}'.format(key), type=deftype, help=comment,
                        nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, str) else '?',
                        default=defvalue)
        #
        parser.error = self._parse_error
        #
        # because of the complexity of having combined and nested workflows, we cannot know how
        # many parameters section a workflow has and therfore have to assume that the unknown parameters
        # are for other sections.
        if check_unused:
            parsed = parser.parse_args(args)
        else:
            parsed, unknown = parser.parse_known_args(args)
            if unknown:
                env.logger.warning('Unparsed arguments [{}] that might be processed by another combined or nested workflow'
                    .format(' '.join(unknown)))
        #
        arguments.update(vars(parsed))
        # now change the value with passed values
        for k, v in arguments.items():
            env.sos_dict[k] = v
            # protect variables from being further modified
            env.readonly_vars.add(k)

    def load_config(self, config_file=None):
        '''load global, local and user-specified config files'''
        cfg = {}
        sos_config_file = os.path.expanduser('~/.sos/config.yaml')
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        #
        # local config file
        sos_config_file = '.sos/config.yaml'
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse local sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        #
        if config_file is not None:
            if not os.path.isfile(config_file):
                raise RuntimeError('Config file {} not found'.format(config_file))
            try:
                with open(config_file) as config:
                    cfg.update(yaml.safe_load(config))
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
            if not (hasattr(env, 'sos_dict') and '__interactive__' in env.sos_dict and env.sos_dict['__interactive__']):
                env.sos_dict = WorkflowDict()
            if env.run_mode != 'run':
                env.sos_dict['__execute_errors__'] = ExecuteError(self.workflow.name)
            # inject a few things
            env.sos_dict.set('__null_func__', __null_func__)
            env.sos_dict.set('__args__', args)
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
        if os.path.isdir('.sos/report'):
            shutil.rmtree('.sos/report')
        os.makedirs('.sos/report')
        env.sos_dict.set('__transcript__', None)

    def finalize(self):
        # collect reports and write to a file
        if env.run_mode != 'run':
            return
        step_reports = glob.glob('.sos/report/*')
        step_reports.sort(key=natural_keys)
        # merge the files
        if step_reports and self.report:
            with open(self.report, 'w') as combined:
                for step_report in step_reports:
                    with open(step_report, 'r') as md:
                        combined.write(md.read())
            env.logger.info('Report saved to {}'.format(self.report))

    def run(self, args=[], nested=False, cmd_name='', config_file=None,
        run_mode='run', sig_mode='default', verbosity=2):
        env.verbosity = verbosity
        if not self.workflow.sections:
            env.logger.trace('Skip because no section is defined')
            return
        if not nested or run_mode == 'inspect':
            env.run_mode = 'inspect'
            try:
                self.setup(args, nested, cmd_name, config_file)
                if run_mode == 'inspect':
                    env.sos_dict.set('__transcript__', self.transcript)
                self.execute(args, nested, cmd_name, config_file)
            except Exception:
                if verbosity and verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise
        if run_mode in ['prepare', 'run'] and (not nested or run_mode == 'prepare'):
            env.run_mode = 'prepare'
            env.sig_mode = sig_mode
            try:
                self.setup(args, nested, cmd_name, config_file)
                if run_mode == 'prepare':
                    env.sos_dict.set('__transcript__', self.transcript)
                self.execute(args, nested, cmd_name, config_file)
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
                self.execute(args, nested, cmd_name, config_file)
            except Exception:
                if verbosity and verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise
        self.finalize()

class Sequential_Executor(Base_Executor):
    #
    # A SoS workflow with multiple steps
    #
    def __init__(self, workflow, report=None, transcript=None):
        Base_Executor.__init__(self, workflow, report, transcript)

    def execute(self, args=[], nested=False, cmd_name='', config_file=None):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        #
        # process step of the pipelinp
        #
        num_parameters_sections = len([x for x in self.workflow.sections if x.is_parameters])
        if num_parameters_sections == 0 and args:
            raise RuntimeError('Unused parameter {}'.format(' '.join(args)))
        #
        # the steps can be executed in the pool (Not implemented)
        # if nested = true, start a new progress bar
        prog = ProgressBar(self.workflow.name, len(self.workflow.sections), disp=env.verbosity == 1 and env.run_mode == 'run')
        for idx, section in enumerate(self.workflow.sections):
            # global section will not change _step etc
            if section.is_parameters:
                # if there is only one parameters section and no nested workflow, check unused section
                self.parse_args(section, args, num_parameters_sections == 1, cmd_name=cmd_name)
                prog.progress(1)
                continue
            # handle skip, which might have to be evaluated till now.
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
            queue = mp.Queue()
            proc = mp.Process(target=Step_Executor(section).run_with_queue,
                args=(queue,))
            proc.start()
            proc.join()
            res = queue.get()
            # if the job is failed
            if isinstance(res, Exception):
                # error must have been displayed.
                #if env.verbosity > 2 and hasattr(res, 'traces'):
                #    env.logger.error(res.traces)
                raise RuntimeError(res)
            #res = section.run()
            for k, v in res.items():
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
    '''Interactive executor called from by iPython'''
    def __init__(self):
        # we actually do not have our own workflow, everything is passed from ipython
        Base_Executor.__init__(self, None, None, None)

    def parse_command_line(self, command_line):
        parser = argparse.ArgumentParser()
        # ignored
        parser.add_argument('-j', type=int, metavar='JOBS', default=1, dest='__max_jobs__')
        # ignored
        parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE')
        # ignored
        parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE', default='.sos/__interactive_report__')
        parser.add_argument('-t', dest='__transcript__', nargs='?', const='__STDERR__',
            metavar='TRANSCRIPT')
        runmode = parser.add_argument_group(title='Run mode options')
        runmode.add_argument('-i', action='store_true', dest='__inspect__')
        runmode.add_argument('-p', action='store_true', dest='__prepare__')
        runmode.add_argument('-f', action='store_true', dest='__rerun__')
        runmode.add_argument('-F', action='store_true', dest='__construct__')
        parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2)
        #
        args, workflow_args = parser.parse_known_args(shlex.split(command_line))
        return args, workflow_args

    def parse_script(self, code):
        '''Used by the kernel to judge if the code is complete'''
        return SoS_Script(content=code)

    def run_interactive(self, block, command_line=''):
        '''Execute a block of SoS script that is sent by iPython.'''
        # first, we try to parse it
        #
        if not block.strip():
            return
        # will exit if there is parsing error
        try:
            script = SoS_Script(content=block)
        except Exception as e:
            # if a step starts with output: etc, try to add a section head
            if 'not allowed outside' in repr(e):
                script = SoS_Script(content = '[0]\n' + block)
            else:
                raise
        #
        if command_line and not command_line.startswith('-'):
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
            if not script.sections:
                # if there is no section, but some actions that are seen as part of a
                # global definition
                if script.global_def:
                    env.sig_mode = 'default'
                    env.run_mode = 'run'
                    if args.__rerun__:
                        env.sig_mode = 'ignore'
                    if args.__prepare__:
                        env.run_mode = 'prepare'
                    if args.__inspect__:
                        env.run_mode = 'inspect'
                    if args.__construct__:
                        env.sig_mode = 'construct'
                    for key in ['__step_input__', '__step_output__']:
                        env.sos_dict.pop(key, None)
                    #
                    self.load_config(args.__config__)
                    env.sos_dict.set('step_name', '__interactive__')
                    env.sos_dict.set('__transcript__', args.__transcript__)
                    env.sos_dict.set('__step_report__', args.__report__)
                    # some actions requires this
                    env.sos_dict.set('_index', 0)
                    try:
                        return SoS_exec(script.global_def)
                    except Exception as e:
                        env.logger.debug('Failed to execute as global_def: {}'.format(e))
                        # add an empty section to execute it as a script
                        # FIXME: not sure if this a good arrangement...
                        section = SoS_Step(SoS_ScriptContent('interactive'), [('default', None)])
                        section.global_def = script.global_def
                        script.sections.append(section)
                else:
                    return None
            #
            sig_mode = 'default'
            run_mode = 'run'
            if args.__rerun__:
                sig_mode = 'ignore'
            if args.__prepare__:
                run_mode = 'prepare'
            if args.__inspect__:
                run_mode = 'inspect'
            if args.__construct__:
                sig_mode = 'construct'
            #
            workflow = script.workflow(wf_name)
            if args.__report__:
                executor = Sequential_Executor(workflow, report=args.__report__, transcript=args.__transcript__)
            else:
                executor = Sequential_Executor(workflow, report='.sos/ipython.md', transcript=args.__transcript__)
            executor.run(workflow_args, cmd_name='<script> {}'.format(wf_name), config_file=args.__config__,
                run_mode=run_mode, sig_mode=sig_mode, verbosity=args.verbosity)
        finally:
            env.verbosity = 2
            env.run_mode = 'run'
            env.sig_mode = 'default'

