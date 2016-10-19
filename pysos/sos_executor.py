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
import shlex
import argparse
import time
import builtins
import keyword
from collections.abc import Sequence
import multiprocessing as mp
from queue import Empty

from io import StringIO
from ._version import __version__
from .sos_step import Dryrun_Step_Executor, Prepare_Step_Executor, SP_Step_Executor, MP_Step_Executor, RQ_Step_Executor, \
    Celery_Step_Executor, Interactive_Step_Executor, analyze_section
from .utils import env, Error, WorkflowDict, get_traceback, ProgressBar, frozendict, dict_merge, short_repr
from .sos_eval import SoS_exec
from .sos_script import SoS_Script
from .sos_syntax import SOS_SECTION_HEADER, SOS_KEYWORDS
from .dag import SoS_DAG
from .target import BaseTarget, FileTarget, UnknownTarget, RemovedTarget, UnavailableLock, sos_variable, textMD5
from .pattern import extract_pattern

__all__ = []


class ExecuteError(Error):
    """Raised when there are errors in prepare mode. Such errors are not raised
    immediately, but will be collected and raised at the end """

    def __init__(self, workflow):
        Error.__init__(self, 'Failed to execute workflow %s' % workflow)
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
        if short_line in self.errors:
            return
        self.errors.append(short_line)
        self.traces.append(get_traceback())
        if isinstance(error, Exception):
            self.message += '\n[%s] %s:\n\t%s' % (short_line, error.__class__.__name__, error)
        else:
            self.message += '\n[%s]:\n\t%s' % (short_line, error)

def __null_func__(*args, **kwargs):
    '''This function will be passed to SoS's namespace and be executed
    to evaluate functions of input, output, and depends directives.'''
    return args, kwargs

class Base_Executor:
    '''This is the base class of all executor that provides common
    set up and tear functions for all executors.'''
    def __init__(self, workflow=None, args=[], config_file=None, nested=False):
        env.__task_engine__ = None
        self.workflow = workflow
        self.args = args
        self.nested = nested
        self.config_file = config_file
        # interactive mode does not pass workflow
        if self.workflow:
            self.md5 = self.create_signature()

    def create_signature(self):
        with StringIO() as sig:
            sig.write('# Sections\n')
            for step in self.workflow.sections + self.workflow.auxiliary_sections:
                sig.write('{}: {}\n'.format(step.step_name(), step.md5))
            sig.write('# Command line options\n')
            sig.write('{}\n'.format(self.args))
            self.sig_content = sig.getvalue()
        return textMD5(self.sig_content)[:16]

    def reset_dict(self):
        # if creating a new dictionary, set it up with some basic varibles
        # and functions
        if self.nested:
            SoS_exec('import os, sys, glob')
            SoS_exec('from pysos.runtime import *')
            self._base_symbols = set(dir(__builtins__)) | set(env.sos_dict.keys()) | set(SOS_KEYWORDS) | set(keyword.kwlist)
            self._base_symbols -= {'dynamic'}
            return

        env.sos_dict = WorkflowDict()

        # inject a few things
        env.sos_dict.set('__null_func__', __null_func__)
        env.sos_dict.set('__args__', self.args)
        env.sos_dict.set('__unknown_args__', self.args)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        env.sos_dict.set('__step_output__', [])

        # load configuration files
        cfg = {}
        sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yaml')
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        # local config file
        sos_config_file = 'config.yaml'
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse local sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        # user-specified configuration file.
        if self.config_file is not None:
            if not os.path.isfile(self.config_file):
                raise RuntimeError('Config file {} not found'.format(self.config_file))
            try:
                with open(self.config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format? ({})'.format(self.config_file, e))
        # set config to CONFIG
        env.sos_dict.set('CONFIG', frozendict(cfg))

        SoS_exec('import os, sys, glob')
        SoS_exec('from pysos.runtime import *')
        self._base_symbols = set(dir(builtins)) | set(env.sos_dict.keys()) | set(SOS_KEYWORDS) | set(keyword.kwlist)
        self._base_symbols -= {'dynamic'}

    def skip(self, section):
        if section.global_def:
            try:
                SoS_exec(section.global_def)
            except RuntimeError as e:
                if env.verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                    section.global_def, e))
        #
        if 'skip' in section.options:
            val_skip = section.options['skip']
            if val_skip is None or val_skip is True:
                env.logger.info('Step ``{}`` is ``ignored`` due to skip option.'.format(section.step_name()))
                return True
            elif val_skip is not False:
                raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
        return False

    def match(self, target, patterns):
        if isinstance(patterns, (str, BaseTarget)):
            patterns = [patterns]
        elif not isinstance(patterns, Sequence):
            raise RuntimeError('Unknown target to match: {}'.format(patterns))
        #
        for p in patterns:
            # other targets has to match exactly
            if isinstance(target, BaseTarget) or isinstance(p, BaseTarget):
                if target == p:
                    return {}
                else:
                    continue
            # if this is a regular string
            res = extract_pattern(p, [target])
            if res and not any(None in x for x in res.values()):
                return res
            # string match
            elif FileTarget(p) == FileTarget(target):
                return True
        return False

    def resolve_dangling_targets(self, dag, targets=None):
        '''Feed dangling targets with their dependncies from auxiliary steps,
        optionally add other targets'''
        resolved = 0
        while True:
            dangling_targets = dag.dangling(targets)
            if not dangling_targets:
                # if no dangling targets, means all objects COULD be produved by DAG
                break
            env.logger.info('Resolving {} objects from {} nodes'.format(len(dangling_targets), dag.number_of_nodes()))
            # find matching steps
            # check auxiliary steps and see if any steps provides it
            for target in dangling_targets:
                # target might no longer be dangling after a section is added.
                if target not in dag.dangling(targets):
                    continue
                mo = [(x, self.match(target, x.options['provides'])) for x in self.workflow.auxiliary_sections]
                mo = [x for x in mo if x[1] is not False]
                if not mo:
                    for x in self.workflow.auxiliary_sections:
                        env.logger.debug('{}: {}'.format(x.step_name(), x.options['provides']))
                    raise RuntimeError('No step to generate target {}{}'.format(target, dag.steps_depending_on(target, self.workflow)))
                if len(mo) > 1:
                    raise RuntimeError('Multiple steps {} to generate target {}'.format(', '.join(x[0].step_name() for x in mo), target))
                #
                # only one step, we need to process it # execute section with specified input
                #
                # NOTE:  Auxiliary can be called with different output files and matching pattern
                # so we are actually creating a new section each time we need an auxillary step.
                #
                section = mo[0][0]
                if isinstance(mo[0][1], dict):
                    for k,v in mo[0][1].items():
                        env.sos_dict.set(k, v[0])
                #
                # for auxiliary, we need to set input and output, here
                # now, if the step does not provide any alternative (e.g. no variable generated
                # from patten), we should specify all output as output of step. Otherwise the
                # step will be created for multiple outputs. issue #243
                if mo[0][1]:
                    env.sos_dict['__default_output__'] = [target]
                elif isinstance(section.options['provides'], Sequence):
                    env.sos_dict['__default_output__'] = section.options['provides']
                else:
                    env.sos_dict['__default_output__'] = [section.options['provides']]
                # will become input, set to None
                env.sos_dict['__step_output__'] = None
                #
                res = analyze_section(section)
                #
                # build DAG with input and output files of step
                env.logger.info('Adding step {} with output {}'.format(res['step_name'], target))
                if isinstance(mo[0][1], dict):
                    context = mo[0][1]
                else:
                    context = {}
                context['__signature_vars__'] = res['signature_vars']
                context['__environ_vars__'] = res['environ_vars']
                context['__changed_vars__'] = res['changed_vars']
                context['__default_output__'] = env.sos_dict['__default_output__']
                # NOTE: If a step is called multiple times with different targets, it is much better
                # to use different names because pydotplus can be very slow in handling graphs with nodes
                # with identical names.
                dag.add_step(section.uuid, '{} {}'.format(section.step_name(),
                    short_repr(env.sos_dict['__default_output__'])), None, res['step_input'],
                    res['step_depends'], res['step_output'], context=context)
                resolved += 1
        return resolved


    def initialize_dag(self, targets=None):
        '''Create a DAG by analyzing sections statically.'''
        '''Run the script in prepare mode to prepare resources.'''
        # this is for testing only and allows tester to call initialize_dag
        # directly to get a DAG
        if not hasattr(self, '_base_symbols'):
            self.reset_dict()

        dag = SoS_DAG()
        default_input = []
        for idx, section in enumerate(self.workflow.sections):
            if self.skip(section):
                continue
            #
            res = analyze_section(section, default_input)

            environ_vars = res['environ_vars'] - self._base_symbols
            signature_vars = res['signature_vars'] - self._base_symbols
            changed_vars = res['changed_vars']
            # parameters, if used in the step, should be considered environmental
            environ_vars |= env.parameter_vars & signature_vars

            # add shared to targets
            if res['changed_vars']:
                if 'provides' in section.options:
                    if isinstance(section.options['provides'], str):
                        section.options.set('provides', [section.options['provides']])
                else:
                    section.options.set('provides', [])
                #
                section.options.set('provides',
                    section.options['provides'] + [sos_variable(var) for var in changed_vars])

            context={'__signature_vars__': signature_vars,
                    '__environ_vars__': environ_vars,
                    '__changed_vars__': changed_vars }

            # for nested workflow, the input is specified by sos_run, not None.
            if self.nested and idx == 0:
                context['__step_output__'] = env.sos_dict['__step_output__']
            # for regular workflow, the output of the last step has
            # to exist (existence of signature does not count)
            if not self.nested and idx + 1 == len(self.workflow.sections):
                context['__hard_target__'] = True

            # NOTE: if a section has option 'shared', the execution of this step would
            # change dictionary, essentially making all later steps rely on this step.
            dag.add_step(section.uuid,
                section.step_name(),
                idx,
                res['step_input'],
                res['step_depends'],
                res['step_output'],
                context = context)
            default_input = res['step_output']
        #
        # analyze auxiliary steps
        for idx, section in enumerate(self.workflow.auxiliary_sections):
            res = analyze_section(section, default_input)
            environ_vars = res['environ_vars'] - self._base_symbols
            signature_vars = res['signature_vars'] - self._base_symbols
            changed_vars = res['changed_vars']
            # parameters, if used in the step, should be considered environmental
            environ_vars |= env.parameter_vars & signature_vars

            # add shared to targets
            if res['changed_vars']:
                if 'provides' in section.options:
                    if isinstance(section.options['provides'], str):
                        section.options.set('provides', [section.options['provides']])
                else:
                    section.options.set('provides', [])
                #
                section.options.set('provides',
                    section.options['provides'] + [sos_variable(var) for var in changed_vars])
        #
        self.resolve_dangling_targets(dag, targets)
        # now, there should be no dangling targets, let us connect nodes
        dag.build(self.workflow.auxiliary_sections)
        # trim the DAG if targets are specified
        if targets:
            dag = dag.subgraph_from(targets)
        # write DAG for debugging purposes
        dag.write_dot(os.path.join(env.exec_dir, '.sos', '{}.dot'.format(self.workflow.name)))
        # check error
        cycle = dag.circular_dependencies()
        if cycle:
            raise RuntimeError('Circular dependency detected {}. It is likely a later step produces input of a previous step.'.format(cycle))

        return dag

    def save_workflow_signature(self, dag):
        '''Save tracked files in .sos so that untracked files can be cleaned by command
        sos clean.
        '''
        with open(os.path.join(env.exec_dir, '.sos', '{}.sig'.format(self.md5)), 'w') as sigfile:
            sigfile.write(self.sig_content)
            sigfile.write('# input and dependent files\n')
            for target in sorted(x for x in dag._all_dependent_files if isinstance(x, str)):
                sigfile.write('{}\n'.format(target))
            sigfile.write('# output files\n')
            for target in sorted(x for x in dag._all_output_files if isinstance(x, str)):
                sigfile.write('{}\n'.format(target))   
        
    def run(self, targets=None, mode='run'):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        self.reset_dict()
        env.run_mode = mode
        # passing run_mode to SoS dict so that users can execute blocks of
        # python statements in different run modes.
        env.sos_dict.set('run_mode', env.run_mode)
        # process step of the pipelinp
        dag = self.initialize_dag(targets=targets)
        #
        # if targets are specified and there are only signatures for them, we need
        # to remove the signature and really generate them
        if targets:
            for t in targets:
                if not FileTarget(t).exists('target'):
                    FileTarget(t).remove('signature')
        #
        SoS_exec('from pysos.runtime import sos_handle_parameter_')
        #
        prog = ProgressBar(self.workflow.name, dag.num_nodes(), disp=dag.num_nodes() > 1 and env.verbosity == 1)
        self.reset_dict()
        env.sos_dict.set('run_mode', env.run_mode)
        exec_error = ExecuteError(self.workflow.name)
        while True:
            # find any step that can be executed and run it, and update the DAT
            # with status.
            runnable = dag.find_executable()
            if runnable is None:
                # no runnable
                #dag.show_nodes()
                break
            # find the section from runnable
            section = self.workflow.section_by_id(runnable._step_uuid)
            #
            # this is to keep compatibility of dag run with sequential run because
            # in sequential run, we evaluate global section of each step in
            # order to determine values of options such as skip.
            # The consequence is that global definitions are available in
            # SoS namespace.
            try:
                SoS_exec(section.global_def)
            except Exception as e:
                if env.verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                    section.global_def, e))

            # clear existing keys, otherwise the results from some random result
            # might mess with the execution of another step that does not define input
            for k in ['__step_input__', '__default_output__', '__step_output__']:
                if k in env.sos_dict:
                    env.sos_dict.pop(k)
            # if the step has its own context
            env.sos_dict.quick_update(runnable._context)
            # execute section with specified input
            runnable._status = 'running'
            q = mp.Queue()
            if mode == 'run':
                executor = SP_Step_Executor(section, q)
            elif mode == 'prepare':
                executor = Prepare_Step_Executor(section, q)
            else:
                executor = Dryrun_Step_Executor(section, q)
            p = mp.Process(target=executor.run)
            p.start()
            #
            res = q.get()
            # if we does get the result
            p.join()
            # if the step says unknown target .... need to check if the target can
            # be build dynamically.
            if isinstance(res, UnknownTarget):
                runnable._status = None
                target = res.target
                if self.resolve_dangling_targets(dag, [target]) == 0:
                    raise RuntimeError('Failed to resolve {}{}.'
                        .format(target, dag.steps_depending_on(target, self.workflow)))
                # now, there should be no dangling targets, let us connect nodes
                # this can be done more efficiently
                dag.build(self.workflow.auxiliary_sections)
                cycle = dag.circular_dependencies()
                if cycle:
                    raise RuntimeError('Circular dependency detected {}. It is likely a later step produces input of a previous step.'.format(cycle))
            elif isinstance(res, RemovedTarget):
                runnable._status = None
                dag.regenerate_target(res.target)
            elif isinstance(res, UnavailableLock):
                runnable._status = 'pending'
                runnable._signature = (res.output, res.sig_file)
                env.logger.info('Waiting on another process for step {}'.format(section.step_name()))
            # if the job is failed
            elif isinstance(res, Exception):
                runnable._status = 'failed'
                exec_error.append(runnable._node_id, res)
                prog.progress(1)
            else:#
                for k, v in res.items():
                    env.sos_dict.set(k, v)
                #
                # set context to the next logic step.
                for edge in dag.out_edges(runnable):
                    node = edge[1]
                    # if node is the logical next step...
                    if node._node_index is not None and runnable._node_index is not None:
                        #and node._node_index == runnable._node_index + 1:
                        node._context.update(env.sos_dict.clone_selected_vars(
                            node._context['__signature_vars__'] | node._context['__environ_vars__'] \
                            | {'_input', '__step_output__', '__default_output__', '__args__'}))
                runnable._status = 'completed'
                prog.progress(1)
            #env.logger.error('completed')
        prog.done()
        if exec_error.errors:
            failed_steps, pending_steps = dag.pending()
            if failed_steps:
                sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in failed_steps]
                exec_error.append(self.workflow.name,
                    RuntimeError('{} failed step{}: {}'.format(len(sections), 
                        's' if len(sections) > 1 else '', ', '.join(sections))))
            if pending_steps:
                sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in pending_steps]
                exec_error.append(self.workflow.name,
                    RuntimeError('{} pending step{}: {}'.format(len(sections),
                        's' if len(sections) > 1 else '', ', '.join(sections))))
            raise exec_error
        else:
            self.save_workflow_signature(dag)
            env.logger.info('Workflow {} (ID={}) is executed successfully.'.format(self.workflow.name, self.md5))

    def dryrun(self, targets=None):
        '''Execute the script in dryrun mode.'''
        try:
            self.run(targets=targets, mode='dryrun')
        # only runtime errors are ignored
        except RuntimeError as e:
            env.logger.warning('Workflow cannot be completed in dryrun mode: {}'.format(e))

    def prepare(self, targets=None):
        '''Execute the script in prepare mode.'''
        try:
            self.run(targets=targets, mode='prepare')
        # only runtime errors are ignored
        except RuntimeError as e:
            env.logger.warning('Workflow cannot be completed in prepare mode: {}'.format(e))

class MP_Executor(Base_Executor):
    #
    # Execute a workflow sequentially in batch mode
    def __init__(self, workflow, args=[], config_file=None, nested=False):
        Base_Executor.__init__(self, workflow, args, config_file, nested=nested)
        if hasattr(env, 'accessed_vars'):
            delattr(env, 'accessed_vars')

    def step_executor(self, section, queue):
        return MP_Step_Executor(section, queue)

    def run(self, targets=None, mode='run'):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        self.reset_dict()
        env.run_mode = mode
        # passing run_mode to SoS dict so that users can execute blocks of
        # python statements in different run modes.
        env.sos_dict.set('run_mode', env.run_mode)
        # process step of the pipelinp
        dag = self.initialize_dag(targets=targets)
        #
        SoS_exec('from pysos.runtime import sos_handle_parameter_')

        # process step of the pipelinp
        #
        procs = [None for x in range(env.max_jobs)]
        prog = ProgressBar(self.workflow.name, dag.num_nodes(), disp=dag.num_nodes() > 1 and env.verbosity == 1)
        exec_error = ExecuteError(self.workflow.name)
        while True:
            # step 1: check existing jobs and see if they are completed
            for idx, proc in enumerate(procs):
                if proc is None:
                    continue
                (p, q, u) = proc
                try:
                    res = q.get_nowait()
                except Empty:
                    # continue waiting
                    continue
                #
                # if we does get the result
                p.join()

                runnable = dag.node_by_id(u)
                if isinstance(res, UnknownTarget):
                    runnable._status = None
                    target = res.target
                    if self.resolve_dangling_targets(dag, [target]) == 0:
                        raise RuntimeError('Failed to resolve {}{}.'
                            .format(target, dag.steps_depending_on(target, self.workflow)))
                    # now, there should be no dangling targets, let us connect nodes
                    # this can be done more efficiently
                    dag.build(self.workflow.auxiliary_sections)
                    cycle = dag.circular_dependencies()
                    if cycle:
                        raise RuntimeError('Circular dependency detected {}. It is likely a later step produces input of a previous step.'.format(cycle))
                elif isinstance(res, RemovedTarget):
                    runnable._status = None
                    dag.regenerate_target(res.target)
                elif isinstance(res, UnavailableLock):
                    runnable._status = 'pending'
                    runnable._signature = (res.output, res.sig_file)
                    section = self.workflow.section_by_id(runnable._step_uuid)
                    env.logger.info('Waiting on another process for step {}'.format(section.step_name()))
                    # move away to let other tasks to run first
                    procs[idx] = None
                # if the job is failed
                elif isinstance(res, Exception):
                    runnable._status = 'failed'
                    exec_error.append(runnable._node_id, res)
                    prog.progress(1)
                    procs[idx] = None
                else:
                    #
                    for k, v in res.items():
                        env.sos_dict.set(k, v)
                    #
                    # set context to the next logic step.
                    for edge in dag.out_edges(runnable):
                        node = edge[1]
                        # if node is the logical next step...
                        if node._node_index is not None and runnable._node_index is not None:
                            #and node._node_index == runnable._node_index + 1:
                            node._context.update(env.sos_dict.clone_selected_vars(
                                node._context['__signature_vars__'] | node._context['__environ_vars__'] \
                                | {'_input', '__step_output__', '__default_output__', '__args__'}))
                    runnable._status = 'completed'
                    prog.progress(1)
                    procs[idx] = None
                #env.logger.error('completed')
                #dag.show_nodes()
            # step 2: submit new jobs if there are empty slots
            for idx, proc in enumerate(procs):
                # if there is empty slot, submit
                if proc is not None:
                    continue
                # find any step that can be executed and run it, and update the DAT
                # with status.
                runnable = dag.find_executable()
                if runnable is None:
                    # no runnable
                    #dag.show_nodes()
                    break
                # find the section from runnable
                section = self.workflow.section_by_id(runnable._step_uuid)
                #
                # this is to keep compatibility of dag run with sequential run because
                # in sequential run, we evaluate global section of each step in
                # order to determine values of options such as skip.
                # The consequence is that global definitions are available in
                # SoS namespace.
                try:
                    SoS_exec(section.global_def)
                except RuntimeError as e:
                    if env.verbosity > 2:
                        sys.stderr.write(get_traceback())
                    raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                        section.global_def, e))

                # clear existing keys, otherwise the results from some random result
                # might mess with the execution of another step that does not define input
                for k in ['__step_input__', '__default_output__', '__step_output__']:
                    if k in env.sos_dict:
                        env.sos_dict.pop(k)
                # if the step has its own context
                env.sos_dict.quick_update(runnable._context)
                # execute section with specified input
                runnable._status = 'running'
                q = mp.Queue()
                executor = self.step_executor(section, q)
                p = mp.Process(target=executor.run)
                procs[idx] = (p, q, runnable._node_uuid)
                p.start()
                #
                #env.logger.error('started')
                #dag.show_nodes()
            #
            if all(x is None for x in procs):
                break
            else:
                time.sleep(0.1)
        prog.done()
        if exec_error.errors:
            failed_steps, pending_steps = dag.pending()
            if failed_steps:
                sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in failed_steps]
                exec_error.append(self.workflow.name,
                    RuntimeError('{} failed step{}: {}'.format(len(sections), 
                        's' if len(sections) > 1 else '', ', '.join(sections))))
            if pending_steps:
                sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in pending_steps]
                exec_error.append(self.workflow.name,
                    RuntimeError('{} pending step{}: {}'.format(len(sections), 
                        's' if len(sections) > 1 else '', ', '.join(sections))))
            raise exec_error
        else:
            sig = self.save_workflow_signature(dag)
            env.logger.info('Workflow {} (ID={}) is executed successfully.'.format(self.workflow.name, sig))

class RQ_Executor(MP_Executor):
    def __init__(self, workflow, args=[], config_file=None, nested=False):
        MP_Executor.__init__(self, workflow, args, config_file, nested=nested)
        env.__task_engine__ = 'rq'

        from rq import Queue as rqQueue
        from redis import Redis
        import yaml

        connection_file = os.path.join(env.exec_dir, '.sos', 'redis_connection.yaml')
        if os.path.isfile(connection_file):
            try:
                with open(connection_file) as conn:
                    cfg = yaml.safe_load(conn)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to parse redis connection file {}. ({}}'.format(connection_file, e))
                sys.exit(1)
            try:
                redis_conn = Redis(host=cfg['host'], port=cfg.get('port', 6379))
            except Exception as e:
                env.logger.error('Failed to connect to redis server: {}'.format(e))
                sys.exit(1)
        else:
            redis_conn = Redis()
        self.redis_queue = rqQueue(connection=redis_conn)

    def step_executor(self, section, queue):
        return RQ_Step_Executor(section, queue, self.redis_queue)

class Celery_Executor(MP_Executor):
    def __init__(self, workflow, args=[], config_file=None, nested=False):
        MP_Executor.__init__(self, workflow, args, config_file, nested=nested)
        env.__task_engine__ = 'Celery'

    def step_executor(self, section, queue):
        # pass celery_app if needed
        return Celery_Step_Executor(section, queue)

class Interactive_Executor(Base_Executor):
    '''Interactive executor called from by iPython Jupyter or Spyder'''
    def __init__(self):
        # we actually do not have our own workflow, everything is passed from ipython
        # by nested = True we actually mean no new dictionary
        Base_Executor.__init__(self, nested=True)

    def parse_command_line(self, command_line):
        parser = argparse.ArgumentParser()
        # no default workflow so it will execute any workflow if the code piece
        # defines only one workflow
        # 
        # parser.add_argument('-j', type=int, metavar='JOBS', default=1, dest='__max_jobs__')
        parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE')
        #parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE',
        #    default=os.path.join('.sos', '__step_report.md'))
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

    def set_dict(self, args):
        env.sos_dict.set('__null_func__', __null_func__)
        env.sos_dict.set('SOS_VERSION', __version__)

        # load configuration files
        cfg = {}
        sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yaml')
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        # local config file
        sos_config_file = 'config.yaml'
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse local sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        if args.__config__ is not None:
            # user-specified configuration file.
            if not os.path.isfile(args.__config__):
                raise RuntimeError('Config file {} not found'.format(args.__config__))
            try:
                with open(args.__config__) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format? ({})'.format(self.config_file, e))
        # set config to CONFIG
        env.sos_dict.set('CONFIG', frozendict(cfg))

    def run(self, block, command_line=''):
        '''Execute a block of SoS script that is sent by iPython/Jupyer/Spyer
        The code can be simple SoS/Python statements, one SoS step, or more
        or more SoS workflows with multiple steps. This executor,
        1. adds a section header to the script if there is no section head
        2. execute the workflow in interactive mode, which is different from
           batch mode in a number of ways, which most notably without support
           for nested workflow.
        3. Optionally execute the workflow in preparation mode for debugging purposes.
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
            env.sos_dict.set('__args__', workflow_args)
            env.sos_dict.set('__unknown_args__', workflow_args)
            self.set_dict(args)
            self.workflow = script.workflow()

            if args.__rerun__:
                env.sig_mode = 'ignore'
            elif args.__construct__:
                env.sig_mode = 'construct'
            else:
                env.sig_mode = 'default'

            #if os.path.isfile(args.__report__):
            #    os.remove(args.__report__)

            # this is the result returned by the workflow, if the
            # last stement is an expression.
            last_res = None
            #
            # clear __step_input__, __step_output__ etc because there is
            # no concept of passing input/outputs across cells.
            env.sos_dict.set('__step_output__', [])
            for k in ['__step_input__', '__default_output__', 'input', 'output', \
                'depends', '_input', '_output', '_depends']:
                env.sos_dict.pop(k, None)

            for idx, section in enumerate(self.workflow.sections):
                if 'skip' in section.options:
                    val_skip = section.options['skip']
                    if val_skip is None or val_skip is True:
                        continue
                    elif val_skip is not False:
                        raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
                #
                last_res = Interactive_Step_Executor(section).run()
                # if the step is failed
                if isinstance(last_res, Exception):
                    raise RuntimeError(last_res)
            return last_res
        finally:
            env.verbosity = 1
            env.sig_mode = 'default'



