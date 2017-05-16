#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
import time
import keyword
import uuid
from collections.abc import Sequence
import multiprocessing as mp

from tqdm import tqdm as ProgressBar
from io import StringIO
from ._version import __version__
from .sos_step import Step_Executor, analyze_section, PendingTasks
from .utils import env, Error, WorkflowDict, get_traceback, frozendict, short_repr, pickleable, \
    load_config_files, save_var, load_var
from .sos_eval import SoS_exec, get_default_global_sigil
from .dag import SoS_DAG
from .target import BaseTarget, FileTarget, UnknownTarget, RemovedTarget, UnavailableLock, sos_variable, textMD5, sos_step
from .pattern import extract_pattern
from .hosts import Host

__all__ = []

class ExecuteError(Error):
    """An exception to collect exceptions raised during run time so that
    other branches of the DAG would continue if some nodes fail to execute."""
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

class SoS_Worker(mp.Process):
    '''
    Worker process to process SoS step or workflow in separate process.
    '''
    def __init__(self,  pipe, config={}, args=[], **kwargs):
        '''
        cmd_queue: a single direction queue for the master process to push
            items to the worker.

        config:
            values for command line options

            config_file: -c
            output_dag: -d

        args:
            command line argument passed to workflow. However, if a dictionary is passed,
            then it is assumed to be a nested workflow where parameters are made
            immediately available.
        '''
        # the worker process knows configuration file, command line argument etc
        super(SoS_Worker, self).__init__(**kwargs)
        #
        self.pipe = pipe
        self.config = config
        self.args = args

    def reset_dict(self):
        env.sos_dict = WorkflowDict()
        env.parameter_vars.clear()
        env.config.update(self.config)

        env.sos_dict.set('__null_func__', __null_func__)
        env.sos_dict.set('__args__', self.args)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        env.sos_dict.set('__step_output__', [])

        # load configuration files
        cfg = load_config_files(self.config['config_file'])
        env.sos_dict.set('CONFIG', frozendict(cfg))

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        self._base_symbols = set(dir(__builtins__)) | set(env.sos_dict['sos_symbols_']) | set(keyword.kwlist)
        self._base_symbols -= {'dynamic'}

        if isinstance(self.args, dict):
            for key, value in self.args.items():
                if not key.startswith('__'):
                    env.sos_dict.set(key, value)

    def run(self):
        # wait to handle jobs
        while True:
            work = self.pipe.recv()
            if work is None:
                break

            env.logger.debug('Worker {} receives request {}'.format(self.name, work))
            try:
                if work[0] == 'step':
                    # this is a step ...
                    self.run_step(*work[1:])
                else:
                    self.run_workflow(*work[1:])
                env.logger.debug('Worker {} completes request {}'.format(self.name, work))
            except KeyboardInterrupt:
                break

    def run_workflow(self, workflow_id, wf, targets, args, shared, config):
        #
        # The pipe is the way to communicate with the master process.
        #
        # get workflow, args, shared, and config
        self.args = args
        self.reset_dict()
        # we are in a separate process and need to set verbosity from workflow config
        # but some tests do not provide verbosity
        env.verbosity = config.get('verbosity', 2)
        env.logger.debug('Worker {} working on a workflow {} with args {}'.format(self.name, workflow_id, args))
        executer = Base_Executor(wf, args=args, shared=shared, config=config)
        # we send the pipe to subworkflow, which would send
        # everything directly to the master process, so we do not
        # have to collect result here
        executer.run(targets=targets, parent_pipe=self.pipe, my_workflow_id=workflow_id)


    def run_step(self, section, context, shared, args, run_mode, sig_mode, verbosity):
        env.logger.debug('Worker {} working on a step with args {}'.format(self.name, args))
        env.config['run_mode'] = run_mode
        env.config['sig_mode'] = sig_mode
        env.verbosity = verbosity
        #
        self.args = args
        self.reset_dict()
        # this is to keep compatibility of dag run with sequential run because
        # in sequential run, we evaluate global section of each step in
        # order to determine values of options such as skip.
        # The consequence is that global definitions are available in
        # SoS namespace.
        try:
            SoS_exec(section.global_def, section.global_sigil)
        except RuntimeError as e:
            if env.verbosity > 2:
                sys.stderr.write(get_traceback())
            raise

        # clear existing keys, otherwise the results from some random result
        # might mess with the execution of another step that does not define input
        for k in ['__step_input__', '__default_output__', '__step_output__']:
            if k in env.sos_dict:
                env.sos_dict.pop(k)
        # if the step has its own context
        env.sos_dict.quick_update(shared)
        # context should be updated after shared because context would contain the
        # correct __step_output__ of the step, whereas shared might contain
        # __step_output__ from auxiliary steps. #526
        env.sos_dict.quick_update(context)

        executor = Step_Executor(section, self.pipe, mode=env.config['run_mode'])
        executor.run()

class dummy_node:
    # a dummy node object to store information of node passed
    # from nested workflow
    def __init__(self):
        pass

class Base_Executor:
    '''This is the base class of all executor that provides common
    set up and tear functions for all executors.'''
    def __init__(self, workflow=None, args=[], shared={}, config={}):
        self.workflow = workflow
        self.args = args
        self.shared = shared
        self.config = config
        env.config.update(config)
        for key in ('config_file', 'output_dag'):
            if key not in self.config:
                self.config[key] = None
        if self.config['config_file'] is not None:
            self.config['config_file'] = os.path.abspath(os.path.expanduser(self.config['config_file']))
        #
        # if the executor is not called from command line, without sigmode setting
        if env.config['sig_mode'] is None:
            env.config['sig_mode'] = 'default'
        # interactive mode does not pass workflow
        self.md5 = self.create_signature()
        if env.config['sig_mode'] != 'ignore' and self.workflow:
            # remove old workflow file.
            with open(os.path.join(env.exec_dir, '.sos', '{}.sig'.format(self.md5)), 'a') as sig:
                sig.write('# workflow: {}\n'.format(self.workflow.name))
                sig.write('# script: {}\n'.format(self.workflow.content.filename))
                sig.write('# included: {}\n'.format(','.join([x[1] for x in self.workflow.content.included])))
                sig.write('# configuration: {}\n'.format(config.get('config_file', '')))
                sig.write('# start time: {}\n'.format(time.strftime('%a, %d %b %Y %H:%M:%S +0000', time.gmtime())))
                sig.write(self.sig_content)
                sig.write('# runtime signatures\n')
        #
        env.config['resumed_tasks'] = set()
        wf_status = os.path.join(os.path.expanduser('~'), '.sos', self.md5 + '.status')
        if env.config['resume_mode']:
            if os.path.isfile(wf_status):
                with open(wf_status) as status:
                    for line in status:
                        if line.startswith('pending_task'):
                            k, v = load_var(line)
                            env.config['resumed_tasks'].add(v[1])
            else:
                env.logger.info('Workflow {} has been completed.'.format(self.md5))
                sys.exit(0)
        # wait is None or True, and there is task
        elif self.config.get('wait_for_task', None) is not True and self.workflow.has_external_task():
            with open(wf_status, 'w') as wf:
                # overwrite previous file
                for key, val in self.config.items():
                    wf.write(save_var(key, val))

        # if this is a resumed task?
        if hasattr(env, 'accessed_vars'):
            delattr(env, 'accessed_vars')

    def save_dag(self, dag):
        if self.config['output_dag'] is None:
            return
        if not hasattr(self, 'dag_count'):
            self.dag_count = 1
            self.last_dag = None
        #
        out = dag.to_string()
        if self.last_dag == out:
            return
        else:
            self.last_dag = out
        # output file name
        if self.config['output_dag'] == '-':
            sys.stdout.write(out)
        else:
            if self.dag_count == 1:
                dag_name = self.config['output_dag'] if self.config['output_dag'].endswith('.dot') else self.config['output_dag'] + '.dot'
            else:
                dag_name = '{}_{}.dot'.format(self.config['output_dag'][:-4] if self.config['output_dag'].endswith('.dot') else self.config['output_dag'], self.dag_count)
            #
            with open(dag_name, 'w') as dfile:
                dfile.write(out)

    def record_quit_status(self, tasks):
        if not self.md5:
            return
        with open(os.path.join(os.path.expanduser('~'), '.sos', self.md5 + '.status'), 'a') as status:
            for q, t in tasks:
                status.write(save_var('pending_task', [q, t]))

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
        env.sos_dict = WorkflowDict()
        env.parameter_vars.clear()
        env.config.update(self.config)

        # inject a few things
        if self.md5:
            env.sos_dict.set('__workflow_sig__', os.path.join(env.exec_dir, '.sos', '{}.sig'.format(self.md5)))
        env.sos_dict.set('__null_func__', __null_func__)
        env.sos_dict.set('__args__', self.args)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        env.sos_dict.set('__step_output__', [])

        # load configuration files
        cfg = load_config_files(self.config['config_file'])
        env.sos_dict.set('CONFIG', frozendict(cfg))

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        self._base_symbols = set(dir(__builtins__)) | set(env.sos_dict['sos_symbols_']) | set(keyword.kwlist)
        self._base_symbols -= {'dynamic'}

        # excute global definition to get some basic setup
        try:
            SoS_exec(self.workflow.global_def, get_default_global_sigil())
        except Exception as e:
            if env.verbosity > 2:
                sys.stderr.write(get_traceback())
            raise

        env.sos_dict.quick_update(self.shared)

        if isinstance(self.args, dict):
            for key, value in self.args.items():
                if not key.startswith('__'):
                    env.sos_dict.set(key, value)

    def skip(self, section):
        if section.global_def:
            try:
                SoS_exec(section.global_def, section.global_sigil)
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

    def match(self, target, step):
        # for sos_step, we need to match step name
        if isinstance(target, sos_step):
            return step.match(target.name())
        if not 'provides' in step.options:
            return False
        patterns = step.options['provides']
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
            added_node = 0
            dangling_targets, existing_targets = dag.dangling(targets)
            if dangling_targets:
                env.logger.info('Resolving {} objects from {} nodes'.format(len(dangling_targets), dag.number_of_nodes()))
            # find matching steps
            # check auxiliary steps and see if any steps provides it
            for target in dangling_targets:
                # target might no longer be dangling after a section is added.
                if target not in dag.dangling(targets)[0]:
                    continue
                mo = [(x, self.match(target, x)) for x in self.workflow.auxiliary_sections]
                mo = [x for x in mo if x[1] is not False]
                if not mo:
                    for x in self.workflow.auxiliary_sections:
                        env.logger.debug('{}: {}'.format(x.step_name(),
                            x.options['provides'] if 'provides' in x.options else ''))
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
                    env.sos_dict['__default_output__'] = [] if isinstance(target, sos_step) else [target]
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
                env.logger.info('Adding step {} with output {}'.format(res['step_name'], short_repr(res['step_output'])))
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
                    res['step_depends'], res['step_output'],
                    res['step_local_input'], res['step_local_output'], context=context)
                added_node += 1
                resolved += 1

            # for existing targets... we should check if still need to be regenerated
            for target in existing_targets:
                if target not in dag.dangling(targets)[1]:
                    continue
                mo = [(x, self.match(target, x)) for x in self.workflow.auxiliary_sections]
                mo = [x for x in mo if x[1] is not False]
                if not mo:
                    # this is ok, this is just an existing target, no one is designed to
                    # generate it.
                    continue
                if len(mo) > 1:
                    # this is not ok.
                    raise RuntimeError('Multiple steps {} to generate target {}'.format(', '.join(x[0].step_name() for x in mo), target))
                #
                # only one step, we need to process it # execute section with specified input
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
                env.logger.info('Adding step {} with output {}'.format(res['step_name'], short_repr(res['step_output'])))
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
                    res['step_depends'], res['step_output'],
                    res['step_local_input'], res['step_local_output'], context=context)
                #
                added_node += 1
                # this case do not count as resolved
                # resolved += 1
            if added_node == 0:
                break
        return resolved

    def initialize_dag(self, targets=None):
        '''Create a DAG by analyzing sections statically.'''
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
            if idx == 0:
                context['__step_output__'] = env.sos_dict['__step_output__']

            # NOTE: if a section has option 'shared', the execution of this step would
            # change dictionary, essentially making all later steps rely on this step.
            dag.add_step(section.uuid,
                section.step_name(),
                idx,
                res['step_input'],
                res['step_depends'],
                res['step_output'],
                res['step_local_input'],
                res['step_local_output'],
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
        #dag.show_nodes()
        # trim the DAG if targets are specified
        if targets:
            dag = dag.subgraph_from(targets)
        # write DAG for debugging purposes
        #dag.write_dot(os.path.join(env.exec_dir, '.sos', '{}.dot'.format(self.workflow.name)))
        # check error
        cycle = dag.circular_dependencies()
        if cycle:
            raise RuntimeError('Circular dependency detected {}. It is likely a later step produces input of a previous step.'.format(cycle))

        self.save_dag(dag)
        return dag

    def save_workflow_signature(self, dag):
        '''Save tracked files in .sos so that untracked files can be cleaned by command
        sos clean.
        '''
        if '__workflow_sig__' in env.sos_dict:
            with open(env.sos_dict['__workflow_sig__'], 'a') as sigfile:
                sigfile.write('# end time: {}\n'.format(time.strftime('%a, %d %b %Y %H:%M:%S +0000', time.gmtime())))
                sigfile.write('# input and dependent files\n')

    def run(self, targets=None, parent_pipe=None, my_workflow_id=None, mode='run'):
        '''Execute a workflow with specified command line args. If sub is True, this
        workflow is a nested workflow and be treated slightly differently.
        '''
        #
        # There are threee cases
        #
        # parent_pipe = None: this is the master workflow executor
        # parent_pipe != None, my_workflow_id != None: this is a nested workflow inside a master workflow
        #   executor and needs to pass tasks etc to master
        # parent_pipe != None, my_workflow_id == None: this is a nested workflow inside a task and needs to
        #   handle its own tasks.
        #
        nested = parent_pipe is not None and my_workflow_id is not None

        def i_am():
            return 'Nested' if nested else 'Master'

        #
        # if the exexcutor is started from sos_run, these should also be passed
        if 'sig_mode' in self.config:
            env.config['sig_mode'] = self.config['sig_mode']
        if 'verbosity' in self.config:
            env.verbosity = self.config['verbosity']

        self.reset_dict()
        env.config['run_mode'] = mode
        # passing run_mode to SoS dict so that users can execute blocks of
        # python statements in different run modes.
        env.sos_dict.set('run_mode', env.config['run_mode'])
        # process step of the pipelinp
        dag = self.initialize_dag(targets=targets)

        # if targets are specified and there are only signatures for them, we need
        # to remove the signature and really generate them
        if targets:
            for t in targets:
                if not FileTarget(t).exists('target') and FileTarget(t).exists('signature'):
                    env.logger.info('Re-generating {}'.format(t))
                    FileTarget(t).remove('signature')
                else:
                    env.logger.info('Target {} already exists'.format(t))
        # process step of the pipelinp
        #
        # running processes. It consisists of
        #
        # [ [proc, queue], pipe, node]
        #
        # where:
        #   proc, queue: process, which is None for the nested workflow.
        #   pipe: pipe to get information from workers
        #   node: node that is being executed, which is a dummy node
        #       created on the fly for steps passed from nested workflow
        #
        procs = []
        #
        # process pool that is used to pool temporarily unused processed.
        pool = []
        #
        wf_result = {'__workflow_id__': my_workflow_id, 'shared': {}}
        #
        # steps sent and queued from the nested workflow
        # they will be executed in random but at a higher priority than the steps
        # on the master process.
        self.step_queue = {}
        try:
            prog = ProgressBar(desc=self.workflow.name, total=dag.num_nodes(), disable=dag.num_nodes() <= 1 or env.verbosity != 1)
            exec_error = ExecuteError(self.workflow.name)
            while True:
                # step 1: check existing jobs and see if they are completed
                for idx, proc in enumerate(procs):
                    if proc is None:
                        continue

                    [p, q, runnable] = proc
                    # echck if there is any message from the pipe
                    if not q.poll():
                        continue

                    # receieve something from the pipe
                    res = q.recv()
                    #
                    # if this is NOT a result, rather some request for task, step, workflow etc
                    if isinstance(res, str):
                        if nested:
                            raise RuntimeError('Nested workflow is not supposed to receive task, workflow, or step requests. {} received.'.format(res))
                        if res.startswith('task'):
                            env.logger.debug('{} receives task request {}'.format(i_am(), res))
                            host = res.split(' ')[1]
                            if host == '__default__':
                                if 'default_queue' in env.config:
                                    host = env.config['default_queue']
                                else:
                                    host = 'localhost'
                            runnable._host = Host(host)
                            runnable._pending_tasks = res.split(' ')[2:]
                            for task in runnable._pending_tasks:
                                runnable._host.submit_task(task)
                            runnable._status = 'task_pending'
                            env.logger.trace('Step becomes task_pending')
                            continue
                        elif res.startswith('step'):
                            # step sent from nested workflow
                            step_id = res.split(' ')[1]
                            step_params = q.recv()
                            env.logger.debug('{} receives step request {} with args {}'.format(i_am(), step_id, step_params[3]))
                            self.step_queue[step_id] = step_params
                            continue
                            #
                        elif res.startswith('workflow'):
                            workflow_id = res.split(' ')[1]
                            # receive the real definition
                            env.logger.debug('{} receives workflow request {}'.format(i_am(), workflow_id))
                            # (wf, args, shared, config)
                            wf, targets, args, shared, config = q.recv()
                            # a workflow needs to be executed immediately because otherwise if all workflows
                            # occupies all workers, no real step could be executed.

                            # if pool is empty, create a new process
                            if not pool:
                                q1, q2 = mp.Pipe()
                                worker = SoS_Worker(pipe=q2, config=config, args=args)
                                worker.start()
                            else:
                                # get worker, q and runnable is not needed any more
                                worker, q1, _ = pool.pop(0)

                            # now we would like to find a worker and
                            runnable._pending_workflow = workflow_id
                            runnable._status = 'workflow_pending'

                            wfrunnable = dummy_node()
                            wfrunnable._node_id = workflow_id
                            wfrunnable._status = 'workflow_running_pending'
                            wfrunnable._pending_workflow = workflow_id
                            #
                            q1.send(('workflow', workflow_id, wf, targets, args, shared, config))
                            procs.append([worker, q1, wfrunnable])
                            #
                            continue
                        else:
                            raise RuntimeError('Unexpected value from step {}'.format(res))

                    # if we does get the result, we send the process to pool
                    pool.append(procs[idx])
                    procs[idx] = None

                    env.logger.debug('{} receive a result'.format(i_am()))
                    if hasattr(runnable, '_from_nested'):
                        # if the runnable is from nested, we will need to send the result back to the workflow
                        env.logger.debug('{} send res to nested'.format(i_am()))
                        runnable._status = 'completed'
                        runnable._child_pipe.send(res)
                    elif isinstance(res, (UnknownTarget, RemovedTarget)):
                        runnable._status = None
                        target = res.target
                        if dag.regenerate_target(target):
                            #runnable._depends_targets.append(target)
                            #dag._all_dependent_files[target].append(runnable)
                            dag.build(self.workflow.auxiliary_sections)
                            #
                            cycle = dag.circular_dependencies()
                            if cycle:
                                raise RuntimeError('Circular dependency detected {} after regeneration. It is likely a later step produces input of a previous step.'.format(cycle))

                        else:
                            if self.resolve_dangling_targets(dag, [target]) == 0:
                                raise RuntimeError('Failed to regenerate or resolve {}{}.'
                                    .format(target, dag.steps_depending_on(target, self.workflow)))
                            runnable._depends_targets.append(target)
                            dag._all_dependent_files[target].append(runnable)
                            dag.build(self.workflow.auxiliary_sections)
                            #
                            cycle = dag.circular_dependencies()
                            if cycle:
                                raise RuntimeError('Circular dependency detected {}. It is likely a later step produces input of a previous step.'.format(cycle))
                        self.save_dag(dag)
                    elif isinstance(res, UnavailableLock):
                        runnable._status = 'signature_pending'
                        runnable._signature = (res.output, res.sig_file)
                        section = self.workflow.section_by_id(runnable._step_uuid)
                        env.logger.info('Waiting on another process for step {}'.format(section.step_name()))
                    # if the job is failed
                    elif isinstance(res, Exception):
                        env.logger.debug('{} received an exception'.format(i_am()))
                        runnable._status = 'failed'
                        exec_error.append(runnable._node_id, res)
                        # if this is a node for a running workflow, need to mark it as failed as well
                        #                        for proc in procs:
                        if isinstance(runnable, dummy_node) and hasattr(runnable, '_pending_workflow'):
                            for proc in procs:
                                if proc is None:
                                    continue
                                if proc[2]._status.endswith('_pending') and proc[2]._pending_workflow == runnable._pending_workflow:
                                    proc[2]._status = 'failed'
                        prog.update(1)
                    elif '__step_name__' in res:
                        env.logger.debug('{} receive step result '.format(i_am()))
                        # if the result of the result of a step
                        svar = {}
                        for k, v in res.items():
                            if k == '__shared__':
                                svar = v
                                env.sos_dict.update(v)
                            else:
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
                            node._context.update(svar)
                            node._context['__completed__'].append(res['__step_name__'])
                        dag.update_step(runnable, env.sos_dict['__step_input__'],
                            env.sos_dict['__step_output__'],
                            env.sos_dict['__step_depends__'],
                            env.sos_dict['__step_local_input__'],
                            env.sos_dict['__step_local_output__'])
                        runnable._status = 'completed'
                        prog.update(1)
                    elif '__workflow_id__' in res:
                        # result from a workflow
                        # the worker process has been returned to the pool, now we need to
                        # notify the step that is waiting for the result
                        env.logger.debug('{} receive workflow result'.format(i_am()))
                        for proc in procs:
                            if proc[2]._status == 'workflow_pending' and proc[2]._pending_workflow == res['__workflow_id__']:
                                proc[1].send(res)
                                proc[2]._status = 'running'
                                break
                    else:
                        raise RuntimeError('Unrecognized response from a step: {}'.format(res))

                # remove None
                procs = [x for x in procs if x is not None]

                # step 2: check if some jobs are done
                for proc_idx, proc in enumerate(procs):
                    # if a job is pending, check if it is done.
                    if proc[2]._status == 'task_pending':
                        res = proc[2]._host.check_status(proc[2]._pending_tasks)
                        #env.logger.warning(res)
                        if any(x in ('aborted', 'failed', 'result-mismatch') for x in res):
                            for t, s in zip(proc[2]._pending_tasks, res):
                                if s in ('aborted', 'failed', 'result-mismatch') and not (hasattr(proc[2], '_killed_tasks') and t in proc[2]._killed_tasks):
                                    env.logger.warning('{} ``{}``'.format(t, s))
                                    if not hasattr(proc[2], '_killed_tasks'):
                                        proc[2]._killed_tasks = {t}
                                    else:
                                        proc[2]._killed_tasks.add(t)
                            if all(x in ('completed', 'aborted', 'failed', 'result-mismatch') for x in res):
                                # we try to get .err .out etc even when jobs are failed.
                                task_status = proc[2]._host.retrieve_results(proc[2]._pending_tasks)
                                proc[1].send(task_status)
                                proc[2]._status == 'failed'
                                raise RuntimeError('Proc {}: {} completed, {} failed, {} aborted, {} mismatch'.format(proc_idx,
                                    len([x for x in res if x=='completed']), len([x for x in res if x=='failed']),
                                    len([x for x in res if x=='aborted']), len([x for x in res if x=='result-mismatch']) ))
                        if any(x in ('pending', 'submitted', 'running') for x in res):
                            continue
                        elif all(x == 'completed' for x in res):
                            env.logger.debug('Proc {} puts results for {} from step {}'.format(proc_idx, ' '.join(proc[2]._pending_tasks), proc[2]._node_id))
                            res = proc[2]._host.retrieve_results(proc[2]._pending_tasks)
                            proc[1].send(res)
                            proc[2]._status = 'running'
                        else:
                            raise RuntimeError('Job returned with status {}'.format(res))

                # step 3: check if there is room and need for another job
                while True:
                    #env.logger.error('{} {}'.format(i_am(), [x[2]._status for x in procs]))
                    num_running = len([x for x in procs if not x[2]._status.endswith('_pending')])
                    if num_running >= env.config['max_procs']:
                        break
                    #
                    # if steps from child nested workflow?
                    if self.step_queue:
                        step_id, step_param = self.step_queue.popitem()
                        section, context, shared, args, run_mode, sig_mode, verbosity, pipe = step_param
                        # run it!
                        # if pool is empty, create a new process
                        if not pool:
                            q1, q2 = mp.Pipe()
                            worker = SoS_Worker(pipe=q2, config=self.config, args=self.args)
                            worker.start()
                        else:
                            # get worker, q and runnable is not needed any more
                            worker, q1, _ = pool.pop(0)

                        runnable = dummy_node()
                        runnable._node_id = step_id
                        runnable._status = 'running'
                        runnable._from_nested = True
                        runnable._child_pipe = pipe

                        env.logger.debug('{} sends {} from step queue with args {}'.format(i_am(), step_id, args))
                        q1.send(('step', section, context, shared, args, run_mode, sig_mode, verbosity))
                        procs.append( [worker, q1, runnable])
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
                    # execute section with specified input
                    runnable._status = 'running'

                    # workflow shared variables
                    shared = {x: env.sos_dict[x] for x in self.shared.keys() if x in env.sos_dict and pickleable(env.sos_dict[x], x)}
                    if 'shared' in section.options:
                        if isinstance(section.options['shared'], str):
                            svars = [section.options['shared']]
                        elif isinstance(section.options['shared'], dict):
                            svars = section.options['shared'].keys()
                        elif isinstance(section.options['shared'], Sequence):
                            svars = []
                            for x in section.options['shared']:
                                if isinstance(x, str):
                                    svars.append(x)
                                elif isinstance(x, dict):
                                    svars.extend(x.keys())
                                else:
                                    raise ValueError('Unacceptable value for parameter shared: {}'.format(section.options['shared']))
                        else:
                            raise ValueError('Unacceptable value for parameter shared: {}'.format(section.options['shared']))
                        shared.update({x: env.sos_dict[x] for x in svars if x in env.sos_dict and pickleable(env.sos_dict[x], x)})

                    if '__workflow_sig__' in env.sos_dict:
                        runnable._context['__workflow_sig__'] = env.sos_dict['__workflow_sig__']

                    if not nested:

                        # if pool is empty, create a new process
                        if not pool:
                            q1, q2 = mp.Pipe()
                            worker = SoS_Worker(pipe=q2, config=self.config, args=self.args)
                            worker.start()
                        else:
                            # get worker, q and runnable is not needed any more
                            worker, q1, _ = pool.pop(0)

                        env.logger.debug('{} execute {} from DAG'.format(i_am(), section.md5))
                        q1.send(('step', section, runnable._context, shared, self.args, env.config['run_mode'], env.config['sig_mode'], env.verbosity))
                        procs.append( [worker, q1, runnable])
                    else:
                        # send the step to the parent
                        step_id = uuid.uuid4()
                        env.logger.debug('{} send step {} to master with args {}'.format(i_am(), step_id, self.args))
                        parent_pipe.send('step {}'.format(step_id))
                        q = mp.Pipe()
                        parent_pipe.send((section, runnable._context, shared, self.args, env.config['run_mode'], env.config['sig_mode'], env.verbosity, q[1]))
                        # this is a real step
                        runnable._status = 'step_pending'
                        procs.append([None, q[0], runnable])

                #
                num_running = len([x for x in procs if not x[2]._status.endswith('_pending')])

                if not procs or all(x[2]._status == 'failed' for x in procs):
                    break
                # if -W is specified, or all task queues are not wait
                elif all(x[2]._status == 'task_pending' for x in procs) and \
                        (env.config['wait_for_task'] is False or \
                        (env.config['wait_for_task'] is None and Host.not_wait_for_tasks())):
                    # if all jobs are pending, let us check if all jbos have been submitted.
                    pending_tasks = []
                    running_tasks = []
                    for n in [x[2] for x in procs]:
                        p, r = n._host._task_engine.get_tasks()
                        pending_tasks.extend(p)
                        running_tasks.extend([(n._host.alias, x) for x in r])
                    if not pending_tasks and running_tasks:
                        env.logger.trace('Exit with {} running tasks: '.format(len(running_tasks), running_tasks))
                        raise PendingTasks(running_tasks)
                else:
                    time.sleep(0.1)
        except PendingTasks as e:
            self.record_quit_status(e.tasks)
            wf_result['pending_tasks'] = [x[1] for x in running_tasks]
            env.logger.info('Workflow {} (ID={}) exits with {} running tasks'.format(self.workflow.name, self.md5, len(e.tasks)))
            for task in e.tasks:
                env.logger.info(task[1])
            # close all processes
        except Exception as e:
            for p, _, _ in procs + pool:
                p.terminate()
            raise e
        finally:
            if not nested:
                for _, p, _ in procs + pool:
                    p.send(None)
                time.sleep(0.1)
                for w, _, _ in procs + pool:
                    if w.is_alive():
                        w.terminate()
                        w.join()
            prog.close()
        #
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
            if parent_pipe is not None:
                parent_pipe.send(exec_error)
                return wf_result
            else:
                raise exec_error
        elif 'pending_tasks' not in wf_result or not wf_result['pending_tasks']:
            # remove task pending status if the workflow is completed normally
            try:
                wf_status = os.path.join(os.path.expanduser('~'), '.sos', self.md5 + '.status')
                if os.path.isfile(wf_status):
                    os.remove(wf_status)
            except Exception as e:
                env.logger.warning('Failed to clear workflow status file: {}'.format(e))
            self.save_workflow_signature(dag)
            env.logger.info('Workflow {} (ID={}) is executed successfully.'.format(self.workflow.name, self.md5))
        else:
            # exit with pending tasks
            pass
        wf_result['shared'] = {x:env.sos_dict[x] for x in self.shared.keys() if x in env.sos_dict}
        if parent_pipe:
            parent_pipe.send(wf_result)
        else:
            return wf_result
