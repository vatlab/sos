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
from .step_executor import Step_Executor, analyze_section, PendingTasks
from .utils import env, Error, WorkflowDict, get_traceback, short_repr, pickleable, \
    load_config_files, save_var, load_var, SlotManager
from .eval import SoS_exec
from .dag import SoS_DAG
from .targets import BaseTarget, file_target, UnknownTarget, RemovedTarget, UnavailableLock, sos_variable, textMD5, sos_step, Undetermined
from .pattern import extract_pattern
from .hosts import Host

__all__ = []

try:
    # https://github.com/pytest-dev/pytest-cov/issues/139
    from pytest_cov.embed import cleanup_on_sigterm
    cleanup_on_sigterm()
except:
    pass

class ExecuteError(Error):
    """An exception to collect exceptions raised during run time so that
    other branches of the DAG would continue if some nodes fail to execute."""
    def __init__(self, workflow):
        Error.__init__(self)
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
        newline = '\n' if self.message else ''
        self.message += f'{newline}[{short_line}]: {error}'

def __null_func__(*args, **kwargs):
    '''This function will be passed to SoS's namespace and be executed
    to evaluate functions of input, output, and depends directives.'''
    def _flatten(x):
        if isinstance(x, str):
            return [x]
        elif isinstance(x, Sequence):
            return sum((_flatten(k) for k in x), [])
        elif hasattr(x, '__flattenable__'):
            return _flatten(x.flatten())
        else:
            return [x]

    return _flatten(args), kwargs

class SoS_Worker(mp.Process):
    '''
    Worker process to process SoS step or workflow in separate process.
    '''
    def __init__(self,  pipe, config=None, args=None, **kwargs):
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
        self.config = {} if config is None else config
        self.args = [] if args is None else args

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
        load_config_files(self.config['config_file'])

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        self._base_symbols = set(dir(__builtins__)) | set(env.sos_dict['sos_symbols_']) | set(keyword.kwlist)
        # if users use sos_run, the "scope" of the step goes beyong names in this step
        # so we cannot save signatures for it.
        self._base_symbols -= {'dynamic', 'sos_run'}

        if isinstance(self.args, dict):
            for key, value in self.args.items():
                if not key.startswith('__'):
                    env.sos_dict.set(key, value)

    def run(self):
        # wait to handle jobs
        while True:
            try:
                work = self.pipe.recv()
                if work is None:
                    break
                env.logger.debug(f'Worker {self.name} receives request {short_repr(work)}')
                if work[0] == 'step':
                    # this is a step ...
                    self.run_step(*work[1:])
                else:
                    self.run_workflow(*work[1:])
                env.logger.debug(f'Worker {self.name} completes request {short_repr(work)}')
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
        env.logger.debug(f'Worker {self.name} working on a workflow {workflow_id} with args {args}')
        executer = Base_Executor(wf, args=args, shared=shared, config=config)
        # we send the pipe to subworkflow, which would send
        # everything directly to the master process, so we do not
        # have to collect result here
        try:
            executer.run(targets=targets, parent_pipe=self.pipe, my_workflow_id=workflow_id)
        except Exception as e:
            self.pipe.send(e)


    def run_step(self, section, context, shared, args, run_mode, sig_mode, verbosity):
        env.logger.debug(f'Worker {self.name} working on {section.step_name()} with args {args}')
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
            SoS_exec(section.global_def)
        except RuntimeError:
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

class ProcInfo(object):
    def __init__(self, worker, pipe, step):
        self.worker = worker
        self.pipe = pipe
        self.step = step

    def set_status(self, status):
        self.step._status = status

    def in_status(self, status):
        return self.step._status == status

    def status(self):
        return self.step._status

    def is_pending(self):
        return self.step._status.endswith('_pending')

class ExecutionManager(object):
    # this class managers workers and their status ...
    def __init__(self, max_workers, master=True):
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
        self.procs = []

        # process pool that is used to pool temporarily unused processed.
        self.pool = []

        self.slot_manager = SlotManager(reset=master)
        self.last_num_procs = None

        self.max_workers = max_workers

    def execute(self, runnable, config, args, spec):
        if not self.pool:
            q1, q2 = mp.Pipe()
            worker = SoS_Worker(pipe=q2, config=config, args=args)
            worker.start()
        else:
            # get worker, q and runnable is not needed any more
            pi = self.pool.pop(0)
            worker = pi.worker
            q1 = pi.pipe

        q1.send(spec)
        self.procs.append(ProcInfo(worker=worker, pipe=q1, step=runnable))

    def add_placeholder_worker(self, runnable, pipe):
        runnable._status = 'step_pending'
        self.procs.append(ProcInfo(worker=None, pipe=pipe, step=runnable))

    def all_busy(self):
        n = len([x for x in self.procs if x and not x.is_pending()])
        if self.last_num_procs is None:
            if n > 0:
                self.slot_manager.acquire(n, self.max_workers)
            self.last_num_procs = n
        elif n != self.last_num_procs:
            if self.last_num_procs > n:
                self.slot_manager.release(self.last_num_procs - n)
            else:
                # we force the increase of numbers because the increment is observed
                self.slot_manager.acquire(n - self.last_num_procs, self.max_workers, force=True)
            self.last_num_procs = n
        return n >= self.max_workers

    def all_done_or_failed(self):
        return not self.procs or all(x.in_status('failed') for x in self.procs)
    
    def mark_idle(self, idx):
        self.pool.append(self.procs[idx])
        self.procs[idx] = None

    def cleanup(self):
        self.procs = [x for x in self.procs if x is not None]

    def terminate(self, brutal=False):
        self.cleanup()
        if not brutal:
            for proc in self.procs + self.pool:
                proc.pipe.send(None)
            time.sleep(0.1)
            for proc in self.procs + self.pool:
                if proc.worker and proc.worker.is_alive():
                    proc.worker.terminate()
                    proc.worker.join()
        else:
            for proc in self.procs + self.pool:
                # proc can be fake if from a nested workflow
                if proc.worker:
                    proc.worker.terminate()


class Base_Executor:
    '''This is the base class of all executor that provides common
    set up and tear functions for all executors.'''
    def __init__(self, workflow=None, args=None, shared=None, config=None):
        self.workflow = workflow
        self.args = [] if args is None else args
        if '__args__' not in self.args:
            # if there is __args__, this is a nested workflow and we do not test this.
            for idx,arg in enumerate(self.args):
                wf_pars = self.workflow.parameters().keys()
                if isinstance(arg, str) and arg.startswith('--'):
                    if not wf_pars:
                        raise ValueError(
                            f'Undefined parameter {arg[2:]} for command line argument "{" ".join(args[idx:])}".')
                    pars = [arg[2:], arg[2:].replace('-', '_').split('=')[0]]
                    if arg[2:].startswith('no-'):
                        pars.extend([arg[5:], arg[5:].replace('-', '_').split('=')[0]])
                    if not any(x in wf_pars for x in pars):
                        raise ValueError(
                            f'Undefined parameter {arg[2:]} for command line argument "{" ".join(args[idx:])}". Acceptable parameters are: {", ".join(wf_pars)}')

        self.shared = {} if shared is None else shared
        self.config = {} if config is None else config
        env.config.update(self.config)
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
        #
        # the md5 of the master workflow would be passed from master workflow...
        if 'master_md5' not in self.config:
            self.config['master_md5'] = self.md5
        if env.config['sig_mode'] != 'ignore' and self.workflow:
            # remove old workflow file.
            with open(os.path.join(env.exec_dir, '.sos', f'{self.md5}.sig'), 'a') as sig:
                sig.write(f'# workflow: {self.workflow.name}\n')
                sig.write(f'# script: {self.workflow.content.filename}\n')
                sig.write(f'# included: {",".join([x[1] for x in self.workflow.content.included])}\n')
                sig.write(f'# configuration: {self.config.get("config_file", "")}\n')
                sig.write(f'# start time: {time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())}\n')
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
                            _, v = load_var(line)
                            env.config['resumed_tasks'].add(v[1])
            else:
                env.logger.info(f'Workflow {self.md5} has been completed.')
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
        if not self.config['output_dag']:
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
                dag_name = f'{self.config["output_dag"][:-4] if self.config["output_dag"].endswith(".dot") else self.config["output_dag"]}_{self.dag_count}.dot'
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
                sig.write(f'{step.step_name()}: {step.md5}\n')
            sig.write('# Command line options\n')
            sig.write(f'{self.args}\n')
            self.sig_content = sig.getvalue()
        return textMD5(self.sig_content)[:16]

    def reset_dict(self):
        env.sos_dict = WorkflowDict()
        env.parameter_vars.clear()
        env.config.update(self.config)

        # inject a few things
        if self.md5:
            env.sos_dict.set('__workflow_sig__', os.path.join(env.exec_dir, '.sos', f'{self.md5}.sig'))
        env.sos_dict.set('__null_func__', __null_func__)
        env.sos_dict.set('__args__', self.args)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        env.sos_dict.set('__step_output__', [])

        # load configuration files
        load_config_files(self.config['config_file'])
        # if check_readonly is set to True, allow checking readonly vars
        #if cfg.get('sos', {}).get('change_all_cap_vars', None) is not None:
        #    if cfg['sos']['change_all_cap_vars'] not in ('warning', 'error'):
        #        env.logger.error(
        #            f'Configuration sos.change_all_cap_vars can only be warning or error: {cfg["sos"]["change_all_cap_vars"]} provided')
        #    else:
        #        env.sos_dict._change_all_cap_vars = cfg['sos']['change_all_cap_vars']

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        self._base_symbols = set(dir(__builtins__)) | set(env.sos_dict['sos_symbols_']) | set(keyword.kwlist)
        self._base_symbols -= {'dynamic'}

        # excute global definition to get some basic setup
        try:
            SoS_exec(self.workflow.global_def)
        except Exception:
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
                SoS_exec(section.global_def)
            except RuntimeError as e:
                if env.verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise RuntimeError(f'Failed to execute statements\n"{section.global_def}"\n{e}')
        #
        if 'skip' in section.options:
            val_skip = section.options['skip']
            if val_skip is None or val_skip is True:
                env.logger.info(f'``{section.step_name(True)}`` is ``ignored`` due to skip option.')
                return True
            elif val_skip is not False:
                raise RuntimeError(
                    f'The value of section option skip can only be None, True or False, {val_skip} provided')
        return False

    def match(self, target, step):
        # for sos_step, we need to match step name
        if isinstance(target, sos_step):
            return step.match(target.target_name())
        if not 'provides' in step.options and not 'autoprovides' in step.options:
            return False
        patterns = step.options['provides'] if 'provides' in step.options else step.options['autoprovides']
        if isinstance(patterns, (str, BaseTarget)):
            patterns = [patterns]
        elif not isinstance(patterns, Sequence):
            raise RuntimeError(f'Unknown target to match: {patterns}')
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
                return {x:y[0] for x,y in res.items()}
            # string match
            elif file_target(p) == file_target(target):
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
                env.logger.debug(f'Resolving {dangling_targets} objects from {dag.number_of_nodes()} nodes')
            # find matching steps
            # check auxiliary steps and see if any steps provides it
            for target in dangling_targets:
                # target might no longer be dangling after a section is added.
                if target not in dag.dangling(targets)[0]:
                    continue
                mo = [(x, self.match(target, x)) for x in self.workflow.auxiliary_sections]
                mo = [x for x in mo if x[1] is not False]
                if not mo:
                    #
                    # if no step produces the target, it is possible that it is an indexed step
                    # so the execution of its previous steps would solves the dependency
                    #
                    # find all the nodes that depends on target
                    nodes = dag._all_dependent_files[target]
                    for node in nodes:
                        # if this is an index step... simply let it depends on previous steps
                        if node._node_index is not None:
                            indexed = [x for x in dag.nodes() if x._node_index is not None and x._node_index < node._node_index and isinstance(x._output_targets, Undetermined)]
                            indexed.sort(key = lambda x: x._node_index)
                            if not indexed:
                                raise RuntimeError(
                                    f'No step to generate target {target}{dag.steps_depending_on(target, self.workflow)}')
                            if isinstance(target, sos_step) and not any(self.workflow.section_by_id(x._step_uuid).match(target.target_name()) for x in indexed):
                                raise RuntimeError(
                                    f'No step to generate target {target}{dag.steps_depending_on(target, self.workflow)}')
                            # now, if it is not a sos_step, but its previous steps have already been executed and still
                            # could not satisfy the requirement..., we should generate an error
                            if not any(x._status is None or x._status.endswith('pending') for x in indexed):
                                # all previous status has been failed or completed...
                                raise RuntimeError(
                                    f'Previous step{" has" if len(indexed) == 1 else "s have"} not generated target {target}{dag.steps_depending_on(target, self.workflow)}')
                            if not isinstance(node._input_targets, Undetermined):
                                node._input_targets = Undetermined('')
                            if not isinstance(node._depends_targets, Undetermined):
                                node._depends_targets = Undetermined('')
                        else:
                            raise RuntimeError(
                                f'No step to generate target {target}{dag.steps_depending_on(target, self.workflow)}')
                    if nodes:
                        resolved += 1
                    continue
                if len(mo) > 1:
                    raise RuntimeError(
                        f'Multiple steps {", ".join(x[0].step_name() for x in mo)} to generate target {target}')
                #
                # only one step, we need to process it # execute section with specified input
                #
                # NOTE:  Auxiliary can be called with different output files and matching pattern
                # so we are actually creating a new section each time we need an auxillary step.
                #
                section = mo[0][0]
                if isinstance(mo[0][1], dict):
                    for k,v in mo[0][1].items():
                        env.sos_dict.set(k, v)
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
                env.logger.debug(
                    f'Adding step {res["step_name"]} with output {short_repr(res["step_output"])} to resolve target {target}')
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
                dag.add_step(section.uuid, f'{section.step_name()} {short_repr(env.sos_dict["__default_output__"])}', None, res['step_input'].targets(),
                             res['step_depends'].targets(), res['step_output'].targets(), context=context)
                added_node += 1
                resolved += 1

            # for existing targets... we should check if it actually exists. If
            # not it would still need to be regenerated
            for target in existing_targets:
                if target not in dag.dangling(targets)[1]:
                    continue
                if file_target(target).target_exists('target') if isinstance(target, str) else target.target_exists('target'):
                    continue
                mo = [(x, self.match(target, x)) for x in self.workflow.auxiliary_sections]
                mo = [x for x in mo if x[1] is not False]
                if not mo:
                    # this is ok, this is just an existing target, no one is designed to
                    # generate it.
                    continue
                if len(mo) > 1:
                    # this is not ok.
                    raise RuntimeError(
                        f'Multiple steps {", ".join(x[0].step_name() for x in mo)} to generate target {target}')
                #
                # only one step, we need to process it # execute section with specified input
                #
                section = mo[0][0]
                if isinstance(mo[0][1], dict):
                    for k,v in mo[0][1].items():
                        env.sos_dict.set(k, v)
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
                env.logger.debug(
                    f'Adding step {res["step_name"]} with output {short_repr(res["step_output"])} to resolve target {target}')
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
                dag.add_step(section.uuid, f'{section.step_name()} {short_repr(env.sos_dict["__default_output__"])}',
                             None, res['step_input'],
                             res['step_depends'], res['step_output'], context=context)
                #
                added_node += 1
                # this case do not count as resolved
                # resolved += 1
            if added_node == 0:
                break
        #dag.show_nodes()
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
                    '__changed_vars__': changed_vars}

            # for nested workflow, the input is specified by sos_run, not None.
            if idx == 0:
                context['__step_output__'] = env.sos_dict['__step_output__']

            # NOTE: if a section has option 'shared', the execution of this step would
            # change dictionary, essentially making all later steps rely on this step.
            dag.add_step(section.uuid,
                section.step_name(),
                idx,
                res['step_input'].targets(),
                res['step_depends'].targets(),
                res['step_output'].targets(),
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
        if self.resolve_dangling_targets(dag, targets) == 0:
            if targets:
                raise RuntimeError(f'No auxiliary step to generate target {targets}.')
        # now, there should be no dangling targets, let us connect nodes
        dag.build(self.workflow.auxiliary_sections)
        #dag.show_nodes()
        # trim the DAG if targets are specified
        if targets:
            dag = dag.subgraph_from(targets)
        # check error
        cycle = dag.circular_dependencies()
        if cycle:
            raise RuntimeError(
                f'Circular dependency detected {cycle}. It is likely a later step produces input of a previous step.')

        self.save_dag(dag)
        return dag

    def save_workflow_signature(self, dag):
        '''Save tracked files in .sos so that untracked files can be cleaned by command
        sos clean.
        '''
        if '__workflow_sig__' in env.sos_dict:
            with open(env.sos_dict['__workflow_sig__'], 'a') as sigfile:
                sigfile.write(f'# end time: {time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime())}\n')
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

        # if targets are specified and there are only signatures for them, we need
        # to remove the signature and really generate them
        if targets:
            for t in targets:
                if file_target(t).target_exists('target'):
                    env.logger.info(f'Target {t} already exists')
                elif file_target(t).target_exists('signature'):
                    env.logger.info(f'Re-generating {t}')
                    file_target(t).remove('signature')
            targets = [x for x in targets if not file_target(x).target_exists('target')]
            if not targets:
                raise RuntimeError('All targets already exists.')

        # process step of the pipelinp
        dag = self.initialize_dag(targets=targets)
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
        manager = ExecutionManager(env.config['max_procs'], master=not nested)
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
                for idx, proc in enumerate(manager.procs):
                    if proc is None:
                        continue

                    runnable = proc.step
                    # echck if there is any message from the pipe
                    if not proc.pipe.poll():
                        continue

                    # receieve something from the pipe
                    res = proc.pipe.recv()
                    #
                    # if this is NOT a result, rather some request for task, step, workflow etc
                    if isinstance(res, str):
                        if nested:
                            raise RuntimeError(
                                f'Nested workflow is not supposed to receive task, workflow, or step requests. {res} received.')
                        if res.startswith('task'):
                            env.logger.debug(f'{i_am()} receives task request {res}')
                            host = res.split(' ')[1]
                            if host == '__default__':
                                if 'default_queue' in env.config:
                                    host = env.config['default_queue']
                                else:
                                    host = 'localhost'
                            runnable._host = Host(host)
                            new_tasks = res.split(' ')[2:]
                            if hasattr(runnable, '_pending_tasks'):
                                runnable._pending_tasks.extend(new_tasks)
                            else:
                                runnable._pending_tasks = new_tasks
                            for task in new_tasks:
                                runnable._host.submit_task(task)
                            runnable._status = 'task_pending'
                            env.logger.trace('Step becomes task_pending')
                            continue
                        elif res.startswith('step'):
                            # step sent from nested workflow
                            step_id = res.split(' ')[1]
                            step_params = proc.pipe.recv()
                            env.logger.debug(f'{i_am()} receives step request {step_id} with args {step_params[3]}')
                            self.step_queue[step_id] = step_params
                            continue
                            #
                        elif res.startswith('workflow'):
                            workflow_id = res.split(' ')[1]
                            # receive the real definition
                            env.logger.debug(f'{i_am()} receives workflow request {workflow_id}')
                            # (wf, args, shared, config)
                            wf, targets, args, shared, config = proc.pipe.recv()
                            # a workflow needs to be executed immediately because otherwise if all workflows
                            # occupies all workers, no real step could be executed.


                            # now we would like to find a worker and
                            runnable._pending_workflow = workflow_id
                            runnable._status = 'workflow_pending'

                            wfrunnable = dummy_node()
                            wfrunnable._node_id = workflow_id
                            wfrunnable._status = 'workflow_running_pending'
                            wfrunnable._pending_workflow = workflow_id
                            #
                            manager.execute(wfrunnable, config=config, args=args,
                                    spec=('workflow', workflow_id, wf, targets, args, shared, config))
                            #
                            continue
                        else:
                            raise RuntimeError(f'Unexpected value from step {short_repr(res)}')

                    # if we does get the result, we send the process to pool
                    manager.mark_idle(idx)

                    env.logger.debug(f'{i_am()} receive a result {short_repr(res)}')
                    if hasattr(runnable, '_from_nested'):
                        # if the runnable is from nested, we will need to send the result back to the workflow
                        env.logger.debug(f'{i_am()} send res to nested')
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
                                raise RuntimeError(
                                    f'Circular dependency detected {cycle} after regeneration. It is likely a later step produces input of a previous step.')

                        else:
                            if self.resolve_dangling_targets(dag, [target]) == 0:
                                raise RuntimeError(
                                    f'Failed to regenerate or resolve {target}{dag.steps_depending_on(target, self.workflow)}.')
                            if not isinstance(runnable._depends_targets, Undetermined):
                                runnable._depends_targets.append(target)
                            if runnable not in dag._all_dependent_files[target]:
                                dag._all_dependent_files[target].append(runnable)                            
                            dag.build(self.workflow.auxiliary_sections)
                            #
                            cycle = dag.circular_dependencies()
                            if cycle:
                                raise RuntimeError(
                                    f'Circular dependency detected {cycle}. It is likely a later step produces input of a previous step.')
                        self.save_dag(dag)
                    elif isinstance(res, UnavailableLock):
                        runnable._status = 'signature_pending'
                        runnable._signature = (res.output, res.sig_file)
                        section = self.workflow.section_by_id(runnable._step_uuid)
                        env.logger.info(f'Waiting on another process for step {section.step_name(True)}')
                    # if the job is failed
                    elif isinstance(res, Exception):
                        env.logger.debug(f'{i_am()} received an exception')
                        runnable._status = 'failed'
                        exec_error.append(runnable._node_id, res)
                        # if this is a node for a running workflow, need to mark it as failed as well
                        #                        for proc in procs:
                        if isinstance(runnable, dummy_node) and hasattr(runnable, '_pending_workflow'):
                            for proc in manager.procs:
                                if proc is None:
                                    continue
                                if proc.is_pending() and hasattr(proc.step, '_pending_workflow') \
                                    and proc.step._pending_workflow == runnable._pending_workflow:
                                    proc.set_status('failed')
                        prog.update(1)
                    elif '__step_name__' in res:
                        env.logger.debug(f'{i_am()} receive step result ')
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
                        dag.update_step(runnable, env.sos_dict['__step_input__'].targets(),
                            env.sos_dict['__step_output__'].targets(),
                            env.sos_dict['__step_depends__'].targets())
                        runnable._status = 'completed'
                        prog.update(1)
                    elif '__workflow_id__' in res:
                        # result from a workflow
                        # the worker process has been returned to the pool, now we need to
                        # notify the step that is waiting for the result
                        env.logger.debug(f'{i_am()} receive workflow result')
                        for proc in manager.procs:
                            if proc is None:
                                continue
                            if proc.in_status('workflow_pending') and proc.step._pending_workflow == res['__workflow_id__']:
                                proc.pipe.send(res)
                                proc.set_status('running')
                                break
                    else:
                        raise RuntimeError(f'Unrecognized response from a step: {res}')

                # remove None
                manager.cleanup()

                # step 2: check if some jobs are done
                for proc_idx, proc in enumerate(manager.procs):
                    # if a job is pending, check if it is done.
                    if proc.in_status('task_pending'):
                        res = proc.step._host.check_status(proc.step._pending_tasks)
                        #env.logger.warning(res)
                        if any(x in ('aborted', 'failed', 'signature-mismatch') for x in res):
                            for t, s in zip(proc.step._pending_tasks, res):
                                if s in ('aborted', 'failed', 'signature-mismatch') and not (hasattr(proc.step, '_killed_tasks') and t in proc.step._killed_tasks):
                                    #env.logger.warning(f'{t} ``{s}``')
                                    if not hasattr(proc.step, '_killed_tasks'):
                                        proc.step._killed_tasks = {t}
                                    else:
                                        proc.step._killed_tasks.add(t)
                            if all(x in ('completed', 'aborted', 'failed', 'signature-mismatch') for x in res):
                                # we try to get .err .out etc even when jobs are failed.
                                task_status = proc.step._host.retrieve_results(proc.step._pending_tasks)
                                proc.pipe.send(task_status)
                                proc.set_status('failed')
                                status = [('completed', len([x for x in res if x=='completed'])),
                                    ('failed', len([x for x in res if x=='failed'])),
                                    ('aborted', len([x for x in res if x=='aborted'])),
                                    ('result mismatch', len([x for x in res if x=='signature-mismatch']))]
                                raise RuntimeError(', '.join([f'{y} job{"s" if y > 1 else ""} {x}' for x, y in status if y > 0]))
                        if any(x in ('pending', 'submitted', 'running') for x in res):
                            continue
                        elif all(x == 'completed' for x in res):
                            env.logger.debug(
                                f'Proc {proc_idx} puts results for {" ".join(proc.step._pending_tasks)} from step {proc.step._node_id}')
                            res = proc.step._host.retrieve_results(proc.step._pending_tasks)
                            proc.pipe.send(res)
                            proc.step._pending_tasks = []
                            proc.set_status('running')
                        else:
                            raise RuntimeError(f'Job returned with status {res}')

                # step 3: check if there is room and need for another job
                while True:
                    #env.logger.error('{} {}'.format(i_am(), [x.status() for x in procs]))
                    if manager.all_busy():
                        break
                    #
                    # if steps from child nested workflow?
                    if self.step_queue:
                        step_id, step_param = self.step_queue.popitem()
                        section, context, shared, args, run_mode, sig_mode, verbosity, pipe = step_param
                        # run it!
                        runnable = dummy_node()
                        runnable._node_id = step_id
                        runnable._status = 'running'
                        runnable._from_nested = True
                        runnable._child_pipe = pipe

                        env.logger.debug(
                            f'{i_am()} sends {section.step_name()} from step queue with args {args} and context {context}')

                        manager.execute(runnable, config=self.config, args=self.args,
                                spec = ('step', section, context, shared, args, run_mode, sig_mode, verbosity))
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
                                    raise ValueError(
                                        f'Unacceptable value for parameter shared: {section.options["shared"]}')
                        else:
                            raise ValueError(f'Unacceptable value for parameter shared: {section.options["shared"]}')
                        shared.update({x: env.sos_dict[x] for x in svars if x in env.sos_dict and pickleable(env.sos_dict[x], x)})

                    if '__workflow_sig__' in env.sos_dict:
                        runnable._context['__workflow_sig__'] = env.sos_dict['__workflow_sig__']

                    if not nested:
                        env.logger.debug(f'{i_am()} execute {section.md5} from DAG')
                        manager.execute(runnable, config=self.config, args=self.args,
                                spec=('step', section, runnable._context, shared, self.args,
                                    env.config['run_mode'], env.config['sig_mode'], env.verbosity))
                    else:
                        # send the step to the parent
                        step_id = uuid.uuid4()
                        env.logger.debug(
                            f'{i_am()} send step {section.step_name()} to master with args {self.args} and context {runnable._context}')
                        parent_pipe.send(f'step {step_id}')
                        q = mp.Pipe()
                        parent_pipe.send((section, runnable._context, shared, self.args, env.config['run_mode'], env.config['sig_mode'], env.verbosity, q[1]))
                        # this is a real step
                        manager.add_placeholder_worker(runnable, q[0])

                if manager.all_done_or_failed():
                    break

                # if -W is specified, or all task queues are not wait
                elif all(x.in_status('task_pending') for x in manager.procs) and \
                        (env.config['wait_for_task'] is False or \
                        (env.config['wait_for_task'] is None and Host.not_wait_for_tasks())):
                    # if all jobs are pending, let us check if all jbos have been submitted.
                    pending_tasks = []
                    running_tasks = []
                    for n in [x.step for x in manager.procs]:
                        p, r = n._host._task_engine.get_tasks()
                        pending_tasks.extend(p)
                        running_tasks.extend([(n._host.alias, x) for x in r])
                    if not pending_tasks and running_tasks:
                        env.logger.trace(f'Exit with {len(running_tasks)} running tasks: {running_tasks}')
                        raise PendingTasks(running_tasks)
                else:
                    time.sleep(0.1)
        except KeyboardInterrupt:
            if exec_error.errors:
                failed_steps, pending_steps = dag.pending()
                if pending_steps:
                    sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in pending_steps]
                    exec_error.append(self.workflow.name,
                                      RuntimeError(
                                          f'{len(sections)} pending step{"s" if len(sections) > 1 else ""}: {", ".join(sections)}'))
                    raise exec_error
            else:
                raise
        except PendingTasks as e:
            self.record_quit_status(e.tasks)
            wf_result['pending_tasks'] = [x[1] for x in running_tasks]
            env.logger.info(f'Workflow {self.workflow.name} (ID={self.md5}) exits with {len(e.tasks)} running tasks')
            for task in e.tasks:
                env.logger.info(task[1])
            # close all processes
        except Exception as e:
            manager.terminate(brutal=True)
            raise e
        finally:
            if not nested:
                manager.terminate()
            prog.close()
        #
        if exec_error.errors:
            failed_steps, pending_steps = dag.pending()
            #if failed_steps:
                #sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in failed_steps]
                #exec_error.append(self.workflow.name,
                #    RuntimeError('{} failed step{}: {}'.format(len(sections),
                #        's' if len(sections) > 1 else '', ', '.join(sections))))
            if pending_steps:
                sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in pending_steps]
                exec_error.append(self.workflow.name,
                                  RuntimeError(
                                      f'{len(sections)} pending step{"s" if len(sections) > 1 else ""}: {", ".join(sections)}'))
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
                env.logger.warning(f'Failed to clear workflow status file: {e}')
            self.save_workflow_signature(dag)
            env.logger.info(f'Workflow {self.workflow.name} (ID={self.md5}) is executed successfully.')
        else:
            # exit with pending tasks
            pass
        wf_result['shared'] = {x:env.sos_dict[x] for x in self.shared.keys() if x in env.sos_dict}
        if parent_pipe:
            parent_pipe.send(wf_result)
        else:
            return wf_result
