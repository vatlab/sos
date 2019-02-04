#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import base64
import os
import subprocess
import sys
import time
import uuid
import zmq

from collections import defaultdict, Sequence
from io import StringIO
from typing import Any, Dict, List, Optional, Tuple, Union
from threading import Event

from ._version import __version__
from .dag import SoS_DAG, SoS_Node
from .eval import SoS_exec
from .hosts import Host
from .parser import SoS_Step, SoS_Workflow
from .pattern import extract_pattern
from .workflow_report import render_report
from .controller import Controller, connect_controllers, disconnect_controllers
from .section_analyzer import analyze_section
from .targets import (BaseTarget, RemovedTarget, UnavailableLock,
                      UnknownTarget, file_target, path, paths,
                      sos_step, sos_targets, sos_variable, textMD5,
                      named_output)
from .utils import (Error, WorkflowDict, env, get_traceback,
                    load_config_files, pickleable, short_repr)
from .workers import SoS_Worker

__all__ = []

# try:
#     # https://github.com/pytest-dev/pytest-cov/issues/139
#     from pytest_cov.embed import cleanup_on_sigterm
#     cleanup_on_sigterm()
# except:
#     pass


class ExecuteError(Error):
    """An exception to collect exceptions raised during run time so that
    other branches of the DAG would continue if some nodes fail to execute."""

    def __init__(self, workflow: str) -> None:
        Error.__init__(self)
        self.workflow = workflow
        self.errors = []
        self.traces = []
        self.args = (workflow, )

    def append(self, line: str, error: Exception) -> None:
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


class dummy_node:
    # a dummy node object to store information of node passed
    # from nested workflow
    def __init__(self) -> None:
        pass


class ProcInfo(object):
    def __init__(self, worker: SoS_Worker, socket, step: Union[SoS_Node, dummy_node]) -> None:
        self.worker = worker
        self.socket = socket
        self.step = step

    def set_status(self, status: str) -> None:
        self.step._status = status

    def in_status(self, status: str) -> bool:
        return self.step._status == status

    def status(self):
        return self.step._status

    def is_pending(self) -> bool:
        return self.step._status.endswith('_pending')


class ExecutionManager(object):
    # this class managers workers and their status ...
    def __init__(self, max_workers: int) -> None:
                #
        # running processes. It consisists of
        #
        # [ [proc, queue], socket, node]
        #
        # where:
        #   proc, queue: process, which is None for the nested workflow.
        #   socket: socket to get information from workers
        #   node: node that is being executed, which is a dummy node
        #       created on the fly for steps passed from nested workflow
        #
        self.procs = []

        # process pool that is used to pool temporarily unused processed.
        self.pool = []
        self.max_workers = max_workers

    def execute(self, runnable: Union[SoS_Node, dummy_node], config: Dict[str, Any], args: Any, spec: Any) -> None:
        if not self.pool:
            socket = env.zmq_context.socket(zmq.PAIR)
            port = socket.bind_to_random_port('tcp://127.0.0.1')
            worker = SoS_Worker(port=port, config=config, args=args)
            worker.start()
        else:
            # get worker, q and runnable is not needed any more
            pi = self.pool.pop(0)
            worker = pi.worker
            socket = pi.socket

        # we need to report number of active works, plus master process itself
        env.controller_push_socket.send_pyobj(
            ['nprocs', self.num_active() + 1])

        socket.send_pyobj(spec)
        self.procs.append(
            ProcInfo(worker=worker, socket=socket, step=runnable))

    def add_placeholder_worker(self, runnable, socket):
        runnable._status = 'step_pending'
        self.procs.append(ProcInfo(worker=None, socket=socket, step=runnable))

    def num_active(self) -> int:
        return len([x for x in self.procs if x and not x.is_pending()
                    and not x.in_status('failed')])

    def all_busy(self) -> bool:
        return self.num_active() >= self.max_workers

    def all_done(self) -> bool:
        return not self.procs or all(x is None for x in self.procs)

    # def all_failed(self) -> bool:
    #     return all(x.in_status('failed') for x in self.procs)

    def mark_idle(self, idx: int) -> None:
        self.pool.append(self.procs[idx])
        self.procs[idx] = None

    def cleanup(self) -> None:
        self.procs = [x for x in self.procs if x is not None]

    def terminate(self, brutal: bool = False) -> None:
        self.cleanup()
        if not brutal:
            for proc in self.procs + self.pool:
                proc.socket.send_pyobj(None)
                proc.socket.LINGER = 0
                proc.socket.close()
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

    def __init__(self, workflow: Optional[SoS_Workflow] = None, args: Optional[Any] = None, shared: None = None,
                 config: Optional[Dict[str, Any]] = {}) -> None:
        self.workflow = workflow
        self.args = [] if args is None else args
        if '__args__' not in self.args:
            # if there is __args__, this is a nested workflow and we do not test this.
            for idx, arg in enumerate(self.args):
                wf_pars = self.workflow.parameters().keys()
                if isinstance(arg, str) and arg.startswith('--'):
                    if not wf_pars:
                        raise ValueError(
                            f'Undefined parameter {arg[2:]} for command line argument "{" ".join(args[idx:])}".')
                    pars = [arg[2:], arg[2:].replace('-', '_').split('=')[0]]
                    if arg[2:].startswith('no-'):
                        pars.extend(
                            [arg[5:], arg[5:].replace('-', '_').split('=')[0]])
                    if not any(x in wf_pars for x in pars):
                        raise ValueError(
                            f'Undefined parameter {arg[2:]} for command line argument "{" ".join(args[idx:])}". Acceptable parameters are: {", ".join(wf_pars)}')

        self.shared = {} if shared is None else shared
        env.config.update(config)
        if env.config['config_file'] is not None:
            env.config['config_file'] = os.path.abspath(
                os.path.expanduser(env.config['config_file']))
        #
        # if the executor is not called from command line, without sigmode setting
        if env.config['sig_mode'] is None:
            env.config['sig_mode'] = 'default'
        # interactive mode does not pass workflow
        self.md5 = self.calculate_md5() if self.workflow else '0'

        env.config['workflow_id'] = self.md5
        env.sos_dict.set('workflow_id', self.md5)

    def write_workflow_info(self):
        # if this is the outter most workflow, master)id should have =
        # not been set so we set it for all other workflows
        workflow_info = {
            'name': self.workflow.name,
            'start_time': time.time(),
        }
        workflow_info['command_line'] = subprocess.list2cmdline(
            [os.path.basename(sys.argv[0])] + sys.argv[1:])
        workflow_info['project_dir'] = os.getcwd()
        workflow_info['script'] = base64.b64encode(
            self.workflow.content.text().encode()).decode('ascii')
        workflow_info['master_id'] = env.config['master_id']
        env.signature_req_socket.send_pyobj(['workflow', 'clear'])
        env.signature_req_socket.recv_pyobj()
        env.signature_push_socket.send_pyobj(
            ['workflow', 'workflow', self.md5, repr(workflow_info)])
        if env.config['exec_mode'] == 'slave':
            env.tapping_listener_socket.send_pyobj(
                {'msg_type': 'workflow_status',
                 'data': {
                     'cell_id': env.config['slave_id'],
                     'workflow_id': self.md5,
                     'workflow_name': self.workflow.name,
                     'start_time': time.time(),
                     'status': 'running'
                 }
                 }
            )

    def run(self, targets: Optional[List[str]]=None, mode=None) -> Dict[str, Any]:
        #
        env.zmq_context = zmq.Context()

        # if this is the executor for the master workflow, start controller
        env.config['master_id'] = self.md5
        #
        # control panel in a separate thread, connected by zmq socket
        ready = Event()
        self.controller = Controller(ready)
        self.controller.start()
        # wait for the thread to start with a signature_req saved to env.config
        ready.wait()

        connect_controllers(env.zmq_context)

        try:
            succ = True
            return self.run_as_master(targets=targets, mode=mode)
        except:
            succ = False
            raise
        finally:
            # end progress bar when the master workflow stops
            env.logger.trace(f'Stop controller from {os.getpid()}')
            env.controller_req_socket.send_pyobj(['done', succ])
            env.controller_req_socket.recv()
            env.logger.trace('disconntecting master')
            # if the process is failed, some workers might be killed, resulting
            # in nonresponseness from the master, and the socket context cannot
            # be killed in this case.
            disconnect_controllers(env.zmq_context if succ else None)
            self.controller.join()
            # when the run() function is called again, the controller
            # thread will be start again.
            env.config['master_id'] = None

    def calculate_md5(self) -> str:
        with StringIO() as sig:
            for step in self.workflow.sections + self.workflow.auxiliary_sections:
                sig.write(f'{step.step_name()}: {step.md5}\n')
            sig.write(f'{self.args}\n')
            return textMD5(sig.getvalue())[:16]

    def reset_dict(self) -> None:
        env.sos_dict = WorkflowDict()
        self.init_dict()

    def init_dict(self) -> None:
        env.parameter_vars.clear()

        env.sos_dict.set('workflow_id', self.md5)
        env.sos_dict.set('master_id', env.config['master_id'])
        env.sos_dict.set('__args__', self.args)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        env.sos_dict.set('__step_output__', sos_targets([]))

        # load configuration files
        load_config_files(env.config['config_file'])

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)

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

    def analyze_auxiliary_step(self, section):
        res = analyze_section(section, vars_and_output_only=True)
        environ_vars = res['environ_vars']
        signature_vars = res['signature_vars']
        changed_vars = res['changed_vars']
        # parameters, if used in the step, should be considered environmental
        environ_vars |= env.parameter_vars & signature_vars

        # add shared to targets
        if res['changed_vars']:
            if 'provides' in section.options:
                if isinstance(section.options['provides'], str):
                    section.options.set(
                        'provides', [section.options['provides']])
            else:
                section.options.set('provides', [])
            #
            section.options.set('provides',
                                section.options['provides'] + [sos_variable(var) for var in changed_vars])


        if not res['step_output'].unspecified():
             section._autoprovides = res['step_output']

    def _build_target_map(self, sections):
        self._target_map = defaultdict(list)
        self._target_patterns = defaultdict(list)

        for step in sections:
            if not hasattr(step, '_analyzed'):
                ret = self.analyze_auxiliary_step(step)
                step._analyzed = True
            # a step first provides sos_step
            for name, index, _ in step.names:
                self._target_map[sos_step(name)].append(step)
                if index is not None:
                    self._target_map[sos_step(f'{name}_{index}')].append(step)
            # named output
            if 'namedprovides' in step.options:
                for x in step.options['namedprovides']:
                    self._target_map[named_output(x)].append(step)
            #
            if hasattr(step, '_autoprovides'):
                for x in step._autoprovides:
                    # x must be BaseTarget (file_target or others)
                    self._target_map[x].append(step)
            #
            if 'provides' in step.options:
                patterns = step.options['provides']
                # str, BaseTarget, path
                # sos_targets, Sequence, paths
                if isinstance(patterns, str):
                    if '{' in patterns and '}' in patterns:
                        self._target_patterns[patterns].append(step)
                    else:
                        self._target_map[file_target(patterns)].append(step)
                elif isinstance(patterns, path):
                    self._target_map[file_target(patterns)].append((step, True))
                elif isinstance(patterns, sos_targets):
                    for x in patterns:
                        self._target_map[x].append(step)
                elif isinstance(patterns, paths):
                    for x in patterns:
                        self._target_map[file_target(x)].append(step)
                elif isinstance(patterns, BaseTarget):
                    self._target_map[patterns].append(step)
                elif isinstance(patterns, Sequence):
                    for pattern in patterns:
                        if isinstance(pattern, str):
                            if '{' in pattern and '}' in pattern:
                                self._target_patterns[pattern].append(step)
                            else:
                                self._target_map[file_target(pattern)].append(step)
                        elif isinstance(pattern, path):
                            self._target_map[file_target(pattern)].append(step)
                        elif isinstance(pattern, sos_targets):
                            for x in pattern:
                                self._target_map[x].append(step)
                        elif isinstance(pattern, paths):
                            for x in pattern:
                                self._target_map[file_target(x)].append(step)
                        elif isinstance(pattern, BaseTarget):
                            self._target_map[pattern].append(step)
                        else:
                            raise ValueError(f'Unacceptable value for option pattern {patterns}')
                else:
                    raise ValueError(f'Unacceptable value for option pattern {patterns}')

    def match(self, target: BaseTarget) -> Union[Dict[str, str], bool]:
        if not hasattr(self, '_target_map'):
            self._build_target_map(self.workflow.auxiliary_sections)
        if target in self._target_map:
            res = self._target_map[target]
            if len(res) > 1:
                self._target_map[target] = list(set(res))
                res = self._target_map[target]
            return res
        if not isinstance(target, file_target) or not self._target_patterns:
            return False
        # try pattern?
        for pattern, steps in self._target_patterns.items():
            if len(steps) > 1:
                raise RuntimeErrors(f'Multiple steps providing the same pattern.')
            # if this is a regular string
            res = extract_pattern(pattern, [str(target)])
            if res and not any(None in x for x in res.values()):
                return [(steps[0], {x: y[0] for x, y in res.items()})]
        return False

    def resolve_dangling_targets(self, dag: SoS_DAG, targets: Optional[sos_targets]=None) -> int:
        '''Feed dangling targets with their dependncies from auxiliary steps,
        optionally add other targets'''
        resolved = 0
        while True:
            added_node = 0
            # first resolve missing
            dangling_targets = dag.dangling(targets)[0]
            if dangling_targets:
                env.logger.debug(
                    f'Resolving {dangling_targets} objects from {dag.number_of_nodes()} nodes')
            # find matching steps
            # check auxiliary steps and see if any steps provides it
            remaining_targets = dangling_targets
            while remaining_targets:
                target = remaining_targets.pop()
                mo = self.match(target)
                if not mo:
                    #
                    # if no step produces the target, it is possible that it is an indexed step
                    # so the execution of its previous steps would solves the dependency
                    #
                    # find all the nodes that depends on target
                    nodes = dag._all_depends_files[target]
                    for node in nodes:
                        # if this is an index step... simply let it depends on previous steps
                        if node._node_index is not None:
                            indexed = [x for x in dag.nodes() if x._node_index is not None and x._node_index <
                                       node._node_index and not x._output_targets.valid()]
                            indexed.sort(key=lambda x: x._node_index)
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
                            if node._input_targets.valid():
                                node._input_targets = sos_targets([])
                            if node._depends_targets.valid():
                                node._depends_targets = sos_targets([])
                        else:
                            raise RuntimeError(
                                f'No step to generate target {target}{dag.steps_depending_on(target, self.workflow)}')
                    if nodes:
                        resolved += 1
                    # cannot resolve, but the step is numerically indexed so we just wait
                    continue
                if len(mo) > 1:
                    # sos_step('a') could match to step a_1, a_2, etc, in this case we are adding a subworkflow
                    if isinstance(target, sos_step):
                        # create a new forward_workflow that is different from the master one
                        # get the step names
                        sections = sorted(mo, key=lambda x: x.step_name())
                        # this is only useful for executing auxiliary steps and
                        # might interfere with the step analysis
                        env.sos_dict.pop('__default_output__', None)
                        env.logger.debug(
                                f'Adding {len(sections)} steps to resolve target {target}')

                        n_added = self.add_forward_workflow(dag, sections, satisfies=target)

                        added_node += n_added
                        resolved += 1
                        # dag.show_nodes()
                        continue
                    else:
                        raise RuntimeError(
                            f'Multiple steps {", ".join(x.step_name() for x in mo)} to generate target {target}')
                #
                # m0[0] can be a tuple
                #   section, context
                # or just a section
                #   section
                if isinstance(mo[0], tuple):
                    # for auxiliary, we need to set input and output, here
                    # now, if the step does not provide any alternative (e.g. no variable generated
                    # from patten), we should specify all output as output of step. Otherwise the
                    # step will be created for multiple outputs. issue #243
                    if mo[0][1]:
                        env.sos_dict['__default_output__'] = sos_targets(target)
                    else:
                        env.sos_dict['__default_output__'] = sos_targets(
                            section.options['provides'])
                    n_added, _ = self.add_backward_step(dag, section=mo[0][0],
                        context = mo[0][1] if isinstance(mo[0][1], dict) else {},
                        target=target)
                else:
                    env.sos_dict['__default_output__'] = sos_targets(target)
                    n_added, _ = self.add_backward_step(dag, section=mo[0], context = {},
                    target=target)

                added_node += n_added
                resolved += n_added
                # this could be made more efficient by removing step output directly
                # but it is inefficient to remove elements from list
                remaining_targets = dag.dangling(remaining_targets)[0]

            # for existing targets that are not in DAG
            traced = [x for x in targets if x.traced]
            if not env.config['trace_existing'] and not traced:
                if added_node == 0:
                    break
                else:
                    continue

            node_added = False
            existing_targets = set(dag.dangling(targets)[1]) if env.config['trace_existing'] else traced

            remaining_targets = existing_targets
            while remaining_targets:
                target = remaining_targets.pop()
                # now we need to build DAG for existing...
                mo = self.match(target)
                if not mo:
                    # this is ok, this is just an existing target, no one is designed to
                    # generate it.
                    continue
                if len(mo) > 1:
                    # this is not ok.
                    raise RuntimeError(
                        f'Multiple steps {", ".join(x.step_name() for x in mo)} to generate target {target}')

                if isinstance(mo[0], tuple):
                    # for auxiliary, we need to set input and output, here
                    # now, if the step does not provide any alternative (e.g. no variable generated
                    # from patten), we should specify all output as output of step. Otherwise the
                    # step will be created for multiple outputs. issue #243
                    if mo[0][1]:
                        env.sos_dict['__default_output__'] = sos_targets(target)
                    else:
                        env.sos_dict['__default_output__'] = sos_targets(
                            section.options['provides'])
                    n_added, resolved_output = self.add_backward_step(dag, section=mo[0][0],
                        context = mo[0][1] if isinstance(mo[0][1], dict) else {},
                        target=target)
                else:
                    env.sos_dict['__default_output__'] = sos_targets(target)
                    n_added, resolved_output = self.add_backward_step(dag, section=mo[0], context = {},
                    target=target)

                node_added = True
                added_node += n_added
                resolved += n_added

                # adding one step can resolve many targets #1199
                if len(resolved_output) > 1:
                    remaining_targets -= set(resolved_output.targets)
                    # set(dag.dangling(remaining_targets)[1])

            if added_node == 0:
                break
        return resolved

    def add_forward_workflow(self, dag, sections, satisfies=None):
        '''Add a forward-workflow, return number of nodes added
        '''
        dag.new_forward_workflow()

        env.logger.debug(f'Adding mini-workflow with {len(sections)} to resolve {satisfies}')
        default_input: sos_targets = sos_targets([])
        for idx, section in enumerate(sections):
            #
            res = analyze_section(section, default_input)

            environ_vars = res['environ_vars']
            signature_vars = res['signature_vars']
            changed_vars = res['changed_vars']
            # parameters, if used in the step, should be considered environmental
            environ_vars |= env.parameter_vars & signature_vars

            # add shared to targets
            if res['changed_vars']:
                if 'provides' in section.options:
                    if isinstance(section.options['provides'], str):
                        section.options.set(
                            'provides', [section.options['provides']])
                else:
                    section.options.set('provides', [])
                #
                section.options.set('provides',
                                    section.options['provides'] + [sos_variable(var) for var in changed_vars])

            context = {'__signature_vars__': signature_vars,
                       '__environ_vars__': environ_vars,
                       '__changed_vars__': changed_vars,
                       '__dynamic_depends__': res['dynamic_depends'],
                       '__dynamic_input__': res['dynamic_input']
                       }

            # for nested workflow, the input is specified by sos_run, not None.
            if idx == 0:
                context['__step_output__'] = env.sos_dict['__step_output__']
            # can be the only step
            if idx == len(sections) - 1 and satisfies is not None:
                res['step_output'].extend(satisfies)

            dag.add_step(section.uuid,
                         section.step_name(),
                         idx,
                         res['step_input'],
                         res['step_depends'],
                         res['step_output'],
                         context=context)
            default_input = res['step_output']
        return len(sections)

    def add_backward_step(self, dag, section, context, target):
        # only one step, we need to process it # execute section with specified input
        #
        for k, v in context.items():
            env.sos_dict.set(k, v)

        # will become input, set to None
        env.sos_dict['__step_output__'] = sos_targets()
        #
        res = analyze_section(section, default_output=env.sos_dict['__default_output__'])
        if isinstance(target, sos_step) and target.target_name() != section.step_name():
            # sos_step target "name" can be matched to "name_10" etc so we will have to
            # ensure that the target is outputted from the "name_10" step.
            # This has been done in a more advanced case when an entire workflow is
            # added
            res['step_output'].extend(target)
        elif isinstance(target, named_output):
            # when a named_output is matched, we add all names from the step
            # to avoid the step to be added multiple times for different named_steps
            # from the same step. #1166
            res['step_output'].extend([named_output(x) for x in section.options['namedprovides']])

        # now, if we are adding a step that is part of a forward-style workflow
        # (with index), and the input of the step is unspecified, the it is possible
        # that the step depends on its previous steps and we will have to add these steps
        # as well. #1206
        if res['step_input'].unspecified() and section.index is not None:
            # there should be one step with that name
            prev_steps = [x for x in self.workflow.auxiliary_sections
                if x.name == section.name and x.index is not None and x.index <= section.index]
            self.add_forward_workflow(dag, prev_steps, target)
            return len(prev_steps), res['step_output']

        # add a single step
        # build DAG with input and output files of step
        env.logger.debug(
            f'Adding step {res["step_name"]} with output {short_repr(res["step_output"])} to resolve target {target}')

        context['__signature_vars__'] = res['signature_vars']
        context['__environ_vars__'] = res['environ_vars']
        context['__changed_vars__'] = res['changed_vars']
        context['__default_output__'] = env.sos_dict['__default_output__']
        context['__dynamic_depends__'] = res['dynamic_depends']
        context['__dynamic_input__'] = res['dynamic_input']

        # NOTE: If a step is called multiple times with different targets, it is much better
        # to use different names because pydotplus can be very slow in handling graphs with nodes
        # with identical names.
        node_name = section.step_name()
        if env.sos_dict["__default_output__"]:
            node_name += f' ({short_repr(env.sos_dict["__default_output__"])})'
        dag.add_step(section.uuid,
                        node_name, None,
                        res['step_input'],
                        res['step_depends'],
                        res['step_output'],
                        context=context)
        return 1, res['step_output']

    def initialize_dag(self, targets: Optional[List[str]] = [], nested: bool = False) -> SoS_DAG:
        '''Create a DAG by analyzing sections statically.'''
        self.reset_dict()

        dag = SoS_DAG(name=self.md5)
        targets = sos_targets(targets)

        self.add_forward_workflow(dag, self.workflow.sections)
        #
        if self.resolve_dangling_targets(dag, targets) == 0:
            if targets:
                raise UnknownTarget(
                    f'No step to generate target {targets}.')
        # now, there should be no dangling targets, let us connect nodes
        dag.build()

        # dag.show_nodes()
        # trim the DAG if targets are specified
        if targets:
            dag = dag.subgraph_from(targets)
        # check error
        cycle = dag.circular_dependencies()
        if cycle:
            raise RuntimeError(
                f'Circular dependency detected {cycle}. It is likely a later step produces input of a previous step.')
        dag.save(env.config['output_dag'])
        return dag

    def describe_completed(self):
        # return a string to summarize completed and skipped steps, substeps, and tasks
        res = []
        # if '__subworkflow_completed__' in self.completed and self.completed['__subworkflow_completed__']:
        #    res.append(f"{self.completed['__subworkflow_completed__']} completed subworkflow{'s' if self.completed['__subworkflow_completed__'] > 1 else ''}")
        # if '__subworkflow_skipped__' in self.completed and self.completed['__subworkflow_skipped__']:
        #    res.append(f"{self.completed['__subworkflow_skipped__']} skipped subworkflow{'s' if self.completed['__subworkflow_skipped__'] > 1 else ''}")
        if '__step_completed__' in self.completed and self.completed['__step_completed__']:
            res.append(
                f"{round(self.completed['__step_completed__'], 1)} completed step{'s' if self.completed['__step_completed__'] > 1 else ''}")
            if '__substep_completed__' in self.completed and self.completed['__substep_completed__'] and self.completed['__substep_completed__'] != self.completed['__step_completed__']:
                res.append(
                    f"{self.completed['__substep_completed__']} completed substep{'s' if self.completed['__substep_completed__'] > 1 else ''}")
        if '__step_skipped__' in self.completed and self.completed['__step_skipped__']:
            res.append(
                f"{round(self.completed['__step_skipped__'], 1)} ignored step{'s' if self.completed['__step_skipped__'] > 1 else ''}")
            if '__substep_skipped__' in self.completed and self.completed['__substep_skipped__'] and self.completed['__substep_skipped__'] != self.completed['__step_skipped__']:
                res.append(
                    f"{self.completed['__substep_skipped__']} ignored substep{'s' if self.completed['__substep_skipped__'] > 1 else ''}")
        if '__task_completed__' in self.completed and self.completed['__task_completed__']:
            res.append(
                f"{self.completed['__task_completed__']} completed task{'s' if self.completed['__task_completed__'] > 1 else ''}")
        if '__task_skipped__' in self.completed and self.completed['__task_skipped__']:
            res.append(
                f"{self.completed['__task_skipped__']} ignored task{'s' if self.completed['__task_skipped__'] > 1 else ''}")
        if len(res) > 1:
            return ', '.join(res[:-1]) + ' and ' + res[-1]
        elif len(res) == 1:
            return res[0]
        else:
            return 'no step executed'

    def step_completed(self, res, dag, runnable):
        # this function can be called by both master step and nested step
        for k, v in res['__completed__'].items():
            self.completed[k] += v
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
                # and node._node_index == runnable._node_index + 1:
                node._context.update(env.sos_dict.clone_selected_vars(
                    node._context['__signature_vars__'] | node._context['__environ_vars__']
                    | {'_input', '__step_output__', '__args__'}))
            node._context.update(svar)
            if node._status == 'target_pending':
                if all(x.target_exists('target') for x in node._pending_targets):
                    # in a master node, this _socket points to the step
                    # in a nested node, this _socket points to the parent socket
                    node._socket.send(b'target_resolved')
                    node._status = 'running'
        dag.update_step(runnable,
                        input_targets = env.sos_dict['__step_input__'],
                        output_targets = env.sos_dict['__step_output__'],
                        depends_targets = env.sos_dict['__step_depends__'])
        runnable._status = 'completed'
        dag.save(env.config['output_dag'])

    def handle_dependent_target(self, dag, targets, runnable) -> int:
        total_added = 0
        resolved = 0
        while True:
            added_node = 0
            # for depending targets... they already exist but we will add
            # nodes that generates them if available.
            node_added = False
            depending_targets = set(dag.dangling(targets)[1])
            for target in depending_targets:
                if node_added:
                    depending_targets = set(dag.dangling(targets)[1])
                    node_added = False
                if target not in depending_targets:
                    continue
                mo = self.match(target)
                if not mo:
                    # this is ok, this is just an existing target, no one is designed to
                    # generate it.
                    env.logger.info(f'{target} already exists')
                    continue
                if len(mo) > 1:
                    # this is not ok.
                    raise RuntimeError(
                        f'Multiple steps {", ".join(x.step_name() for x in mo)} to generate target {target}')

                if isinstance(mo[0], tuple):
                    if mo[0][1]:
                        env.sos_dict['__default_output__'] = sos_targets(target)
                    else:
                        env.sos_dict['__default_output__'] = sos_targets(
                            section.options['provides'])
                    n_added, _ = self.add_backward_step(dag, section=mo[0][0],
                        context = mo[0][1] if isinstance(mo[0][1], dict) else {},
                        target=target)
                else:
                    env.sos_dict['__default_output__'] = sos_targets(target)
                    n_added, _ = self.add_backward_step(dag, section=mo[0], context = {},
                    target=target)

                node_added = True
                added_node += n_added
                resolved += n_added

            total_added += added_node
            if added_node == 0:
                break
        #
        if total_added:
            if runnable._depends_targets.valid():
                runnable._depends_targets.extend(targets)
            for taget in targets:
                if runnable not in dag._all_depends_files[target]:
                    dag._all_depends_files[target].append(
                        runnable)
            dag.build()
            #
            cycle = dag.circular_dependencies()
            if cycle:
                raise RuntimeError(
                    f'Circular dependency detected {cycle}. It is likely a later step produces input of a previous step.')
            dag.save(env.config['output_dag'])
        return total_added

    def handle_unknown_target(self, target, dag, runnable):
        runnable._status = None
        dag.save(env.config['output_dag'])

        if isinstance(target, file_target) and (target + '.zapped').exists():
            (target + '.zapped').unlink()

        if dag.regenerate_target(target):
            # runnable._depends_targets.append(target)
            # dag._all_depends_files[target].append(runnable)
            dag.build()
            #
            cycle = dag.circular_dependencies()
            if cycle:
                raise RuntimeError(
                    f'Circular dependency detected {cycle} after regeneration. It is likely a later step produces input of a previous step.')

        else:
            if self.resolve_dangling_targets(dag, sos_targets(target)) == 0:
                raise RuntimeError(
                    f'Failed to regenerate or resolve {target}{dag.steps_depending_on(target, self.workflow)}.')
            if runnable._depends_targets.valid():
                runnable._depends_targets.extend(target)
            if runnable not in dag._all_depends_files[target]:
                dag._all_depends_files[target].append(
                    runnable)
            dag.build()
            #
            cycle = dag.circular_dependencies()
            if cycle:
                raise RuntimeError(
                    f'Circular dependency detected {cycle}. It is likely a later step produces input of a previous step.')
        dag.save(env.config['output_dag'])

    def handle_unavailable_lock(self, res, dag, runnable):
        runnable._status = 'signature_pending'
        dag.save(env.config['output_dag'])
        runnable._signature = (res.output, res.sig_file)
        section = self.workflow.section_by_id(
            runnable._step_uuid)
        env.logger.info(
            f'Waiting on another process for step {section.step_name(True)}')

    def get_shared_vars(self, option):
        if isinstance(option, str):
            svars = [option]
        elif isinstance(option, dict):
            svars = option.keys()
        elif isinstance(option, Sequence):
            svars = []
            for x in option:
                if isinstance(x, str):
                    svars.append(x)
                elif isinstance(x, dict):
                    svars.extend(x.keys())
                else:
                    raise ValueError(
                        f'Unacceptable value for parameter shared: {option}')
        else:
            raise ValueError(
                f'Unacceptable value for parameter shared: {option}')
        return {x: env.sos_dict[x] for x in svars if x in env.sos_dict and pickleable(env.sos_dict[x], x)}

    def finalize_and_report(self):
        # remove task pending status if the workflow is completed normally
        if self.workflow.name != 'scratch':
            if self.completed["__step_completed__"] == 0:
                sts = 'ignored'
            elif env.config["run_mode"] == 'dryrun':
                sts = 'tested successfully'
            else:
                sts = 'executed successfully'
            env.logger.info(
                f'Workflow {self.workflow.name} (ID={self.md5}) is {sts} with {self.describe_completed()}.')
        if env.config['output_dag']:
            env.logger.info(
                f"Workflow DAG saved to {env.config['output_dag']}")
        workflow_info = {
            'end_time': time.time(),
            'stat': dict(self.completed),
        }
        if env.config['output_dag'] and env.config['master_id'] == self.md5:
            workflow_info['dag'] = env.config['output_dag']
        env.signature_push_socket.send_pyobj(
            ['workflow', 'workflow', self.md5, repr(workflow_info)])
        if env.config['master_id'] == env.sos_dict['workflow_id'] and env.config['output_report']:
            # if this is the outter most workflow
            render_report(env.config['output_report'],
                          env.sos_dict['workflow_id'])
        if env.config['run_mode'] == 'dryrun':
            env.signature_req_socket.send_pyobj(
                ['workflow', 'placeholders', env.sos_dict['workflow_id']])
            for filename in env.signature_req_socket.recv_pyobj():
                try:
                    if os.path.getsize(file_target(filename)) == 0:
                        file_target(filename).unlink()
                        env.logger.debug(f'Remove placeholder {filename}')
                except Exception as e:
                    env.logger.trace(
                        f'Failed to remove placeholder {filename}: {e}')

    def run_as_master(self, targets=None, mode=None) -> Dict[str, Any]:
        self.completed = defaultdict(int)

        self.write_workflow_info()

        def i_am():
            return 'Master'

        self.reset_dict()

        env.config['run_mode'] = env.config.get(
            'run_mode', 'run') if mode is None else mode

        # passing run_mode to SoS dict so that users can execute blocks of
        # python statements in different run modes.
        env.sos_dict.set('run_mode', env.config['run_mode'])

        wf_result = {'__workflow_id__': self.md5, 'shared': {}}

        try:
            if env.config['output_dag'] and os.path.isfile(env.config['output_dag']):
                os.unlink(env.config['output_dag'])
        except Exception as e:
            env.logger.warning(
                f'Failed to remove existing DAG file {env.config["output_dag"]}: {e}')

        # process step of the pipelinp
        try:
            dag = self.initialize_dag(targets=targets)
        except UnknownTarget as e:
            # if the names cannot be found, try to see if they are named_output
            try:
                named_targets = [named_output(x) for x in targets]
                dag = self.initialize_dag(targets=named_targets)
            except UnknownTarget as e:
                raise RuntimeError(f'No step to generate target {targets}')

        # manager of processes
        manager = ExecutionManager(env.config['max_procs'])
        #
        # steps sent and queued from the nested workflow
        # they will be executed in random but at a higher priority than the steps
        # on the master process.
        self.step_queue = {}
        try:
            exec_error = ExecuteError(self.workflow.name)
            while True:
                # step 1: check existing jobs and see if they are completed
                for idx, proc in enumerate(manager.procs):
                    if proc is None:
                        continue

                    # echck if there is any message from the socket
                    if not proc.socket.poll(0):
                        continue

                    # receieve something from the pipe
                    res = proc.socket.recv_pyobj()
                    runnable = proc.step
                    # if this is NOT a result, rather some request for task, step, workflow etc
                    if isinstance(res, list):
                        if res[0] == 'tasks':
                            env.logger.debug(
                                f'{i_am()} receives task request {res}')
                            host = res[1]
                            if host == '__default__':
                                if 'default_queue' in env.config:
                                    host = env.config['default_queue']
                                else:
                                    host = 'localhost'
                            try:
                                new_tasks = res[2:]
                                runnable._host = Host(host)
                                if hasattr(runnable, '_pending_tasks'):
                                    runnable._pending_tasks.extend(new_tasks)
                                else:
                                    runnable._pending_tasks = new_tasks
                                for task in new_tasks:
                                    runnable._host.submit_task(task)
                                runnable._status = 'task_pending'
                                dag.save(env.config['output_dag'])
                                env.logger.trace('Step becomes task_pending')
                            except Exception as e:
                                proc.socket.send_pyobj(
                                    {x: {'ret_code': 1, 'task': x, 'output': {}, 'exception': e} for x in new_tasks})
                                env.logger.error(e)
                                proc.set_status('failed')
                        elif res[0] == 'missing_target':
                            # the target that is missing from the running step
                            missed = res[1]
                            if hasattr(runnable, '_from_nested'):
                                # if the step is from a subworkflow, then the missing target
                                # should be resolved by the nested workflow
                                runnable._child_socket.send_pyobj(res)
                                reply = runnable._child_socket.recv_pyobj()
                                if reply: # if the target is resolvable in nested workflow
                                    runnable._status = 'target_pending'
                                    runnable._pending_targets = [missed]
                                    # tell the step that the target is resolved and it can continue
                                    runnable._socket = proc.socket
                                else:
                                    # otherwise say the target cannot be resolved
                                    proc.socket.send(b'')
                                    proc.set_status('failed')
                            else:
                                # if the missing target is from master, resolve from here
                                try:
                                    self.handle_unknown_target(missed, dag, runnable)
                                    runnable._status = 'target_pending'
                                    runnable._pending_targets = [missed]
                                    runnable._socket = proc.socket
                                except Exception as e:
                                    env.logger.error(e)
                                    proc.socket.send(b'')
                                    proc.set_status('failed')
                        elif res[0] == 'dependent_target':
                            # The target might be dependent on other steps and we
                            # are trying to extend the DAG to verify the target
                            # if possible. It does not matter if the DAG cannot be
                            # extended.
                            if hasattr(runnable, '_from_nested'):
                                # if the step is from a subworkflow, then the missing target
                                # should be resolved by the nested workflow
                                runnable._child_socket.send_pyobj(res)
                                reply = runnable._child_socket.recv_pyobj()
                                if reply:
                                    # if there are dependent steps, the current step
                                    # has to wait
                                    runnable._status = 'target_pending'
                                    runnable._pending_targets = res[1:]
                                    # tell the step that the target is resolved and it can continue
                                    runnable._socket = proc.socket
                                else:
                                    # otherwise there is no target to verify
                                    # and we just continue
                                    proc.socket.send(b'target_resolved')
                            else:
                                # if the missing target is from master, resolve from here
                                reply = self.handle_dependent_target(dag, sos_targets(res[1:]), runnable)
                                if reply:
                                    runnable._status = 'target_pending'
                                    runnable._pending_targets = res[1:]
                                    runnable._socket = proc.socket
                                else:
                                    proc.socket.send(b'target_resolved')
                        elif res[0] == 'step':
                            # step sent from nested workflow
                            step_id = res[1]
                            step_params = res[2:]
                            env.logger.debug(
                                f'{i_am()} receives step request {step_id} with args {step_params[3]}')
                            self.step_queue[step_id] = step_params
                        elif res[0] == 'workflow':
                            workflow_ids, wfs, targets, args, shared, config = res[1:]
                            # receive the real definition
                            env.logger.debug(
                                f'{i_am()} receives workflow request {workflow_ids}')

                            # now we would like to find a worker and
                            if hasattr(runnable, '_pending_workflows'):
                                runnable._pending_workflows.extend(workflow_ids)
                            else:
                                runnable._pending_workflows = workflow_ids
                            runnable._status = 'workflow_pending'
                            dag.save(env.config['output_dag'])

                            for wid, wf in zip(workflow_ids, wfs):
                                # for each subworkflow, create a dummy node to track
                                # its status. Note that this dummy_node does not have
                                # the _from_nested flag.
                                wfrunnable = dummy_node()
                                wfrunnable._node_id = wid
                                wfrunnable._status = 'workflow_running_pending'
                                dag.save(env.config['output_dag'])
                                wfrunnable._pending_workflows = [wid]
                                #
                                manager.execute(wfrunnable, config=config, args=args,
                                                spec=('workflow', wid, wf, targets, args, shared, config))
                        else:
                            raise RuntimeError(
                                f'Unexpected value from step {short_repr(res)}')
                        continue
                    elif isinstance(res, str):
                        raise RuntimeError(
                            f'Unexpected value from step {short_repr(res)}')

                    # if we does get the result, we send the process to pool

                    manager.mark_idle(idx)

                    env.logger.debug(
                        f'{i_am()} receive a result {short_repr(res)}')
                    if hasattr(runnable, '_from_nested'):
                        # if the runnable is from nested, we will need to send the result back
                        # to the nested workflow
                        env.logger.debug(f'{i_am()} send res to nested')
                        runnable._status = 'completed'
                        dag.save(env.config['output_dag'])
                        runnable._child_socket.send_pyobj(res)
                        # this is a onetime use socket that passes results from
                        # nested workflow to master
                        runnable._child_socket.LINGER = 0
                        runnable._child_socket.close()
                    elif isinstance(res, UnavailableLock):
                        self.handle_unavailable_lock(res, dag, runnable)
                    elif isinstance(res, RemovedTarget):
                        self.handle_unknown_target(res.target, dag, runnable)
                    # if the job is failed
                    elif isinstance(res, Exception):
                        env.logger.debug(f'{i_am()} received an exception')
                        # env.logger.error(res)
                        runnable._status = 'failed'
                        dag.save(env.config['output_dag'])
                        exec_error.append(runnable._node_id, res)
                        raise exec_error
                    elif '__step_name__' in res:
                        env.logger.debug(f'{i_am()} receive step result ')
                        self.step_completed(res, dag, runnable)
                    elif '__workflow_id__' in res:
                        # result from a workflow
                        # the worker process has been returned to the pool, now we need to
                        # notify the step that is waiting for the result
                        env.logger.debug(f'{i_am()} receive workflow result')
                        # aggregate steps etc with subworkflows
                        for k, v in res['__completed__'].items():
                            self.completed[k] += v
                        # if res['__completed__']['__step_completed__'] == 0:
                        #    self.completed['__subworkflow_skipped__'] += 1
                        # else:
                        #    self.completed['__subworkflow_completed__'] += 1
                        for proc in manager.procs:
                            # do not care about dummy processes
                            if proc is None:
                                continue
                            if proc.in_status('workflow_pending') and res['__workflow_id__'] in proc.step._pending_workflows:
                                proc.step._pending_workflows.remove(res['__workflow_id__'])
                                if not proc.step._pending_workflows:
                                    proc.set_status('running')
                                proc.socket.send_pyobj(res)
                                break
                        dag.save(env.config['output_dag'])
                    else:
                        raise RuntimeError(
                            f'Unrecognized response from a step: {res}')

                # remove None
                manager.cleanup()

                # step 2: check if some jobs are done
                for proc_idx, proc in enumerate(manager.procs):
                    # if a job is pending, check if it is done.
                    if proc.in_status('task_pending'):
                        res = proc.step._host.check_status(
                            proc.step._pending_tasks)
                        # env.logger.warning(res)
                        if all(x in ('completed', 'aborted', 'failed') for x in res):
                            env.logger.debug(
                                f'Proc {proc_idx} puts results for {" ".join(proc.step._pending_tasks)} from step {proc.step._node_id}')
                            res = proc.step._host.retrieve_results(
                                proc.step._pending_tasks)
                            proc.socket.send_pyobj(res)
                            proc.step._pending_tasks = []
                            proc.set_status('running')
                            #proc.set_status('failed')
                        elif any(x in ('new', 'pending', 'submitted', 'running') for x in res):
                            continue
                        elif 'missing' in res:
                            raise RuntimeError(
                                f'Task no longer exists: {" ".join(x for x,y in zip(proc.step._pending_tasks, res) if y == "missing")}')
                        else:
                            raise RuntimeError(
                                f'Task {" ".join(proc.step._pending_tasks)} returned with status {" ".join(res)}')
                    elif proc.in_status('target_pending') and hasattr(proc.step, '_from_nested') and proc.step._child_socket.poll(0):
                        # see if the child node has sent something
                        res = proc.step._child_socket.recv()
                        if res == b'target_resolved':
                            # this _socket is the socket to the step
                            proc.step._socket.send_pyobj(res)
                            proc.step._status = 'running'
                        else:
                            raise RuntimeError(f'Unrecognized response from child process {res}')

                # step 3: check if there is room and need for another job
                while True:
                    if manager.all_busy():
                        break
                    #
                    # if steps from child nested workflow?
                    if self.step_queue:
                        step_id, step_param = self.step_queue.popitem()
                        section, context, shared, args, config, verbosity, port = step_param
                        # run it!
                        runnable = dummy_node()
                        runnable._node_id = step_id
                        runnable._status = 'running'
                        dag.save(env.config['output_dag'])
                        runnable._from_nested = True
                        runnable._child_socket = env.zmq_context.socket(zmq.PAIR)
                        runnable._child_socket.connect(f'tcp://127.0.0.1:{port}')

                        env.logger.debug(
                            f'{i_am()} sends {section.step_name()} from step queue with args {args} and context {context}')

                        manager.execute(runnable, config=env.config, args=self.args,
                                        spec=('step', section, context, shared, args, config, verbosity))
                        continue

                    # find any step that can be executed and run it, and update the DAT
                    # with status.
                    runnable = dag.find_executable()
                    if runnable is None:
                        # no runnable
                        # dag.show_nodes()
                        break

                    # find the section from runnable
                    section = self.workflow.section_by_id(runnable._step_uuid)
                    # execute section with specified input
                    runnable._status = 'running'
                    dag.save(env.config['output_dag'])

                    # workflow shared variables
                    shared = {x: env.sos_dict[x] for x in self.shared.keys(
                    ) if x in env.sos_dict and pickleable(env.sos_dict[x], x)}
                    if 'shared' in section.options:
                        shared.update(self.get_shared_vars(
                            section.options['shared']))

                    if 'workflow_id' in env.sos_dict:
                        runnable._context['workflow_id'] = env.sos_dict['workflow_id']

                    env.logger.debug(
                        f'{i_am()} execute {section.md5} from DAG')
                    manager.execute(runnable, config=env.config, args=self.args,
                                    spec=('step', section, runnable._context, shared, self.args,
                                          env.config, env.verbosity))

                if manager.all_done():
                    break
                # if manager.all_failed():
                #     steps = list(dict.fromkeys(
                #         [str(x.step) for x in manager.procs]))
                #     raise RuntimeError(
                #         f'Workflow exited due to failed step{"s" if len(steps) > 1 else ""} {", ".join(steps)}.')
                else:
                    time.sleep(0.1)
        except KeyboardInterrupt:
            if exec_error.errors:
                failed_steps, pending_steps = dag.pending()
                if pending_steps:
                    sections = [self.workflow.section_by_id(
                        x._step_uuid).step_name() for x in pending_steps]
                    exec_error.append(self.workflow.name,
                                      RuntimeError(
                                          f'{len(sections)} pending step{"s" if len(sections) > 1 else ""}: {", ".join(sections)}'))
                    raise exec_error
            else:
                raise
        # close all processes
        except Exception as e:
            if not isinstance(e, ExecuteError):
                exec_error.append(self.workflow.name, e)
            manager.terminate(brutal=True)
        finally:
            manager.terminate()
        #

        if exec_error.errors:
            failed_steps, pending_steps = dag.pending()
            # if failed_steps:
            # sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in failed_steps]
            # exec_error.append(self.workflow.name,
            #    RuntimeError('{} failed step{}: {}'.format(len(sections),
            #        's' if len(sections) > 1 else '', ', '.join(sections))))
            if pending_steps:
                sections = [self.workflow.section_by_id(
                    x._step_uuid).step_name() for x in pending_steps]
                exec_error.append(self.workflow.name,
                                  RuntimeError(
                                      f'{len(sections)} pending step{"s" if len(sections) > 1 else ""}: {", ".join(sections)}'))
            raise exec_error
        elif 'pending_tasks' not in wf_result or not wf_result['pending_tasks']:
            self.finalize_and_report()
        else:
            # exit with pending tasks
            pass
        wf_result['shared'] = {x: env.sos_dict[x]
                               for x in self.shared.keys() if x in env.sos_dict}
        wf_result['__completed__'] = self.completed
        return wf_result

    def run_as_nested(self, parent_socket, targets=None, my_workflow_id='', mode=None) -> Dict[str, Any]:
        #
        # run a nested workflow, it simply send all steps and tasks to the master to execute
        #
        # this function still uses a manager, but never calls manager.execute()
        #
        self.completed = defaultdict(int)

        def i_am():
            return 'Nested'

        self.reset_dict()
        env.config['run_mode'] = env.config.get(
            'run_mode', 'run') if mode is None else mode
        env.sos_dict.set('run_mode', env.config['run_mode'])

        wf_result = {'__workflow_id__': my_workflow_id, 'shared': {}}

        # this is the initial targets specified by subworkflow, users
        # should specify named_output directly if needed.
        dag = self.initialize_dag(targets=targets)

        # the mansger will have all fake executors
        manager = ExecutionManager(env.config['max_procs'])
        #
        # steps sent and queued from the nested workflow
        # they will be executed in random but at a higher priority than the steps
        # on the master process.
        self.step_queue = {}
        try:
            exec_error = ExecuteError(self.workflow.name)
            while True:
                # step 1: check existing jobs and see if they are completed
                for idx, proc in enumerate(manager.procs):
                    if proc is None:
                        continue

                    # echck if there is any message from the socket
                    if not proc.socket.poll(0):
                        continue

                    # receieve something from the pipe
                    res = proc.socket.recv_pyobj()
                    runnable = proc.step

                    if isinstance(res, list):
                        # missing target from nested workflow
                        if res[0] == 'missing_target':
                            missed = res[1]
                            try:
                                self.handle_unknown_target(missed, dag, runnable)
                                # tell the master that the nested can resolve the target
                                proc.socket.send_pyobj(True)
                                runnable._status = 'target_pending'
                                runnable._pending_targets = [missed]
                                # when the target is resolved, tell the parent that
                                # the target is resolved and the step can continue
                                runnable._socket = proc.socket
                            except Exception as e:
                                env.logger.error(e)
                                # tell the master that nested cannot resolve the
                                # target so the workflow should stop
                                proc.socket.send_pyobj(False)
                            continue
                        elif res[0] == 'dependent_target':
                            reply = self.handle_dependent_target(dag, sos_targets(res[1:]), runnable)
                            proc.socket.send_pyobj(reply)
                            if reply:
                                # tell the master that the nested can resolve the target
                                runnable._status = 'target_pending'
                                runnable._pending_targets = res[1:]
                                # when the target is resolved, tell the parent that
                                # the target is resolved and the step can continue
                                runnable._socket = proc.socket
                            continue
                        else:
                            raise RuntimeError(
                                f'Unexpected value from step {short_repr(res)}')

                    manager.mark_idle(idx)
                    env.logger.debug(
                        f'{i_am()} receive a result {short_repr(res)}')
                    if isinstance(res, UnavailableLock):
                        self.handle_unavailable_lock(res, dag, runnable)
                    elif isinstance(res, RemovedTarget):
                        # RemovedTarget is usually hanled at the step level
                        # by sending a missing-target message here. However,
                        # if the exception is raised from substep workers, it is
                        # difficult for a substep to rerun particular substep
                        # so we handle it here by rerunning the entire step
                        self.handle_unknown_target(res.target, dag, runnable)
                    # if the job is failed
                    elif isinstance(res, Exception):
                        env.logger.debug(f'{i_am()} received an exception')
                        runnable._status = 'failed'
                        dag.save(env.config['output_dag'])
                        exec_error.append(runnable._node_id, res)
                        raise exec_error
                    elif '__step_name__' in res:
                        env.logger.debug(f'{i_am()} receive step result ')
                        self.step_completed(res, dag, runnable)
                    else:
                        raise RuntimeError(
                            f'Nested wokflow received an unrecognized response: {res}')

                # remove None
                manager.cleanup()
                # step 3: check if there is room and need for another job
                while True:
                    # with status.
                    runnable = dag.find_executable()
                    if runnable is None:
                        break

                    # find the section from runnable
                    section = self.workflow.section_by_id(runnable._step_uuid)
                    # execute section with specified input
                    runnable._status = 'running'
                    dag.save(env.config['output_dag'])

                    # workflow shared variables
                    shared = {x: env.sos_dict[x] for x in self.shared.keys(
                    ) if x in env.sos_dict and pickleable(env.sos_dict[x], x)}
                    if 'shared' in section.options:
                        shared.update(self.get_shared_vars(
                            section.options['shared']))

                    if 'workflow_id' in env.sos_dict:
                        runnable._context['workflow_id'] = env.sos_dict['workflow_id']

                    # send the step to the parent
                    step_id = uuid.uuid4()
                    env.logger.debug(
                        f'{i_am()} send step {section.step_name()} to master with args {self.args} and context {runnable._context}')

                    socket = env.zmq_context.socket(zmq.PAIR)
                    port = socket.bind_to_random_port('tcp://127.0.0.1')
                    parent_socket.send_pyobj(['step', step_id, section, runnable._context, shared, self.args,
                                              env.config, env.verbosity, port])
                    # the nested workflow also needs a step to receive result
                    manager.add_placeholder_worker(runnable, socket)

                if manager.all_done():
                    break
                else:
                    time.sleep(0.01)
        except KeyboardInterrupt:
            if exec_error.errors:
                failed_steps, pending_steps = dag.pending()
                if pending_steps:
                    sections = [self.workflow.section_by_id(
                        x._step_uuid).step_name() for x in pending_steps]
                    exec_error.append(self.workflow.name,
                                      RuntimeError(
                                          f'{len(sections)} pending step{"s" if len(sections) > 1 else ""}: {", ".join(sections)}'))
                    raise exec_error
            else:
                raise
        except Exception as e:
            if not isinstance(e, ExecuteError):
                exec_error.append(self.workflow.name, e)
            manager.terminate(brutal=True)

        if exec_error.errors:
            failed_steps, pending_steps = dag.pending()
            # if failed_steps:
            # sections = [self.workflow.section_by_id(x._step_uuid).step_name() for x in failed_steps]
            # exec_error.append(self.workflow.name,
            #    RuntimeError('{} failed step{}: {}'.format(len(sections),
            #        's' if len(sections) > 1 else '', ', '.join(sections))))
            if pending_steps:
                sections = [self.workflow.section_by_id(
                    x._step_uuid).step_name() for x in pending_steps]
                exec_error.append(self.workflow.name,
                                  RuntimeError(
                                      f'{len(sections)} pending step{"s" if len(sections) > 1 else ""}: {", ".join(sections)}'))
            parent_socket.send_pyobj(exec_error)
        else:
            wf_result['shared'] = {x: env.sos_dict[x]
                                   for x in self.shared.keys() if x in env.sos_dict}
            wf_result['__completed__'] = self.completed
            parent_socket.send_pyobj(wf_result)
