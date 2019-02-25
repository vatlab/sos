#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import multiprocessing as mp
import os
import signal
import subprocess
import sys
import time
from typing import Any, Dict, Optional

import zmq

from ._version import __version__
from .controller import (close_socket, connect_controllers, create_socket,
                         disconnect_controllers)
from .eval import SoS_exec
from .executor_utils import kill_all_subprocesses
from .targets import sos_targets
from .utils import (WorkflowDict, env, get_traceback, load_config_files,
                    short_repr, ProcessKilled)

def signal_handler(*args, **kwargs):
    raise ProcessKilled()

class SoS_Worker(mp.Process):
    '''
    Worker process to process SoS step or workflow in separate process.
    '''

    def __init__(self, config: Optional[Dict[str, Any]] = None, args: Optional[Any] = None,
            **kwargs) -> None:
        '''

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
        self.config = config

        self.args = [] if args is None else args

        # there can be multiple jobs for this worker, each using their own port and socket
        self._master_sockets = []
        self._master_ports = []
        self._stack_idx = 0

    def reset_dict(self):
        env.sos_dict = WorkflowDict()
        env.parameter_vars.clear()

        env.sos_dict.set('__args__', self.args)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        env.sos_dict.set('__step_output__', sos_targets())

        # load configuration files
        load_config_files(env.config['config_file'])

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)

        if isinstance(self.args, dict):
            for key, value in self.args.items():
                if not key.startswith('__'):
                    env.sos_dict.set(key, value)

    def run(self):
        # env.logger.warning(f'Worker created {os.getpid()}')
        env.config.update(self.config)
        env.zmq_context = connect_controllers()

        # create controller socket
        env.ctrl_socket = create_socket(env.zmq_context, zmq.REQ, 'worker backend')
        env.ctrl_socket.connect(f'tcp://127.0.0.1:{self.config["sockets"]["worker_backend"]}')

        signal.signal(signal.SIGTERM, signal_handler)

        # create at last one master socket
        env.master_socket = create_socket(env.zmq_context, zmq.PAIR)
        port = env.master_socket.bind_to_random_port('tcp://127.0.0.1')
        self._master_sockets.append(env.master_socket)
        self._master_ports.append(port)

        # result socket used by substeps
        env.result_socket = None
        env.result_socket_port = None

        # wait to handle jobs
        while True:
            try:
                if not self._process_job():
                    break
            except ProcessKilled:
                # in theory, this will not be executed because the exception
                # will be caught by the step executor, and then sent to the master
                # process, which will then trigger terminate() and send a None here.
                break
            except KeyboardInterrupt:
                break
        # Finished
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        kill_all_subprocesses(os.getpid())

        close_socket(env.result_socket, 'substep result', now=True)

        for socket in self._master_sockets:
            close_socket(socket, 'worker master', now=True)
        close_socket(env.ctrl_socket, now=True)
        disconnect_controllers(env.zmq_context)

    def push_env(self):
        self._stack_idx += 1
        env.switch(self._stack_idx)
        if len(self._master_sockets) > self._stack_idx:
            # if current stack is ok
            env.master_socket = self._master_sockets[self._stack_idx]
        else:
            # a new socket is needed
            env.master_socket = create_socket(env.zmq_context, zmq.PAIR)
            port = env.master_socket.bind_to_random_port('tcp://127.0.0.1')
            self._master_sockets.append(env.master_socket)
            self._master_ports.append(port)
            env.logger.trace(f'WORKER {self.name} ({os.getpid()}) creates ports {self._master_ports}')

    def pop_env(self):
        self._stack_idx -= 1
        env.switch(self._stack_idx)
        env.master_socket = self._master_sockets[self._stack_idx]

    def _type_of_work(self, work):
        if 'section' in work:
            return 'step'
        elif 'wf' in work:
            return 'workflow'
        else:
            return 'substep'

    def _name_of_work(self, work):
        if 'section' in work:
            return work['section'].step_name()
        elif 'wf' in work:
            return work['workflow_id']
        else:
            return 'substep'

    def _process_job(self):
        # send all the available sockets...
        env.ctrl_socket.send_pyobj([self._stack_idx] + self._master_ports[self._stack_idx:])
        work = env.ctrl_socket.recv_pyobj()

        if work is None:
            if self._stack_idx != 0:
                env.logger.error(f'WORKER terminates with pending tasks. sos might not be termianting properly.')
            env.logger.trace(f'WORKER {self.name} ({os.getpid()}) quits after receiving None.')
            return False
        elif not work: # an empty task {}
            time.sleep(0.1)
            return True

        env.logger.trace(
            f'WORKER {self.name} ({os.getpid()}, level {self._stack_idx}) receives {self._type_of_work(work)} request {self._name_of_work(work)} with master port {self._master_ports[self._stack_idx]}')

        if 'task' in work:
            self.run_substep(work)
            env.logger.trace(
                f'WORKER {self.name} ({os.getpid()}) completes substep {self._name_of_work(work)}')
            return True

        master_port = work['config']['sockets']['master_port']
        if master_port != self._master_ports[self._stack_idx]:
            idx = self._master_ports.index(master_port)
            self._master_sockets[idx], self._master_sockets[self._stack_idx] = self._master_sockets[self._stack_idx], self._master_sockets[idx]
            self._master_ports[idx], self._master_ports[self._stack_idx] = self._master_ports[self._stack_idx], self._master_ports[idx]
            env.master_socket = self._master_sockets[self._stack_idx]

        # step and workflow can yield
        runner = self.run_step(**work) if 'section' in work else self.run_workflow(**work)
        try:
            poller = next(runner)
            while True:
                # if request is None, it is a normal "break" and
                # we do not need to jump off
                if poller is None:
                    poller = runner.send(None)
                    continue

                while True:
                    if poller.poll(200):
                        poller = runner.send(None)
                        break
                    # now let us ask if the master has something else for us
                    self.push_env()
                    self._process_job()
                    self.pop_env()
        except StopIteration as e:
            pass
        env.logger.trace(
            f'WORKER {self.name} ({os.getpid()}) completes request {self._type_of_work(work)} request {self._name_of_work(work)}')
        return True

    def run_workflow(self, workflow_id, wf, targets, args, shared, config, **kwargs):
        #
        #
        # get workflow, args, shared, and config
        from .workflow_executor import Base_Executor

        self.args = args
        env.config.update(config)
        self.reset_dict()
        # we are in a separate process and need to set verbosity from workflow config
        # but some tests do not provide verbosity
        env.verbosity = config.get('verbosity', 2)
        env.logger.debug(
            f'Worker {self.name} working on a workflow {workflow_id} with args {args}')
        executer = Base_Executor(wf, args=args, shared=shared, config=config)
        # we send the socket to subworkflow, which would send
        # everything directly to the master process, so we do not
        # have to collect result here
        try:
            runner = executer.run_as_nested(targets=targets, parent_socket=env.master_socket,
                         my_workflow_id=workflow_id)
            try:
                yreq = next(runner)
                while True:
                    yres = yield yreq
                    yreq = runner.send(yres)
            except StopIteration:
                pass

        except Exception as e:
            env.master_socket.send_pyobj(e)

    def run_step(self, section, context, shared, args, config, verbosity):
        from .step_executor import Step_Executor

        env.logger.debug(
            f'Worker {self.name} working on {section.step_name()} with args {args}')
        env.config.update(config)
        env.verbosity = verbosity
        #
        self.args = args
        self.reset_dict()

        # Execute global namespace. The reason why this is executed outside of
        # step is that the content of the dictioary might be overridden by context
        # variables.
        try:
            SoS_exec(section.global_def)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(e.stderr)
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

        executor = Step_Executor(
            section, env.master_socket, mode=env.config['run_mode'])

        runner = executor.run()
        try:
            yreq = next(runner)
            while True:
                yres = yield yreq
                yreq = runner.send(yres)
        except StopIteration:
            pass

    def run_substep(self, work):
        from .substep_executor import execute_substep
        execute_substep(**work)


class WorkerManager(object):
    # manager worker processes

    def __init__(self, max_workers, backend_socket):
        self._max_workers = max_workers

        self._workers = []
        self._num_workers = 0
        self._n_requested = 0
        self._n_processed = 0

        self._worker_alive_time = time.time()
        self._last_avail_time = time.time()

        self._substep_requests = []
        self._step_requests = {}

        self._worker_backend_socket = backend_socket

        # ports of workers working for blocking workflow
        self._blocking_ports = set()

        self._available_ports = set()
        self._claimed_ports = set()

        # start a worker
        self.start()

    def report(self, msg):
        return
        env.logger.trace(f'{msg.upper()}: {self._num_workers} workers (of which {len(self._blocking_ports)} is blocking), {self._n_requested} requested, {self._n_processed} processed')

    def add_request(self, msg_type, msg):
        self._n_requested += 1
        if msg_type == 'substep':
            self._substep_requests.insert(0, msg)
            self.report(f'Substep requested')
        else:
            port = msg['config']['sockets']['master_port']
            self._step_requests[port] = msg
            self.report(f'Step {port} requested')

        # start a worker is necessary (max_procs could be incorrectly set to be 0 or less)
        # if we are just starting, so do not start two workers
        if self._n_processed > 0 and not self._available_ports and self._num_workers < self._max_workers:
            self.start()

    def worker_available(self, blocking):
        if self._available_ports:
            claimed = self._available_ports.pop()
            self._claimed_ports.add(claimed)
            return claimed

        if not blocking:
            # no available port, can we start a new worker?
            if self._num_workers < self._max_workers:
                self.start()
            return None

        # we start a worker right now.
        self.start()
        while True:
            if not self._worker_backend_socket.poll(5000):
                raise RuntimeError('No worker is started after 5 seconds')
            msg = self._worker_backend_socket.recv_pyobj()
            port = self.process_request(msg[0], msg[1:], request_blocking=True)
            if port is None:
                continue
            self._claimed_ports.add(port)
            self._max_workers += 1
            self._blocking_ports.add(port)
            env.logger.debug(f'Increasing maximum number of workers to {self._max_workers} to accommodate a blocking subworkflow.')
            return port

    def process_request(self, level, ports, request_blocking=False):
        '''port is the open port at the worker, level is the level of stack.
        A non-zero level means that the worker is pending on something while
        looking for new job, so the worker should not be killed.
        '''
        if any(port in self._step_requests for port in ports):
            # if the port is available
            port = [x for x in ports if x in self._step_requests][0]
            self._worker_backend_socket.send_pyobj(self._step_requests.pop(port))
            self._last_avail_time = time.time()
            self._n_processed += 1
            self.report(f'Step {port} processed')
            # port should be in claimed ports
            self._claimed_ports.remove(port)
        elif any(port in self._claimed_ports for port in ports):
            # the port is claimed, but the real message is not yet available
            self._worker_backend_socket.send_pyobj({})
            self.report(f'pending with claimed {ports}')
        elif any(port in self._blocking_ports for port in ports):
            # in block list but appear to be idle, kill it
            self._max_workers -= 1
            env.logger.debug(f'Reduce maximum number of workers to {self._max_workers} after completion of a blocking subworkflow.')
            for port in ports:
                if port in self._blocking_ports:
                    self._blocking_ports.remove(port)
                if port in self._available_ports:
                    self._available_ports.remove(port)
            self._worker_backend_socket.send_pyobj(None)
            self._num_workers -= 1
            self.report(f'Blocking worker {ports} killed')
        elif self._substep_requests:
            # port is not claimed, free to use for substep worker
            msg = self._substep_requests.pop()
            self._worker_backend_socket.send_pyobj(msg)
            self._last_avail_time = time.time()
            self._n_processed += 1
            self.report(f'Substep processed with {ports[0]}')
            # port can however be in available ports
            for port in ports:
                if port in self._available_ports:
                    self._available_ports.remove(port)
        elif request_blocking:
            self._worker_backend_socket.send_pyobj({})
            return ports[0]
        else:
            # the port will be available for others to use
            self._available_ports.add(ports[0])
            self._worker_backend_socket.send_pyobj({})
            self.report(f'pending with port {ports}')

    def start(self):
        worker = SoS_Worker(env.config)
        worker.start()
        self._workers.append(worker)
        self._num_workers += 1
        self.report('start worker')

    def check_workers(self):
        '''Kill workers that have been pending for a while and check if all workers
        are alive. '''
        if time.time() - self._worker_alive_time > 5:
            self._worker_alive_time = time.time()
            self._workers = [worker for worker in self._workers if worker.is_alive()]
            if len(self._workers) < self._num_workers:
                raise ProcessKilled('One of the workers has been killed.')
        # if there is at least one request has been processed in 5 seconds
        if time.time() - self._last_avail_time < 5:
            return
        # we keep at least one worker
        attempts = self._num_workers - 1
        while attempts > 0:
            attempts -= 1
            if not self._worker_backend_socket.poll(100):
                continue
            msg = self._worker_backend_socket.recv_pyobj()
            if any(port in self._claimed_ports for port in msg[1:]) or msg[0] > 0:
                self._worker_backend_socket.send_pyobj({})
                continue
            for port in msg[1:]:
                if port in self._available_ports:
                    self._available_ports.remove(port)
            self._worker_backend_socket.send_pyobj(None)
            self._num_workers -= 1
            self.report(f'Kill standing {msg[1:]}')

    def kill_all(self):
        '''Kill all workers'''
        while self._num_workers > 0 and self._worker_backend_socket.poll(1000):
            msg = self._worker_backend_socket.recv_pyobj()
            self._worker_backend_socket.send_pyobj(None)
            self._num_workers -= 1
            self.report('Kill {msg[1:]}')