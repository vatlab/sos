#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import multiprocessing as mp
import os
import signal
import time
from typing import Any, Dict, Optional

import zmq

from .controller import (close_socket, connect_controllers, create_socket,
                         disconnect_controllers)
from .executor_utils import kill_all_subprocesses, prepare_env
from .utils import env, ProcessKilled, short_repr


def signal_handler(*args, **kwargs):
    raise ProcessKilled()


class Runner(object):
    '''
    This runner class takea a generator function and run it.
    1. When the generator returns None, continue to run without yielding.
    2. When the generator returns a poller, continue is it receives any message in 0.2s.
    3. Otherwise return False.
    4. The the generator completes, return True.

    So in summary, the Runner returns
    1. True if all completed, return value from generator is ignored.
    2. `self` if waiting.
    '''

    def __init__(self, runner, name):
        self._runner = runner
        self._poller = 0
        self._name = name

    def __repr__(self):
        return self._name

    def run_until_waiting(self):
        try:
            # if poller is not initialized, run the runner to get the first poller
            if self._poller == 0:
                self._poller = next(self._runner)

            while True:
                if self._poller is None:
                    self._poller = self._runner.send(None)
                    continue

                if self._poller.poll(200):
                    self._poller = self._runner.send(None)
                    continue

                # the poller is not ready, let us break
                return self
        except StopIteration:
            return True

    def can_proceed(self):
        return self._poller.poll(0)


class SoS_Worker(mp.Process):
    '''
    Worker process to process SoS step or workflow in separate process.
    '''

    def __init__(self, config: Optional[Dict[str, Any]] = None,
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

        # there can be multiple jobs for this worker, each using their own port and socket
        self._master_sockets = []
        self._master_ports = []
        # current runner, which can be a runner or True if the runner has completed
        self._runners = []
        # env index, which contains sos_dict for each runner
        self._env_idx = []

    def waiting_runners(self):
        # check if
        return [
            idx for idx, runner in enumerate(self._runners)
            if isinstance(runner, Runner) and runner.can_proceed()
        ]

    def completed_runners(self):
        return [
            idx for idx, runner in enumerate(self._runners) if runner is True
        ]

    def available_ports(self):
        return [
            port for port, runner in zip(self._master_ports, self._runners)
            if runner is True
        ]

    def num_pending(self):
        return len(
            [runner for runner in self._runners if isinstance(runner, Runner)])

    def switch_to(self, idx):
        if len(self._master_sockets) > idx:
            # if current stack is ok
            env.master_socket = self._master_sockets[idx]
            env.switch(self._env_idx[idx])
        else:
            assert idx == len(self._master_ports)
            # a new socket is needed
            env.master_socket = create_socket(env.zmq_context, zmq.PAIR)
            port = env.master_socket.bind_to_random_port('tcp://127.0.0.1')
            # switch to a new env_idx and returns new_idx, old_idx
            self._env_idx.append(env.request_new()[0])
            self._master_sockets.append(env.master_socket)
            self._master_ports.append(port)
            self._runners.append(True)
            if 'WORKER' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                env.log_to_file(
                    'WORKER',
                    f'WORKER {self.name} ({os.getpid()}) creates ports {self._master_ports}'
                )

    def __repr__(self):
        return self.name + ' ' + ' '.join(
            str(x) if isinstance(x, Runner) else str(idx)
            for idx, x in enumerate(self._runners))

    def run(self):
        # env.logger.warning(f'Worker created {os.getpid()}')
        env.config.update(self.config)
        env.zmq_context = connect_controllers()

        # create controller socket
        env.ctrl_socket = create_socket(env.zmq_context, zmq.REQ,
                                        'worker backend')
        env.ctrl_socket.connect(
            f'tcp://127.0.0.1:{self.config["sockets"]["worker_backend"]}')

        signal.signal(signal.SIGTERM, signal_handler)
        # result socket used by substeps
        env.result_socket = None
        env.result_socket_port = None

        # wait to handle jobs
        while True:
            try:
                wr = self.waiting_runners()
                if wr:
                    for idx in wr:
                        self.switch_to(idx)
                        # it can be True for completion and Runner itself for continue
                        self._runners[idx] = self._runners[
                            idx].run_until_waiting()
                    continue

                cr = self.completed_runners()
                # using an completed slot or create a new one
                new_idx = len(self._runners) if not cr else cr[0]
                self.switch_to(new_idx)

                # although we have chosen one port, but we hae advertised multiple ports
                # and the executor might choose another one. We therefore need to send all
                # avilable ports to the controller. We also need to send a flag to let the
                # controller know if we have any pending job, and the controller might decide
                # to kill this worker.
                env.ctrl_socket.send_pyobj([self.num_pending()] +
                                           self.available_ports())
                reply = env.ctrl_socket.recv_pyobj()

                if reply is None:
                    if len(wr) != 0:
                        env.logger.error(
                            f'WORKER terminates with pending tasks. sos might not be termianting properly.'
                        )
                    if 'WORKER' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                        env.log_to_file(
                            'WORKER',
                            f'WORKER {self.name} ({os.getpid()}) quits after receiving None.'
                        )
                    break
                if not reply:  # if an empty job is returned
                    time.sleep(0.1)
                    continue

                #
                # if a real job is returned, run it. _process_job will either return True
                # or a runner in case it is interrupted.
                env.log_to_file(
                    'WORKER',
                    f'WORKER {self.name} ({os.getpid()}, {self.num_pending()} pending) receives {self._type_of_work(reply)} request {self._name_of_work(reply)} with master port {self._master_ports[new_idx]}'
                )

                if 'task' in reply:
                    self.run_substep(reply)
                    env.log_to_file(
                        'WORKER',
                        f'WORKER {self.name} ({os.getpid()}) completes substep {self._name_of_work(reply)}'
                    )
                    self._runners[new_idx] = True
                    continue

                master_port = reply['config']['sockets']['master_port']
                if master_port != self._master_ports[new_idx]:
                    new_idx = self._master_ports.index(master_port)
                    self.switch_to(new_idx)

                # step and workflow can yield. Here we call run_until_waiting directly because we know the Runner can proceed.
                self._runners[new_idx] = Runner(
                    self.run_step(**reply)
                    if 'section' in reply else self.run_workflow(**reply),
                    name=self._name_of_work(reply)).run_until_waiting()
                if 'WORKER' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                    env.log_to_file(
                        'WORKER',
                        'STATUS ' + self._name_of_work(reply) + str(self))
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

    def run_workflow(self, workflow_id, wf, targets, args, shared, config,
                     **kwargs):
        #
        #
        # get workflow, args, shared, and config
        from .workflow_executor import Base_Executor
        env.config.update(config)
        # we are in a separate process and need to set verbosity from workflow config
        # but some tests do not provide verbosity
        env.verbosity = config.get('verbosity', 2)
        env.log_to_file(
            'WORKER',
            f'Worker {self.name} working on a workflow {workflow_id} with args {short_repr(args)}'
        )

        executer = Base_Executor(wf, args=args, shared=shared, config=config)
        # we send the socket to subworkflow, which would send
        # everything directly to the master process, so we do not
        # have to collect result here
        try:
            runner = executer.run_as_nested(
                targets=targets,
                parent_socket=env.master_socket,
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

        env.log_to_file(
            'WORKER',
            f'Worker {self.name} working on {section.step_name()} with args {short_repr(args)}'
        )
        env.config.update(config)
        env.verbosity = verbosity
        #
        # Execute global namespace. The reason why this is executed outside of
        # step is that the content of the dictioary might be overridden by context
        # variables.
        prepare_env(section.global_def, section.global_vars,
                    env.config['workflow_vars'])

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
        self._last_pending_time = {}

        self._substep_requests = []
        self._step_requests = {}

        self._worker_backend_socket = backend_socket

        # ports of workers working for blocking workflow
        self._blocking_ports = set()

        self._available_ports = set()
        self._claimed_ports = set()

        self._last_pending_msg = {}

        # start a worker
        self.start()

    def report(self, msg):
        if 'WORKER' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
            env.log_to_file(
                'WORKER',
                f'{msg.upper()}: {self._num_workers} workers (of which {len(self._blocking_ports)} is blocking), {self._n_requested} requested, {self._n_processed} processed'
            )

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

    def worker_available(self, blocking, excluded):
        if self._available_ports:
            usable = [x for x in self._available_ports if x not in excluded]
            if usable:
                claimed = usable[0]
                self._available_ports.remove(usable[0])
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
            if port is None or port in excluded:
                continue
            self._claimed_ports.add(port)
            self._max_workers += 1
            self._blocking_ports.add(port)
            env.logger.debug(
                f'Increasing maximum number of workers to {self._max_workers} to accommodate a blocking subworkflow.'
            )
            return port

    def process_request(self, num_pending, ports, request_blocking=False):
        '''port is the open port at the worker, num_pending is the num_pending of stack.
        A non-zero num_pending means that the worker is pending on something while
        looking for new job, so the worker should not be killed.
        '''
        if any(port in self._step_requests for port in ports):
            # if the port is available
            port = [x for x in ports if x in self._step_requests][0]
            self._worker_backend_socket.send_pyobj(
                self._step_requests.pop(port))
            self._n_processed += 1
            self.report(f'Step {port} processed')
            # port should be in claimed ports
            self._claimed_ports.remove(port)
            if ports[0] in self._last_pending_time:
                self._last_pending_time.pop(ports[0])
        elif any(port in self._claimed_ports for port in ports):
            # the port is claimed, but the real message is not yet available
            self._worker_backend_socket.send_pyobj({})
            self.report(f'pending with claimed {ports}')
        elif any(port in self._blocking_ports for port in ports):
            # in block list but appear to be idle, kill it
            self._max_workers -= 1
            env.logger.debug(
                f'Reduce maximum number of workers to {self._max_workers} after completion of a blocking subworkflow.'
            )
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
            self._n_processed += 1
            self.report(f'Substep processed with {ports[0]}')
            # port can however be in available ports
            for port in ports:
                if port in self._available_ports:
                    self._available_ports.remove(port)
                if port in self._last_pending_time:
                    self._last_pending_time.pop(port)
        elif request_blocking:
            self._worker_backend_socket.send_pyobj({})
            return ports[0]
        elif num_pending == 0 and ports[
                0] in self._last_pending_time and time.time(
                ) - self._last_pending_time[ports[0]] > 5:
            # kill the worker
            for port in ports:
                if port in self._available_ports:
                    self._available_ports.remove(port)
            self._worker_backend_socket.send_pyobj(None)
            self._num_workers -= 1
            self.report(f'Kill standing {ports}')
            self._last_pending_time.pop(ports[0])
        else:
            if num_pending == 0 and ports[0] not in self._last_pending_time:
                self._last_pending_time[ports[0]] = time.time()
            self._available_ports.add(ports[0])
            self._worker_backend_socket.send_pyobj({})
            ports = tuple(ports)
            if (ports, num_pending) not in self._last_pending_msg or time.time(
            ) - self._last_pending_msg[(ports, num_pending)] > 1.0:
                self.report(
                    f'pending with port {ports} at num_pending {num_pending}')
                self._last_pending_msg[(ports, num_pending)] = time.time()

    def start(self):
        worker = SoS_Worker(env.config)
        worker.start()
        self._worker_alive_time = time.time()
        self._workers.append(worker)
        self._num_workers += 1
        self.report('start worker')

    def check_workers(self):
        '''Kill workers that have been pending for a while and check if all workers
        are alive. '''
        if time.time() - self._worker_alive_time > 5:
            self._worker_alive_time = time.time()
            # join processes if they are now gone, it should not do anything bad
            # if the process is still running
            [worker.join() for worker in self._workers if not worker.is_alive()]
            self._workers = [
                worker for worker in self._workers if worker.is_alive()
            ]
            if len(self._workers) < self._num_workers:
                raise ProcessKilled('One of the workers has been killed.')

    def kill_all(self):
        '''Kill all workers'''
        while self._num_workers > 0 and self._worker_backend_socket.poll(1000):
            msg = self._worker_backend_socket.recv_pyobj()
            self._worker_backend_socket.send_pyobj(None)
            self._num_workers -= 1
            self.report(f'Kill {msg[1:]}')
        # join all processes
        [worker.join() for worker in self._workers]
