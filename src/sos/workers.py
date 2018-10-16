#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import subprocess
import sys
import zmq
import multiprocessing as mp

from typing import Any, Dict, Optional

from ._version import __version__
from .eval import SoS_exec
from .controller import connect_controllers, disconnect_controllers
from .targets import sos_targets
from .utils import WorkflowDict, env, get_traceback, load_config_files, short_repr
from .executor_utils import  __null_func__

class SoS_Worker(mp.Process):
    '''
    Worker process to process SoS step or workflow in separate process.
    '''

    def __init__(self, port: int, config: Optional[Dict[str, Any]] = None, args: Optional[Any] = None, **kwargs) -> None:
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
        self.port = port
        self.config = config

        self.args = [] if args is None else args

    def reset_dict(self):
        env.sos_dict = WorkflowDict()
        self.init_dict()

    def init_dict(self):
        env.parameter_vars.clear()

        env.sos_dict.set('__null_func__', __null_func__)
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
        env.config.update(self.config)
        env.zmq_context = connect_controllers()

        env.master_socket = env.zmq_context.socket(zmq.PAIR)
        env.master_socket.connect(f'tcp://127.0.0.1:{self.port}')

        # wait to handle jobs
        while True:
            try:
                work = env.master_socket.recv_pyobj()
                if work is None:
                    break
                env.logger.debug(
                    f'Worker {self.name} receives request {short_repr(work)}')
                if work[0] == 'step':
                    # this is a step ...
                    self.run_step(*work[1:])
                else:
                    self.run_workflow(*work[1:])
                env.logger.debug(
                    f'Worker {self.name} completes request {short_repr(work)}')
            except KeyboardInterrupt:
                break
        # Finished
        env.master_socket.LINGER = 0
        env.master_socket.close()
        disconnect_controllers(env.zmq_context)


    def run_workflow(self, workflow_id, wf, targets, args, shared, config):
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
            executer.run_as_nested(targets=targets, parent_socket=env.master_socket,
                         my_workflow_id=workflow_id)
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
        # this is to keep compatibility of dag run with sequential run because
        # in sequential run, we evaluate global section of each step in
        # order to determine values of options such as skip.
        # The consequence is that global definitions are available in
        # SoS namespace.
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
        executor.run()


class SoS_SubStep_Worker(mp.Process):
    '''
    Worker process to process SoS step or workflow in separate process.
    '''
    LRU_READY = b"\x01"

    def __init__(self, config={}, **kwargs) -> None:
        # the worker process knows configuration file, command line argument etc
        super(SoS_SubStep_Worker, self).__init__(**kwargs)
        self.config = config
        self.daemon = True

    def run(self):
        env.config.update(self.config)
        env.zmq_context = connect_controllers()
        from .substep_executor import execute_substep
        env.master_socket = env.zmq_context.socket(zmq.REQ)
        env.master_socket.connect(f'tcp://127.0.0.1:{self.config["sockets"]["substep_backend"]}')
        env.logger.trace(f'Substep worker {os.getpid()} started')

        while True:
            env.master_socket.send(self.LRU_READY)
            msg = env.master_socket.recv_pyobj()
            if not msg:
                env.logger.debug(f'stop substep worker {os.getpid()}')
                break

            env.logger.debug(f'Substep worker {os.getpid()} receives request {short_repr(msg)}')
            execute_substep(**msg)

        env.master_socket.LINGER = 0
        env.master_socket.close()
        disconnect_controllers(env.zmq_context)
