#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import zmq

import threading
from .utils import env
from .signatures import TargetSignatures, StepSignatures, WorkflowSignatures
# from zmq.utils.monitor import recv_monitor_message

def connect_controllers(context=None):
    if not context:
        env.zmq_context = zmq.Context()
    env.signature_push_socket = env.zmq_context.socket(zmq.PUSH)
    env.signature_push_socket.connect(f'tcp://127.0.0.1:{env.config["sockets"]["signature_push"]}')
    env.signature_req_socket = env.zmq_context.socket(zmq.REQ)
    env.signature_req_socket.connect(f'tcp://127.0.0.1:{env.config["sockets"]["signature_req"]}')

    env.controller_push_socket = env.zmq_context.socket(zmq.PUSH)
    env.controller_push_socket.connect(f'tcp://127.0.0.1:{env.config["sockets"]["controller_push"]}')
    env.controller_req_socket = env.zmq_context.socket(zmq.REQ)
    env.controller_req_socket.connect(f'tcp://127.0.0.1:{env.config["sockets"]["controller_req"]}')

class Controller(threading.Thread):
    def __init__(self, ready):
        threading.Thread.__init__(self)
        self.daemon = True

        self.target_signatures = TargetSignatures()
        self.step_signatures = StepSignatures()
        self.workflow_signatures = WorkflowSignatures()

        self.ready = ready

        # number of active running master processes
        self._nprocs = 0

        # self.event_map = {}
        # for name in dir(zmq):
        #     if name.startswith('EVENT_'):
        #         value = getattr(zmq, name)
        #          self.event_map[value] = name

    def run(self):
        # there are two sockets
        #
        # signature_push is used to write signatures. It is a single push operation with no reply.
        # signature_req is used to query information. The sender would need to get an response.
        env.zmq_context = zmq.Context()

        sig_push_socket = env.zmq_context.socket(zmq.PULL)
        env.config['sockets']['signature_push'] = sig_push_socket.bind_to_random_port('tcp://127.0.0.1')
        sig_req_socket = env.zmq_context.socket(zmq.REP)
        env.config['sockets']['signature_req'] = sig_req_socket.bind_to_random_port('tcp://127.0.0.1')

        ctl_push_socket = env.zmq_context.socket(zmq.PULL)
        env.config['sockets']['controller_push'] = ctl_push_socket.bind_to_random_port('tcp://127.0.0.1')
        ctl_req_socket = env.zmq_context.socket(zmq.REP)
        env.config['sockets']['controller_req'] = ctl_req_socket.bind_to_random_port('tcp://127.0.0.1')

        #monitor_socket = sig_req_socket.get_monitor_socket()
        # tell others that the sockets are ready
        self.ready.set()
        # Process messages from receiver and controller
        poller = zmq.Poller()
        poller.register(sig_push_socket, zmq.POLLIN)
        poller.register(sig_req_socket, zmq.POLLIN)
        poller.register(ctl_push_socket, zmq.POLLIN)
        poller.register(ctl_req_socket, zmq.POLLIN)
        #poller.register(monitor_socket, zmq.POLLIN)

        while True:
            try:
                socks = dict(poller.poll())

                if sig_push_socket in socks:
                    msg = sig_push_socket.recv_pyobj()
                    try:
                        if msg[0] == 'workflow':
                            self.workflow_signatures.write(*msg[1:])
                        elif msg[0] == 'target':
                            self.target_signatures.set(*msg[1:])
                        elif msg[0] == 'step':
                            self.step_signatures.set(*msg[1:])
                        else:
                            env.logger.warning(f'Unknown message passed {msg}')
                    except Exception as e:
                        env.logger.warning(f'Failed to push signature {msg}: {e}')

                if sig_req_socket in socks:
                    msg = sig_req_socket.recv_pyobj()
                    try:
                        if msg[0] == 'workflow':
                            if msg[1] == 'clear':
                                self.workflow_signatures.clear()
                                sig_req_socket.send_pyobj('ok')
                            elif msg[1] == 'placeholders':
                                sig_req_socket.send_pyobj(self.workflow_signatures.placeholders(msg[2]))
                            else:
                                env.logger.warning(f'Unknown signature request {msg}')
                        elif msg[0] == 'target':
                            if msg[1] == 'get':
                                sig_req_socket.send_pyobj(self.target_signatures.get(msg[2]))
                            else:
                                env.logger.warning(f'Unknown signature request {msg}')
                        elif msg[0] == 'step':
                            if msg[1] == 'get':
                                sig_req_socket.send_pyobj(self.step_signatures.get(*msg[2:]))
                            else:
                                env.logger.warning(f'Unknown signature request {msg}')
                        else:
                            raise RuntimeError(f'Unrecognized signature request {msg}')
                    except Exception as e:
                        env.logger.warning(f'Failed to respond to signature request {msg}: {e}')
                        sig_req_socket.send_pyobj(None)

                if ctl_push_socket in socks:
                    msg = ctl_push_socket.recv_pyobj()
                    try:
                        if msg[0] == 'nprocs':
                            env.logger.trace(f'Active running process set to {msg[1]}')
                            self._nprocs = msg[1]
                        else:
                            raise RuntimeError(f'Unrecognized request {msg}')
                    except Exception as e:
                        env.logger.warning(f'Failed to push controller {msg}: {e}')

                if ctl_req_socket in socks:
                    msg = ctl_req_socket.recv_pyobj()
                    try:
                        if msg[0] == 'nprocs':
                            ctl_req_socket.send_pyobj(self._nprocs)
                        else:
                            raise RuntimeError(f'Unrecognized request {msg}')
                    except Exception as e:
                        env.logger.warning(f'Failed to respond controller {msg}: {e}')
                # if monitor_socket in socks:
                #     evt = recv_monitor_message(monitor_socket)
                #     if evt['event'] == zmq.EVENT_ACCEPTED:
                #         self._num_clients += 1
                #     elif evt['event'] == zmq.EVENT_DISCONNECTED:
                #         self._num_clients -= 1
            except KeyboardInterrupt:
                break
