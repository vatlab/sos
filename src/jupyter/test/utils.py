#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
#
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

#
# NOTE: for some namespace reason, this test can only be tested using
# nose. 
#
# % nosetests test_kernel.py
#
#
from contextlib import contextmanager
from ipykernel.tests.utils import start_new_kernel

import atexit
from queue import Empty
KM = None
KC = None


@contextmanager
def sos_kernel():
    """Context manager for the global kernel instance
    Should be used for most kernel tests
    Returns
    -------
    kernel_client: connected KernelClient instance
    """
    yield start_sos_kernel()


def flush_channels(kc=None):
    """flush any messages waiting on the queue"""

    if kc is None:
        kc = KC
    for channel in (kc.shell_channel, kc.iopub_channel):
        while True:
            try:
                channel.get_msg(block=True, timeout=0.1)
            except Empty:
                break
            # do not validate message because SoS has special sos_comm
            #else:
            #    validate_message(msg)

def start_sos_kernel():
    """start the global kernel (if it isn't running) and return its client"""
    global KM, KC
    if KM is None:
        KM, KC = start_new_kernel(kernel_name='sos')
        atexit.register(stop_sos_kernel)
    else:
        flush_channels(KC)
    return KC

def stop_sos_kernel():
    """Stop the global shared kernel instance, if it exists"""
    global KM, KC
    KC.stop_channels()
    KC = None
    if KM is None:
        return
    KM.shutdown_kernel(now=False)
    KM = None

def get_result(iopub):
    """retrieve result from an execution"""
    result = None
    while True:
        msg = iopub.get_msg(block=True, timeout=1)
        msg_type = msg['msg_type']
        content = msg['content']
        if msg_type == 'status' and content['execution_state'] == 'idle':
            # idle message signals end of output
            break
        elif msg['msg_type'] == 'execute_result':
            result = content['data']
        elif msg['msg_type'] == 'display_data':
            result = content['data']
        else:
            # other output, ignored
            pass
    # text/plain can have fronzen dict, this is ok,
    from sos.utils import frozendict
    from numpy import array, matrix
    # suppress pyflakes warning
    frozendict
    array
    matrix
    # it can also have dict_keys, we will have to redefine it
    def dict_keys(args):
        return args
    if result is None:
        return None
    else:
        return eval(result['text/plain'])

def get_display_data(iopub):
    """retrieve display_data from an execution from subkernel
    because subkernel (for example irkernel) does not return
    execution_result
    """
    result = None
    while True:
        msg = iopub.get_msg(block=True, timeout=1)
        msg_type = msg['msg_type']
        content = msg['content']
        if msg_type == 'status' and content['execution_state'] == 'idle':
            # idle message signals end of output
            break
        elif msg['msg_type'] == 'display_data':
            result = content['data']['text/plain']
        # some early version of IRKernel still passes execute_result
        elif msg['msg_type'] == 'execute_result':
            result = content['data']['text/plain']
        else:
            # other output, ignored
            pass
    return result

def clear_channels(iopub):
    """assemble stdout/err from an execution"""
    while True:
        msg = iopub.get_msg(block=True, timeout=1)
        msg_type = msg['msg_type']
        content = msg['content']
        if msg_type == 'status' and content['execution_state'] == 'idle':
            # idle message signals end of output
            break
        else:
            # other output, ignored
            pass

