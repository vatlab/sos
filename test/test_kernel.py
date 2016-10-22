#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
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
import unittest
from contextlib import contextmanager
from pysos.kernel import SoS_Kernel
from jupyter_client import manager

from subprocess import STDOUT
import os
import nose as nt
import atexit

KM = None
KC = None

@contextmanager
def kernel():
    """Context manager for the global kernel instance
    Should be used for most kernel tests
    Returns
    -------
    kernel_client: connected KernelClient instance
    """
    yield start_global_kernel()


def start_new_kernel(**kwargs):
    """start a new kernel, and return its Manager and Client
    Integrates with our output capturing for tests.
    """
    try:
        stdout = nt.iptest_stdstreams_fileno()
    except AttributeError:
        stdout = open(os.devnull)
    kwargs.update(dict(stdout=stdout, stderr=STDOUT))
    return manager.start_new_kernel(startup_timeout=10, **kwargs)


def start_global_kernel():
    """start the global kernel (if it isn't running) and return its client"""
    global KM, KC
    if KM is None:
        KM, KC = start_new_kernel()
        atexit.register(stop_global_kernel)
    else:
        flush_channels(KC)
    return KC

def uses_kernel(test_f):
    """Decorator for tests that use the global kernel"""
    def wrapped_test():
        with kernel() as kc:
            test_f(kc)
    wrapped_test.__doc__ = test_f.__doc__
    wrapped_test.__name__ = test_f.__name__
    return wrapped_test

def stop_global_kernel():
    """Stop the global shared kernel instance, if it exists"""
    global KM, KC
    KC.stop_channels()
    KC = None
    if KM is None:
        return
    KM.shutdown_kernel(now=True)
    KM = None

def execute(code='', kc=None, **kwargs):
    """wrapper for doing common steps for validating an execution request"""
    #from .test_message_spec import validate_message
    if kc is None:
        kc = KC
    msg_id = kc.execute(code=code, **kwargs)
    reply = kc.get_shell_msg(timeout=10)
    #validate_message(reply, 'execute_reply', msg_id)
    busy = kc.get_iopub_msg(timeout=10)
    #validate_message(busy, 'status', msg_id)
    #nt.assert_equal(busy['content']['execution_state'], 'busy')

    if not kwargs.get('silent'):
        execute_input = kc.get_iopub_msg(timeout=10)
        # validate_message(execute_input, 'execute_input', msg_id)
        #nt.assert_equal(execute_input['content']['code'], code)

    return msg_id, reply['content']


class TestKernel(unittest.TestCase):
    def testDict(self):
        with kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="print ('hi')")
            #stdout, stderr = assemble_output(iopub)
            #nt.assert_equal(stdout, 'hi\n')
            #nt.assert_equal(stderr, '')
            #_check_master(kc, expected=True)

if __name__ == '__main__':
    unittest.main()
