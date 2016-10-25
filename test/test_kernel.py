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
from ipykernel.tests.utils import assemble_output, start_new_kernel, flush_channels, stop_global_kernel, execute, wait_for_idle

import os
import atexit

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
    KM.shutdown_kernel(now=True)
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
        else:
            # other output, ignored
            pass
    # text/plain can have fronzen dict, this is ok,
    from pysos.utils import frozendict
    # it can also have dict_keys, we will have to redefine it
    def dict_keys(args):
        return args
    if result is None:
        return None
    else:
        try:
            return eval(result['text/plain'])
        except:
            # result from R etc are not python syntax
            return result['text/plain']

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

class TestKernel(unittest.TestCase):
    def testInterpolation(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="print('a=${100+11}')")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stdout, 'a=111\n')
            self.assertEqual(stderr, '')

    def testMagicDict(self):
        '''Test %dict magic'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="a=12345")
            clear_channels(iopub)
            msg_id, content = execute(kc=kc, code="%dict")
            self.assertEqual(get_result(iopub)['a'], 12345)
            #
            msg_id, content = execute(kc=kc, code="%dict keys")
            self.assertTrue('a' in get_result(iopub))
            #
            msg_id, content = execute(kc=kc, code="%dict reset")
            self.assertTrue('a' not in get_result(iopub))
            #
            msg_id, content = execute(kc=kc, code="%dict keys all")
            res = get_result(iopub)
            for key in ('run', 'sh', 'tcsh', 'expand_pattern'):
                self.assertTrue(key in res)

    def testShell(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="!ls test_kernel.py")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stdout, 'test_kernel.py\n')
            self.assertEqual(stderr, '')

    def testCD(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="!cd ..")
            clear_channels(iopub)
            msg_id, content = execute(kc=kc, code="print(os.getcwd())")
            stdout, stderr = assemble_output(iopub)
            self.assertTrue(stdout.strip().upper().endswith('SOS'))
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="!cd test")
        
    def testSubKernel(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a <- 1024")
            stdout, stderr = assemble_output(iopub)
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_result(iopub)
            self.assertEqual(res, '[1] 1024')
    
    def testMagicPut(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a <- 1024")
            stdout, stderr = assemble_output(iopub)
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put a")
            clear_channels(iopub)
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_result(iopub)
            self.assertEqual(res, 1024)

    def testMagicGet(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="a = 1025")
            clear_channels(iopub)
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%get a")
            clear_channels(iopub)
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_result(iopub)
            self.assertEqual(res, '[1] 1025')

    def testPutPythonDataFrameToR(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='''
import pandas as pd
import numpy as np
arr = np.random.randn(1000) 
arr[::10] = np.nan
df = pd.DataFrame({'column_{0}'.format(i): arr for i in range(10)})
''')
            clear_channels(iopub)
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%get df")
            clear_channels(iopub)
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="dim(df)")
            res = get_result(iopub)
            self.assertEqual(res, '[1] 1000   10')

    def testPutRDataFrameToPython(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='''
%use R
''')
            clear_channels(iopub)
            msg_id, content = execute(kc=kc, code="%put mtcars")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use sos")
            clear_channels(iopub)
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="mtcars.shape")
            res = get_result(iopub)
            self.assertEqual(res, '32, 11')

if __name__ == '__main__':
    unittest.main()
