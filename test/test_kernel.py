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
from ipykernel.tests.utils import assemble_output, start_new_kernel,\
    flush_channels, stop_global_kernel, execute, wait_for_idle
from pysos.utils import frozendict

from numpy import array
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
    from numpy import array
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
            self.assertTrue(stdout.endswith('a=111\n'))
            self.assertEqual(stderr, '')

    def testMagicDict(self):
        '''Test %dict magic'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="a=12345")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict a")
            self.assertEqual(get_result(iopub)['a'], 12345)
            msg_id, content = execute(kc=kc, code="%dict --keys")
            self.assertTrue('a' in get_result(iopub))
            msg_id, content = execute(kc=kc, code="%dict --reset")
            self.assertTrue('a' not in get_result(iopub))
            msg_id, content = execute(kc=kc, code="%dict --keys --all")
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
            wait_for_idle(kc)
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
            msg_id, content = execute(kc=kc, code="a <- 1024")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_result(iopub)
            self.assertEqual(res, '[1] 1024')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
    
    def testMagicPut(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="a <- 1024")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put a")
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
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="%get a")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_result(iopub)
            self.assertEqual(res, '[1] 1025')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)

    def testGetPythonDataFrameFromR(self):
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
            msg_id, content = execute(kc=kc, code="%get df")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="dim(df)")
            res = get_result(iopub)
            self.assertEqual(res, '[1] 1000   10')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)

    def testGetPythonDataFromR(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="null_var = None")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="num_var = 123")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="import numpy\nnum_arr_var = numpy.array([1, 2, 3])")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="logic_var = True")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="logic_arr_var = [True, False, True]")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="char_var = '123'")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="char_arr_var = ['1', '2', '3']")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="list_var = [1, 2, '3']")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="dict_var = dict(a=1, b=2, c='3')")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="set_var = {1, 2, '3'}")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="mat_var = numpy.matrix([[1,2],[3,4]])")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use R")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%get null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var set_var list_var dict_var")
            wait_for_idle(kc)
            # need to test passed values
            # but let us cheat by passing data back
            msg_id, content = execute(kc=kc, code="%dict -r")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var set_var list_var dict_var")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var set_var list_var dict_var")
            res = get_result(iopub)
            self.assertEqual(res['null_var'], None)
            self.assertEqual(res['num_var'], 123)
            self.assertEqual(res['num_arr_var'], [1,2,3])
            self.assertEqual(res['logic_var'], True)
            self.assertEqual(res['logic_arr_var'], [True, False, True])
            self.assertEqual(res['char_var'], '123')
            self.assertEqual(res['char_arr_var'], ['1', '2', '3'])
            self.assertEqual(res['list_var'], ['1','2','3'])
            self.assertEqual(res['dict_var'], {'a': 1, 'b': 2, 'c': '3'})
            self.assertEqual(res['mat_var'].shape, (2,2))
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)

    def testPutRDataFrameToPython(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='%use R')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put mtcars")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="mtcars.shape")
            res = get_result(iopub)
            self.assertEqual(res, (32, 11))
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)

    def testPutRDataToPython(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='%use R')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="null_var = NULL")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="num_var = 123")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="num_arr_var = c(1, 2, 3)")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="logic_var = TRUE")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="logic_arr_var = c(TRUE, FALSE, TRUE)")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="char_var = '123'")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="char_arr_var = c(1, 2, '3')")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="list_var = list(1, 2, '3')")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="named_list_var = list(a=1, b=2, c='3')")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="mat_var = matrix(c(1,2,3,4), nrow=2)")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var list_var named_list_var")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var list_var named_list_var")
            res = get_result(iopub)
            self.assertEqual(res['null_var'], None)
            self.assertEqual(res['num_var'], 123)
            self.assertEqual(res['num_arr_var'], [1,2,3])
            self.assertEqual(res['logic_var'], True)
            self.assertEqual(res['logic_arr_var'], [True, False, True])
            self.assertEqual(res['char_var'], '123')
            self.assertEqual(res['char_arr_var'], ['1', '2', '3'])
            self.assertEqual(res['list_var'], [1,2,'3'])
            self.assertEqual(res['named_list_var'], {'a': 1, 'b': 2, 'c': '3'})
            self.assertEqual(res['mat_var'].shape, (2,2))
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)

if __name__ == '__main__':
    unittest.main()
