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

###############################################################################
#                                                                             #
# To test support for a language:                                             #
#                                                                             #
# 1. Import needed functions from ipykernel.tests.utils and                   #
#    sos.jupyter.test.utils.                                                  #
# 2. Define a test case with functions to get and put variables               #
#                                                                             #
# The following is a placeholder class for a new unittest. Please refer to    #
# tests for other languages for details.                                      #
#                                                                             #
###############################################################################

# import unittest
# from ipykernel.tests.utils import execute, wait_for_idle
# from sos.jupyter.test_utils import sos_kernel, get_result
#
# class TestKernel(unittest.TestCase):
#     def setUp(self):
#         self.olddir = os.getcwd()
#         if os.path.dirname(__file__):
#             os.chdir(os.path.dirname(__file__))
#
#     def tearDown(self):
#         os.chdir(self.olddir)
#
#     def testGetData(self):
#         with sos_kernel() as kc:
#             iopub = kc.iopub_channel
#             msg_id, content = execute(kc=kc, code='')
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code="%use LAN")
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code="%get VAR")
#             wait_for_idle(kc)
#
#     def testPutData(self):
#         with sos_kernel() as kc:
#             iopub = kc.iopub_channel
#             msg_id, content = execute(kc=kc, code="%use LAN")
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code='')
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code="%put VAR")
#             wait_for_idle(kc)
#
#
# if __name__ == '__main__':
#     unittest.main()




import os
import unittest
from ipykernel.tests.utils import assemble_output, execute, wait_for_idle
from sos.jupyter.test_utils import sos_kernel, get_result, get_display_data, \
    clear_channels

class TestOctaveKernel(unittest.TestCase):
    #
    # Beacuse these tests would be called from sos/test, we
    # should switch to this directory so that some location
    # dependent tests could run successfully
    #
    def setUp(self):
        self.olddir = os.getcwd()
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.olddir)
    
    # Fixme
    def testGetPythonDataFrameFromOctave(self):
        # Python -> Matlab/Octave
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            execute(kc=kc, code='''
import pandas as pd
import numpy as np
import scipy.io as sio
arr = np.random.randn(1000)
arr[::10] = np.nan
df = pd.DataFrame({'column_{0}'.format(i): arr for i in range(10)})
''')
            clear_channels(iopub)
            execute(kc=kc, code="%use Octave")
            _, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            execute(kc=kc, code="%get df")
            wait_for_idle(kc)
            execute(kc=kc, code="dim(df)")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 1000   10')
            execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
    #

    def testGetPythonDataFromOctave(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
null_var = None
num_var = 123
import numpy
import scipy.io as sio
num_arr_var = numpy.array([1, 2, 3])
logic_var = True
logic_arr_var = [True, False, True]
char_var = '1"23'
char_arr_var = ['1', '2', '3']
list_var = [1, 2, '3']
dict_var = dict(a=1, b=2, c='3')
set_var = {1, 2, '3'}
mat_var = numpy.matrix([[1,2],[3,4]])
recursive_var = {'a': {'b': 123}, 'c': True}
''')
            wait_for_idle(kc)
            execute(kc=kc, code='''
%use Octave
%get null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var set_var list_var dict_var recursive_var
%dict -r
%put null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var set_var list_var dict_var recursive_var
%use sos
%dict null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var set_var list_var dict_var recursive_var
                ''')
            res = get_result(iopub)
            self.assertEqual(res['null_var'], None)
            self.assertEqual(res['num_var'], 123)
            self.assertEqual(res['num_arr_var'], [1,2,3])
            self.assertEqual(res['logic_var'], True)
            self.assertEqual(res['logic_arr_var'], [True, False, True])
            self.assertEqual(res['char_var'], '1"23')
            self.assertEqual(res['char_arr_var'], ['1', '2', '3'])
            self.assertEqual(res['list_var'], [1,2,'3'])
            self.assertEqual(res['dict_var'], {'a': 1, 'b': 2, 'c': '3'})
            self.assertEqual(res['set_var'], {1, 2, '3'})
            self.assertEqual(res['mat_var'].shape, (2,2))
            self.assertEqual(res['recursive_var'],  {'a': {'b': 123}, 'c': True})

#def testPutOctaveDataFrameToPython(self):
    # Octave -> Python
    # Fixme
    #with sos_kernel() as kc:
    #iopub = kc.iopub_channel
            # create a data frame
            # execute(kc=kc, code='%use Octave')
            # wait_for_idle(kc)
            #  execute(kc=kc, code="%put mtcars")
            #  assemble_output(iopub)
            # the message can contain "Loading required package feathre"
            #self.assertEqual(stderr, '')
            #  execute(kc=kc, code="%use sos")
            # wait_for_idle(kc)
            # execute(kc=kc, code="mtcars.shape")
            #res = get_result(iopub)
            #self.assertEqual(res, (32, 11))
            #execute(kc=kc, code="mtcars.index[0]")
            #res = get_result(iopub)
#self.assertEqual(res, 'Mazda RX4')

    def testPutOctaveDataToPython(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            execute(kc=kc, code='%use Octave')
            wait_for_idle(kc)
            execute(kc=kc, code="null_var = NaN")
            wait_for_idle(kc)
            execute(kc=kc, code="num_var = 123")
            wait_for_idle(kc)
            execute(kc=kc, code="num_arr_var = [1, 2, 3]")
            wait_for_idle(kc)
            execute(kc=kc, code="logic_var = true")
            wait_for_idle(kc)
            execute(kc=kc, code="logic_arr_var = [true, false, true]")
            wait_for_idle(kc)
            execute(kc=kc, code="char_var = '1\"23'")
            wait_for_idle(kc)
            execute(kc=kc, code="char_arr_var = [1, 2, '3']")
            wait_for_idle(kc)
            execute(kc=kc, code="list_var = list(1, 2, '3')")
            wait_for_idle(kc)
            # Not in Octave
            #execute(kc=kc, code="named_list_var = list(a=1, b=2, c='3')")
            #wait_for_idle(kc)
            execute(kc=kc, code="mat_var = [1:3; 2:4]")
            wait_for_idle(kc)
            # Not in Octave
            #execute(kc=kc, code="recursive_var = list(a=1, b=list(c=3, d='whatever'))")
            #wait_for_idle(kc)
            execute(kc=kc, code="%put null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var list_var")
            wait_for_idle(kc)
            execute(kc=kc, code="%dict null_var num_var num_arr_var logic_var logic_arr_var char_var char_arr_var mat_var list_var")
            res = get_result(iopub)
            self.assertEqual(res['null_var'], None)
            self.assertEqual(res['num_var'], 123)
            self.assertEqual(res['num_arr_var'], [1,2,3])
            self.assertEqual(res['logic_var'], True)
            self.assertEqual(res['logic_arr_var'], [True, False, True])
            self.assertEqual(res['char_var'], '1"23')
            self.assertEqual(res['char_arr_var'], ['1', '2', '3'])
            self.assertEqual(res['list_var'], [1,2,'3'])
            #self.assertEqual(res['named_list_var'], {'a': 1, 'b': 2, 'c': '3'})
            self.assertEqual(res['mat_var'].shape, (2,2))
            #self.assertEqual(res['recursive_var'], {'a': 1, 'b': {'c': 3, 'd': 'whatever'}})
            execute(kc=kc, code="%use sos")
            wait_for_idle(kc)


if __name__ == '__main__':
    unittest.main()

