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
import os
import unittest
from ipykernel.tests.utils import assemble_output, execute, wait_for_idle
from sos.jupyter.test_utils import sos_kernel, get_result, get_display_data, \
    clear_channels

class TestJuliaKernel(unittest.TestCase):
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

    def testGetPythonDataFrameFromJulia(self):
        # Python -> Julia
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            execute(kc=kc, code='''
import pandas as pd
import numpy as np
arr = np.random.randn(1000)
arr[::10] = np.nan
df = pd.DataFrame({'column_{0}'.format(i): arr for i in range(10)})
''')
            clear_channels(iopub)
            execute(kc=kc, code="%use Julia")
            wait_for_idle(kc)
            execute(kc=kc, code="%get df")
            wait_for_idle(kc)
            execute(kc=kc, code="size(df)")
            res = get_display_data(iopub)
            self.assertEqual(res, '(1000, 10)')
            execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            #

    def testGetPythonDataFromJulia(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
num_var = 123
import numpy
import pandas
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
comp_var = 1+2j
seri_var = pandas.Series([1,2,3,3,3,3])
''')
            wait_for_idle(kc)
            execute(kc=kc, code='''
%use Julia
%get num_var num_arr_var logic_var logic_arr_var char_var char_arr_var set_var list_var dict_var recursive_var comp_var seri_var
%dict -r
%put num_var num_arr_var logic_var logic_arr_var char_var char_arr_var set_var list_var dict_var recursive_var comp_var seri_var
%use sos
seri_var = list(seri_var)
''')
            wait_for_idle(kc)
            execute(kc=kc, code='''
%dict num_var num_arr_var logic_var logic_arr_var char_var char_arr_var set_var list_var dict_var recursive_var comp_var seri_var
''')
            res = get_result(iopub)
            #self.assertEqual(res['null_var'], None)
            self.assertEqual(res['num_var'], 123)
            self.assertEqual(res['num_arr_var'], [1,2,3])
            self.assertEqual(res['logic_var'], True)
            self.assertEqual(res['logic_arr_var'], [True, False, True])
            self.assertEqual(res['char_var'], '1"23')
            self.assertEqual(res['char_arr_var'], ['1', '2', '3'])
            self.assertEqual(res['set_var'], {1, 2, '3'})
            self.assertEqual(res['list_var'], [1,2,'3'])
            self.assertEqual(res['dict_var'], {'a': 1, 'b': 2, 'c': '3'})
            #self.assertEqual(res['mat_var'].shape, (2, 2))
            self.assertEqual(res['recursive_var'],  {'a': {'b': 123}, 'c': True})
            self.assertEqual(res['comp_var'], 1+2j)
            self.assertEqual(res['seri_var'], [1,2,3,3,3,3])

#dataframe

    def testPutJuliaDataToPython(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            execute(kc=kc, code="%use Julia")
            wait_for_idle(kc)
            #execute(kc=kc, code="null_var = NaN")
            #wait_for_idle(kc)
            execute(kc=kc, code="num_var = 123")
            wait_for_idle(kc)
            execute(kc=kc, code="num_arr_var = [1, 2, 3]")
            wait_for_idle(kc)
            execute(kc=kc, code="logic_var = true")
            wait_for_idle(kc)
            execute(kc=kc, code="logic_arr_var = [true, true, false]")
            wait_for_idle(kc)
            execute(kc=kc, code='''char_var = "1\"23"''')
            wait_for_idle(kc)
            execute(kc=kc, code='''char_arr_var = [1, 2, "3"]''')
            wait_for_idle(kc)
            #execute(kc=kc, code='''named_list_var = NamedArray([1,2,3],(["a","b","c"],))''')
            #wait_for_idle(kc)
            execute(kc=kc, code="mat_var = [1 2; 3 4]")
            wait_for_idle(kc)
            execute(kc=kc, code='''recursive_var = Dict("a" => 1, "b" => Dict("c" => 3),"d" => "whatever")''')
            wait_for_idle(kc)
            execute(kc=kc, code="comp_var = 1+2im")
            wait_for_idle(kc)
            #execute(kc=kc, code="seri_var = setNames(c(1,2,3,3,3,3),c(0:5))")
            #wait_for_idle(kc)
            execute(kc=kc, code="%put num_var num_arr_var logic_var logic_arr_var char_var char_arr_var recursive_var comp_var")
            wait_for_idle(kc)
            execute(kc=kc, code='''
%use sos
seri_var = list(seri_var)
''')
            wait_for_idle(kc)
#           execute(kc=kc, code='''
#named_list_var = list(named_list_var)
#''')
#            wait_for_idle(kc)
            execute(kc=kc, code="%dict num_var num_arr_var logic_var logic_arr_var char_var char_arr_var recursive_var comp_var")
            res = get_result(iopub)
            #self.assertEqual(res['null_var'], None)
            self.assertEqual(res['num_var'], 123)
            self.assertEqual(res['num_arr_var'], [1,2,3])
            self.assertEqual(res['logic_var'], True)
            self.assertEqual(res['logic_arr_var'], [True, True, False])
            self.assertEqual(res['char_var'], '1"23')
            self.assertEqual(res['char_arr_var'], [1, 2, '3'])
            #self.assertEqual(res['named_list_var'], [1, 2, 3])
            #self.assertEqual(res['mat_var'].shape, (2, 2))
            self.assertEqual(res['recursive_var'], {'a': 1, 'b': {'c': 3}, 'd': 'whatever'})
            self.assertEqual(res['comp_var'], 1+2j)
            #self.assertEqual(res['seri_var'], [1,2,3,3,3,3])
            execute(kc=kc, code="%use sos")
            wait_for_idle(kc)


if __name__ == '__main__':
    unittest.main()
