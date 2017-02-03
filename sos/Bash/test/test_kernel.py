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
from ipykernel.tests.utils import execute, wait_for_idle
from sos.jupyter.test.utils import sos_kernel, get_result

try:
    import feather
    feather
    with_feather = True
except ImportError:
    with_feather = False

class TestKernel(unittest.TestCase):
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

    @unittest.skipIf(not with_feather, 'Skip test because of no feather module')
    def testGetDataFromBash(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code='''
null_var = None
num_var = 123
logic_var = True
char_var = '123'
char_arr_var = ['1', '2', '3']
list_var = [1, 2, '3']
dict_var = dict(a=1, b=2, c='3')
set_var = {1, 2, '3'}
''')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use Bash")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%get null_var num_var logic_var char_var char_arr_var list_var dict_var set_var")
            wait_for_idle(kc)
            # need to test passed values
            # but let us cheat by passing data back
            msg_id, content = execute(kc=kc, code="%dict -r")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put null_var num_var logic_var char_var char_arr_var list_var dict_var set_var")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict null_var num_var logic_var char_var char_arr_var list_var dict_var set_var")
            res = get_result(iopub)
            self.assertEqual(res['null_var'], '')
            self.assertEqual(res['num_var'], '123')
            self.assertEqual(res['logic_var'], 'TRUE')
            self.assertEqual(res['char_var'], '123')
            self.assertEqual(res['char_arr_var'], '1 2 3')
            self.assertEqual(res['list_var'], '1 2 3')
            self.assertEqual(sorted(res['dict_var']), [' ', ' ', 'a', 'b', 'c'])
            self.assertEqual(sorted(res['set_var']), [' ', ' ', '1', '2', '3'])

    @unittest.skipIf(not with_feather, 'Skip test because of no feather module')
    def testPutBashDataToPython(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='%use Bash')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="export null_var=")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="export num_var=123")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put null_var num_var")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict null_var num_var")
            res = get_result(iopub)
            self.assertEqual(res['null_var'], '')
            self.assertEqual(res['num_var'], '123')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)


if __name__ == '__main__':
    unittest.main()
