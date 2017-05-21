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
from sos.jupyter.test.utils import sos_kernel, get_result, get_display_data

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

    def testInterpolation(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code='print("a=${100+11}")')
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
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict --keys --all")
            res = get_result(iopub)
            self.assertTrue('a' not in res)
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
            msg_id, content = execute(kc=kc, code="%cd ..")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="print(os.getcwd())")
            stdout, stderr = assemble_output(iopub)
            self.assertFalse(stdout.strip().endswith('test'))
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="%cd test")
        
    def testMagicUse(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="%use R1 -l sos.R.kernel:sos_R -c #CCCCCC")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use R2 -k ir")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a <- 1024")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 1024')
            msg_id, content = execute(kc=kc, code="%use R3 -k ir -l R")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a <- 233")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 233')
            msg_id, content = execute(kc=kc, code="%use R2 -c red")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 1024')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
      
    def testSubKernel(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="a <- 1024")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_display_data(iopub)
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
            #msg_id, content = execute(kc=kc, code="%put __k_k")
            #wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_result(iopub)
            self.assertEqual(res, 1024)
            # strange name
            msg_id, content = execute(kc=kc, code="%use R")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code=".a.b <- 22")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put .a.b")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="_a_b")
            res = get_result(iopub)
            self.assertEqual(res, 22)
            #
            # test to yet another kernel
            #
            msg_id, content = execute(kc=kc, code="%put --to Python3 _a_b")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use Python3")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="_a_b")
            res = get_result(iopub)
            self.assertEqual(res, 22)
            #
            msg_id, content = execute(kc=kc, code="kkk = 'ast'")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put --to R kkk")
            res = get_result(iopub)
            msg_id, content = execute(kc=kc, code="%use R")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="kkk <- paste0(kkk, '1')")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%put --to Python3 kkk")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use Python3")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="kkk")
            res = get_result(iopub)
            self.assertEqual(res, 'ast1')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)


    def testMagicGet(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="a = 1025")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="_b_a = 22")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use R")
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            msg_id, content = execute(kc=kc, code="%get a")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="a")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 1025')
            msg_id, content = execute(kc=kc, code="b <- 122\nc<-555")
            wait_for_idle(kc)
            #
            msg_id, content = execute(kc=kc, code="%get _b_a")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code=".b_a")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 22')
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            #
            # get from a sub kernel
            msg_id, content = execute(kc=kc, code="%get --from R b")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="b")
            res = get_result(iopub)
            self.assertEqual(res, 122)
            # get from a third kernel
            msg_id, content = execute(kc=kc, code="%use Python3")
            wait_for_idle(kc) 
            msg_id, content = execute(kc=kc, code="%get --from R c")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="c")
            res = get_result(iopub)
            self.assertEqual(res, 555)
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)


    def testAutoSharedVars(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            msg_id, content = execute(kc=kc, code="sos_null = None")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="sos_num = 123")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use R")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="sos_num")
            res = get_display_data(iopub)
            self.assertEqual(res, '[1] 123')
            msg_id, content = execute(kc=kc, code="sos_num = sos_num + 10")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%use sos")
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="sos_num")
            res = get_display_data(iopub)
            self.assertEqual(res, '133')

    def testWith(self):
        '''Test magic with'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='var = [1, 2, 3, 4]')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code='%with R -i var -o m\nm=mean(var)')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict m")
            res = get_result(iopub)
            self.assertEqual(res['m'], 2.5)

    def testSetSigil(self):
        '''Test set_options of sigil'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # create a data frame
            msg_id, content = execute(kc=kc, code='%set_options sigil=None')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code='a="${}".format(100)')
            wait_for_idle(kc)
            msg_id, content = execute(kc=kc, code="%dict a")
            res = get_result(iopub)
            self.assertEqual(res['a'], "$100")

if __name__ == '__main__':
    unittest.main()
