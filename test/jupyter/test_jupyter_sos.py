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
import subprocess
from ipykernel.tests.utils import execute, wait_for_idle, assemble_output
from sos.jupyter.test_utils import sos_kernel, get_result
from sos.target import FileTarget
from sos.utils import env

class TestJupyterSoS(unittest.TestCase):
    #
    # Beacuse these tests would be called from sos/test, we
    # should switch to this directory so that some location
    # dependent tests could run successfully
    #
    def setUp(self):
        self.olddir = os.getcwd()
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))
        subprocess.call('sos remove -s', shell=True)

    def tearDown(self):
        os.chdir(self.olddir)

    def testInterpolation(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%run
[a]
b=10

[default]
sos_run('a')
''')
            wait_for_idle(kc)
            execute(kc=kc, code="b")
            res = get_result(iopub)
            self.assertEqual(res, 10)

    def testReadonlyALLCAPVars(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%run
A=10

[default]
A=20
''')
            _, stderr = assemble_output(iopub)
            self.assertTrue('A' in stderr, 'Expect an error {}'.format(stderr))

    def testRerun(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%run
parameter: a=10

[default]
b = a
''')
            wait_for_idle(kc)
            #
            execute(kc=kc, code='''
%rerun --a 20
''')
            wait_for_idle(kc)
            execute(kc=kc, code="b")
            res = get_result(iopub)
            self.assertEqual(res, 20)


    def testDAG(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%run
[a]
b=10

[default]
sos_run('a')
''')
            wait_for_idle(kc)
            execute(kc=kc, code="b")
            res = get_result(iopub)
            self.assertEqual(res, 10)

    def testTarget(self):
        for f in ['A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
            FileTarget(f).remove('both')
        #
        #  A1 <- B1 <- B2 <- B3
        #   |
        #   |
        #  \/
        #  A2 <- B2 <- C1 <- C2 <- C4
        #                    C3
        #
        script = '''\
%run -t B1.txt -s force
[A_1]
input: 'B1.txt'
output: 'A1.txt'
run:
    touch A1.txt

[A_2]
depends:  'B2.txt'
run:
    touch A2.txt

[B1: provides='B1.txt']
depends: 'B2.txt'
run:
    touch B1.txt

[B2: provides='B2.txt']
depends: 'B3.txt', 'C1.txt'
run:
    touch B2.txt

[B3: provides='B3.txt']
run:
    touch B3.txt

[C1: provides='C1.txt']
depends: 'C2.txt', 'C3.txt'
run:
    touch C1.txt

[C2: provides='C2.txt']
depends: 'C4.txt'
run:
    touch C2.txt

[C3: provides='C3.txt']
depends: 'C4.txt'
run:
    touch C3.txt

[C4: provides='C4.txt']
run:
    touch C4.txt

        '''
        script2 = '''
import os
fail = 0
for f in ['A1.txt', 'A2.txt']:
    fail += os.path.exists(f)
for f in ['C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
    fail += not os.path.exists(f)
fail
'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code=script)
            wait_for_idle(kc)
            execute(kc=kc, code=script2)
            res = get_result(iopub)
            self.assertEqual(res, 0)


    def testReverseSharedVariable(self):
        '''Test shared variables defined in auxiliary steps'''
        FileTarget('a.txt').remove('both')
        script = r'''
%run B
[A: shared='b', provides='a.txt']
b = 1
run:
    touch a.txt

[B_1]
depends: 'a.txt'

[B_2]
print(b)

'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code=script)
            wait_for_idle(kc)
            execute(kc=kc, code="b")
            res = get_result(iopub)
            self.assertEqual(res, 1)

if __name__ == '__main__':
    unittest.main()
