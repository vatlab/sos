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
from sos.jupyter.test_utils import sos_kernel, get_result, get_display_data

class TestKernel(unittest.TestCase):
    #
    def testMagicSave(self):
        with sos_kernel() as kc:
            if os.path.isfile('test.txt'):
                os.remove('test.txt')
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%save test.txt
print("a=${100+11}")
''')
            stdout, stderr = assemble_output(iopub)
            self.assertTrue(stdout.endswith('a=111\n'))
            self.assertEqual(stderr, '')
            with open('test.txt') as tt:
                self.assertEqual(tt.read(), 'print("a=${100+11}")\n')

if __name__ == '__main__':
    unittest.main()
