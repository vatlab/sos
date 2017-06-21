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

class TestSoSMagics(unittest.TestCase):
    #
    def testMagicShutdown(self):
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%use R
a = 100
cat(a)
''')
            stdout, stderr = assemble_output(iopub)
            self.assertTrue(stdout.endswith('100'), 'Should have output {}'.format(stdout))
            self.assertEqual(stderr, '')
            # now let us restart
            execute(kc=kc, code='''
%shutdown --restart R
%use R
cat(a)
''')
            stdout, stderr = assemble_output(iopub)
            # not sure what is going on
            # we should have error message, right?
            self.assertEqual(stdout, '')
            self.assertEqual(stderr, '')
            execute(kc=kc, code='%use SoS')
            wait_for_idle(kc)

    def testMagicSave(self):
        with sos_kernel() as kc:
            if os.path.isfile('test.txt'):
                os.remove('test.txt')
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%save ~/test.txt
a=1
''')
            wait_for_idle(kc)
            with open(os.path.join(os.path.expanduser('~'), 'test.txt')) as tt:
                self.assertEqual(tt.read(), 'a=1\n')

if __name__ == '__main__':
    unittest.main()
