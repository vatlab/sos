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
import unittest
import os
import subprocess
import tempfile

class TestWorkflow(unittest.TestCase):
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

    def testSoSSave(self):
        for f in ('test_wf.sos', 'test.out'):
            if os.path.isfile(f):
                os.remove(f)
        #
        subprocess.call('sos run test.ipynb test', shell=True)
        #
        self.assertTrue(os.path.isfile('test.out'))
        with open('test.out') as to:
            self.assertEqual(to.read(), 'test line \n')

        #
        subprocess.call('sos convert test.ipynb test_wf.sos', shell=True)
        self.assertTrue(os.path.isfile('test_wf.sos'))
        with open('test_wf.sos') as to:
            self.assertEqual(to.read(), '''\
#!/usr/bin/env sos-runner
#fileformat=SOS1.0

#! ## Notebook for testing purpose

#! ## Section 1

[10]
a = 100

[test]
output: 'test.out'
sh:
   echo "test line " >> ${output}

''')
        for f in ('test_wf.sos', 'test.out'):
            if os.path.isfile(f):
                os.remove(f)




if __name__ == '__main__':
    unittest.main()
