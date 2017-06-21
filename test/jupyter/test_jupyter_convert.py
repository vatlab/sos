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

import os
import sys
import unittest
import shutil
import subprocess

from sos.utils import env
from sos.jupyter.converter import script_to_notebook, notebook_to_script

file_dir = os.path.split(__file__)[0]
if not file_dir:
    file_dir = '.'

class TestConvert(unittest.TestCase):
    def setUp(self):
        env.reset()
        if not os.path.isdir('temp'):
            os.mkdir('temp')
        with open('temp/script1.sos', 'w') as script:
            script.write('''
[0]
seq = range(3)
input: for_each='seq'
output: 'test${_seq}.txt'
print(output)
''')
        with open('temp/script2.sos', 'w') as script:
            # with tab after run:
            script.write('''
[0]
seq = range(3)
input: for_each='seq'
output: 'test${_seq}.txt'
run:			concurrent=True
    echo 'this is test script'
[10]
report('this is action report')
''')
        self.scripts = ['temp/script1.sos', 'temp/script2.sos']

    def tearDown(self):
        shutil.rmtree('temp')

    def testScriptToAndFromNotebook(self):
        '''Test sos show script --notebook'''
        for script_file in self.scripts:
            script_to_notebook(script_file, script_file + '.ipynb')
            notebook_to_script(script_file + '.ipynb', script_file)

    def testConvertAll(self):
        olddir = os.getcwd()
        os.chdir(file_dir)
        subprocess.call('sos convert test.ipynb test_wf.sos --all', shell=True)
        self.assertTrue(os.path.isfile('test_wf.sos'))
        subprocess.call('sos convert test_wf.sos test2.ipynb', shell=True)
        self.assertTrue(os.path.isfile('test2.ipynb'))
        # --execute does not yet work
        os.remove('test_wf.sos')
        os.remove('test2.ipynb')
        os.chdir(olddir)

    def testConvertHTML(self):
        olddir = os.getcwd()
        os.chdir(file_dir)
        subprocess.call('sos convert test.ipynb test_wf.html', shell=True)
        self.assertTrue(os.path.isfile('test_wf.html'))
        os.chdir(olddir)

    @unittest.skipIf(sys.platform == 'win32', 'No XeLatex under windows to compile pdf')
    def testConvertPDF(self):
        olddir = os.getcwd()
        os.chdir(file_dir)
        subprocess.call('sos convert test.ipynb test_wf.pdf', shell=True)
        self.assertTrue(os.path.isfile('test_wf.pdf'))
        os.chdir(olddir)

    def testConvertMD(self):
        olddir = os.getcwd()
        os.chdir(file_dir)
        subprocess.call('sos convert test.ipynb test_wf.md', shell=True)
        self.assertTrue(os.path.isfile('test_wf.md'))
        os.chdir(olddir)

if __name__ == '__main__':
    #suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestConvert)
    #unittest.TextTestRunner().run(suite)
    unittest.main()
