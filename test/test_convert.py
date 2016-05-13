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

import os
import unittest
import shutil
import subprocess

from pysos.utils import env
from pysos.sos_script import SoS_Script, ParsingError
from pysos.converter import script_to_html, script_to_markdown, script_to_term, script_to_notebook, \
    workflow_to_html, workflow_to_markdown, workflow_to_term, workflow_to_notebook, notebook_to_script


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
! This is a report
report('this is action report')
''')
        self.scripts = ['temp/script1.sos', 'temp/script2.sos']

    def tearDown(self):
        shutil.rmtree('temp')

    def testScriptToHtml(self):
        '''Test sos show script --html'''
        for script_file in self.scripts:
            with open('temp/test.transcript', 'w') as transcript:
                script = SoS_Script(filename=script_file, transcript=transcript)
            script_to_html('temp/test.transcript', script_file, script_file + '.html', [])
    
    def testScriptToMarkdown(self):
        '''Test sos show script --markdown'''
        for script_file in self.scripts:
            with open('temp/test.transcript', 'w') as transcript:
                script = SoS_Script(filename=script_file, transcript=transcript)
            script_to_markdown('temp/test.transcript', script_file, script_file + '.md')

    def testScriptToAndFromNotebook(self):
        '''Test sos show script --notebook'''
        for script_file in self.scripts:
            with open('temp/test.transcript', 'w') as transcript:
                script = SoS_Script(filename=script_file, transcript=transcript)
            script_to_notebook('temp/test.transcript', script_file, script_file + '.ipynb')
            #
            notebook_to_script(script_file + '.ipynb', script_file) 

    def testWorkflowToHtml(self):
        '''Test sos show script --html'''
        for script_file in self.scripts:
            script = SoS_Script(filename=script_file)
            wf = script.workflow()
            workflow_to_html(wf, script_file, script_file + '.html', [])
    
    def testWorkflowToMarkdown(self):
        '''Test sos show script --markdown'''
        for script_file in self.scripts:
            script = SoS_Script(filename=script_file)
            wf = script.workflow()
            workflow_to_markdown(wf, script_file, script_file + '.md')

    def testWorkflowToNotebook(self):
        '''Test converting sos script to notebook'''
        for script_file in self.scripts:
            script = SoS_Script(filename=script_file)
            wf = script.workflow()
            workflow_to_notebook(wf, script_file, script_file + '.ipynb')

if __name__ == '__main__':
    unittest.main()

