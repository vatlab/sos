#!/usr/bin/env python
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

import unittest
import time
import subprocess
import shutil

from sos.sos_script import SoS_Script
from sos.utils import env
from sos.sos_executor import Base_Executor
from sos.target import FileTarget

class TestTarget(unittest.TestCase):
    def setUp(self):
        env.reset()
        subprocess.call('sos remove -s', shell=True)
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            FileTarget(f).remove('both')

    def testRLibrary(self):
        '''Test target R_Library'''
        script = SoS_Script('''
[default]
depends: R_library("dplyr")
R:
    library('dplyr')
''')
        wf = script.workflow()
        Base_Executor(wf).run()


    def testDependsRLibrary(self):
        '''Testing depending on R_library'''
        # first remove xtable package
        if not shutil.which('R'):
            return 
        subprocess.call('R CMD REMOVE xtable', shell=True)
        script = SoS_Script('''
[0]

depends: R_library('xtable')
R:
    library('xtable')
    ## Demonstrate data.frame
    tli.table <- xtable(cars)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    def testReexecution(self):
        '''Test re-execution of steps with R_library'''
        script = SoS_Script('''
[1]
depends: R_library("ggplot2")
output: '1.txt'
run:
    sleep 5
    touch ${output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        #st = time.time()
        Base_Executor(wf).run()
        #self.assertLess(time.time() - st, 4)

if __name__ == '__main__':
    unittest.main()
