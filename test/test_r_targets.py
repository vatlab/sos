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
import subprocess
import shutil

from sos.parser import SoS_Script
from sos.utils import env
from sos.workflow_executor import Base_Executor
from sos.targets import file_target

class TestTarget(unittest.TestCase):
    def setUp(self):
        env.reset()
        subprocess.call('sos remove -s', shell=True)
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            file_target(f).remove('both')

    @unittest.skipIf(not shutil.which('octave'), 'Octave not installed')
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


    @unittest.skipIf(not shutil.which('Rscript'), 'R not installed')
    def testDependsRLibrary(self):
        '''Testing depending on R_library'''
        # first remove xtable package
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

    @unittest.skipIf(not shutil.which('Rscript'), 'R not installed')
    def testReexecution(self):
        '''Test re-execution of steps with R_library'''
        script = SoS_Script('''
[1]
depends: R_library("ggplot2", "2.2+")
output: '1.txt'
run: expand=True
    sleep 5
    touch {_output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        Base_Executor(wf).run()

if __name__ == '__main__':
    unittest.main()
