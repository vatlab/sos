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
import time
import glob
import unittest
import shutil

from pysos import SoS_Script
from pysos.dag import SoS_DAG
from pysos.utils import env
from pysos.sos_eval import Undetermined
from pysos.sos_executor import Sequential_Executor, Interactive_Executor, ExecuteError
from pysos.sos_script import ParsingError


import matplotlib.pyplot as plt


class TestDAG(unittest.TestCase):
    def testLinearDag(self):
        '''Test DAG with linear dependency'''
        script = SoS_Script('''
[1]
input: 'a.txt'
output: 'b.txt'

[2]
input: 'b.txt'
output: 'c.txt'

[3]
input: 'c.txt'
output: 'd.txt'

[4]
input: 'd.txt'
output: 'e.txt'

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        dag.write_dot('linear.dot')

if __name__ == '__main__':
    unittest.main()
