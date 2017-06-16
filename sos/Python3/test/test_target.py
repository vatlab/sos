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

from sos.utils import env
from sos.sos_script import SoS_Script
from sos.sos_executor import Base_Executor
from sos.target import FileTarget

class TestTarget(unittest.TestCase):
    def setUp(self):
        env.reset()
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            FileTarget(f).remove('both')

    def testPy_Module(self):
        '''Test target Py_Module'''
        FileTarget('report.md').remove('both')
        script = SoS_Script(r'''
[10]
depends: Py_Module('tabulate')
from tabulate import tabulate
table = [["Sun",696000,1989100000],["Earth",6371,5973.6],
    ["Moon",1737,73.5],["Mars",3390,641.85]]
report: output='report.md'
    ${tabulate(table)}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(FileTarget('report.md').exists('target'))
        with open('report.md') as rep:
            self.assertEqual(rep.read().strip(), '''
-----  ------  -------------
Sun    696000     1.9891e+09
Earth    6371  5973.6
Moon     1737    73.5
Mars     3390   641.85
-----  ------  -------------
'''.strip())


if __name__ == '__main__':
    unittest.main()
