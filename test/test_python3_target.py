#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import unittest

from sos.parser import SoS_Script
from sos.targets import file_target
from sos.utils import env
from sos.workflow_executor import Base_Executor


class TestTarget(unittest.TestCase):
    def setUp(self):
        env.reset()
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            file_target(f).remove('both')

    def testPy_Module(self):
        '''Test target Py_Module'''
        file_target('report.md').remove('both')
        script = SoS_Script(r'''
[10]
depends: Py_Module('tabulate')
from tabulate import tabulate
table = [["Sun",696000,1989100000],["Earth",6371,5973.6],
    ["Moon",1737,73.5],["Mars",3390,641.85]]
report: output='report.md', expand=True
    {tabulate(table)}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(file_target('report.md').target_exists('target'))
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
