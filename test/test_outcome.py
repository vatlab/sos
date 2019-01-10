#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import glob
import os
import sys
import shutil
import subprocess
import unittest

from sos._version import __version__
from sos.parser import SoS_Script
from sos.targets import file_target
from sos.utils import env
# if the test is imported under sos/test, test interacive executor
from sos.workflow_executor import Base_Executor



class TestExecute(unittest.TestCase):
    def setUp(self):
        env.reset()
        subprocess.call('sos remove -s', shell=True)
        # self.resetDir('~/.sos')
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            if file_target(f).exists():
                file_target(f).unlink()

    def touch(self, files):
        '''create temporary files'''
        if isinstance(files, str):
            files = [files]
        #
        for f in files:
            with open(f, 'w') as tmp:
                tmp.write('test')
        #
        self.temp_files.extend(files)

    def resetDir(self, dirname):
        if os.path.isdir(os.path.expanduser(dirname)):
            shutil.rmtree(os.path.expanduser(dirname))
        os.mkdir(os.path.expanduser(dirname))

    
    def removeIfExists(self, targets):
        if isinstance(targets, str):
            targets = [targets]
        for target in targets:
            if os.path.isfile(target):
                os.remove(target)

    def testPlainTarget(self):
        '''Test sos run -t filename'''
        self.removeIfExists('t_a.txt')
        script = SoS_Script('''
[A]
output: 't_a.txt'
_output.touch()
''')
        wf = script.workflow()
        Base_Executor(wf, config={'target': 't_a.txt'}).run()
        self.assertTrue(os.path.isfile('t_a.txt'))

    def testTargetInNamed(self):
        '''Test sos run -t filename'''
        self.removeIfExists('t_na.txt')
        script = SoS_Script('''
[A]
output: res='t_na.txt'
_output.touch()
''')
        wf = script.workflow()
        Base_Executor(wf, config={'target': 't_na.txt'}).run()
        self.assertTrue(os.path.isfile('t_na.txt'))


if __name__ == '__main__':
    unittest.main()
