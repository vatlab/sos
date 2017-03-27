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
import time
import unittest
import shutil
import glob

from sos.sos_script import SoS_Script, ParsingError
from sos.utils import env
from sos.sos_executor import Base_Executor
from sos.target import FileTarget
from sos.hosts import Host
import subprocess

try:
    subprocess.check_output('docker ps | grep test_sshd', shell=True).decode()
except subprocess.CalledProcessError:
    subprocess.call('sh build_test_docker.sh', shell=True)
    try:
        subprocess.check_output('docker ps | grep test_sshd', shell=True).decode()
    except subprocess.CalledProcessError:
        sys.exit('Failed to set up a docker machine with sos')

class TestRemote(unittest.TestCase):
    def setUp(self):
        env.reset()
        #self.resetDir('~/.sos')
        self.temp_files = []
        Host.reset()

    def tearDown(self):
        for f in self.temp_files:
            FileTarget(f).remove('both')


    def testRemoteExecute(self):
        script = SoS_Script('''
[10]
output: 'result.txt'
task:

run:
  echo 'a' > 'result.txt'

''')
        wf = script.workflow()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                'wait_for_task': True,
                'default_queue': 'docker',
                }).run()
        self.assertTrue(FileTarget('result.txt').exists())
        with open('result.txt') as res:
            self.assertEqual(res.read(), 'a\n')

if __name__ == '__main__':
    unittest.main()
