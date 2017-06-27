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

#from sos.sos_script import SoS_Script
from sos.target import FileTarget
from sos.sos_script import SoS_Script
from sos.sos_executor import Base_Executor

import unittest
import subprocess


has_docker = True
try:
    subprocess.check_output('docker ps | grep test_sos', shell=True).decode()
except subprocess.CalledProcessError:
    subprocess.call('sh ../build_test_docker.sh', shell=True)
    try:
        subprocess.check_output('docker ps | grep test_sos', shell=True).decode()
    except subprocess.CalledProcessError:
        print('Failed to set up a docker machine with sos')
        has_docker = False

class TestRQQueue(unittest.TestCase):

    @unittest.skipIf(not has_docker, "Docker container not usable")
    def testRemoteExecute(self):
        script = SoS_Script('''
[10]
output: 'result.txt'
task:

run:
  echo 'rq' > 'result.txt'

''')
        wf = script.workflow()
        Base_Executor(wf, config={
                'config_file': '~/docker.yml',
                'wait_for_task': True,
                'default_queue': 'docker_rq',
                'sig_mode': 'force',
                }).run()
        self.assertTrue(FileTarget('result.txt').exists())
        with open('result.txt') as res:
            self.assertEqual(res.read(), 'rq\n')

if __name__ == '__main__':
    unittest.main()
