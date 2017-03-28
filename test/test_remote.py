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

    def testStatus(self):
        script = SoS_Script('''
[10]
input: for_each={'i': range(5)}
task:

run:
    echo I am ${i}
    sleep ${i*2}
''')
        wf = script.workflow()
        res = Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': False,
                'default_queue': 'docker',
                'sig_mode': 'force',
                }).run()
        import time
        # we should be able to get status
        tasks = ' '.join(res['pending_tasks'])
        time.sleep(2)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertGreater(out.count('running'), 1)
        # local should all be pending but we might have results from before
        subprocess.call('cd ~/.sos/tasks; rm -f {}'.format(' '.join(x+'.res' for x in res['pending_tasks'])), shell=True)
        # wait another 20 seconds?
        time.sleep(10)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']))
        # now, local status should still be pending
        out = subprocess.check_output('sos status {} -c docker.yml'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('pending'), len(res['pending_tasks']))
        # until we run the workflow again
        st = time.time()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker',
                }).run()
        # should finish relatively fast?
        #self.assertLess(time.time() - st, 5)
        out = subprocess.check_output('sos status {} -c docker.yml'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']))

    def testSendSymbolicLink(self):
        '''Test to_host symbolic link or directories that contain symbolic link. #508'''
        # create a symbloc link
        subprocess.call('ln -s scripts ll', shell=True)
        script = SoS_Script('''
[10]
task: to_host='ll'
files = os.listdir('ll')
''')
        wf = script.workflow()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker',
                }).run()
        os.remove('ll')

    @unittest.skipIf(not os.path.isfile(os.path.expanduser('~')),
            'Skip test for case sensitive file system')
    def testCaseInsensitiveLocalPath(self):
        '''Test path_map from a case insensitive file system.'''
        FileTarget('test_remote.py.bak').remove('both')
        script = SoS_Script('''
[10]
output: 'test_remote.py.bak'
task: to_host='{}'
run:
    cat test_remote.py > ${{output}}
'''.format(os.path.join(os.path.abspath('.').upper(), 'test_remote.py')))
        wf = script.workflow()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker',
                'sig_mode': 'force',
                }).run()
        self.assertTrue(FileTarget('test_remote.py.bak').exists('target'))
        # the files should be the same
        with open('test_remote.py') as ori, open('test_remote.py.bak') as bak:
            self.assertEqual(ori.read(), bak.read())
        FileTarget('test_remote.py.bak').remove('both')



if __name__ == '__main__':
    unittest.main()
