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
import unittest

from sos.sos_script import SoS_Script
from sos.utils import env
from sos.sos_executor import Base_Executor
from sos.target import FileTarget
from sos.hosts import Host
import subprocess

try:
    subprocess.check_output('docker ps | grep test_sos', shell=True).decode()
except subprocess.CalledProcessError:
    subprocess.call('sh build_test_docker.sh', shell=True)
    try:
        subprocess.check_output('docker ps | grep test_sos', shell=True).decode()
    except subprocess.CalledProcessError:
        sys.exit('Failed to set up a docker machine with sos')

class TestRemote(unittest.TestCase):
    def setUp(self):
        env.reset()
        #self.resetDir('~/.sos')
        self.temp_files = []
        Host.reset()
        # remove .status file left by failed workflows.
        subprocess.call('rm -f ~/.sos/*.status', shell=True)

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
                'sig_mode': 'force',
                }).run()
        self.assertTrue(FileTarget('result.txt').exists())
        with open('result.txt') as res:
            self.assertEqual(res.read(), 'a\n')

    def testRemoteExecution(self):
        subprocess.check_output('cd ~/.sos/tasks; rm -f *.res *.sh *.pulse', shell=True).decode()
        script = SoS_Script('''
[10]
input: for_each={'i': range(5)}
task:

run:
    echo I am ${i}
    sleep ${5+i}
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
        # wait another 15 seconds?
        time.sleep(15)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']), 'Expect all completed jobs: ' + out)

        Host.reset()
        # until we run the workflow again
        #st = time.time()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker',
                'resume_mode': True,
                }).run()
        # should finish relatively fast?
        #self.assertLess(time.time() - st, 5)
        out = subprocess.check_output('sos status {} -c docker.yml'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']), 'Expect all completed jobs: ' + out)

    def testTaskSpooler(self):
        subprocess.check_output('cd ~/.sos/tasks; rm -f *.res *.sh *.pulse', shell=True).decode()
        script = SoS_Script('''
[10]
input: for_each={'i': range(3)}
task:

run:
    echo I am task spooler ${i}
    sleep ${5+i*2}
''')
        wf = script.workflow()
        res = Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': False,
                'default_queue': 'ts',
                'sig_mode': 'force',
                }).run()
        import time
        # we should be able to get status
        tasks = ' '.join(res['pending_tasks'])
        time.sleep(2)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertGreaterEqual(out.count('running'), 1, 'Expect at least one running job: ' + out)
        # wait another 20 seconds?
        time.sleep(15)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']))
        # until we run the workflow again
        #st = time.time()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'ts',
                'resume_mode': True,
                }).run()
        # should finish relatively fast?
        #self.assertLess(time.time() - st, 5)
        out = subprocess.check_output('sos status {} -c docker.yml'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']), 'Expect all completed jobs: ' + out)

    def testTaskSpoolerWithForceSigMode(self):
        subprocess.check_output('cd ~/.sos/tasks; rm -f *.res *.sh *.pulse', shell=True).decode()
        script = SoS_Script('''
[10]
input: for_each={'i': range(3)}
task:

run:
    echo I am spooler with force ${i}
    sleep ${10 + i*2}
''')
        wf = script.workflow()
        res = Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': False,
                'default_queue': 'ts',
                'sig_mode': 'force',
                }).run()
        import time
        # we should be able to get status
        tasks = ' '.join(res['pending_tasks'])
        time.sleep(2)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertGreaterEqual(out.count('running'), 1, 'Expect at least one running jobs: ' + out)
        # wait another 20 seconds?
        time.sleep(20)
        out = subprocess.check_output('sos status {} -c docker.yml -q docker'.format(tasks), shell=True).decode()
        self.assertEqual(out.count('completed'), len(res['pending_tasks']), 'Expect all completed jobs: ' + out)
        # until we run the workflow again
        st = time.time()
        Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'ts',
                #
                # This is the only difference, because running with -s force would still 
                # skip some of the completed task.
                'sig_mode': 'force',
                'resume_mode': True,
                }).run()
        # should finish relatively fast?
        self.assertLess(time.time() - st, 9)
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

    @unittest.skipIf(not os.path.exists(os.path.expanduser('~').upper()),
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

    def testMaxMem(self):
        '''Test server restriction max_mem'''
        script = SoS_Script('''
[10]
task: mem='2G'
print('a')
'''.format(os.path.join(os.path.abspath('.').upper(), 'test_remote.py')))
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker_limited',
                'sig_mode': 'force',
                }).run)

    def testRuntimeMaxWalltime(self):
        '''Test server max_walltime option'''
        script = SoS_Script('''
[10]
task:
import time
time.sleep(15)
'''.format(os.path.join(os.path.abspath('.').upper(), 'test_remote.py')))
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker_limited',
                'sig_mode': 'force',
                }).run)

    def testMaxProcs(self):
        '''Test server restriction max_procs'''
        script = SoS_Script('''
[10]
task: procs=8
print('a')
'''.format(os.path.join(os.path.abspath('.').upper(), 'test_remote.py')))
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker_limited',
                'sig_mode': 'force',
                }).run)

    def testMaxWalltime(self):
        '''Test server restriction max_walltime'''
        script = SoS_Script('''
[10]
task: walltime='1:00:00'
print('a')
'''.format(os.path.join(os.path.abspath('.').upper(), 'test_remote.py')))
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf, config={
                'config_file': 'docker.yml',
                # do not wait for jobs
                'wait_for_task': True,
                'default_queue': 'docker_limited',
                'sig_mode': 'force',
                }).run)


if __name__ == '__main__':
    unittest.main()
