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
import time
import unittest
import shutil
import glob

from sos.sos_script import SoS_Script, ParsingError
from sos.utils import env
from sos.sos_executor import Base_Executor
from sos.target import FileTarget
import subprocess

class TestTask(unittest.TestCase):
    def setUp(self):
        env.reset()
        subprocess.call('sos remove -s', shell=True)
        #self.resetDir('~/.sos')
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            FileTarget(f).remove('both')

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


    def testWorkdir(self):
        '''Test workdir option for runtime environment'''
        script =  SoS_Script(r"""
[0]
output: 'result.txt'
task: workdir='..'

with open('test/result.txt', 'w') as res:
   for file in os.listdir('test'):
       res.write(file + '\n')
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        with open('result.txt') as res:
            content = [x.strip() for x in res.readlines()]
            self.assertTrue('test_execute.py' in content)
        os.remove('result.txt')

    def testConcurrency(self):
        '''Test concurrency option for runtime environment'''
        env.max_jobs = 5
        env.sig_mode = 'force'
        script =  SoS_Script(r"""
[0]

repeat = range(4)
input: for_each='repeat'

task: concurrent=False

import time
time.sleep(_repeat + 1)
print('I am {}, waited {} seconds'.format(_index, _repeat + 1))
""")
        wf = script.workflow()
        start = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - start, 11)
        #
        #
        script =  SoS_Script(r"""
[0]

repeat = range(4)
input: for_each='repeat'

task: concurrent=True

if run_mode == 'run':
    import time
    time.sleep(_repeat + 1)
    print('I am {}, waited {} seconds'.format(_index, _repeat + 1))
""")
        wf = script.workflow()
        start = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - start, 9)

    def testPrependPath(self):
        '''Test prepend path'''
        import stat
        if not os.path.isdir('temp'):
            os.mkdir('temp')
        with open('temp/temp_cmd', 'w') as tc:
            tc.write('echo "a"')
        os.chmod('temp/temp_cmd', stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)
        #
        script = SoS_Script(r"""
[1]
task:
sh:
    temp_cmd
""")
        wf = script.workflow()
        env.sig_mode = 'force'
        self.assertRaises(Exception, Base_Executor(wf).run)
        # use option env
        script = SoS_Script(r"""
[1]
task: env={'PATH': 'temp' + os.pathsep + os.environ['PATH']}
sh:
    temp_cmd
""")
        wf = script.workflow()
        env.sig_mode = 'force'
        Base_Executor(wf).run()
        #
        #
        script = SoS_Script(r"""
[1]
task: prepend_path='temp'
sh:
    temp_cmd
""")
        wf = script.workflow()
        Base_Executor(wf).run()


    def testActiveActionOption(self):
        '''Test the active option of actions'''
        # disallow
        self.assertRaises(ParsingError, SoS_Script, '''
[1]
rep = range(5)
input: for_each = 'rep'
# ff should change and be usable inside run
ff = "${_rep}.txt"
run:  active=1,2
echo ${ff}
touch temp/${ff}
''')
        #
        for active, result in [
            ('0', ['temp/0.txt']),
            ('-1', ['temp/4.txt']),
            ('(1,2)', ['temp/1.txt', 'temp/2.txt']),
            ('[2,3]', ['temp/2.txt', 'temp/3.txt']),
            ('(0,2,4)', ['temp/0.txt', 'temp/2.txt', 'temp/4.txt']),
            ('slice(1,None)', ['temp/1.txt', 'temp/2.txt', 'temp/3.txt', 'temp/4.txt']),
            ('slice(1,-2)', ['temp/1.txt', 'temp/2.txt']),
            ('slice(None,None,2)', ['temp/0.txt', 'temp/2.txt', 'temp/4.txt']),
            ]:
            if os.path.isdir('temp'):
                shutil.rmtree('temp')
            os.mkdir('temp')
            # test first iteration
            script = SoS_Script('''
[1]
rep = range(5)
input: for_each = 'rep'
# ff should change and be usable inside run
ff = "${_rep}.txt"
run:  active=%s
echo ${ff}
touch temp/${ff}
''' % active)
            wf = script.workflow()
            env.sig_mode = 'force'
            Base_Executor(wf).run()
            files = list(glob.glob('temp/*.txt'))
            self.assertEqual(files, result)
            #
            # test last iteration
            shutil.rmtree('temp')
            #
            # test active option for task
            os.mkdir('temp')
            script = SoS_Script('''
[1]
rep = range(5)
input: for_each = 'rep'
# ff should change and be usable inside run
ff = "${_rep}.txt"
task:  active=%s
run:
echo ${ff}
touch temp/${ff}
''' % active)
            wf = script.workflow()
            env.sig_mode = 'force'
            Base_Executor(wf).run()
            files = list(glob.glob('temp/*.txt'))
            self.assertEqual(files, result)
            #
            # test last iteration
            shutil.rmtree('temp')


    def testNestedWorkdir(self):
        '''Test nested runtime option for work directory'''
        if os.path.isdir('tmp'):
            shutil.rmtree('tmp')
        script = SoS_Script('''
[step]
task: workdir='tmp'
bash:
    touch 'a.txt'

[default]
task: workdir='tmp'
sos_run('step')
''')
        wf = script.workflow()
        # this should be ok.
        Base_Executor(wf).run()
        os.path.isfile('tmp/tmp/a.txt')
        shutil.rmtree('tmp')



    def testPassingVarToTask(self):
        '''Test passing used variable to tasks'''
        for i in range(10, 13):
            FileTarget('myfile_{}.txt'.format(i)).remove('both')
        #
        env.sig_mode = 'force'
        script = SoS_Script(r'''
parameter: gvar = 10

[10]
# generate a file
tt = range(gvar, gvar + 3)
input: for_each='tt'
output: "myfile_${_tt}.txt"
# additional comment

# _tt should be used in task
task: concurrent=True
python:
    # ${gvar}
    with open(${_output!r}, 'w') as tmp:
        tmp.write('${_tt}_${_index}')

''')
        wf = script.workflow()
        env.max_jobs = 4
        Base_Executor(wf).run()
        for t in range(10, 13):
            with open('myfile_{}.txt'.format(t)) as tmp:
                self.assertEqual(tmp.read(), str(t) + '_' + str(t-10))
            FileTarget('myfile_{}.txt'.format(t)).remove('both')

    def testMaxJobs(self):
        '''Test default max number of jobs'''
        script = SoS_Script(r'''

[10]
input: for_each=[{'a': range(10)}, {'b': range(3)}]

task: concurrent=True
run:
    echo "a = ${a}, b = ${b}"
    sleep ${a + b}
''')
        wf = script.workflow()
        Base_Executor(wf).run()

if __name__ == '__main__':
    unittest.main()
