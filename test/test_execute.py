#!/usr/bin/env python
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

# passing string as unicode to python 2 version of SoS
# to ensure compatibility
from __future__ import unicode_literals

import os
import unittest

from pysos import *
from pysos import __version__
import subprocess

class TestRun(unittest.TestCase):
    def testCommandLine(self):
        '''Test command line arguments'''
        result = subprocess.check_output('sos --version', stderr=subprocess.STDOUT, shell=True).decode()
        self.assertTrue(result.startswith('SoS {}'.format(__version__)))
        if hasattr(subprocess, 'DEVNULL'):
            devnull = subprocess.DEVNULL
        else:
            devnull = open(os.devnull, 'w')
        self.assertEqual(subprocess.call('sos -h', stderr=devnull, stdout=devnull, shell=True), 0)
        self.assertEqual(subprocess.call('sos run -h', stderr=devnull, stdout=devnull, shell=True), 0)
        self.assertEqual(subprocess.call('sos dryrun -h', stderr=devnull, stdout=devnull, shell=True), 0)
        self.assertEqual(subprocess.call('sos dryrun scripts/master.sos', stderr=devnull, stdout=devnull, shell=True), 1)
        self.assertEqual(subprocess.call('sos dryrun scripts/master.sos L', stderr=devnull, stdout=devnull, shell=True), 0)
        self.assertEqual(subprocess.call('sos show -h', stderr=devnull, stdout=devnull, shell=True), 0)

    def testInterpolation(self):
        '''Test string interpolation during execution'''
        script = SoS_Script(r"""
[0]
res = ''
b = 200
res += '${b}'
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['res'], '200')
        #
        script = SoS_Script(r"""
[0]
res = ''
for b in range(5):
    res += '${b}'
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['res'], '01234')
        
    def testGlobalVars(self):
        '''Test SoS defined variables'''
        script = SoS_Script(r"""
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['SOS_VERSION'], __version__)

    def testSignature(self):
        '''Test recognizing the format of SoS script'''
        env.run_mode = 'run'
        script = SoS_Script(r"""
[*_0]
output: 'temp/a.txt', 'temp/b.txt'

run('''echo "a.txt" > 'temp/a.txt' ''')
run('''echo "b.txt" > 'temp/b.txt' ''')

[1: alias='oa']
dest = ['temp/c.txt', 'temp/d.txt']
input: group_by='single', paired_with='dest'
output: dest

run(''' cp ${_input} ${_dest} ''')
""")
        wf = script.workflow('default_0')
        wf.run()
        # not the default value of 1.0
        self.assertTrue(os.path.isfile('temp/a.txt'))
        self.assertTrue(os.path.isfile('temp/b.txt'))
        with open('temp/a.txt') as ta:
            self.assertTrue(ta.read(), 'a.txt')
        with open('temp/b.txt') as tb:
            self.assertTrue(tb.read(), 'b.txt')
        #
        wf = script.workflow()
        wf.run()
        # not the default value of 1.0
        self.assertTrue(os.path.isfile('temp/c.txt'))
        self.assertTrue(os.path.isfile('temp/d.txt'))
        with open('temp/c.txt') as tc:
            self.assertTrue(tc.read(), 'a.txt')
        with open('temp/d.txt') as td:
            self.assertTrue(td.read(), 'b.txt')
        self.assertEqual(env.locals['oa'].output, ['temp/c.txt', 'temp/d.txt'])
        

    def testInput(self):
        '''Test input specification'''
        env.run_mode = 'dryrun'
        script = SoS_Script(r"""
[0]
input: '*.py'
output: _input
""")
        wf = script.workflow()
        wf.run()
        self.assertTrue('test_execute.py' in env.locals['_step'].output)

    def testForEach(self):
        '''Test for_each option of input'''
        env.run_mode = 'dryrun'
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
names = ['a', 'b', 'c']
c = ['1', '2']
counter = 0
all_names = ''
all_loop = ''

input: 'a.pdf', files, group_by='single', paired_with='names', for_each='c'

all_names += _names[0] + " "
all_loop += _c + " "

counter = counter + 1
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['counter'], 6)
        self.assertEqual(env.locals['all_names'], "a b c a b c ")
        self.assertEqual(env.locals['all_loop'], "1 1 1 2 2 2 ")
        #
        # test same-level for loop and parameter with nested list
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
par = [(1, 2), (1, 3), (2, 3)]
res = ['p1.txt', 'p2.txt', 'p3.txt']
processed = []

input: files, for_each='par,res'
output: res

processed.append((_par, _res))
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['processed'], [((1, 2), 'p1.txt'), ((1, 3), 'p2.txt'), ((2, 3), 'p3.txt')])


    def testPairedWith(self):
        '''Test option paired_with '''
        pass

    def testInputPattern(self):
        '''Test option pattern of step input '''
        env.run_mode = 'dryrun'
        #env.verbosity = 4
        script = SoS_Script(r"""
[0]

files = ['a-20.txt', 'b-10.txt']
input: files, pattern=['{name}-{par}.txt', '{base}.txt']
output: ['{}-{}-{}.txt'.format(x,y,z) for x,y,z in zip(_base, _name, _par)]

""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['base'], ["a-20", 'b-10'])
        self.assertEqual(env.locals['name'], ["a", 'b'])
        self.assertEqual(env.locals['par'], ["20", '10'])
        self.assertEqual(env.locals['_output'], ["a-20-a-20.txt", 'b-10-b-10.txt'])

    def testOutputPattern(self):
        '''Test option pattern of step output'''
        env.run_mode = 'dryrun'
        #env.verbosity = 4
        script = SoS_Script(r"""
[0]

files = ['a-20.txt', 'b-10.txt']
input: files, pattern=['{name}-{par}.txt', '{base}.txt']
output: pattern=['{base}-{name}-{par}.txt', '{par}.txt']

""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['base'], ["a-20", 'b-10'])
        self.assertEqual(env.locals['name'], ["a", 'b'])
        self.assertEqual(env.locals['par'], ["20", '10'])
        self.assertEqual(env.locals['_output'], ['a-20-a-20.txt', 'b-10-b-10.txt', '20.txt', '10.txt'])


    def testAlias(self):
        '''Test option alias'''
        env.run_mode = 'dryrun'
        script = SoS_Script(r"""
[0: alias='oa']
files = ['a.txt', 'b.txt']
names = ['a', 'b', 'c']
c = ['1', '2']
counter = "0"

input: 'a.pdf', files, group_by='single', paired_with='names', for_each='c'

counter = str(int(counter) + 1)
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['oa'].input, ["a.pdf", 'a.txt', 'b.txt'])

    def testFileType(self):
        '''Test input option filetype'''
        env.run_mode = 'dryrun'
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', files, filetype='*.txt', group_by='single'

output: _input

""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['_step'].output, ['a.txt', 'b.txt'])
        #
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', 'b.html', files, filetype=('*.txt', '*.pdf'), group_by='single'

counter += 1
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['counter'], 3)
        #
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', 'b.html', files, filetype=lambda x: 'a' in x, group_by='single'

counter += 1
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['counter'], 2)

    def testSkip(self):
        '''Test input option skip'''
        env.run_mode = 'dryrun'
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', 'b.html', files, skip=counter == 0

counter += 1
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['counter'], 0)

    def testOutputFromInput(self):
        '''Test deriving output files from input files'''
        env.run_mode = 'dryrun'
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
counter = 0

input: files, group_by='single'
output: _input[0] + '.bak'

counter += 1
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['counter'], 2)
        self.assertEqual(env.locals['_step'].output, ['a.txt.bak', 'b.txt.bak'])

    def testWorkdir(self):
        '''Test workdir option for runtime environment'''
        script =  SoS_Script(r"""
[0]

process: workdir='..'

with open('test/result.txt', 'w') as res:
   for file in os.listdir('test'):
       res.write(file + '\n')
""")
        wf = script.workflow()
        wf.run()
        content = [x.strip() for x in open('result.txt').readlines()]
        self.assertTrue('test_execute.py' in content)

    def testRunmode(self):
        '''Test the runmode decoration'''
        script = SoS_Script(r"""
from pysos import SoS_Action

@SoS_Action(run_mode='run')
def fail():
    return 1

[0]
a = fail()
""")
        wf = script.workflow()
        env.run_mode = 'dryrun'
        wf.run()
        # should return 0 in dryrun mode
        self.assertEqual(env.locals['a'], 0)
        #

        env.run_mode = 'run'
        wf.run()
        # shoulw return 1 in run mode
        self.assertEqual(env.locals['a'], 1)

    def testReadonlyVarsInGalobal(self):
        '''Test vars defined in global section are readonly'''
        script = SoS_Script(r"""

a = 10

[0]

a += 1

""")
        wf = script.workflow()
        env.run_mode = 'dryrun'
        self.assertRaises(RuntimeError, wf.run)

    def testReadonlyVarsInParameters(self):
        '''Test vars defined in global section are readonly'''
        script = SoS_Script(r"""

a = 10
import random

[parameters]

b = random.randint(0, 100000)

[a_1]

[default=a+a]

""")
        wf = script.workflow()
        env.run_mode = 'dryrun'
        self.assertRaises(RuntimeError, wf.run)

    def testExecuteInParallel(self):
        '''Test execution of step process in parallel'''
        script = SoS_Script(r"""

import time
import random

[1]
repeat=range(5)
input: for_each='repeat'

wait = random.randint(1,5)
time.sleep(wait)
print('I am {} after {} seconds'.format(_index, wait))

""")
        env.max_jobs = 5
        wf = script.workflow()
        wf.run()
        #
        # execute in parallel mode


if __name__ == '__main__':
    unittest.main()
