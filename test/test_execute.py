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

from sos.sos_script import SoS_Script
from sos._version import __version__
from sos.utils import env
from sos.sos_eval import Undetermined
from sos.sos_executor import Base_Executor, MP_Executor, ExecuteError
from sos.sos_script import ParsingError
from sos.target import FileTarget
import subprocess

class TestExecute(unittest.TestCase):
    def setUp(self):
        env.reset()
        self.resetDir('.sos')
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

    def resetDir(self, dirname):
        if os.path.isdir(dirname):
            shutil.rmtree(dirname)
        os.mkdir(dirname)

    def testCommandLine(self):
        '''Test command line arguments'''
        result = subprocess.check_output('sos --version', stderr=subprocess.STDOUT, shell=True).decode()
        self.assertTrue(result.startswith('sos {}'.format(__version__)))
        self.assertEqual(subprocess.call('sos', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos -h', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos run -h', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos-runner -h', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos dryrun -h', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos dryrun scripts/master', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        self.assertEqual(subprocess.call('sos dryrun scripts/master.sos', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        self.assertEqual(subprocess.call('sos-runner scripts/master', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        self.assertEqual(subprocess.call('sos dryrun file://{}/scripts/master.sos'.format(os.getcwd()), stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        self.assertEqual(subprocess.call('sos dryrun scripts/master.sos L', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos-runner scripts/master L', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos convert -h', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        #
        self.assertEqual(subprocess.call('sos config -h', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos config -g --get', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos config --get', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos config --set a 5', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos config --get a', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        self.assertEqual(subprocess.call('sos config --unset a', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        # a redirect bug related to blessing, not sure why the test fails
        #self.assertEqual(subprocess.call('sos run scripts/slave1.sos -v1 > /dev/null', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)


    def testInterpolation(self):
        '''Test string interpolation during execution'''
        self.touch(['a_1.txt', 'b_2.txt', 'c_2.txt'])
        script = SoS_Script(r"""
[0: shared='res']
res = ''
b = 200
res += "${b}"
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], '200')
        #
        script = SoS_Script(r"""
[0: shared='res']
res = ''
for b in range(5):
    res += "${b}"
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], '01234')
        #
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
output: ['{}_{}_processed.txt'.format(x,y) for x,y in zip(name, model)]

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'],  ['a_1_processed.txt', 'b_2_processed.txt', 'c_2_processed.txt'])
        #
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
output: ["${x}_${y}_process.txt" for x,y in zip(name, model)]

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'],  ['a_1_process.txt', 'b_2_process.txt', 'c_2_process.txt'])
        #
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
def add_a(x):
    return ['a'+_x for _x in x]

input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
output: add_a(["${x}_${y}_process.txt" for x,y in zip(name, model)])

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'],  ['aa_1_process.txt', 'ab_2_process.txt', 'ac_2_process.txt'])

    def testGlobalVars(self):
        '''Test SoS defined variables'''
        script = SoS_Script(r"""
[0]
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['SOS_VERSION'], __version__)

    def testFuncDef(self):
        '''Test defintion of function that can be used by other steps'''
        self.touch(['aa.txt', 'ab.txt'])
        script = SoS_Script(r"""
def myfunc(a):
    sum(range(5))
    return ['a' + x for x in a]

[0: shared={'test':'input'}]
input: myfunc(['a.txt', 'b.txt'])
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['test'], ['aa.txt', 'ab.txt'])
        # in nested workflow?
        script = SoS_Script(r"""
def myfunc(a):
    return ['a' + x for x in a]

[mse: shared={'test':'output'}]
input: myfunc(['a.txt', 'b.txt'])

[1]
sos_run('mse')
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        #
        # Names defined in subworkflow is not returned to the master dict
        self.assertTrue('test' not in env.sos_dict)

    def testInput(self):
        '''Test input specification'''
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
input: '*.py'
output: [x + '.res' for x in _input]
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertTrue('test_execute.py.res' in env.sos_dict['res'])

    def testForEach(self):
        '''Test for_each option of input'''
        self.touch(['a.txt', 'b.txt', 'a.pdf'])
        script = SoS_Script(r"""
[0: shared=['counter', 'all_names', 'all_loop']]
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
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['counter'], 6)
        self.assertEqual(env.sos_dict['all_names'], "a b c a b c ")
        self.assertEqual(env.sos_dict['all_loop'], "1 1 1 2 2 2 ")
        #
        # test same-level for loop and parameter with nested list
        script = SoS_Script(r"""
[0: shared=['processed']]
files = ['a.txt', 'b.txt']
par = [(1, 2), (1, 3), (2, 3)]
res = ['p1.txt', 'p2.txt', 'p3.txt']
processed = []

input: files, for_each='par,res'
output: res

processed.append((_par, _res))
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['processed'], [((1, 2), 'p1.txt'), ((1, 3), 'p2.txt'), ((2, 3), 'p3.txt')])
        #
        # test for each for pandas dataframe
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
import pandas as pd
data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

input: for_each='data'
output: "${_data['A']}_${_data['B']}_${_data['C']}.txt"
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], ['1_2_Hello.txt', '2_4_World.txt'])

    def testPairedWith(self):
        '''Test option paired_with '''
        pass

    def testInputPattern(self):
        '''Test option pattern of step input '''
        #env.verbosity = 4
        self.touch(['a-20.txt', 'b-10.txt'])
        script = SoS_Script(r"""
[0: shared=['base', 'name', 'par', '_output']]

files = ['a-20.txt', 'b-10.txt']
input: files, pattern=['{name}-{par}.txt', '{base}.txt']
output: ['{}-{}-{}.txt'.format(x,y,z) for x,y,z in zip(_base, _name, _par)]

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['base'], ["a-20", 'b-10'])
        self.assertEqual(env.sos_dict['name'], ["a", 'b'])
        self.assertEqual(env.sos_dict['par'], ["20", '10'])
        self.assertEqual(env.sos_dict['_output'], ["a-20-a-20.txt", 'b-10-b-10.txt'])

    def testOutputPattern(self):
        '''Test option pattern of step output'''
        #env.verbosity = 4
        self.touch(['a-20.txt', 'b-10.txt'])
        script = SoS_Script(r"""
[0: shared=['base', 'name', 'par', '_output']]

files = ['a-20.txt', 'b-10.txt']
input: files, pattern=['{name}-{par}.txt', '{base}.txt']
output: expand_pattern('{base}-{name}-{par}.txt'), expand_pattern('{par}.txt')

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['base'], ["a-20", 'b-10'])
        self.assertEqual(env.sos_dict['name'], ["a", 'b'])
        self.assertEqual(env.sos_dict['par'], ["20", '10'])
        self.assertEqual(env.sos_dict['_output'], ['a-20-a-20.txt', 'b-10-b-10.txt', '20.txt', '10.txt'])

    def testShared(self):
        '''Test option shared'''
        script = SoS_Script(r"""
parameter: res = 1

[0]
res = 2

[1]
res = 3
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], 1)
        #
        script = SoS_Script(r"""
parameter: res = 1

[0: shared='res']
res = 2

[1]
res = 3
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], 2)
        #
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[0: shared='a']
res = 2

[1: shared='res']
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], 3)
        self.assertEqual(env.sos_dict['a'], 30)
        # test multiple vars
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=('res', 'a')]
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], 3)
        self.assertEqual(env.sos_dict['a'], 5)
        #
        # test expression
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared={'res': 'res + 6', 'c': 'a'}]
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], 9)
        self.assertEqual(env.sos_dict['c'], 5)
        # test mixed vars and mapping
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=['res', {'c': 'a'}]]
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], 3)
        self.assertEqual(env.sos_dict['c'], 5)

#    def testSectionOptionWorkdir(self):
#        '''Test section option workdir'''
#        script = SoS_Script(r"""
#
#[1: workdir='tmp']
#run:
#    touch 'a.txt'
#""")
#        wf = script.workflow()
#        Base_Executor(wf).run()
#        self.assertTrue(os.path.isfile('tmp/a.txt'))
#        shutil.rmtree('tmp')

    def testFileType(self):
        '''Test input option filetype'''
        self.touch(['a.txt', 'b.txt', 'a.pdf', 'b.html'])
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', files, filetype='*.txt', group_by='single'

output: "${_input}.res"

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['res'], ['a.txt.res', 'b.txt.res'])
        #
        script = SoS_Script(r"""
[0: shared='counter']
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', 'b.html', files, filetype=('*.txt', '*.pdf'), group_by='single'

counter += 1
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['counter'], 3)
        #
        script = SoS_Script(r"""
[0: shared='counter']
files = ['a.txt', 'b.txt']
counter = 0

input: 'a.pdf', 'b.html', files, filetype=lambda x: 'a' in x, group_by='single'

counter += 1
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['counter'], 2)

    def testOutputFromInput(self):
        '''Test deriving output files from input files'''
        self.touch(['a.txt', 'b.txt'])
        script = SoS_Script(r"""
[0: shared={'counter':'counter', 'step':'output'}]
files = ['a.txt', 'b.txt']
counter = 0

input: files, group_by='single'
output: _input[0] + '.bak'

counter += 1
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['counter'], 2)
        self.assertEqual(env.sos_dict['step'], ['a.txt.bak', 'b.txt.bak'])

    def testWorkdir(self):
        '''Test workdir option for runtime environment'''
        script =  SoS_Script(r"""
[0]

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
        '''Test workdir option for runtime environment'''
        env.max_jobs = 5
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
        MP_Executor(wf).run()
        self.assertGreater(time.time() - start, 9)
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
        MP_Executor(wf).run()
        self.assertLess(time.time() - start, 6)

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
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        # use option env
        script = SoS_Script(r"""
[1]
task: env={'PATH': 'temp' + os.pathsep + os.environ['PATH']}
sh:
    temp_cmd
""")
        wf = script.workflow()
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


    def testRunmode(self):
        '''Test the runmode decoration'''
        script = SoS_Script(r"""
from sos.actions import SoS_Action

@SoS_Action(run_mode='run')
def fail():
    return 1

[0: shared='a']
a = fail()
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        # should return 0 in dryrun mode
        self.assertTrue(isinstance(env.sos_dict['a'], Undetermined))
        #
        Base_Executor(wf).run()
        # shoulw return 1 in run mode
        self.assertEqual(env.sos_dict['a'], 1)

    def testReadonlyVarsInGalobal(self):
        '''Test vars defined in global section are readonly'''
        script = SoS_Script(r"""

a = 10

[0]

a += 1

""")
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testPassingVarsToNestedWorkflow(self):
        '''Test if variables can be passed to nested workflows'''
        script = SoS_Script(r"""
%set_options sigil='[ ]'
import time
import random

[nested]
print('I am nested [nested] with seed [seed]')

[0]
reps = range(5)
input: for_each='reps'
task: concurrent=True
nested = _reps
seed = random.randint(1, 1000)
print('Passing [seed] to [nested]')
sos_run('nested')

""")
        wf = script.workflow()
        Base_Executor(wf).run()

    def testUserDefinedFunc(self):
        '''Test the use of user-defined functions in SoS script'''
        script = SoS_Script(r"""

def myfunc():
  return 'a'

[1: shared={'test':'output'}]
output: myfunc()

myfunc()

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['test'], ['a'])
        # User defined function should also work under nested workflows
        # This is difficult because the 'local namespace' is usually
        # not seen inside function definition. The solution now is to
        # use a single workspace.
        script = SoS_Script(r"""

def myfunc():
    # test if builtin functions (sum and range) can be used here.
    return 'a' + str(sum(range(10)))

[1: shared={'test':'output'}]
output: [myfunc() for i in range(10)][0]

myfunc()

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['test'], ['a45'])

    def testReadOnlyStepVars(self):
        '''Test if the step variables can be changed.'''
        #
        script = SoS_Script(r"""
[1: shared={'test':'output'}]
output: 'a.txt'

[2]
test.output=['ab.txt']

""")
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testReadOnlyInputOutputVars(self):
        '''Test readonly input output vars'''
        script = SoS_Script(r"""
[1: shared='test']
output: 'a.txt'
_output = ['b.txt']

""")
        wf = script.workflow()
        env.run_mode = 'dryrun'
        # I would like to disallow setting _output directly, but this is
        # not the case now.
        self.assertRaises(ExecuteError, Base_Executor(wf).dryrun)

    def testLocalNamespace(self):
        '''Test if steps are well separated.'''
        self.touch('a.txt')
        script = SoS_Script(r"""
[1]
a = 1

[2]
# this should fail because a is defined in another step
print(a)

""")
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        # however, alias should be sent back
        script = SoS_Script(r"""
[1: shared={'shared': 'output'}]
input: 'a.txt'
output: 'b.txt'

[2: shared={'tt':'output'}]
print(shared)

output: [x + '.res' for x in shared]

""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['shared'], ['b.txt'])
        self.assertEqual(env.sos_dict['tt'], ['b.txt.res'])
        #
        # this include other variables set in the step
        script = SoS_Script(r"""
[1: shared={'shared':'c', 'd':'d'}]
input: 'a.txt'
output: 'b.txt'

c = 'c.txt'
d = 1

[2: shared={'d': 'e'}]
# this should fail because a is defined in another step
print(shared)

output: shared

e = d + 1

""")
        wf = script.workflow()
        # I would like to disallow accessing variables defined
        # in other cases.
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['shared'], 'c.txt')
        self.assertEqual(env.sos_dict['d'], 2)

    def testSklearnImportFailure(self):
        '''Test problem with Sklean when using Celery/multiprocessing'''
        script = SoS_Script('''
import sklearn

[run]
print(0)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

#
#    def testCollectionOfErrors(self):
#        '''Test collection of errors when running in dryrun mode.'''
#        script = SoS_Script('''
#[0]
#depends: executable('a1')
#[1: skip='blah']
#depends: executable('a2')
#
#[2]
#input: None
#depends: executable('a3')
#[3]
#depends: executable('a4')
#
#''')
#        wf = script.workflow()
#        # we should see a single error with 2 messages.
#        # because 2 being on a separate branch will be executed but
#        # the later steps will not be executed
#        try:
#            Base_Executor(wf).dryrun()
#        except Exception as e:
#            self.assertEqual(len(e.errors), 3)

    def testSearchPath(self):
        '''Test if any action should exit in five seconds in dryrun mode'''
        sos_config_file = 'config.yml'
        #
        with open(sos_config_file, 'w') as sos_config:
            sos_config.write('''
#
# global sos configuration file
#
{{
    "sos_path": ["{0}/crazy_path", "{0}/crazy_path/more_crazy/"]
}}
'''.format(os.getcwd()))
        #
        if not os.path.isdir('crazy_path'):
            os.mkdir('crazy_path')
            os.mkdir('crazy_path/more_crazy')
        with open('crazy_path/crazy_master.sos', 'w') as crazy:
            crazy.write('''
[0]
sos_run('cc', source='crazy_slave.sos')

''')
        with open('crazy_path/more_crazy/crazy_slave.sos', 'w') as crazy:
            crazy.write('''
[cc_0]
print('hay, I am crazy')
''')

        script = SoS_Script(filename='crazy_master.sos')
        script.workflow()
        #
        shutil.rmtree('crazy_path')
        os.remove(sos_config_file)

    def testDynamicOutput(self):
        '''Testing dynamic output'''
        #
        if not os.path.isdir('temp'):
            os.mkdir('temp')
        #
        script = SoS_Script('''
[10: shared={'test':'output'}]
ofiles = []
output: dynamic(ofiles)

for i in range(4):
    ff = 'temp/something{}.html'.format(i)
    ofiles.append(ff)
    with open(ff, 'w') as h:
       h.write('a')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['test'], ['temp/something{}.html'.format(x) for x in range(4)])
        #
        shutil.rmtree('temp')

    def testDynamicInput(self):
        '''Testing dynamic input'''
        #
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        #
        script = SoS_Script('''
%set_options sigil='< >'
[1]

for i in range(5):
    run("touch temp/test_<i>.txt")


[10: shared={'test':'output'}]
input: dynamic('temp/*.txt'), group_by='single'
output: dynamic('temp/*.txt.bak')

run:
touch <_input>.bak
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['test'], ['temp/test_{}.txt.bak'.format(x) for x in range(5)])
        # this time we use th existing signature
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['test'], ['temp/test_{}.txt.bak'.format(x) for x in range(5)])
        #
        shutil.rmtree('temp')


    def testAssignmentAfterInput(self):
        '''Testing assignment after input should be usable inside step process.'''
        #
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        #
        env.sig_mode = 'ignore'
        script = SoS_Script('''
%set_options sigil='%( )'
[1]
rep = range(5)
input:  for_each='rep'
output: "temp/%(_rep).txt"

# ff should change and be usable inside run
ff = "%(_rep).txt"
run:
echo %(ff)
touch temp/%(ff)
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        #
        shutil.rmtree('temp')

    def testUseOfRunmode(self):
        '''Test the use of run_mode variable in SoS script'''
        #
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        env.sig_mode = 'ignore'
        script = SoS_Script('''
[1: shared={'res':'output'}]
import random
if run_mode == 'run':
    for i in range(3):
        with open("temp/test_${random.randint(1, 100000)}.txt", 'w') as res:
            res.write(str(i))

''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # we should have 9 files
        files = glob.glob('temp/*.txt')
        self.assertEqual(len(files), 3)
        shutil.rmtree('temp')

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
            Base_Executor(wf).run()
            files = list(glob.glob('temp/*.txt'))
            self.assertEqual(files, result)
            #
            # test last iteration
            shutil.rmtree('temp')

    def testActionBeforeInput(self):
        '''Testing the execution of actions before input directive
        (variables such as _index should be made available). '''
        script = SoS_Script('''
[0]
bash('echo "A"')
input: 
''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()

    def testDuplicateIOFiles(self):
        '''Test interpretation of duplicate input/output/depends'''
        self.resetDir('temp')
        # Test duplicate input
        os.system('touch temp/1.txt')
        script = SoS_Script('''
[1]
input: ['temp/1.txt' for x in range(5)]
python:
with open('temp/{}.input'.format(len([${input!r,}])), 'w') as f: f.write('')
        ''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('temp/5.input'))
        # Test duplicate output
        script = SoS_Script('''
[1]
output: ['temp/2.txt' for x in range(5)]
python:
with open('temp/2.txt', 'w') as f: f.write('')
with open('temp/{}.output'.format(len([${output!r,}])), 'w') as f: f.write('')
        ''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('temp/5.output'))
        # Test duplicate depends
        script = SoS_Script('''
[1]
input: 'temp/1.txt'
depends: ['temp/2.txt' for x in range(5)]
output: 'temp/3.txt'
python:
with open('temp/3.txt', 'w') as f: f.write('')
with open('temp/{}.depends'.format(len([${depends!r,}])), 'w') as f: f.write('')
        ''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('temp/5.depends'))
        shutil.rmtree('temp')

    def testOutputInLoop(self):
        '''Test behavior of ${output} when used in loop'''
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        env.sig_mode = 'ignore'
        script = SoS_Script('''
[default]
s = [x for x in range(5)]
output_files = ['temp/{}.txt'.format(x) for x in range(5)]
input: for_each = ['s']
output: output_files[_index]
run: active = 0
rm -f temp/out.log
run:
echo ${output} >> temp/out.log
touch ${output}
        ''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # output should have 1, 2, 3, 4, 5, respectively, and
        # the total record files would be 1+2+3+4+5=15
        with open('temp/out.log') as out:
            self.assertEqual(len(out.read().split()), 15)
        shutil.rmtree('temp')
        #
        os.mkdir('temp')
        script = SoS_Script('''
[default]
s = [x for x in range(5)]
output_files = ['temp/{}.txt'.format(x) for x in range(5)]
input: for_each = ['s']
output: output_files[_index]
run: active = 0
rm -f temp/out.log
task:
run:
echo ${output} >> temp/out.log
touch ${output}
        ''')
        wf = script.workflow()
        env.sig_mode = 'ignore'
        Base_Executor(wf).run()
        with open('temp/out.log') as out:
            self.assertEqual(len(out.read().split()), 15)
        shutil.rmtree('temp')

    def testSignature(self):
        self._testSignature(r"""
import time
[*_0]
output: 'temp/a.txt', 'temp/b.txt'
task:
if run_mode == 'run':
   time.sleep(1)
   run('''echo "a.txt" > 'temp/a.txt' ''')
   run('''echo "b.txt" > 'temp/b.txt' ''')

[1: shared={'oa':'output'}]
dest = ['temp/c.txt', 'temp/d.txt']
input: group_by='single', paired_with='dest'
output: _dest

task:
if run_mode == 'run':
    time.sleep(0.5)
    run(" cp ${_input} ${_dest} ")
""")
        #
        env.max_jobs = 4
        self._testSignature(r"""
import time
[*_0]
output: 'temp/a.txt', 'temp/b.txt'

task:
if run_mode == 'run':
    time.sleep(1)
    run('''echo "a.txt" > 'temp/a.txt' ''')
    run('''echo "b.txt" > 'temp/b.txt' ''')

[1: shared={'oa':'output'}]
dest = ['temp/c.txt', 'temp/d.txt']
input: group_by='single', paired_with='dest'
output: _dest

task:
if run_mode == 'run':
   time.sleep(0.5)
   run(" cp ${_input} ${_dest} ")
""")
        # script format
        env.max_jobs = 4
        self._testSignature(r"""
import time
[*_0]
output: 'temp/a.txt', 'temp/b.txt'

run:
sleep 1
echo "a.txt" > 'temp/a.txt'

run:

echo "b.txt" > 'temp/b.txt'

[1: shared={'oa':'output'}]
dest = ['temp/c.txt', 'temp/d.txt']
input: group_by='single', paired_with='dest'
output: _dest

task:
if run_mode == 'run':
    time.sleep(0.5)
run:
cp ${_input} ${_dest}
""")
        # reset env mode
        env.sig_mode = 'default'
        shutil.rmtree('temp')

    def testSignatureWithSharedVariable(self):
        '''Test restoration of signature from variables.'''
        FileTarget('a.txt').remove('both')
        # shared 
        script = SoS_Script(r"""
import time
[0: shared='a']
output: 'a.txt'
run:
   sleep 3
   touch a.txt

a= 5

[1]
print(a)

""")
        # alias should also be recovered.
        wf = script.workflow('default')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 3)
        # rerun
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 1)
        FileTarget('a.txt').remove('both')

    def testSignatureWithoutOutput(self):
        # signature without output file
        self._testSignature(r"""
import time
[*_0]
output: []

run:
sleep 1
[ -d temp ] || mkdir temp
echo "a.txt" > 'temp/a.txt'

run:

echo "b.txt" > 'temp/b.txt'

[1: shared={'oa':'output'}]
dest = ['temp/c.txt', 'temp/d.txt']
input: 'temp/a.txt', 'temp/b.txt', group_by='single', paired_with='dest'
output: _dest

run:
sleep 0.5
cp ${_input} ${_dest}
""")
        # reset env mode
        env.sig_mode = 'default'
        shutil.rmtree('temp')

    def _testSignature(self, text):
        '''Test recognizing the format of SoS script'''
        script = SoS_Script(text)
        for f in ['temp/a.txt', 'temp/b.txt']:
            FileTarget(f).remove('both')
        #
        # only the first step
        wf = script.workflow('default:0')
        start = time.time()
        env.sig_mode = 'default'
        if env.max_jobs == 1:
            Base_Executor(wf).run()
        else:
            MP_Executor(wf).run()
        self.assertGreater(time.time() - start, 1)
        self.assertTrue(os.path.isfile('temp/a.txt'))
        self.assertTrue(os.path.isfile('temp/b.txt'))
        with open('temp/a.txt') as ta:
            self.assertTrue(ta.read(), 'a.txt')
        with open('temp/b.txt') as tb:
            self.assertTrue(tb.read(), 'b.txt')
        env.sig_mode = 'assert'
        if env.max_jobs == 1:
            Base_Executor(wf).run()
        else:
            MP_Executor(wf).run()
        #
        wf = script.workflow()
        start = time.time()
        env.sig_mode = 'build'
        if env.max_jobs == 1:
            Base_Executor(wf).run()
        else:
            MP_Executor(wf).run()

        self.assertLess(time.time() - start, 1.5)
        #
        self.assertTrue(os.path.isfile('temp/c.txt'))
        self.assertTrue(os.path.isfile('temp/d.txt'))
        with open('temp/c.txt') as tc:
            self.assertTrue(tc.read(), 'a.txt')
        with open('temp/d.txt') as td:
            self.assertTrue(td.read(), 'b.txt')
        self.assertEqual(env.sos_dict['oa'], ['temp/c.txt', 'temp/d.txt'])
        #
        # now in assert mode, the signature should be there
        env.sig_mode = 'assert'
        if env.max_jobs == 1:
            Base_Executor(wf).run()
        else:
            MP_Executor(wf).run()

        #
        start = time.time()
        env.sig_mode = 'default'
        if env.max_jobs == 1:
            Base_Executor(wf).run()
        else:
            MP_Executor(wf).run()
        
        self.assertLess(time.time() - start, 1.5)
        #
        # change script a little bit
        script = SoS_Script('# comment\n' + text)
        wf = script.workflow()
        env.sig_mode = 'assert'
        if env.max_jobs == 1:
            Base_Executor(wf).run()
        else:
            MP_Executor(wf).run()

        # add some other variable?
        #script = SoS_Script('comment = 1\n' + text)
        #wf = script.workflow()
        #env.sig_mode = 'assert'
        #self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testReexecution(self):
        '''Test -f option of sos run'''
        script = SoS_Script('''
import time

[0]
output: 'a.txt'
task:
time.sleep(3)
run("touch ${output}")
''')
        wf = script.workflow()
        try:
            # remove existing output if exists
            FileTarget('a.txt').remove('both')
        except:
            pass
        start = time.time()
        Base_Executor(wf).run()
        # regularly take more than 5 seconds to execute
        self.assertGreater(time.time() - start, 2)
        # now, rerun should be much faster
        start = time.time()
        Base_Executor(wf).run()
        # rerun takes less than 1 second
        self.assertLess(time.time() - start, 1)
        #
        # force rerun mode
        start = time.time()
        env.sig_mode = 'ignore'
        Base_Executor(wf).run()
        # regularly take more than 5 seconds to execute
        self.assertGreater(time.time() - start, 2)
        try:
            # remove existing output if exists
            os.remove('a.txt')
        except:
            pass

    def testDependsExecutable(self):
        '''Testing target executable.'''
        script = SoS_Script('''
[0]
depends: executable('ls')
sh:
    touch a.txt
''')
        wf = script.workflow()
        FileTarget('a.txt').remove('both')
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('a.txt'))
        FileTarget('a.txt').remove('both')
        
    def testOutputExecutable(self):
        '''Testing target executable.'''
        # change $PATH so that lls can be found at the current
        # directory.
        os.environ['PATH'] += os.pathsep + '.'
        script = SoS_Script('''
[0]
output: executable('lls')
sh:
    touch lls
    sleep 3
    chmod +x lls
''')
        wf = script.workflow()
        FileTarget('lls').remove('both')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2)
        # test validation
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 2)
        FileTarget('lls').remove('both')

    def testDependsRLibrary(self):
        '''Testing depending on R_library'''
        # first remove xtable package
        if not shutil.which('R'):
            return 
        subprocess.call('R CMD REMOVE xtable', shell=True)
        script = SoS_Script('''
[0]

depends: R_library('xtable')
R:
    library('xtable')
    ## Demonstrate data.frame
    tli.table <- xtable(cars)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    def testDependsEnvVariable(self):
        '''Testing target env_variable.'''
        FileTarget('a.txt').remove('both')
        script = SoS_Script('''
[0]
depends: env_variable('AA')
output:  'a.txt'
sh:
    sleep 2
    echo $AA > a.txt
''')
        wf = script.workflow()
        os.environ['AA'] = 'A1'
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 1.5)
        with open('a.txt') as at:
            self.assertEqual(at.read(), 'A1\n')
        # test validation
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 1)
        # now if we change var, it should be rerun
        os.environ['AA'] = 'A2'
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 1.5)
        with open('a.txt') as at:
            self.assertEqual(at.read(), 'A2\n')
        FileTarget('a.txt').remove('both')

    def testProvidesExecutable(self):
        '''Testing provides executable target.'''
        # change $PATH so that lls can be found at the current
        # directory.
        os.environ['PATH'] += os.pathsep + '.'
        FileTarget('lls').remove('both')
        script = SoS_Script('''
[lls: provides=executable('lkls')]
sh:
    touch lkls
    sleep 3
    chmod +x lkls

[c]
depends: executable('lkls')

''')
        wf = script.workflow('c')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2)
        FileTarget('lkls').remove('both')


    def testSignatureAfterRemovalOfFiles(self):
        '''test action shrink'''
        if os.path.isfile('largefile.txt'):
            os.remove('largefile.txt')
        script = SoS_Script(r'''
[10]

# generate a file
output: 'largefile.txt'

python:
    import time
    time.sleep(3)
    with open("${output}", 'w') as out:
        for i in range(1000):
            out.write('{}\n'.format(i))

''')
        wf = script.workflow()
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 3)
        # rerun, because this is the final target, it has to be
        # re-generated
        os.remove('largefile.txt')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 3)
        # 
        self.assertTrue(os.path.isfile('largefile.txt'))
        # we discard just the signature, the step will be ignored
        # as long as the file is not touched.
        st = time.time()
        FileTarget('largefile.txt').remove('signature')
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 0.5)
        #
        # now if we touch the file, it needs to be regenerated
        st = time.time()
        with open('largefile.txt', 'a') as lf:
            lf.write('something')
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 3)
        FileTarget('largefile.txt').remove('both')

    def testSignatureWithParameter(self):
        '''Test signature'''
        FileTarget('myfile.txt').remove('both')
        #
        script = SoS_Script(r'''
parameter: gvar = 10

[10]
# generate a file
output: 'myfile.txt'
# additional comment
python:
    import time
    time.sleep(3)
    with open(${output!r}, 'w') as tmp:
        tmp.write('${gvar}')

''')
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2.5)
        with open('myfile.txt') as tmp:
            self.assertEqual(tmp.read(), '10')
        #
        # now if we change parameter, the step should be rerun
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf, args=['--gvar', '20']).run()
        with open('myfile.txt') as tmp:
            self.assertEqual(tmp.read(), '20')
        self.assertGreater(time.time() - st, 2.5)
        #
        # do it again, signature should be effective
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf, args=['--gvar', '20']).run()
        with open('myfile.txt') as tmp:
            self.assertEqual(tmp.read(), '20')
        self.assertLess(time.time() - st, 2.5)

        #
        script = SoS_Script(r'''
parameter: gvar = 10

[10]
# generate a file
output: 'myfile.txt'
# additional comment
task:
python:
    import time
    time.sleep(3)
    with open(${output!r}, 'w') as tmp:
        tmp.write('${gvar}')

''')
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2.5)
        with open('myfile.txt') as tmp:
            self.assertEqual(tmp.read(), '10')
        #
        # now if we change parameter, the step should be rerun
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf, args=['--gvar', '20']).run()
        with open('myfile.txt') as tmp:
            self.assertEqual(tmp.read(), '20')
        self.assertGreater(time.time() - st, 2.5)
        #
        # do it again, signature should be effective
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf, args=['--gvar', '20']).run()
        with open('myfile.txt') as tmp:
            self.assertEqual(tmp.read(), '20')
        self.assertLess(time.time() - st, 2.5)
        FileTarget('myfile.txt').remove('both')

    def testPassingVarToTask(self):
        '''Test passing used variable to tasks'''
        for i in range(10, 13):
            FileTarget('myfile_{}.txt'.format(i)).remove('both')
        #
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
        MP_Executor(wf).run()
        for t in range(10, 13):
            with open('myfile_{}.txt'.format(t)) as tmp:
                self.assertEqual(tmp.read(), str(t) + '_' + str(t-10))
            FileTarget('myfile_{}.txt'.format(t)).remove('both')


    def testLoopWiseSignature(self):
        '''Test partial signature'''
        for i in range(10, 12):
            FileTarget('myfile_{}.txt'.format(i)).remove('both')
        #
        script = SoS_Script(r'''
parameter: gvar = 10

[10]
tt = [gvar]
input: for_each='tt'
output: "myfile_${_tt}.txt"
python:
    import time
    time.sleep(3)
    print("DO ${_tt}")
    with open(${_output!r}, 'w') as tmp:
        tmp.write('${_tt}')
''')
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2.5)
        # now we modify the script 
        script = SoS_Script(r'''
parameter: gvar = 10

[10]
tt = [gvar, gvar + 1]
input: for_each='tt'
output: "myfile_${_tt}.txt"
python:
    import time
    time.sleep(3)
    print("DO ${_tt}")
    with open(${_output!r}, 'w') as tmp:
        tmp.write('${_tt}')
''')
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2.5)
        self.assertLess(time.time() - st, 5)
        #
        # run it again, neither needs to be rerun
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 2)
        #
        # change again, the second one is already there.
        script = SoS_Script(r'''
parameter: gvar = 10

[10]
tt = [gvar + 1]
input: for_each='tt'
output: "myfile_${_tt}.txt"
python:
    import time
    time.sleep(3)
    print("DO ${_tt}")
    with open(${_output!r}, 'w') as tmp:
        tmp.write('${_tt}')
''')
        st = time.time()
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 2)
        #
        for t in range(10, 12):
            with open('myfile_{}.txt'.format(t)) as tmp:
                self.assertEqual(tmp.read(), str(t))
            FileTarget('myfile_{}.txt'.format(t)).remove('both')


    def testExecutionLock(self):
        '''Test execution lock of two processes'''
        with open('lock.sos', 'w') as lock:
            lock.write(r'''
import time
[A_1]
output: 'a.txt'
time.sleep(3)
with open('a.txt', 'w') as txt:
    txt.write('A1\n')

# A1 and A2 are independent
[A_2]
input: None
output: 'b.txt'
time.sleep(3)
with open('b.txt', 'w') as txt:
    txt.write('A2\n')
        ''')
        st = time.time()
        ret1 = subprocess.Popen('sos run lock -j1', shell=True)
        ret2 = subprocess.Popen('sos run lock -j1', shell=True)
        ret1.wait()
        ret2.wait()
        # two processes execute A_1 and A_2 separately, usually
        # takes less than 5 seconds
        self.assertLess(time.time() - st, 7)
        FileTarget('lock.sos').remove('both')


    def testOutputFromSignature(self):
        'Test restoration of output from signature'''
        self.touch(['1.txt', '2.txt'])
        script = SoS_Script('''
parameter: K = [2,3]

[work_1]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run:
  touch ${_output}

[work_2]

input: group_by = 'single', pattern = '{name}.{ext}', paired_with = ['K']
output: expand_pattern('{_name}.{_K}.out')
run: 
  touch ${_output}
    ''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # for the second run, output should be correctly constructed
        Base_Executor(wf).run()
        for file in ['1.out', '2.out', '1.2.out', '2.3.out']:
            FileTarget(file).remove('both')


    def testRemovedIntermediateFiles(self):
        '''Test behavior of workflow with removed internediate files'''
        FileTarget('a.txt').remove('both')
        FileTarget('aa.txt').remove('both')
        script = SoS_Script('''
[10]
output: 'a.txt'
sh:
    sleep 2
    echo "a" > a.txt

[20]
output: 'aa.txt'
sh:
    sleep 2
    cat ${input} > ${output}
''')
        wf = script.workflow()
        st = time.time()
        Base_Executor(wf).run()
        self.assertTrue(FileTarget('aa.txt').exists())
        self.assertGreater(time.time() - st, 3.5)
        # rerun should be faster
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 1)
        # if we remove the middle result, it should not matter
        os.remove('a.txt')
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, 1)
        #
        # if we remove the final result, it will be rebuilt
        os.remove('aa.txt')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 1.5)
        #
        # now we request the generation of target
        FileTarget('a.txt').remove('target')
        FileTarget('aa.txt').remove('both')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 3.5)
        #
        FileTarget('a.txt').remove('both')
        FileTarget('aa.txt').remove('both')

    def testSharedVarInPairedWith(self):
        self.touch(['1.txt', '2.txt'])
        script = SoS_Script('''
[work_1: shared = {'data': 'output'}]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run:
  touch ${_output}

[work_2]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}', paired_with = ['data']
output: expand_pattern('{_name}.out2')
run:
  touch ${_data} ${_output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for file in ('1.out', '2.out', '1.out2', '2.out2'):
            FileTarget(file).remove('both')

    def testSharedVarInForEach(self):
        self.touch(['1.txt', '2.txt'])
        script = SoS_Script('''
[work_1: shared = {'data': 'output'}]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run:
  touch ${_output}

[work_2]
input: "1.txt", "2.txt", group_by = 'single', for_each = 'data, data',  pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out2')
run:
  touch ${_data} ${_output}

''')
        wf = script.workflow()
        Base_Executor(wf).run()

    def testRemovedDepends(self):
        '''Test a case where a dependent file has signature, but
        gets removed before the next run.'''
        script = SoS_Script('''
[tet: provides='a.txt']
run:
    echo "something" > a.txt

[20]
depends: 'a.txt'
output: 'b.txt'
run:
    cat b.txt > b.txt
''')
        wf = script.workflow()
        # this should be ok.
        Base_Executor(wf).run()
        # now let us remove a.txt (but the signature is still there)
        os.remove('a.txt')
        os.remove('b.txt')
        Base_Executor(wf).run()

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

if __name__ == '__main__':
    unittest.main()
