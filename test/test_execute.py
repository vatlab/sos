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
import glob
import unittest
import shutil

from sos.sos_script import SoS_Script
from sos._version import __version__
from sos.utils import env
from sos.sos_eval import Undetermined
from sos.sos_executor import Base_Executor, ExecuteError
from sos.target import FileTarget
import subprocess

def multi_attempts(fn):
    def wrapper(*args, **kwargs):
        for n in range(4):
            try:
                fn(*args, **kwargs)
                break
            except:
                if n > 1:
                    raise
    return wrapper

class TestExecute(unittest.TestCase):
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

    def resetDir(self, dirname):
        if os.path.isdir(os.path.expanduser(dirname)):
            shutil.rmtree(os.path.expanduser(dirname))
        os.mkdir(os.path.expanduser(dirname))

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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['res'],  ['a_1_processed.txt', 'b_2_processed.txt', 'c_2_processed.txt'])
        #
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
output: ["${x}_${y}_process.txt" for x,y in zip(name, model)]

""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['res'],  ['aa_1_process.txt', 'ab_2_process.txt', 'ac_2_process.txt'])

    def testGlobalVars(self):
        '''Test SoS defined variables'''
        script = SoS_Script(r"""
[0]
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['test'], ['aa.txt', 'ab.txt'])

    def testInput(self):
        '''Test input specification'''
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
input: '*.py'
output: [x + '.res' for x in _input]
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['res'], ['1_2_Hello.txt', '2_4_World.txt'])

        # test dictionary format of for_each
        self.touch(['a.txt', 'b.txt', 'a.pdf'])
        script = SoS_Script(r"""
[0: shared=['counter', 'all_names', 'all_loop']]
files = ['a.txt', 'b.txt']
names = ['a', 'b', 'c']
counter = 0
all_names = ''
all_loop = ''

input: 'a.pdf', files, group_by='single', paired_with='names', for_each={'c':  ['1', '2']}

all_names += _names[0] + " "
all_loop += c + " "

counter = counter + 1
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['counter'], 6)
        self.assertEqual(env.sos_dict['all_names'], "a b c a b c ")
        self.assertEqual(env.sos_dict['all_loop'], "1 1 1 2 2 2 ")
        #
        # test multi-key dictionary format of for_each
        self.touch(['a.txt'])
        script = SoS_Script(r"""
import itertools
[0: shared=['counter', 'all_names', 'all_loop']]
parameter: n = [300, 100]
parameter: p = [50, 200, 100]
parameter: outfile = ['1', '2', '3', '4', '5', '6']
counter = 0
all_names = ''
all_loop = ''
input: 'a.txt', group_by='single', for_each={'_n,_p': [(_n,_p) for _n,_p in itertools.product(n,p) if _n > _p]}

all_names += outfile[_index] + " "
all_loop += '{} {} '.format(_n, _p)
counter = counter + 1
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['counter'], 4)
        self.assertEqual(env.sos_dict['all_names'], "1 2 3 4 ")
        self.assertEqual(env.sos_dict['all_loop'], "300 50 300 200 300 100 100 50 ")
        #
        # test same-level for loop and parameter with nested list
        script = SoS_Script(r"""
[0: shared=['processed']]
files = ['a.txt', 'b.txt']
processed = []

input: files, for_each={'par':[(1, 2), (1, 3), (2, 3)], 'res': ['p1.txt', 'p2.txt', 'p3.txt']}
output: res

processed.append((par, res))
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['processed'], [((1, 2), 'p1.txt'), ((1, 3), 'p2.txt'), ((2, 3), 'p3.txt')])
        #
        # test for each for pandas dataframe
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
import pandas as pd
input: for_each={'data': pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])}
output: "${data['A']}_${data['B']}_${data['C']}.txt"
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['res'], ['1_2_Hello.txt', '2_4_World.txt'])
        # 
        # support for pands Series and Index types
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
import pandas as pd
data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])
input: for_each={'A': data['A']}
output: "a_${A}.txt"
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['res'], ['a_1.txt', 'a_2.txt'])
        #
        script = SoS_Script(r"""
[0: shared={'res':'output'}]
import pandas as pd
data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])
data.set_index('C', inplace=True)
input: for_each={'A': data.index}
output: "${A}.txt"
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['res'], ['Hello.txt', 'World.txt'])



    def testPairedWith(self):
        '''Test option paired_with '''
        self.touch(['a.txt', 'b.txt'])
        for ofile in ['a.txt1', 'b.txt2']:
            FileTarget(ofile).remove('both')
        #
        # string input
        script = SoS_Script(r'''
[0]
files = ['a.txt', 'b.txt']
vars = [1, 2]

input: files, paired_with='vars', group_by=1
output: "${_input}${_vars}"
run:
    touch ${output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for ofile in ['a.txt1', 'b.txt2']:
            self.assertTrue(FileTarget(ofile).exists('target'))
            FileTarget(ofile).remove('both')
        #
        # list input
        script = SoS_Script(r'''
[0]
files = ['a.txt', 'b.txt']
vars = [1, 2]
vars2 = ['a', 'b']

input: files, paired_with=('vars', 'vars2'), group_by=1
output: "${_input}${_vars}"
run:
    touch ${output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for ofile in ['a.txt1', 'b.txt2']:
            self.assertTrue(FileTarget(ofile).exists('target'))
            FileTarget(ofile).remove('both')
        #
        # dict input
        script = SoS_Script(r'''
[0]
files = ['a.txt', 'b.txt']
input: files, paired_with={'var': [1,2], 'var2': ['a', 'b']}, group_by=1
output: "${_input}${var}"
run:
    touch ${output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for ofile in ['a.txt1', 'b.txt2']:
            self.assertTrue(FileTarget(ofile).exists('target'))
            FileTarget(ofile).remove('both')



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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['base'], ["a-20", 'b-10'])
        self.assertEqual(env.sos_dict['name'], ["a", 'b'])
        self.assertEqual(env.sos_dict['par'], ["20", '10'])
        self.assertEqual(env.sos_dict['_output'], ['a-20-a-20.txt', 'b-10-b-10.txt', '20.txt', '10.txt'])


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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['counter'], 2)
        self.assertEqual(env.sos_dict['step'], ['a.txt.bak', 'b.txt.bak'])


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
        Base_Executor(wf).run(mode='dryrun')
        # should return 0 in dryrun mode
        self.assertTrue(isinstance(env.sos_dict['a'], Undetermined))
        #
        Base_Executor(wf).run()
        # shoulw return 1 in run mode
        self.assertEqual(env.sos_dict['a'], 1)

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
        env.config['run_mode'] = 'dryrun'
        # I would like to disallow setting _output directly, but this is
        # not the case now.
        self.assertRaises(Exception, Base_Executor(wf).run, mode="dryrun")

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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['shared'], 'c.txt')
        self.assertEqual(env.sos_dict['d'], 2)

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
#            Base_Executor(wf).run(mode='dryrun')
#        except Exception as e:
#            self.assertEqual(len(e.errors), 3)



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
        self.assertEqual(env.sos_dict['test'], [os.path.join('temp', 'test_{}.txt.bak'.format(x)) for x in range(5)])
        # this time we use th existing signature
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['test'], [os.path.join('temp', 'test_{}.txt.bak'.format(x)) for x in range(5)])
        #
        shutil.rmtree('temp')


    def testAssignmentAfterInput(self):
        '''Testing assignment after input should be usable inside step process.'''
        #
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        #
        env.config['sig_mode'] = 'ignore'
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
        env.config['sig_mode'] = 'ignore'
        script = SoS_Script('''
[1: shared={'res':'output'}]
import random
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

    def testActionBeforeInput(self):
        '''Testing the execution of actions before input directive
        (variables such as _index should be made available). '''
        script = SoS_Script('''
[0]
bash('echo "A"')
input: 
''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')

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
        self.assertTrue(os.path.isfile(os.path.join('temp', '5.depends')))
        shutil.rmtree('temp')

    def testOutputInLoop(self):
        '''Test behavior of ${output} when used in loop'''
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        env.config['sig_mode'] = 'ignore'
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
run:
echo ${output} >> temp/out.log
touch ${output}
        ''')
        wf = script.workflow()
        env.config['sig_mode'] = 'ignore'
        Base_Executor(wf).run()
        with open('temp/out.log') as out:
            self.assertEqual(len(out.read().split()), 15)
        shutil.rmtree('temp')

    @multi_attempts
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
        #self.assertLess(time.time() - st, 11)
        FileTarget('lock.sos').remove('both')


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
        elapsed = time.time() - st
        # rerun should be faster
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, elapsed)
        # if we remove the middle result, it should not matter
        os.remove('a.txt')
        st = time.time()
        Base_Executor(wf).run()
        self.assertLess(time.time() - st, elapsed)
        #
        # if we remove the final result, it will be rebuilt
        os.remove('aa.txt')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, 2.5)
        #
        # now we request the generation of target
        FileTarget('a.txt').remove('target')
        FileTarget('aa.txt').remove('both')
        st = time.time()
        Base_Executor(wf).run()
        self.assertGreater(time.time() - st, elapsed - 1)
        #
        FileTarget('a.txt').remove('both')
        FileTarget('aa.txt').remove('both')

    def testReexecutionOfAuxiliaryStep(self):
        '''Test re-execution of auxiliary step if the step has been changed.'''
        script = SoS_Script('''
[process: provides='a.txt']
sh:
    echo 'a.txt' > a.txt

[default]
depends: 'a.txt'
output: 'a.txt.bak'
sh:
    cp a.txt ${output}
''')
        wf = script.workflow()
        #
        Base_Executor(wf).run()
        for f in ['a.txt', 'a.txt.bak']:
            self.assertTrue(FileTarget(f).exists())
            with open(f) as ifile:
                self.assertEqual(ifile.read(), 'a.txt\n')
        # now let us change how a.txt should be generated
        script = SoS_Script('''
[process: provides='a.txt']
sh:
    echo 'aa.txt' > a.txt

[default]
depends: 'a.txt'
output: 'a.txt.bak'
sh:
    cp a.txt ${output}
''')
        wf = script.workflow()
        #
        Base_Executor(wf).run()
        for f in ['a.txt', 'a.txt.bak']:
            self.assertTrue(FileTarget(f).exists())
            with open(f) as ifile:
                self.assertEqual(ifile.read(), 'aa.txt\n')


    def testStoppedOutput(self):
        '''test output with stopped step'''
        for file in ["${a}.txt" for a in range(10)]:
            FileTarget(file).remove('both')

        script = SoS_Script('''
[test_1]
input: for_each={'a': range(10)}
output: "${a}.txt"

stop_if(a % 2 == 0)
run:
    touch ${_output}

[test_2]
assert(len(input) == 5)
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for idx in range(10):
            if idx % 2 == 0:
                self.assertFalse(FileTarget("{}.txt".format(idx)).exists())
            else:
                self.assertTrue(FileTarget("{}.txt".format(idx)).exists())
                FileTarget("${idx}.txt").remove('both')

    def testAllowError(self):
        '''Test option allow error'''
        FileTarget('a.txt').remove('all')
        script = SoS_Script('''
[test]
sh:  allow_error=True
    something_wrong

sh:
    touch a.txt
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(FileTarget('a.txt').exists())
        FileTarget('a.txt').remove('all')

    def testConcurrentWorker(self):
        '''Test the starting of multiple workers #493 '''
        with open('test_script.sos', 'w') as script:
            script.write('''
[10]
input: for_each={'i': range(1)}

[20]
input: for_each={'i': range(2)}
''')
        subprocess.call('sos run test_script', shell=True)
        os.remove('test_script.sos')

if __name__ == '__main__':
    unittest.main()
