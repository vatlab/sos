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
import unittest
import shutil
import subprocess

from sos.utils import env
from sos.parser import SoS_Script
from sos.workflow_executor import Base_Executor, ExecuteError
from sos.targets import file_target

class TestNested(unittest.TestCase):
    def setUp(self):
        env.reset()
        subprocess.call('sos remove -s', shell=True)
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            file_target(f).remove('both')

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

    def testProgressBar(self):
        # progress bar with nested workflow
        script = SoS_Script('''
import time
time.sleep(0)
[sub_1]
[sub_2]
[sub_3]
[sub_4]
[a_1]
[a_2]
[a_3]
sos_run('sub')
[a_4]
[a_5]
''')
        env.verbosity = 1
        wf = script.workflow('a')
        Base_Executor(wf).run()

    def testNestedWorkflow(self):
        '''Test the creation and execution of combined workfow'''
        self.touch(['a.txt', 'b.txt', 'b.begin'])
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
if 'inputs' not in locals():
    inputs = []

[a_1: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(_input)
[a_2: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(_input)
[a_3: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(_input)
[a_4: shared=['executed', 'inputs']]
executed.append(step_name)
output: 'a.done'
inputs.append(_input)
run: expand=True
    touch {_output}
[b_1: shared=['executed', 'inputs']]
executed.append(step_name)
input: 'b.begin'
inputs.append(_input)
[b_2: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(_input)
[b_3: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(_input)
[b_4: shared=['executed', 'inputs']]
executed.append(step_name)
output: 'b.txt'
inputs.append(_input)
[c: shared=['executed', 'inputs']]
executed.append(step_name)
input: 'a.txt'
output: 'b.txt'
inputs.append(_input)
sos_run('a+b', shared=['executed', 'inputs'])
''')
        env.config['sig_mode'] = 'ignore'
        wf = script.workflow('c')
        Base_Executor(wf).run()
        # order of execution is not guaranteed
        self.assertEqual(sorted(env.sos_dict['executed']), sorted(['c_0', 'a_1', 'a_2', 'a_3', 'a_4',
            'b_1', 'b_2', 'b_3', 'b_4']))
        # step will be looped
        self.touch(['a.txt', 'b.txt'])
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
if 'inputs' not in locals():
    inputs = []

[a_1:shared=['executed', 'inputs']]
executed.append(step_name)
output: _input[0] + '.a1'
inputs.append(_input)
run: expand=True
    touch {_output}
[a_2:shared=['executed', 'inputs']]
executed.append(step_name)
output: _input[0] + '.a2'
inputs.append(_input)
run: expand=True
    touch {_output}
[c:shared=['executed', 'inputs']]
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
inputs.append(_input)
sos_run('a', shared=['executed', 'inputs'])
''')
        wf = script.workflow('c')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'a_1', 'a_2', 'a_1', 'a_2'])
        #self.assertEqual(env.sos_dict['inputs'], [['a.txt'], ['a.txt'], ['a.txt.a1'], ['b.txt'], ['b.txt'], ['b.txt.a1']])
        for file in ('a.txt.a1', 'a.txt.a1.a2', 'b.txt.a1', 'b.txt.a1.a2'):
            file_target(file).remove('both')
        #
        # allow specifying a single step
        # step will be looped
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
[a_1:shared='executed']
executed.append(step_name)
[a_2:shared='executed']
executed.append(step_name)
[c_0:shared='executed']
executed.append(step_name)
[c_1:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a:2', shared='executed')
''')
        wf = script.workflow('c')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'c_1', 'a_2',  'a_2'])
        # allow specifying a single step
        # step will be looped
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
[a_1:shared='executed']
executed.append(step_name)
[a_2:shared='executed']
executed.append(step_name)
[c_0:shared='executed']
executed.append(step_name)
[c_1:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a:2', shared='executed')
''')
        wf = script.workflow('c')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'c_1', 'a_2', 'a_2'])
        #
        # recursive subworkflow not allowed
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
[a_1:shared='executed']
executed.append(step_name)
[a_2:shared='executed']
executed.append(step_name)
[c_0:shared='executed']
executed.append(step_name)
[c_1:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a_2+c', shared='executed')
''')
        wf = script.workflow('c')
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        # nested subworkflow is allowed
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
[a_1:shared='executed']
executed.append(step_name)
[a_2:shared='executed']
executed.append(step_name)
[a_3:shared='executed']
executed.append(step_name)
[b_1:shared='executed']
executed.append(step_name)
[b_2:shared='executed']
executed.append(step_name)
sos_run('a:1-2', shared='executed')
[c_0:shared='executed']
executed.append(step_name)
[c_1:shared='executed']
executed.append(step_name)
input: 'a.txt'
sos_run('a+b', shared='executed')
''')
        wf = script.workflow('c')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'c_1', 'a_1', 'a_2', 'a_3',
            'b_1', 'b_2', 'a_1', 'a_2'])
        #
        #
        # nested subworkflow with step option and others
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
[a_1:shared='executed']
executed.append(step_name)
[a_2:shared='executed']
executed.append(step_name)
[a_3:shared='executed']
executed.append(step_name)
[b:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a:3+a:1', shared='executed')
[d:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a:2', shared='executed')
[e2_2:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
''')
        wf = script.workflow('b')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['b_0', 'a_3', 'a_1', 'a_3', 'a_1'])
        wf = script.workflow('d')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['d_0', 'a_2', 'a_2'])
        wf = script.workflow('e2')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['e2_2'])
        #
        # clean up
        file_target('a.done').remove('both')

    def testDynamicNestedWorkflow(self):
        '''Test nested workflow controlled by command line option'''
        script = SoS_Script('''
if 'executed' not in locals():
    executed = []
parameter: wf='a'

[a_1:shared='executed']
executed.append(step_name)
[a_2:shared='executed']
executed.append(step_name)
[a_3:shared='executed']
executed.append(step_name)
[b_1:shared='executed']
executed.append(step_name)
[b_2:shared='executed']
executed.append(step_name)
[b_3:shared='executed']
executed.append(step_name)

[default:shared='executed']
executed.append(step_name)
sos_run(wf, shared='executed')
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--wf', 'b']).run()
        self.assertEqual(env.sos_dict['executed'], [ 'default_0', 'b_1', 'b_2', 'b_3'])
        #
        Base_Executor(wf, args=['--wf', 'a']).run()
        self.assertEqual(env.sos_dict['executed'], ['default_0', 'a_1', 'a_2', 'a_3'])

    def testIncludedNestedWorkFlow(self):
        '''Test the source option of sos_run'''
        # nested subworkflow with step option and others
        self.touch(['a.txt', 'b.txt'])
        #
        with open('inc.sos', 'w') as sos:
            sos.write('''
# test sos script

# global definition
GLB = 5
parameter: parB = 10

[A_1: shared='executed']
executed.append('t.' + step_name)
output: _input[0] + '.a1'
run: expand=True
    touch {_output}

[A_2: shared='executed']
executed.append('t.' + step_name)
output: _input[0] + '.a2'
run: expand=True
    touch {_output}
''')
        script = SoS_Script('''
%from inc include *

if 'executed' not in locals():
    executed = []

[b_1: skip=False, shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('A', shared='executed')
''')
        wf = script.workflow('b')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['GLB'], 5)
        self.assertEqual(env.sos_dict['parB'], 10)
        self.assertEqual(env.sos_dict['executed'], ['b_1', 't.A_1', 't.A_2', 't.A_1', 't.A_2'])
        #
        subprocess.call('sos remove -s', shell=True)
        for file in ('a.txt.a1', 'a.txt.a1.a2', 'b.txt.a1', 'b.txt.a1.a2'):
            file_target(file).remove('both')
        script = SoS_Script('''
%include inc as k

if 'executed' not in locals():
    executed = []

[b_1: skip=False, shared='executed']
executed.append('g.' + step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('k.A', shared='executed')
''')
        wf = script.workflow('b')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['k'].GLB, 5)
        self.assertEqual(env.sos_dict['k'].parB, 10)
        self.assertEqual(env.sos_dict['executed'], ['g.b_1', 't.k.A_1', 't.k.A_2', 't.k.A_1', 't.k.A_2'])
        #
        os.remove('inc.sos')


    def testSoSRun(self):
        '''Test action sos_run with keyword parameters'''
        for f in ['0.txt', '1.txt']:
            file_target(f).remove('both')
        script = SoS_Script(r'''
[A]
parameter: num=5
run: expand=True
    touch {num}.txt

[batch]
for k in range(2):
    sos_run('A', num=k)
''')
        wf = script.workflow('batch')
        Base_Executor(wf).run()
        for f in ['0.txt', '1.txt']:
            self.assertTrue(file_target(f).target_exists())
            file_target(f).remove('both')
        #
        # if we do not pass num, parameter would not change
        for f in ['0.txt', '1.txt']:
            file_target(f).remove('both')
        script = SoS_Script(r'''
[A]
parameter: num=5
run: expand=True
    touch {num}.txt

[batch]
for num in range(2):
    sos_run('A')
''')
        wf = script.workflow('batch')
        Base_Executor(wf).run()
        for f in ['0.txt', '1.txt']:
            self.assertFalse(file_target(f).target_exists())
        self.assertTrue(file_target('5.txt').target_exists())
        file_target('5.txt').remove('both')
        #
        # test parameter shared to send and return vars
        #
        script = SoS_Script(r'''
[A: shared='k']
k += 10

[batch]
for k in range(2):
    sos_run('A', shared='k')
    run(f"touch {k}.txt")
''')
        wf = script.workflow('batch')
        Base_Executor(wf).run()
        for f in ['10.txt', '11.txt']:
            self.assertTrue(file_target(f).target_exists())
            file_target(f).remove('both')

    def testDAGofDynamicNestedWorkflow(self):
        #
        # Because we are not sure which workflows would be executed
        # until run time, the DAG should not contain nested workflow
        # until runtime.
        #
        for f in ['B0.txt', 'B0.txt.p', 'B1.txt', 'B1.txt.p', 'B2.txt', 'B2.txt.p']:
            file_target(f).remove('both')
        #
        #  A1 <- P <- B
        #  A1 <- P <- B
        #  A2
        #
        #  ALL calls A and B with parameter
        #
        script = SoS_Script('''
[A_1]
parameter: num = 2
input: f"B{num}.txt.p"

[B: provides='B{num}.txt']
run: expand=True
    touch 'B{num[0]}.txt'

[P: provides='{filename}.p']
input: filename
run: expand=True
    touch {_output}

[ALL]

for i in range(3):
    sos_run('A', num=i)


''')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow('ALL')
        Base_Executor(wf).run()
        for f in ['B0.txt', 'B0.txt.p', 'B1.txt', 'B1.txt.p', 'B2.txt', 'B2.txt.p']:
            self.assertTrue(file_target(f).target_exists())
            file_target(f).remove('both')


    def testPassingVarsToNestedWorkflow(self):
        '''Test if variables can be passed to nested workflows'''
        script = SoS_Script(r"""
import time
import random

[nested]
print(f'I am nested {nested} with seed {seed}')

[0]
reps = range(5)
input: for_each='reps'
import random
nested = _reps
seed = random.randint(1, 1000)
print(f'Passing {seed} to {nested}')
sos_run('nested', nested=nested, seed=seed)

""")
        wf = script.workflow()
        Base_Executor(wf).run()

    def testUserDefinedFunc(self):
        '''Test the use of user-defined functions in SoS script'''
        script = SoS_Script(r"""

def myfunc():
  return 'a'

[1: shared={'test':'_output'}]
output: myfunc()

myfunc()

""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['test'], ['a'])
        # User defined function should also work under nested workflows
        # This is difficult because the 'local namespace' is usually
        # not seen inside function definition. The solution now is to
        # use a single workspace.
        script = SoS_Script(r"""

def myfunc():
    # test if builtin functions (sum and range) can be used here.
    return 'a' + str(sum(range(10)))

[1: shared={'test':'_output'}]
output: [myfunc() for i in range(10)][0]

myfunc()

""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['test'], ['a45'])


    def testConfigFileOfNestedWorkflow(self):
        '''Test passing of configurationg to nested workflow'''
        script = SoS_Script('''
[test_1]
parameter: key = None
print(CONFIG[key])

[default_1]
sos_run('test:1', key = '1')
    ''')
        with open('test.conf', 'w') as conf:
            conf.write("""{'1':'hi'}""")
        wf = script.workflow()
        Base_Executor(wf, config={'config_file': 'test.conf'}).run()


    def testErrorFromSubworkflow(self):
        '''Test if error from subworkflow is passed to master (#396)'''
        script = SoS_Script('''
[test_1]
R:
  set.seed(xxx)

[default]
sos_run('test')
''')
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)

    def testFunDef(self):
        '''Test defintion of function that can be used by other steps'''
        self.touch(['aa.txt', 'ab.txt'])
        # in nested workflow?
        script = SoS_Script(r"""
def myfunc(a):
    return ['a' + x for x in a]

[mse: shared={'test':'_output'}]
input: myfunc(['a.txt', 'b.txt'])

[1]
sos_run('mse')
""")
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        #
        # Names defined in subworkflow is not returned to the master dict
        self.assertTrue('test' not in env.sos_dict)

    def testSearchPath(self):
        '''Test if any action should exit in five seconds in dryrun mode'''
        sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yml')
        shutil.copy(sos_config_file, 'test.yml')
        #
        subprocess.call('sos config --set sos_path {0}/crazy_path {0}/crazy_path/more_crazy/'.format(os.getcwd()), shell=True)
        #
        if not os.path.isdir('crazy_path'):
            os.mkdir('crazy_path')
            os.mkdir(os.path.join('crazy_path', 'more_crazy'))
        with open(os.path.join('crazy_path', 'crazy_master.sos'), 'w') as crazy:
            crazy.write('''
[0]
sos_run('cc', source='crazy_slave.sos')

''')
        with open(os.path.join('crazy_path', 'more_crazy', 'crazy_slave.sos'), 'w') as crazy:
            crazy.write('''
[cc_0]
print('hay, I am crazy')
''')

        script = SoS_Script(filename='crazy_master.sos')
        script.workflow()
        #
        shutil.rmtree('crazy_path')
        os.remove(sos_config_file)
        shutil.copy('test.yml', sos_config_file)


    def testNestedWorkdir(self):
        '''Test nested runtime option for work directory'''
        if os.path.isdir('tmp'):
            shutil.rmtree('tmp')
        script = SoS_Script('''
[step]
task: workdir='tmp'
run:
    touch 'a.txt'

[default]
sos_run('step', workdir='tmp')
''')
        wf = script.workflow()
        # this should be ok.
        Base_Executor(wf).run()
        os.path.isfile('tmp/tmp/a.txt')
        shutil.rmtree('tmp')


    def testFailureOfNestedWorkflow(self):
        '''Test failure of nested workflow #838'''
        if os.path.isfile('a.txt'):
            os.remove('a.txt')
        script = SoS_Script('''
[something]
input: 'a.txt'

[default]
sos_run('something')
''')
        wf = script.workflow()
        # this should be ok.
        self.assertRaises(Exception, Base_Executor(wf).run)


    def testNestedFromAnotherFile(self):
        '''Test nested runtime option for work directory'''
        if os.path.isdir('a.txt'):
            os.remove('a.txt')
        with open('another.sos', 'w') as another:
            another.write('''
[whatever]
run:
    touch 'a.txt'

''')
        script = SoS_Script('''
[default]
sos_run('whatever', source='another.sos')
''')
        wf = script.workflow()
        # this should be ok.
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('a.txt'), 'a.txt should have been created by nested workflow from another file')


if __name__ == '__main__':
    #suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestParser)
    #unittest.TextTestRunner(, suite).run()
    unittest.main()
