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
import unittest
import shutil
import subprocess

from pysos.utils import env
from pysos.sos_script import SoS_Script, ParsingError
from pysos.sos_executor import Sequential_Executor, ArgumentError, ExecuteError

class TestParser(unittest.TestCase):
    def setUp(self):
        env.reset()

    def testFileFormat(self):
        '''Test recognizing the format of SoS script'''
        # file format must be 'fileformat=SOSx.x'
        self.assertRaises(ParsingError, SoS_Script,
            '#fileformat=SS2')
        self.assertRaises(ParsingError, SoS_Script,
            '#fileformat=SOS1.0beta')
        #
        # parse a larger script with gormat 1.1
        script = SoS_Script(filename='scripts/section1.sos')
        # not the default value of 1.0
        self.assertEqual(script.format_version, '1.1')

    def testStringLiteral(self):
        '''Test string literals of SoS'''
        env.shared_vars = ['a', 'b', 'c']
        script = SoS_Script(r"""
[0]
a = 'a\n'
b = "a\n"
c = '''a\n

b'''
""")
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['a'], 'a\n')
        self.assertEqual(env.sos_dict['b'], 'a\n')
        # MAYJOR difference
        self.assertEqual(env.sos_dict['c'], 'a\\n\nb')
        script = SoS_Script(r'''
[0]
c = """a\n

b"""
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        # Note the difference between """ and ''' quotes
        self.assertEqual(env.sos_dict['c'], 'a\n\nb')

    def testWorkflows(self):
        '''Test workflows defined in SoS script'''
        script = SoS_Script('''[0]''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[1]''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[*_1]''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[*_1]\n[auxillary:target='*.txt']''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[*_1]\n[human_1]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human'])
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])
        # skip option
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2:skip]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human'])
        # unnamed
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])

    def testSkipStep(self):
        '''Test the skip option to skip certain steps'''
        script = SoS_Script('''
[parameters]
skip = 0

[0: alias='a', skip=skip==0]
var = 0

[1: alias='b', skip=skip==1]
var = 1

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--skip', '0'])
        self.assertEqual(env.sos_dict['b'].var, 1)
        #
        Sequential_Executor(wf).run(args=['--skip', '1'])
        self.assertEqual(env.sos_dict['a'].var, 0)

    def testSections(self):
        '''Test section definitions'''
        # bad names
        for badname in ['56_1', '_a', 'a_', '1x', '*', '?']:
            self.assertRaises(ParsingError, SoS_Script, '[{}]'.format(badname))
        # bad options
        for badoption in ['ss', 'skip a', 'skip:_', 'skip, skip']:
            self.assertRaises(ParsingError, SoS_Script, '[0:{}]'.format(badoption))
        # option value should be a valid python expression
        for badoption in ['sigil=a', 'sigil="[]"', 'sigil="| |"']:
            self.assertRaises(ParsingError, SoS_Script, '[0:{}]'.format(badoption))
        # good options
        for goodoption in ['sigil="[ ]"', 'alias="a"']:
            SoS_Script('[0:{}]'.format(goodoption))
        # allowed names
        for name in ['a5', 'a_5', '*_0', 'a*1_100']:
            SoS_Script('[{}]'.format(name))
        # no directive in global section
        self.assertRaises(ParsingError, SoS_Script,
            '''input: 'filename' ''')
        # duplicate sections
        self.assertRaises(ParsingError, SoS_Script,
            '''[1]\n[1]''')
        self.assertRaises(ParsingError, SoS_Script,
            '''[1]\n[3]\n[2,1]''')
        self.assertRaises(ParsingError, SoS_Script,
            '''[a_1]\n[a_3]\n[*_1]''')
        # no duplicated section header
        SoS_Script('''[a_1]\n[a_3]\n[b*_1]''')

    def testGlobalVariables(self):
        '''Test definition of variables'''
        # global section cannot have directive
        self.assertRaises(ParsingError, SoS_Script,
            '''input: 'filename' ''')
        # allow definition
        SoS_Script('''a = '1' ''')
        SoS_Script('''a = ['a', 'b'] ''')
        # but this one has incorrect syntax
        self.assertRaises(ParsingError, SoS_Script,
            '''a = 'b  ''')
        # This one also does not work because b is not defined.
        delattr(env, 'sos_dict')
        script = SoS_Script('''a = b\n[0] ''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Sequential_Executor(wf).run)
        # multi-line string literal
        SoS_Script('''a = """
this is a multi line
string """
''')
        # multi-line list literal, even with newline in between
        SoS_Script('''a = [
'a',

'b'
]
''')
        #
        script = SoS_Script(filename='scripts/section1.sos')
        # not the default value of 1.0

    def testParameters(self):
        '''Test parameters section'''
        # directive not allowed in parameters
        self.assertRaises(ParsingError, SoS_Script,
            '''
[parameters]
# comment
input: 'filename' 
''')
        self.assertRaises(ParsingError, SoS_Script,
            '''
[parameters]
# comment
func()
''')    
        script = SoS_Script(filename='scripts/section1.sos')
        wf = script.workflow('chapter_0')
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run,
            args=['--not_exist'])
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run,
            args=['--par1', 'a', 'b'])
        # 
        script = SoS_Script('''
[parameters]
a = [1, 2]
[0]
''')    
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['a'], [1,2])
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--a', '3'])
        self.assertEqual(env.sos_dict['a'], [3])
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--a', '3', '5'])
        self.assertEqual(env.sos_dict['a'], [3, 5])
        #
        script = SoS_Script('''
[parameters]
# comment
# comment
a = ['a.txt', 'b.txt']
[0]
''')    
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['a'], ['a.txt', 'b.txt'])
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--a', '3'])
        self.assertEqual(env.sos_dict['a'], ['3'])
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--a', '3', '5'])
        self.assertEqual(env.sos_dict['a'], ['3', '5'])
        #
        # test parameter using global definition
        script = SoS_Script('''
a="100"

[parameters]
# comment
b=str(int(a)+1)
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], '101')
        #
        env.sos_dict.clear()
        script = SoS_Script('''
a=100

[parameters]
b=a+1
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], 101)
        #
        script = SoS_Script('''
a=100

[parameters]
# comment
b=a+1
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run, args=['--b', 'a'])
        #
        script = SoS_Script('''
a=100

[parameters]
b=a+1.
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', '1000'])
        #
        self.assertEqual(env.sos_dict['b'], 1000)
        #
        # test string interpolation of the parameter section
        script = SoS_Script('''
a=100

[parameters]
# comment

b='${a+1}'
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], '101')
        # test alternative sigil
        script = SoS_Script('''
a=100

[parameters: sigil='[ ]']
b='[a+1]'
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], '101')
        #
        # argument has hve a value
        self.assertRaises(ParsingError, SoS_Script, '''
[parameters]
b=
''')
        # if it is a type, must provide value
        script = SoS_Script('''
[parameters]
# comment
b = int
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run)
        #
        script = SoS_Script('''
[parameters]
b = list
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run)
        # also require the type
        script = SoS_Script('''
[parameters]
b = int
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run, args=['--b', 'file'])
        #
        script = SoS_Script('''
[parameters]
b = int
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', '5'])
        self.assertEqual(env.sos_dict['b'], 5)
        # list is ok
        script = SoS_Script('''
[parameters]
b = list
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', '5'])
        self.assertEqual(env.sos_dict['b'], ['5'])
        # bool
        script = SoS_Script('''
[parameters]
# comment
b = bool
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', 't'])
        self.assertEqual(env.sos_dict['b'], True)
        script = SoS_Script('''
[parameters]
b = True
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', 'False'])
        self.assertEqual(env.sos_dict['b'], False)
        script = SoS_Script('''
[parameters]
b = True
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', '1'])
        self.assertEqual(env.sos_dict['b'], True)
        script = SoS_Script('''
[parameters]
b = bool
''')
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--b', 'no'])
        self.assertEqual(env.sos_dict['b'], False)
        #
        # should fail for undefined variables
        script = SoS_Script('''
[parameters]
a = 5
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Sequential_Executor(wf).run,
            args=['--b', 'file'])
        # and for cases without parameter section
        script = SoS_Script('''
[0]
''')
        wf = script.workflow()
        self.assertRaises(RuntimeError, Sequential_Executor(wf).run,
            args=['--b', 'file'])

    def testSectionVariables(self):
        '''Test variables in sections'''
        # directive name cannot be used as variable
        self.assertRaises(ParsingError, SoS_Script,
            '''[0]
input='a.txt' ''')
        self.assertRaises(ParsingError, SoS_Script,
            '''[0]
output='a.txt' ''')
        self.assertRaises(ParsingError, SoS_Script,
            '''[0]
depends='a.txt' ''')

    def testSectionDirectives(self):
        '''Test directives of sections'''
        # cannot be in the global section
        self.assertRaises(ParsingError, SoS_Script,
            '''input: 'filename' ''')
        # multi-line OK
        SoS_Script('''
[0]
input: 'filename',
    'filename1'

''')
        # An abusive case with multi-line OK, from first column ok, blank line ok
        SoS_Script('''
[0]
input: 'filename',
'filename1',

filename4,
opt1=value
output: 
    blah

depends:
'something else'
''')
        # option with expression ok
        SoS_Script('''
[0]
input: 'filename',  'filename2', opt=value==1

''')
        # unrecognized directive, allowed now
        SoS_Script('''
[0]
something: 'filename',  filename2, opt=value==1
''')
        # need commma
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
input: 'filename'  filename2
''')
        # can be after action
        SoS_Script('''
[0]
func()        
input: 'filename',  'filename2', opt=value==1
''')
        # assignments between directives are allowed
        SoS_Script('''
[0]
input: 'filename',  'filename2', opt=value==1
a = 'some text'
output: 'filename',  'filename2', opt=value==1
''')
        # can have action between directives
        SoS_Script('''
[0]
input: 'filename',  'filename2', opt=value==1
abc
output: 'filename',  'filename2', opt=value==1
''')

    def testScriptFormat(self):
        '''Test putting scripts directly under action'''
        script = SoS_Script('''
[0]
input: 'filename',  'filename2', opt=value==1
R:

open('something')
save.put()

''')
        script = SoS_Script('''
[0]
input: 'filename',  'filename2', opt=value==1
R: concurrent = True

open('something')
save.put()

''')
        script = SoS_Script('''
[0]
input: 'filename',  'filename2', opt=value==1
R: concurrent = True,
    workdir = 'someelse else'

open('something')
save.put()

''')
        # test dedent
        script = SoS_Script('''
[0]
python3:
    from pysos import logger
    logger.warning('I am from a dented text')
    if 1:
        logger.warning('Another warning')
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        # with section head in the script, 
        # this will not work even if the embedded
        # python script is perfectly valid.
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
input: 'filename',  'filename2', opt=value==1
python3:

with open('something') as e:
   e.write("""
[section]
""")

''')

    def testInput(self):
        '''Test input directive'''
        script = SoS_Script('''
[0]
files = ['a.txt', 'b.txt']

input: 'a.pdf', files, skip=False

''')
        wf = script.workflow('default')
        Sequential_Executor(wf).run(run_mode='inspect')
        #
        # test input types
        script = SoS_Script('''
[0:alias='test']
files = ('a${i}' for i in range(2))
input: {'a.txt', 'b.txt'}, files
output: ('a${x}' for x in _input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(sorted(env.sos_dict['test'].input), ['a.txt', 'a0', 'a1', 'b.txt'])
        self.assertEqual(sorted(env.sos_dict['test'].output), ['aa.txt', 'aa0', 'aa1', 'ab.txt'])


    def testGroupBy(self):
        '''Test group_by parameter of step input'''
        # group_by = 'all'
        env.shared_vars = ['executed']
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='all'

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a2.txt', 'a3.txt', 'a4.txt']])
        # group_by = 'single'
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='single'

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt'], ['a2.txt'], ['a3.txt'], ['a4.txt']])
        # group_by = 'pairs'
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='pairs'

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a3.txt'], ['a2.txt', 'a4.txt']])
        # group_by = 'pairwise'
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='pairwise'

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a2.txt'], ['a2.txt', 'a3.txt'], ['a3.txt', 'a4.txt']])
        # group_by = 'combinations'
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='combinations'

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a2.txt'], ['a1.txt', 'a3.txt'], 
            ['a1.txt', 'a4.txt'], ['a2.txt', 'a3.txt'], ['a2.txt', 'a4.txt'], ['a3.txt', 'a4.txt']])
        # group_by chunks specified as integers
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=3

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],
                         [['a1.txt', 'a2.txt', 'a3.txt'],
                         ['a4.txt', 'a5.txt', 'a6.txt'],
                         ['a7.txt', 'a8.txt', 'a9.txt']])
        # group_by chunks specified as integer strings
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='3'

executed.append(_input)

''')
        wf = script.workflow()
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'],
                         [['a1.txt', 'a2.txt', 'a3.txt'],
                         ['a4.txt', 'a5.txt', 'a6.txt'],
                         ['a7.txt', 'a8.txt', 'a9.txt']])
        # number of files should be divisible by group_by
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=4

executed.append(_input)

''')
        wf = script.workflow()
        self.assertRaises((ExecuteError, RuntimeError), Sequential_Executor(wf).run, run_mode='inspect')
        # incorrect value causes an exception
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='something'

executed.append(_input)

''')
        wf = script.workflow()
        self.assertRaises((ExecuteError, RuntimeError), Sequential_Executor(wf).run, run_mode='inspect')

    def testSectionActions(self):
        '''Test actions of sections'''
        SoS_Script(
            """
[0]
func('''
multiline 
string''', with_option=1
)
""")
        self.assertRaises(ParsingError, SoS_Script,
            '''
[0]
func(
''')

    def testDescriptions(self):
        '''Test script and workflow descriptions'''
        script = SoS_Script('''# first block

# global
# description

# human
# description of human

# description of human continued

[human_1]

a = '1'

# mouse
# mouse description
#

[mouse_1]
''')
        self.assertEqual(script.description, 'global\ndescription\n\n')
        self.assertEqual(script.workflow('human').description, 'description of human\ndescription of human continued\n')
        self.assertEqual(script.workflow('mouse').description, 'mouse description\n')

    def testLongerCode(self):
        '''Test definition of classes (with intermediate newlines) in step.'''
        script = SoS_Script('''# first block

[0]
class A:
    def __init__(self):
        pass

    # the newline above should be fine because SoS treat this as
    # regular lines
    def __call__(self):
        return 0

b = A()()

''')
        env.shared_vars=['b']
        wf = script.workflow()
        Sequential_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], 0)

    def testCombinedWorkflow(self):
        '''Test the creation and execution of combined workfow'''
        env.shared_vars = ['executed', 'input_b1', 'a']
        script = SoS_Script('''
a0 = 0
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[parameters]
a = a0 + 1
[a_1]
[a_2]
[a_3]
[a_4]
output: 'out_a_4'
[b_1]
input_b1 = input
[b_2]
[b_3]
[b_4]
[c]
[d]
''')
        wf = script.workflow('a+b')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['a_parameters', 'a_1', 'a_2', 'a_3', 'a_4', 'b_parameters', 'b_1', 'b_2', 'b_3', 'b_4'])
        self.assertEqual(env.sos_dict['a'], 1)
        self.assertEqual(env.sos_dict['input_b1'], ['out_a_4'])
        #
        wf = script.workflow('a_ 1-2 + a_4 + b_3-')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['a_parameters', 'a_1', 'a_2', 'a_parameters', 'a_4', 
            'b_parameters', 'b_3', 'b_4'])
        #
        wf = script.workflow('a+c+d')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['a_parameters', 'a_1', 'a_2', 'a_3', 'a_4', 'c_parameters', 'c_0', 'd_parameters', 'd_0'])

    def testNestedWorkflow(self):
        '''Test the creation and execution of combined workfow'''
        env.shared_vars = ['executed', 'inputs']
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
if 'inputs' not in locals():
    inputs = []

[paramters]
[a_1]
inputs.append(input)
[a_2]
inputs.append(input)
[a_3]
inputs.append(input)
[a_4]
output: 'a.done'
inputs.append(input)
[b_1]
input: 'b.begin'
inputs.append(input)
[b_2]
inputs.append(input)
[b_3]
inputs.append(input)
[b_4]
output: 'b.txt'
inputs.append(input)
[c]
input: 'a.txt'
output: 'b.txt'
inputs.append(input)
sos_run('a+b')
''')
        wf = script.workflow('c')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'a_1', 'a_2', 'a_3', 'a_4',
            'b_1', 'b_2', 'b_3', 'b_4'])
        self.assertEqual(env.sos_dict['inputs'], [['a.txt'], ['a.txt'], [], [], [], ['b.begin'], [], [], []])
        # step will be looped
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
if 'inputs' not in locals():
    inputs = []

[a_1]
output: _input[0] + '.a1'
inputs.append(input)
[a_2]
output: _input[0] + '.a2'
inputs.append(input)
[c]
input: 'a.txt', 'b.txt', group_by='single'
inputs.append(_input)
sos_run('a')
''')
        wf = script.workflow('c')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'a_1', 'a_2', 'a_1', 'a_2'])
        #self.assertEqual(env.sos_dict['inputs'], [['a.txt'], ['a.txt'], ['a.txt.a1'], ['b.txt'], ['b.txt'], ['b.txt.a1']])
        #
        # allow specifying a single step
        # step will be looped
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[a_1]
[a_2]
[c_0]
[c_1]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a_2')
''')
        wf = script.workflow('c')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'c_1', 'a_2',  'a_2'])
        # allow specifying a single step
        # step will be looped
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[a_1]
[a_2]
[c_0]
[c_1]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a_2')
''')
        wf = script.workflow('c')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'c_1', 'a_2', 'a_2'])
        #
        # recursive subworkflow not allowed
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[a_1]
[a_2]
[c_0]
[c_1]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a_2+c')
''')
        wf = script.workflow('c')
        self.assertRaises((ExecuteError, RuntimeError), Sequential_Executor(wf).run, run_mode='inspect')
        #
        # nested subworkflow is allowed
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[a_1]
[a_2]
[a_3]
[b_1]
[b_2]
sos_run('a_1+a_2')
[c_0]
[c_1]
input: 'a.txt'
sos_run('a+b')
''')
        wf = script.workflow('c')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'c_1', 'a_1', 'a_2', 'a_3',
            'b_1', 'b_2', 'a_1', 'a_2'])
        #
        #
        # nested subworkflow with step option and others
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[a_1]
[a_2]
[a_3]
[b]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a_3+a_1')
[d]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a_2')
[e2_2]
input: 'a.txt', 'b.txt', group_by='single'
''')
        wf = script.workflow('b')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['b_0', 'a_3', 'a_1', 'a_3', 'a_1'])
        wf = script.workflow('d')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['d_0', 'a_2', 'a_2'])
        wf = script.workflow('e2')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['e2_2'])

    def testDynamicNestedWorkflow(self):
        '''Test nested workflow controlled by command line option'''
        script = SoS_Script('''
if 'executed' in locals():
    executed.append(step_name)
else:
    executed = [step_name]
[parameters]
wf='a'

[a_1]
[a_2]
[a_3]
[b_1]
[b_2]
[b_3]

[default]
sos_run(wf)
''')
        env.shared_vars = ['executed']
        wf = script.workflow()
        Sequential_Executor(wf).run(args=['--wf', 'b'], run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['default_parameters', 'default_0', 'b_parameters', 'b_1', 'b_2', 'b_3'])
        #
        Sequential_Executor(wf).run(args=['--wf', 'a'], run_mode='inspect')
        self.assertEqual(env.sos_dict['executed'], ['default_parameters', 'default_0', 'a_parameters', 'a_1', 'a_2', 'a_3'])

    def testSourceOption(self):
        '''Test the source option of sos_run'''
        # nested subworkflow with step option and others
        env.shared_vars = ['executed', 'GLB', 'parB']
        if not os.path.isdir('temp'):
            os.mkdir('temp')
        with open('temp/test.sos', 'w') as sos:
            sos.write('''
# test sos script

# global definition
GLB = 5
if 'executed' in locals():
    executed.append('t.' + step_name)
else:
    executed = ['t.' + step_name]

[parameters]
parB = 10

[A_1]
output: _input[0] + 'a1'

[A_2]
output: _input[0] + 'a2'

''')
        script = SoS_Script('''
if 'executed' in locals():
    executed.append('g.' + step_name)
else:
    executed = ['g.' + step_name]
[b_1: skip=False]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('A', source='temp/test.sos')
''')
        wf = script.workflow('b')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['GLB'], 5)
        self.assertEqual(env.sos_dict['parB'], 10)
        self.assertEqual(env.sos_dict['executed'], ['g.b_1', 't.A_parameters', 't.A_1', 't.A_2', 't.A_parameters', 't.A_1', 't.A_2'])
        #
        script = SoS_Script('''
if 'executed' in locals():
    executed.append('g.' + step_name)
else:
    executed = ['g.' + step_name]
[b_1: skip=False]
input: 'a.txt', 'b.txt', group_by='single'
sos_run('k.A', source={'k':'temp/test.sos'})
''')
        wf = script.workflow('b')
        Sequential_Executor(wf).run(run_mode='inspect')
        self.assertEqual(env.sos_dict['GLB'], 5)
        self.assertEqual(env.sos_dict['parB'], 10)
        self.assertEqual(env.sos_dict['executed'], ['g.b_1', 't.A_parameters', 't.A_1', 't.A_2', 't.A_parameters', 't.A_1', 't.A_2'])
        #
        shutil.rmtree('temp')


    def testYAMLConfig(self):
        '''Test config file in yaml format'''
        with open('config.yaml', 'w') as config:
            config.write('''
# Lines beginning with # are skipped when the JSON is parsed, so we can
# put comments into our JSON configuration files
{
    StoreOwner : "John Doe",

    # List of items that we sell
    Fruits: [ "apple", "banana", "pear" ],
    Price: 1.05
}            
''')
        with open('config.sos', 'w') as sos:
            sos.write('''
[0]
print(CONFIG['StoreOwner'])
print(CONFIG.get('StoreOwner', 'something'))
print(CONFIG.get('StoreOwnerSpouse', 'someone else'))
print(CONFIG.StoreOwner)
'''
)
        # run the command
        self.assertEqual(subprocess.call('sos run config.sos -c config.yaml', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        # now test the value
        script = SoS_Script(filename='config.sos')
        wf = script.workflow()
        Sequential_Executor(wf).run(config_file='config.yaml')
        self.assertEqual(env.sos_dict['CONFIG']['Price'], 1.05)
        self.assertEqual(env.sos_dict['CONFIG']['StoreOwner'], 'John Doe')
        self.assertEqual(env.sos_dict['CONFIG']['Fruits'], ['apple', 'banana', 'pear'])
        # configuration items should be readonly
        with open('config.sos', 'w') as sos:
            sos.write('''
[0]
CONFIG['a'] = 'b'
'''
)
        # the command would fail with error message
        # ERROR: Failed to process statement CONFIG['a'] = 'b'
        # : Cannot modify a readonly dictionary.
        self.assertEqual(subprocess.call('sos run config.sos -c config.yaml', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        #
        with open('config.sos', 'w') as sos:
            sos.write('''
[0]
CONFIG.a = 'b'
'''
)
        # the command would fail with error message
        # ERROR: Failed to process statement CONFIG['a'] = 'b'
        # : Cannot modify a readonly dictionary.
        self.assertEqual(subprocess.call('sos run config.sos -c config.yaml', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        #
        for filename in ['config.sos', 'config.yaml']:
            os.remove(filename)

    def testReport(self):
        '''Test report lines'''
        script = SoS_Script('''
[0]
! this is report
!
! this is another line 
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        #
        self.assertRaises(ParsingError, SoS_Script, '''
[parameters]
! this is report in parameters section, which is not alloed
[0]
''')
        #

    def testVarOutput(self):
        '''Test early appearance of variable output'''
        script = SoS_Script('''
[0]
seq = range(3)
input: for_each='seq'
output: 'test${_seq}.txt'
print(output)
''')
        wf = script.workflow()
        # this does not work before until we make variable output available sooner
        Sequential_Executor(wf).run(run_mode='inspect')
        # however, output should be the same in task
        script = SoS_Script('''
[0]
seq = range(3)
input: for_each='seq'
output: 'test${_seq}.txt'
assert(len(output), _index + 1)
task:
assert(len(output), 3)
''')
        wf = script.workflow()
        # this does not work before until we make variable output available sooner
        Sequential_Executor(wf).run(run_mode='inspect')

if __name__ == '__main__':
    unittest.main()
