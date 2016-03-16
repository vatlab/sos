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

class TestParser(unittest.TestCase):
    def testFileFormat(self):
        '''Test recognizing the format of SoS script'''
        # file format must be 'fileformat=SOSx.x'
        self.assertRaises(ParsingError, SoS_Script,
            '#fileformat=SS2')
        self.assertRaises(ParsingError, SoS_Script,
            '#fileformat=SOS1.0beta')
        #
        # parse a larger script with gormat 1.1
        script = SoS_Script('scripts/section1.sos')
        # not the default value of 1.0
        self.assertEqual(script.format_version, '1.1')

    def testStringLiteral(self):
        '''Test string literals of SoS'''
        script = SoS_Script(r"""
a = 'a\n'
b = "a\n"
c = '''a\n

b'''
""")
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['a'], 'a\n')
        self.assertEqual(env.locals['b'], 'a\n')
        # MAYJOR difference
        self.assertEqual(env.locals['c'], 'a\\n\nb')
        script = SoS_Script(r'''
c = """a\n

b"""
''')
        wf = script.workflow()
        wf.run()
        # Note the difference between """ and ''' quotes
        self.assertEqual(env.locals['c'], 'a\n\nb')

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
        #
        # skip option
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2:skip]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human'])
        # unnamed
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])

    def testSections(self):
        '''Test section definitions'''
        # bad names
        for badname in ['56_1', '_a', 'a_', '1x', '*', '?']:
            self.assertRaises(ParsingError, SoS_Script, '[{}]'.format(badname))
        # bad options
        for badoption in ['ss', 'skip a', 'skip:_', 'skip, skip']:
            self.assertRaises(ParsingError, SoS_Script, '[0:{}]'.format(badoption))
        # option value should be a valid python expression
        for badoption in ['sigil=a', 'alias=a', 'sigil="[]"', 'sigil="| |"']:
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
        # or unrecognized directive
        self.assertRaises(ParsingError, SoS_Script,
            '''inputs: 'filename' ''')
        # allow definition
        SoS_Script('''a = '1' ''')
        SoS_Script('''a = ['a', 'b'] ''')
        # but this one has incorrect syntax
        self.assertRaises(ParsingError, SoS_Script,
            '''a = 'b  ''')
        # This one also does not work because b is not defined.
        script = SoS_Script('''a = b  ''')
        wf = script.workflow()
        self.assertRaises(RuntimeError, wf.run)
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
        script = SoS_Script('scripts/section1.sos')
        # not the default value of 1.0

    def testParameters(self):
        '''Test parameters section'''
        # directive not allowed in parameters
        self.assertRaises(ParsingError, SoS_Script,
            '''
[parameters]
input: 'filename' 
''')
        self.assertRaises(ParsingError, SoS_Script,
            '''
[parameters]
func()
''')    
        script = SoS_Script('scripts/section1.sos')
        self.assertRaises(ArgumentError, script.workflow('chapter_0').run,
            args=['--not_exist'])
        self.assertRaises(ArgumentError, script.workflow('chapter_0').run,
            args=['--par1', 'a', 'b'])
        # 
        # test parameter using global definition
        script = SoS_Script('''
a="100"

[parameters]
b=str(int(a)+1)
''')
        script.workflow().run()
        self.assertEqual(env.locals['b'], '101')
        script = SoS_Script('''
a=100

[parameters]
b=a+1
''')
        script.workflow().run()
        self.assertEqual(env.locals['b'], 101)
        #
        script = SoS_Script('''
a=100

[parameters]
b=a+1
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, wf.run, args=['--b', 'a'])
        #
        script = SoS_Script('''
a=100

[parameters]
b=a+1.
''')
        wf = script.workflow()
        wf.run(args=['--b', '1000'])
        #
        self.assertEqual(env.locals['b'], 1000)
        #
        # test string interpolation of the parameter section
        script = SoS_Script('''
a=100

[parameters]
b='${a+1}'
''')
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['b'], '101')
        # test alternative sigil
        script = SoS_Script('''
a=100

[parameters: sigil='[ ]']
b='[a+1]'
''')
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['b'], '101')
        #
        # argument has hve a value
        self.assertRaises(ParsingError, SoS_Script, '''
[parameters]
b=
''')
        # if it is a type, must provide value
        script = SoS_Script('''
[parameters]
b = int
''')
        self.assertRaises(ArgumentError, script.workflow().run)
        #
        script = SoS_Script('''
[parameters]
b = list
''')
        self.assertRaises(ArgumentError, script.workflow().run)
        # also require the type
        script = SoS_Script('''
[parameters]
b = int
''')
        self.assertRaises(ArgumentError, script.workflow().run, args=['--b', 'file'])
        #
        script = SoS_Script('''
[parameters]
b = int
''')
        wf = script.workflow()
        wf.run(args=['--b', '5'])
        self.assertEqual(env.locals['b'], 5)
        # list is ok
        script = SoS_Script('''
[parameters]
b = list
''')
        wf = script.workflow()
        wf.run(args=['--b', '5'])
        self.assertEqual(env.locals['b'], ['5'])
        # bool
        script = SoS_Script('''
[parameters]
b = bool
''')
        wf = script.workflow()
        wf.run(args=['--b', 't'])
        self.assertEqual(env.locals['b'], True)
        script = SoS_Script('''
[parameters]
b = True
''')
        wf = script.workflow()
        wf.run(args=['--b', 'False'])
        self.assertEqual(env.locals['b'], False)
        script = SoS_Script('''
[parameters]
b = True
''')
        wf = script.workflow()
        wf.run(args=['--b', '1'])
        self.assertEqual(env.locals['b'], True)
        script = SoS_Script('''
[parameters]
b = bool
''')
        wf = script.workflow()
        wf.run(args=['--b', 'no'])
        self.assertEqual(env.locals['b'], False)
        #
        # should fail for undefined variables
        script = SoS_Script('''
[parameters]
a = 5
''')
        self.assertRaises(ArgumentError, script.workflow().run,
            args=['--b', 'file'])
        # and for cases without parameter section
        script = SoS_Script('''
''')        
        self.assertRaises(ArgumentError, script.workflow().run,
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
        # unrecognized directive
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
something: 'filename',  filename2, opt=value==1
''')
        # need commma
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
input: 'filename'  filename2
''')
        # cannot be after action
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
func()        
input: 'filename',  'filename2', opt=value==1
''')
        # cannot be assigned between directives
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
input: 'filename',  'filename2', opt=value==1
a = 'some text'
output: 'filename',  'filename2', opt=value==1
''')
        # cannot be action between directives
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
input: 'filename',  'filename2', opt=value==1
abc
output: 'filename',  'filename2', opt=value==1
''')

    def testInput(self):
        '''Test input directive'''
        script = SoS_Script('''
[0]
files = ['a.txt', 'b.txt']

input: 'a.pdf', files, skip=False

''')
        env.run_mode = 'dryrun'
        script.workflow('default').run()

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
        wf = script.workflow()
        wf.run()
        self.assertEqual(env.locals['b'], 0)

    def testCombinedWorkflow(self):
        '''Test the creation and execution of combined workfow'''
        script = SoS_Script('''
executed = []
a = 0
[parameters]
a = a + 1
[a_1]
executed.append(_step.name)
[a_2]
executed.append(_step.name)
[a_3]
executed.append(_step.name)
[a_4]
executed.append(_step.name)
[b_1]
executed.append(_step.name)
[b_2]
executed.append(_step.name)
[b_3]
executed.append(_step.name)
[b_4]
executed.append(_step.name)
[c]
executed.append(_step.name)
[d]
executed.append(_step.name)
''')
        wf = script.workflow('a+b')
        env.verbosity = 4
        wf.run()
        self.assertEqual(env.locals['executed'], ['a_1', 'a_2', 'a_3', 'a_4', 'b_1', 'b_2', 'b_3', 'b_4'])
        self.assertEqual(env.locals['a'], 1)
        #
        wf = script.workflow('a_ 1 + a_2 + a_4 + b_3-')
        wf.run()
        self.assertEqual(env.locals['executed'], ['a_1', 'a_2', 'a_4', 'b_3', 'b_4'])
        #
        wf = script.workflow('a+c+d')
        wf.run()
        self.assertEqual(env.locals['executed'], ['a_1', 'a_2', 'a_3', 'a_4', 'c_0', 'd_0'])

    def testNestedWorkflow(self):
        '''Test the creation and execution of combined workfow'''
        script = SoS_Script('''
executed = []
inputs = []
[a_1]
inputs.append(_step.input)
executed.append(_step.name)
[a_2]
inputs.append(_step.input)
executed.append(_step.name)
[a_3]
inputs.append(_step.input)
executed.append(_step.name)
[a_4]
output: 'a.done'
inputs.append(_step.input)
executed.append(_step.name)
[b_1]
input: 'b.begin'
inputs.append(_step.input)
executed.append(_step.name)
[b_2]
inputs.append(_step.input)
executed.append(_step.name)
[b_3]
inputs.append(_step.input)
executed.append(_step.name)
[b_4]
output: 'b.txt'
inputs.append(_step.input)
executed.append(_step.name)
[c=a+b]
input: 'a.txt'
output: 'b.txt'
inputs.append(_step.input)
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('c')
        wf.run()
        self.assertEqual(env.locals['executed'], ['c_0', 'a_1', 'a_2', 'a_3', 'a_4', 'b_1', 'b_2', 'b_3', 'b_4'])
        self.assertEqual(env.locals['inputs'], [['a.txt'], ['a.txt'], [], [], [], ['b.begin'], [], [], []])
        # step will be looped
        script = SoS_Script('''
executed = []
inputs = []
[a_1]
output: _input[0] + '.a1'
inputs.append(_step.input)
executed.append(_step.name)
[a_2]
output: _input[0] + '.a2'
inputs.append(_step.input)
executed.append(_step.name)
[c=a]
input: 'a.txt', 'b.txt', group_by='single'
inputs.append(_input)
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('c')
        wf.run()
        self.assertEqual(env.locals['executed'], ['c_0', 'a_1', 'a_2', 'c_0', 'a_1', 'a_2'])
        self.assertEqual(env.locals['inputs'], [['a.txt'], ['a.txt'], ['a.txt.a1'], ['b.txt'], ['b.txt'], ['b.txt.a1']])
        #
        # allow specifying a single step
        # step will be looped
        script = SoS_Script('''
executed = []
[a_1]
executed.append(_step.name)
[a_2]
executed.append(_step.name)
[c_0]
executed.append(_step.name)
[c_1=a_2]
input: 'a.txt', 'b.txt', group_by='single'
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('c')
        wf.run()
        self.assertEqual(env.locals['executed'], ['c_0', 'c_1', 'a_2', 'c_1', 'a_2'])
        # allow specifying a single step
        # step will be looped
        script = SoS_Script('''
executed = []
[a_1]
executed.append(_step.name)
[a_2]
executed.append(_step.name)
[c_0]
executed.append(_step.name)
[c_1=a_2]
input: 'a.txt', 'b.txt', group_by='single'
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('c')
        wf.run()
        self.assertEqual(env.locals['executed'], ['c_0', 'c_1', 'a_2', 'c_1', 'a_2'])
        #
        # recursive subworkflow not allowed
        script = SoS_Script('''
executed = []
[a_1]
executed.append(_step.name)
[a_2]
executed.append(_step.name)
[c_0]
executed.append(_step.name)
[c_1=a_2+c]
input: 'a.txt', 'b.txt', group_by='single'
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        self.assertRaises(RuntimeError, script.workflow, 'c')
        #
        # nested subworkflow is allowed
        script = SoS_Script('''
executed = []
[a_1]
executed.append(_step.name)
[a_2]
executed.append(_step.name)
[a_3]
executed.append(_step.name)
[b_1]
executed.append(_step.name)
[b_2=a_1+a_2]
executed.append(_step.name)
[c_0]
executed.append(_step.name)
[c_1=a+b]
input: 'a.txt'
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('c')
        wf.run()
        self.assertEqual(env.locals['executed'], ['c_0', 'c_1', 'a_1', 'a_2', 'a_3',
            'b_1', 'b_2', 'a_1', 'a_2'])
        #
        #
        # nested subworkflow with step option and others
        script = SoS_Script('''
executed = []
[a_1]
executed.append(_step.name)
[a_2]
executed.append(_step.name)
[a_3]
executed.append(_step.name)
[b=a_3+a_1, d=a_2, e2_2]
input: 'a.txt', 'b.txt', group_by='single'
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('b')
        wf.run()
        self.assertEqual(env.locals['executed'], ['b_0', 'a_3', 'a_1', 'b_0', 'a_3', 'a_1'])
        wf = script.workflow('d')
        wf.run()
        self.assertEqual(env.locals['executed'], ['d_0', 'a_2', 'd_0', 'a_2'])
        wf = script.workflow('e2')
        wf.run()
        self.assertEqual(env.locals['executed'], ['e2_2', 'e2_2'])

    def testSourceOption(self):
        '''Test the source section option'''
        # nested subworkflow with step option and others
        with open('temp/test.sos', 'w') as sos:
            sos.write('''
# test sos script

# global definition
GLB = 5

[parameters]
parB = 10

[A_1]
output: _input + 'a1'

[A_2]
output: _input + 'a2'

''')
        script = SoS_Script('''
executed = []
[b_1=A : source=['temp/test.sos', 'temp/test.sos'], skip=False]
input: 'a.txt', 'b.txt', group_by='single'
executed.append(_step.name)
''')
        env.run_mode = 'dryrun'
        wf = script.workflow('b')
        print(wf)
        wf.run()
        self.assertEqual(env.locals['GLB'], 5)
        self.assertEqual(env.locals['parB'], 5)
        self.assertEqual(env.locals['executed'], [])
        #
        os.remove('temp/test.sos')


if __name__ == '__main__':
    unittest.main()
