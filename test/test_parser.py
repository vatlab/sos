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
import subprocess
import shutil

from sos.utils import env, ArgumentError
from sos.sos_script import SoS_Script, ParsingError
from sos.sos_executor import Base_Executor, ExecuteError
from sos.target import FileTarget

class TestParser(unittest.TestCase):
    def setUp(self):
        env.reset()
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

    def testSetSigil(self):
        '''Test %set_options sigil'''
        script = SoS_Script('''%set_options sigil='[ ]' ''')
        self.assertEqual(script.global_sigil, '[ ]')

    def testMixedTabAndSpace(self):
        '''Test handling of mixed tab and space'''
        script = SoS_Script('''
[1: shared=['a', 'b', 'c']]
if True:
    a = 1
\tb = 2
\tc= 3
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['a'], 1)
        self.assertEqual(env.sos_dict['b'], 2)
        self.assertEqual(env.sos_dict['c'], 3)

    def testWorkflows(self):
        '''Test workflows defined in SoS script'''
        script = SoS_Script('''[0]''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[1]''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[*_1]''')
        self.assertEqual(sorted(script.workflows), ['default'])
        script = SoS_Script('''[0]\n[*_1]\n[auxiliary:provides='{a}.txt']''')
        self.assertEqual(sorted(script.workflows), ['auxiliary', 'default'])
        script = SoS_Script('''[0]\n[*_1]\n[human_1]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human'])
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])
        # skip option is not effective at parsing time
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse_2:skip]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])
        # unnamed
        script = SoS_Script('''[0]\n[*_1]\n[human_1]\n[mouse]\n[s*_2]''')
        self.assertEqual(sorted(script.workflows), ['default', 'human', 'mouse'])
        #
        # workflow name with -
        script = SoS_Script('''[proc-1]\n[test-case_2]''')
        self.assertEqual(sorted(script.workflows), ['proc-1', 'test-case'])
        script.workflow('proc-1')
        script.workflow('proc-1 + test-case:2')

    def testSkipStep(self):
        '''Test the skip option to skip certain steps'''
        script = SoS_Script('''
parameter: skip = 0

[0: shared={'a':'var'}, skip=skip==0]
var = 0

[1: shared={'b': 'var'}, skip=skip==1]
var = 1

''')
        wf = script.workflow()
        Base_Executor(wf, args=['--skip', '0']).run()
        self.assertEqual(env.sos_dict['b'], 1)
        #
        Base_Executor(wf, args=['--skip', '1']).run()
        self.assertEqual(env.sos_dict['a'], 0)

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
            self.assertRaises((ValueError, ParsingError), SoS_Script, '[0:{}]'.format(badoption))
        # good options
        for goodoption in ['sigil="[ ]"', 'shared="a"']:
            SoS_Script('[0:{}]'.format(goodoption))
        # allowed names
        for name in ['a5', 'a_5', '*_0', 'a*1_100']:
            SoS_Script('[{}]'.format(name))
        # allowed names with alias
        for name in ['a5 (p1)', 'a_5 (something fun)', '*_0 (no way)', 'a*1_100']:
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
        #
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
        #delattr(env, 'sos_dict')
        #script = SoS_Script('''a = b\n[0] ''')
        #wf = script.workflow()
        #dag = Base_Executor(wf).dryrun()
        #self.assertRaises(ValueError, Base_Executor(wf).run, dag)
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
        SoS_Script(filename='scripts/section1.sos')
        # not the default value of 1.0

    def testParameters(self):
        '''Test parameters section'''
        # directive not allowed in parameters
        script = SoS_Script(filename='scripts/section1.sos')
        wf = script.workflow('chapter:0')
        #self.assertRaises(ArgumentError, Base_Executor(wf).run,
        #    args=['--not_exist'])
        #self.assertRaises(ArgumentError, Base_Executor(wf).run,
        #    args=['--par1', 'a', 'b'])
        #
        script = SoS_Script('''
parameter: a = [1, 2]
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['a'], [1,2])
        wf = script.workflow()
        Base_Executor(wf, args=['--a', '3']).run()
        self.assertEqual(env.sos_dict['a'], [3])
        wf = script.workflow()
        Base_Executor(wf, args=['--a', '3', '5']).run()
        self.assertEqual(env.sos_dict['a'], [3, 5])
        #
        script = SoS_Script('''
# comment
# comment
parameter: a = ['a.txt', 'b.txt']
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['a'], ['a.txt', 'b.txt'])
        wf = script.workflow()
        Base_Executor(wf, args=['--a', '3']).run()
        self.assertEqual(env.sos_dict['a'], ['3'])
        wf = script.workflow()
        Base_Executor(wf, args=['--a', '3', '5']).run()
        self.assertEqual(env.sos_dict['a'], ['3', '5'])
        #
        # test parameter using global definition
        script = SoS_Script('''
a="100"

# comment
parameter: b=str(int(a)+1)
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], '101')
        #
        env.sos_dict.clear()
        script = SoS_Script('''
a=100

parameter: b=a+1
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], 101)
        #
        script = SoS_Script('''
a=100

parameter: b=a+1.
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '1000']).run()
        #
        self.assertEqual(env.sos_dict['b'], 1000)
        #
        # test string interpolation of the parameter section
        script = SoS_Script('''
a=100

# comment

parameter: b="${a+1}"
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], '101')
        #
        # argument has hve a value
        self.assertRaises(ParsingError, SoS_Script, '''
[0]
parameter: b=

''')
        # if it is a type, must provide value
        script = SoS_Script('''
# comment
parameter: b = int
[0]
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Base_Executor(wf).dryrun)
        #
        script = SoS_Script('''
parameter: b = list
[0]
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Base_Executor(wf).dryrun)
        # also require the type
        script = SoS_Script('''
parameter: b = int
[0]
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Base_Executor(wf, args=['--b', 'file']).dryrun)
        #
        script = SoS_Script('''
parameter: b = int
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).dryrun()
        self.assertEqual(env.sos_dict['b'], 5)
        # string
        script = SoS_Script('''
parameter: b = str
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).dryrun()
        self.assertEqual(env.sos_dict['b'], '5')
        # list is ok
        script = SoS_Script('''
parameter: b = list
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).dryrun()
        self.assertEqual(env.sos_dict['b'], ['5'])
        # bool required
        script = SoS_Script('''
# comment
parameter: b = bool
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b']).dryrun()
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--no-b']).dryrun()
        self.assertEqual(env.sos_dict['b'], False)
        # bool with default True
        script = SoS_Script('''
parameter: b = True
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=[]).dryrun()
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--b']).dryrun()
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--no-b']).dryrun()
        self.assertEqual(env.sos_dict['b'], False)
        # bool with default False
        script = SoS_Script('''
parameter: b = False
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=[]).dryrun()
        self.assertEqual(env.sos_dict['b'], False)
        Base_Executor(wf, args=['--b']).dryrun()
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--no-b']).dryrun()
        self.assertEqual(env.sos_dict['b'], False)
        #
        # parameters cannot coincide with a readonly global variable
        # are masked by previous definition
        script = SoS_Script('''
a = 4
parameter: a = 5
[0]
''')
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf, args=['--a', 7]).dryrun)
        #self.assertEqual(env.sos_dict['a'], 4)
        #
        # test parameters with dash
        script = SoS_Script('''
parameter: a_b = 5
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a_b', '10']).dryrun()
        self.assertEqual(env.sos_dict['a_b'], 10)
        Base_Executor(wf, args=['--a-b', '10']).dryrun()
        self.assertEqual(env.sos_dict['a_b'], 10)
        #
        #
        script = SoS_Script('''
parameter: a_b = int
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a_b', '10']).dryrun()
        self.assertEqual(env.sos_dict['a_b'], 10)
        Base_Executor(wf, args=['--a-b', '10']).dryrun()
        self.assertEqual(env.sos_dict['a_b'], 10)
        #
        # parameter cannot be any keyword
        for key in ['input', 'output', '_input', 'with']:
            self.assertRaises(Exception, SoS_Script, '''

parameter: {} = int
[0]
'''.format(key))


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
    from sos.runtime import logger
    logger.warning('I am from a dented text')
    if 1:
        logger.warning('Another warning')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
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
        # scripts with section head-like lines
        script = SoS_Script('''
[0]
R:
some::function(param)
''')
        wf = script.workflow()
        #
        # script with first-line indent
        #
        script = SoS_Script('''
[0]
sh:
  echo "a"

sh('echo "b"')
''')
        wf = script.workflow()


    def testInput(self):
        '''Test input directive'''
        self.touch(['a.txt', 'b.txt', 'a.pdf', 'a0', 'a1'])
        script = SoS_Script('''
[0]
files = ['a.txt', 'b.txt']

input: 'a.pdf', files

''')
        wf = script.workflow('default')
        Base_Executor(wf).dryrun()
        #
        # test input types
        script = SoS_Script('''
[0:shared={'i':'input', 'o':'output'}]
files = ("a${i}" for i in range(2))
input: {'a.txt', 'b.txt'}, files
output: ("a${x}" for x in _input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(sorted(env.sos_dict['i']), ['a.txt', 'a0', 'a1', 'b.txt'])
        self.assertEqual(sorted(env.sos_dict['o']), ['aa.txt', 'aa0', 'aa1', 'ab.txt'])


    def testGroupBy(self):
        '''Test group_by parameter of step input'''
        # group_by = 'all'
        self.touch(['a{}.txt'.format(x) for x in range(11)])
        #
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='all'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a2.txt', 'a3.txt', 'a4.txt']])
        # group_by = 'single'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='single'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt'], ['a2.txt'], ['a3.txt'], ['a4.txt']])
        # group_by = 'pairs'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='pairs'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a3.txt'], ['a2.txt', 'a4.txt']])
        # group_by = 'pairwise'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='pairwise'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a2.txt'], ['a2.txt', 'a3.txt'], ['a3.txt', 'a4.txt']])
        # group_by = 'combinations'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='combinations'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],  [['a1.txt', 'a2.txt'], ['a1.txt', 'a3.txt'],
            ['a1.txt', 'a4.txt'], ['a2.txt', 'a3.txt'], ['a2.txt', 'a4.txt'], ['a3.txt', 'a4.txt']])
        # group_by chunks specified as integers
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=3

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],
                         [['a1.txt', 'a2.txt', 'a3.txt'],
                         ['a4.txt', 'a5.txt', 'a6.txt'],
                         ['a7.txt', 'a8.txt', 'a9.txt']])
        # group_by chunks specified as integer strings
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='3'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'],
                         [['a1.txt', 'a2.txt', 'a3.txt'],
                         ['a4.txt', 'a5.txt', 'a6.txt'],
                         ['a7.txt', 'a8.txt', 'a9.txt']])
        # number of files should be divisible by group_by
        self.touch(['a{}.txt'.format(x) for x in range(1, 10)])
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=4

executed.append(_input)

''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).dryrun)
        # incorrect value causes an exception
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='something'

executed.append(_input)

''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).dryrun)

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

    def testLongerCode(self):
        '''Test definition of classes (with intermediate newlines) in step.'''
        script = SoS_Script('''# first block

[0: shared='b']
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
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], 0)

    def testCombinedWorkflow(self):
        '''Test the creation and execution of combined workfow'''
        script = SoS_Script('''
a0 = 0
if 'executed' not in locals():
    executed = []
parameter: a = a0 + 1
[a_1: shared='executed']
executed.append(step_name)
[a_2: shared='executed']
executed.append(step_name)
[a_3: shared='executed']
executed.append(step_name)
[a_4: shared='executed']
executed.append(step_name)
output: 'out_a_4'
[b_1: shared=['executed', 'input_b1']]
executed.append(step_name)
input_b1 = input
[b_2: shared='executed']
executed.append(step_name)
[b_3: shared='executed']
executed.append(step_name)
[b_4: shared='executed']
executed.append(step_name)
[c: shared='executed']
executed.append(step_name)
[d: shared='executed']
executed.append(step_name)
''')
        wf = script.workflow('a+b')
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'], ['a_1', 'a_2', 'a_3', 'a_4', 'b_1', 'b_2', 'b_3', 'b_4'])
        self.assertEqual(env.sos_dict['a'], 1)
        self.assertEqual(env.sos_dict['input_b1'], ['out_a_4'])
        #
        wf = script.workflow('a: 1-2 + a:4 + b:3-')
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'], ['a_1', 'a_2', 'a_4',
            'b_3', 'b_4'])
        #
        wf = script.workflow('a+c+d')
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['executed'], ['a_1', 'a_2', 'a_3', 'a_4', 'c_0', 'd_0'])

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
inputs.append(input)
[a_2: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(input)
[a_3: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(input)
[a_4: shared=['executed', 'inputs']]
executed.append(step_name)
output: 'a.done'
inputs.append(input)
sh:
    touch ${output}
[b_1: shared=['executed', 'inputs']]
executed.append(step_name)
input: 'b.begin'
inputs.append(input)
[b_2: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(input)
[b_3: shared=['executed', 'inputs']]
executed.append(step_name)
inputs.append(input)
[b_4: shared=['executed', 'inputs']]
executed.append(step_name)
output: 'b.txt'
inputs.append(input)
[c: shared=['executed', 'inputs']]
executed.append(step_name)
input: 'a.txt'
output: 'b.txt'
inputs.append(input)
sos_run('a+b')
''')
        env.sig_mode = 'ignore'
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
inputs.append(input)
sh:
    touch ${output}
[a_2:shared=['executed', 'inputs']]
executed.append(step_name)
output: _input[0] + '.a2'
inputs.append(input)
sh:
    touch ${output}
[c:shared=['executed', 'inputs']]
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
inputs.append(_input)
sos_run('a')
''')
        wf = script.workflow('c')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['executed'], ['c_0', 'a_1', 'a_2', 'a_1', 'a_2'])
        #self.assertEqual(env.sos_dict['inputs'], [['a.txt'], ['a.txt'], ['a.txt.a1'], ['b.txt'], ['b.txt'], ['b.txt.a1']])
        for file in ('a.txt.a1', 'a.txt.a1.a2', 'b.txt.a1', 'b.txt.a1.a2'):
            FileTarget(file).remove('both')
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
sos_run('a:2')
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
sos_run('a:2')
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
sos_run('a_2+c')
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
sos_run('a:1-2')
[c_0:shared='executed']
executed.append(step_name)
[c_1:shared='executed']
executed.append(step_name)
input: 'a.txt'
sos_run('a+b')
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
sos_run('a:3+a:1')
[d:shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('a:2')
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
        FileTarget('a.done').remove('both')

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
sos_run(wf)
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
        shutil.rmtree('.sos')
        os.makedirs('.sos/.runtime')
        with open('inc.sos', 'w') as sos:
            sos.write('''
# test sos script

# global definition
GLB = 5
parameter: parB = 10

[A_1: shared='executed']
executed.append('t.' + step_name)
output: _input[0] + '.a1'
sh:
    touch ${output}

[A_2: shared='executed']
executed.append('t.' + step_name)
output: _input[0] + '.a2'
sh:
    touch ${output}
''')
        script = SoS_Script('''
%from inc include *

if 'executed' not in locals():
    executed = []

[b_1: skip=False, shared='executed']
executed.append(step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('A')
''')
        wf = script.workflow('b')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['GLB'], 5)
        self.assertEqual(env.sos_dict['parB'], 10)
        self.assertEqual(env.sos_dict['executed'], ['b_1', 't.A_1', 't.A_2', 't.A_1', 't.A_2'])
        #
        shutil.rmtree('.sos')
        os.makedirs('.sos/.runtime')
        for file in ('a.txt.a1', 'a.txt.a1.a2', 'b.txt.a1', 'b.txt.a1.a2'):
            FileTarget(file).remove('both')
        script = SoS_Script('''
%include inc as k

if 'executed' not in locals():
    executed = []

[b_1: skip=False, shared='executed']
executed.append('g.' + step_name)
input: 'a.txt', 'b.txt', group_by='single'
sos_run('k.A')
''')
        wf = script.workflow('b')
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['k'].GLB, 5)
        self.assertEqual(env.sos_dict['k'].parB, 10)
        self.assertEqual(env.sos_dict['executed'], ['g.b_1', 't.k.A_1', 't.k.A_2', 't.k.A_1', 't.k.A_2'])
        #
        os.remove('inc.sos')

    def testIncludeWithNamespace(self):
        '''Test include a workflow that uses variables from its own global module'''
        self.touch(['a.txt', 'b.txt'])
        #
        with open('inc.sos', 'w') as sos:
            sos.write('''
# test sos script

# global definition
parameter: parB = 10

[A_1]
a = parB + 1

''')
        script = SoS_Script('''
%include inc

''')
        wf = script.workflow('inc.A')
        Base_Executor(wf).dryrun()


    def testYAMLConfig(self):
        '''Test config file in yaml format'''
        with open('myconfig.yml', 'w') as config:
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
        self.assertEqual(subprocess.call('sos run config.sos -c myconfig.yml', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 0)
        # now test the value
        script = SoS_Script(filename='config.sos')
        wf = script.workflow()
        Base_Executor(wf, config={'config_file': 'myconfig.yml'}).run()
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
        self.assertEqual(subprocess.call('sos run config.sos -c myconfig.yml', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
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
        self.assertEqual(subprocess.call('sos run config.sos -c myconfig.yml', stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True), 1)
        #
        for filename in ['config.sos', 'myconfig.yml']:
            os.remove(filename)

    def testVarOutput(self):
        '''Test early appearance of variable output'''
        script = SoS_Script('''
[0]
seq = range(3)
input: for_each='seq'
output: "test${_seq}.txt"
print(output)
''')
        wf = script.workflow()
        # this does not work before until we make variable output available sooner
        Base_Executor(wf).dryrun()
        # however, output should be the same in task
        script = SoS_Script('''
[0]
seq = range(3)
input: for_each='seq'
output: "test${_seq}.txt"
assert(len(output), _index + 1)
task:
assert(len(output), 3)
''')
        wf = script.workflow()
        # this does not work before until we make variable output available sooner
        Base_Executor(wf).dryrun()

    def testInclude(self):
        '''Test include keyword'''
        with open('inc.sos', 'w') as ts:
            ts.write('''
# a slave script to be included
gv = 1
[A_1]
[A_2]
[B]
''')
        script = SoS_Script('''
%include inc
res = inc.gv
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 1)
        #
        # include with alias
        script = SoS_Script('''
%include inc as tt
res1 = tt.gv
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res1'], 1)
        os.remove('inc.sos')

    def testFromInclude(self):
        '''Test include keyword'''
        with open('inc.sos', 'w') as ts:
            ts.write('''
# a slave script to be included
gv = 1
[A_1]
[A_2]
[B]
''')
        script = SoS_Script('''
%from inc include gv
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['gv'], 1)
        #
        # include with alias
        script = SoS_Script('''
%from inc include gv as g
res1 = g
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res1'], 1)

    def testCell(self):
        '''Test ignoring %cell'''
        SoS_Script('''
%cell 1
[step ]
a = 1
''')

    def testIfElse(self):
        '''Test if/elif/else/endif structural directive'''
        # no matching %endif
        self.assertRaises(ParsingError, SoS_Script, '''
%if 1
a = 1
%else
a=2
''')
        # no if for else
        self.assertRaises(ParsingError, SoS_Script, '''
%else
a=2
''')
        # no conditon for if
        self.assertRaises(ParsingError, SoS_Script, '''
%if
a=2
%endif
''')
        # no conditon for elif
        self.assertRaises(ParsingError, SoS_Script, '''
%if 1
%elif
a=2
%endif
[0]
''')
        # test if else
        script = SoS_Script('''
%if 0
a = 1
%else
a = 2
%endif
[0]
''')
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertEqual(env.sos_dict['a'], 2)

    def testOverwriteKeyword(self):
        '''Test overwrite sos keyword with user defined one.'''
        FileTarget('a.txt').remove('both')
        #
        script = SoS_Script('''
def run(script):
    pass

[1]
run:
    touch a.txt
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertFalse(os.path.isfile('a.txt'))
        #
        script = SoS_Script('''
parameter: run = 5

[1]
run:
    touch a.txt
''')
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)

if __name__ == '__main__':
    #suite = unittest.defaultTestLoader.loadTestsFromTestCase(TestParser)
    #unittest.TextTestRunner(, suite).run()
    unittest.main()
