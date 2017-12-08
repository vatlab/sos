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
import subprocess

from sos.utils import env, ArgumentError
from sos.parser import SoS_Script, ParsingError
from sos.workflow_executor import Base_Executor
from sos.targets import file_target, sos_targets

section1_sos = '''
#!/usr/bin/env sos-runner
#fileformat=SOS1.1

#
# this is a test sos script
#
var1='value1'
var2 = 'value2'
var3 = [var1,
  var2]

# section
# description for workflow section
#
[parameters]
# par1 string
par1 = 'var1'

# par2 list
par2 = ['a', 'b', 'c']

# par3 multiline
par3 = ['a',
	'b']

[*_0]
var0 = '0'

[section_10]
#
#step 10
var1 = 'a'

[section_2 : shared='var3']
var2 = 'a'
input: var1
output: var2

var3 = 'a'

[section_3, *_4 : shared='var4', skip ]
output: 
	var2,
	var3

print()
var4= 'value4'

[chapter_5]
var4='5'
 '''

class TestParser(unittest.TestCase):
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

    def testFileFormat(self):
        '''Test recognizing the format of SoS script'''
        # file format must be 'fileformat=SOSx.x'
        self.assertRaises(ParsingError, SoS_Script,
            '#fileformat=SS2')
        self.assertRaises(ParsingError, SoS_Script,
            '#fileformat=SOS1.0beta')
        #
        # parse a larger script with gormat 1.1
        script = SoS_Script(section1_sos)
        # not the default value of 1.0
        #self.assertEqual(script.format_version, '1.1')

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

    def testTypeHint(self):
        '''We should be able to differentiate type hint and sos action
        '''
        script = SoS_Script('''a : list = 5''')
        script = SoS_Script('''a : list''')
        # action
        script = SoS_Script('''a : input='filename' ''')
        # action
        script = SoS_Script('''a : expand='${ }' ''')


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
        # allowed names
        for name in ['a5', 'a_5', '*_0', 'a*1_100']:
            SoS_Script('[{}]'.format(name))
        # allowed names with alias
        for name in ['a5 (p1)', 'a_5 (something fun)', '*_0 (no way)', 'a*1_100']:
            SoS_Script('[{}]'.format(name))
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
        #
        # global section
        self.assertRaises(ParsingError, SoS_Script,
        '''[global: skip]''')
        self.assertRaises(ParsingError, SoS_Script,
        '''[global, step_10]''')

    def testGlobalVariables(self):
        '''Test definition of variables'''
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
        #dag = Base_Executor(wf).run(mode='dryrun')
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
        SoS_Script(section1_sos)
        # not the default value of 1.0
        #
        script = SoS_Script('''
[global]
a = 1

[b]
print(a)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    def testParameters(self):
        '''Test parameters section'''
        # directive not allowed in parameters
        script = SoS_Script(section1_sos)
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
        self.assertEqual(list(wf.parameters().keys()), ['a'])
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
        self.assertEqual(list(wf.parameters().keys()), ['b'])
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
        self.assertRaises(ArgumentError, Base_Executor(wf).run, mode='dryrun')
        #
        script = SoS_Script('''
parameter: b = list
[0]
''')
        self.assertEqual(list(wf.parameters().keys()), ['b'])
        wf = script.workflow()
        self.assertRaises(ArgumentError, Base_Executor(wf).run, mode='dryrun')
        # also require the type
        script = SoS_Script('''
parameter: b = int
[0]
''')
        wf = script.workflow()
        self.assertRaises(ArgumentError, Base_Executor(wf, args=['--b', 'file']).run, mode='dryrun')
        #
        script = SoS_Script('''
parameter: b = int
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], 5)
        # string
        script = SoS_Script('''
parameter: b = str
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], '5')
        # list is ok
        script = SoS_Script('''
parameter: b = list
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], ['5'])
        # bool required
        script = SoS_Script('''
# comment
parameter: b = bool
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--no-b']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], False)
        # bool with default True
        script = SoS_Script('''
parameter: b = True
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=[]).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--b']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--no-b']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], False)
        # bool with default False
        script = SoS_Script('''
parameter: b = False
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=[]).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], False)
        Base_Executor(wf, args=['--b']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], True)
        Base_Executor(wf, args=['--no-b']).run(mode='dryrun')
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
        self.assertRaises(Exception, Base_Executor(wf, args=['--a', 7]).run, mode='dryrun')
        #self.assertEqual(env.sos_dict['a'], 4)
        #
        # test parameters with dash
        script = SoS_Script('''
parameter: a_b = 5
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a_b', '10']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'], 10)
        Base_Executor(wf, args=['--a-b', '10']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'], 10)
        #
        #
        script = SoS_Script('''
parameter: a_b = int
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a_b', '10']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'], 10)
        Base_Executor(wf, args=['--a-b', '10']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'], 10)
        #
        # parameter cannot be any keyword
        for key in ['input', 'output', '_input', 'with']:
            self.assertRaises(Exception, SoS_Script, '''

parameter: {} = int
[0]
'''.format(key))
        # multiple parameters
        script = SoS_Script('''
parameter: a_b = int
[0]
parameter: c_b = list
''')
        wf = script.workflow()
        self.assertEqual(sorted(list(wf.parameters().keys())), ['a_b', 'c_b'])

#    this test is no longer valid because sos has stopped parsing assignments
#
#    def testSectionVariables(self):
#        '''Test variables in sections'''
#        # directive name cannot be used as variable
#        self.assertRaises(ParsingError, SoS_Script,
#            '''[0]
#input='a.txt' ''')
#        self.assertRaises(ParsingError, SoS_Script,
#            '''[0]
#output='a.txt' ''')
#        self.assertRaises(ParsingError, SoS_Script,
#            '''[0]
#depends='a.txt' ''')
        
    def testTypeTraitParameter(self):
        # type trait
        script = SoS_Script('''
parameter: b
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], '5')
        #
        script = SoS_Script('''
parameter: b :str
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], '5')        
        #
        script = SoS_Script('''
parameter: b : list
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], ['5'])  

        #
        script = SoS_Script('''
parameter: b : list = 5
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], 5)
        #
        script = SoS_Script('''
parameter: b : int
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', '5']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'], 5)

    def testInputTarget(self):
        # test input of targets
        script = SoS_Script('''
parameter: b : file_target
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'file_target')
        #
        script = SoS_Script('''
parameter: b = file_target('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'file_target')
        #
        script = SoS_Script('''
parameter: b : sos_targets
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa', 'bbb']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'sos_targets')
        #
        script = SoS_Script('''
parameter: b = sos_targets('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'sos_targets')
        #
        script = SoS_Script('''
parameter: a_b : file_target
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'file_target')
        #
        script = SoS_Script('''
parameter: a_b = file_target('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'file_target')
        #
        script = SoS_Script('''
parameter: a_b : sos_targets
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa', 'bbb']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'sos_targets')
        #
        script = SoS_Script('''
parameter: a_b = sos_targets('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'sos_targets')
        #
        #
        #
        #
        script = SoS_Script('''
parameter: b : path
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'path')
        #
        script = SoS_Script('''
parameter: b = path('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'path')
        #
        script = SoS_Script('''
parameter: b : paths
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa', 'bbb']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'paths')
        #
        script = SoS_Script('''
parameter: b = paths('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['b'].__class__.__name__, 'paths')
        #
        script = SoS_Script('''
parameter: a_b : path
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'path')
        #
        script = SoS_Script('''
parameter: a_b = path('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'path')
        #
        script = SoS_Script('''
parameter: a_b : paths
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa', 'bbb']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'paths')
        #
        script = SoS_Script('''
parameter: a_b = paths('file')
[0]
''')
        wf = script.workflow()
        Base_Executor(wf, args=['--a-b', 'aaa']).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a_b'].__class__.__name__, 'paths')
        #

    def testSectionDirectives(self):
        '''Test directives of sections'''
        # cannot be in the global section
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
        Base_Executor(wf).run(mode='dryrun')
        #
        # test input types
        script = SoS_Script('''
[0:shared={'i':'_input', 'o':'_output'}]
files = (f"a{i}" for i in range(2))
input: {'a.txt', 'b.txt'}, files
output: (f"a{x}" for x in _input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
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
input: [f'a{x}.txt' for x in range(1, 5)], group_by='all'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'],  [sos_targets('a1.txt', 'a2.txt', 'a3.txt', 'a4.txt')])
        # group_by = 'single'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='single'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'],  [sos_targets('a1.txt'), sos_targets('a2.txt'), sos_targets('a3.txt'), sos_targets('a4.txt')])
        # group_by = 'pairs'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: [f'a{x}.txt' for x in range(1, 5)], group_by='pairs'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'],  [sos_targets('a1.txt', 'a3.txt'), sos_targets('a2.txt', 'a4.txt')])
        # group_by = 'pairwise'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='pairwise'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'],  [sos_targets('a1.txt', 'a2.txt'), sos_targets('a2.txt', 'a3.txt'), sos_targets('a3.txt', 'a4.txt')])
        # group_by = 'combinations'
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='combinations'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'],  [sos_targets('a1.txt', 'a2.txt'), sos_targets('a1.txt', 'a3.txt'),
            sos_targets('a1.txt', 'a4.txt'), sos_targets('a2.txt', 'a3.txt'), sos_targets('a2.txt', 'a4.txt'), sos_targets('a3.txt', 'a4.txt')])
        # group_by chunks specified as integers
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=3

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'], [
                        sos_targets('a1.txt', 'a2.txt', 'a3.txt'),
                        sos_targets('a4.txt', 'a5.txt', 'a6.txt'),
                        sos_targets('a7.txt', 'a8.txt', 'a9.txt')])
        # group_by chunks specified as integer strings
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='3'

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'], [
                         sos_targets('a1.txt', 'a2.txt', 'a3.txt'),
                         sos_targets('a4.txt', 'a5.txt', 'a6.txt'),
                         sos_targets('a7.txt', 'a8.txt', 'a9.txt')])
        # number of files should be divisible by group_by
        self.touch(['a{}.txt'.format(x) for x in range(1, 10)])
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=4

executed.append(_input)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode="dryrun")
        # incorrect value causes an exception
        script = SoS_Script('''
[0]

executed = []
input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='something'

executed.append(_input)

''')
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run, mode="dryrun")

    def testOutputGroupBy(self):
        '''Test group_by parameter of step output'''
        # group_by = 'all'
        self.touch(['a{}.txt'.format(x) for x in range(4)])
        #
        script = SoS_Script('''
[0: shared='executed']

executed = []
input: ['a{}.txt'.format(x) for x in range(4)], group_by=2
output: ['a{}.txt.bak'.format(x) for x in range(4)], group_by=2

executed.append(_output)

''')
        wf = script.workflow()
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'], [sos_targets('a0.txt.bak', 'a1.txt.bak'), sos_targets('a2.txt.bak', 'a3.txt.bak')])

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
input_b1 = _input
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'], ['a_1', 'a_2', 'a_3', 'a_4', 'b_1', 'b_2', 'b_3', 'b_4'])
        self.assertEqual(env.sos_dict['a'], 1)
        self.assertEqual(env.sos_dict['input_b1'].targets(), ['out_a_4'])
        #
        wf = script.workflow('a: 1-2 + a:4 + b:3-')
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'], ['a_1', 'a_2', 'a_4',
            'b_3', 'b_4'])
        #
        wf = script.workflow('a+c+d')
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['executed'], ['a_1', 'a_2', 'a_3', 'a_4', 'c_0', 'd_0'])

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
        Base_Executor(wf).run(mode='dryrun')


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
#print(CONFIG.StoreOwner)
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


    def testVarOutput(self):
        '''Test early appearance of variable output'''
        script = SoS_Script('''
[0]
seq = range(3)
input: for_each='seq'
output: "test${_seq}.txt"
print(_output)
''')
        wf = script.workflow()
        # this does not work before until we make variable output available sooner
        Base_Executor(wf).run(mode='dryrun')


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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run(mode='dryrun')
        self.assertEqual(env.sos_dict['a'], 2)

    def testOverwriteKeyword(self):
        '''Test overwrite sos keyword with user defined one.'''
        file_target('a.txt').remove('both')
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
