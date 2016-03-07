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

from sostestcase import SoS_TestCase

from pysos import *


class TestParser(unittest.TestCase): #SoSTestCase):
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

    def testSections(self):
        '''Test section definitions'''
        # bad names
        for badname in ['56_1', '_a', 'a_', '1x']:
            self.assertRaises(ParsingError, SoS_Script, '[{}]'.format(badname))
        # bad options
        for badoption in ['ss', 'skip a', 'skip:_']:
            self.assertRaises(ParsingError, SoS_Script, '[0:{}]'.format(badoption))
        # allowed names
        for name in ['a5', 'a_5', '*_0', '*1_100']:
            SoS_Script('[{}]'.format(name))
        # glo
        self.assertRaises(ParsingError, SoS_Script,
            '''input: 'filename' ''')

        script = SoS_Script('scripts/section1.sos')
        self.assertTrue('section' in script.workflows.keys())
        self.assertTrue('chapter' in script.workflows.keys())

    def testGlobalVariables(self):
        '''Test definition of variables'''
        # global section cannot have directive
        self.assertRaises(ParsingError, SoS_Script,
            '''input: 'filename' ''')
        # or unrecognized directive
        self.assertRaises(ParsingError, SoS_Script,
            '''inputs: 'filename' ''')
        # or unrecoginzied varialbe
        self.assertRaises(ParsingError, SoS_Script,
            '''something ''')
        # or function call
        self.assertRaises(ParsingError, SoS_Script,
            '''somefunc() ''')
        # allow definition
        SoS_Script('''a = '1' ''')
        SoS_Script('''a = b''')
        # but this one has incorrect syntax
        self.assertRaises(ParsingError, SoS_Script,
            '''a = 'b  ''')
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
        self.assertRaises(ArgumentError, SoS_Script,
            'scripts/section1.sos', args=['--not_exist'])
        self.assertRaises(ArgumentError, SoS_Script,
            'scripts/section1.sos', args=['--par1', 'a', 'b'])
        script = SoS_Script('scripts/section1.sos', args=['--par1', 'var2'])
        # need to check if par1 is set to correct value
        self.assertEqual(script.parameter('par1'), "var2")


    def testSectionVariables(self):
        '''Test variables in sections'''
        pass

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

    def testSectionActions(self):
        '''Test actions of sections'''
        self.assertRaises(ParsingError, SoS_Script,
            '''func()''')
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

if __name__ == '__main__':
    unittest.main()
