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

import unittest
from IPython.terminal.embed import InteractiveShellEmbed


class TestIpython(unittest.TestCase):
    def setUp(self):
        self.ipshell = InteractiveShellEmbed()
        self.ipshell.run_cell('%load_ext sos_magic')
        
    def assertDictValue(self, key, value):
        self.ipshell.run_cell('__sos_dict__ = %sosdict')
        self.assertEqual(self.ipshell.user_ns['__sos_dict__'][key], value)
        
    def testSoSDict(self):
        '''Test sos dict magic'''
        self.ipshell.run_cell('%sosdict')
        self.ipshell.run_cell('keys = %sosdict --keys')
        # this make sure the symbols are imported
        for key in ['run', 'R', 'bash', 'python', 'sos_variable']:
            self.assertTrue(key not in self.ipshell.user_ns['keys'])
        #
        self.ipshell.run_cell('keys = %sosdict --keys --all')
        # this make sure the symbols are imported
        for key in ['run', 'R', 'bash', 'python', 'sos_variable']:
            self.assertTrue(key in self.ipshell.user_ns['keys'])
        #
        # add something
        self.ipshell.run_cell('%sos a=1')
        self.assertDictValue('a', 1)
        #
        # reset
        self.ipshell.run_cell('%sosdict --reset')
        self.ipshell.run_cell('__sos_dict__ = %sosdict')
        self.assertTrue('a' not in self.ipshell.user_ns['__sos_dict__'])

    def testSoS(self):
        '''Test magic %sos'''
        self.ipshell.run_cell('''%sos a=10''')
        self.ipshell.run_cell('''%sos b=['file1.txt', 'file2.txt']''')
        self.ipshell.run_cell('''%sos c="${b!r,}" ''')
        self.assertDictValue('a', 10)
        self.assertDictValue('b', ['file1.txt', 'file2.txt'])
        self.assertDictValue('c', "'file1.txt', 'file2.txt'")

    def testSet(self):
        '''test %sosset'''
        self.ipshell.run_cell('''%sosset --rep 3''')
        self.ipshell.run_cell('''%%sos
parameter: rep = 5
''')
        self.assertDictValue('rep', 3)

    def testSoSPut(self):
        '''test %put'''
        self.ipshell.run_cell('a = 12345')
        self.ipshell.run_cell('b = "12345"')
        self.ipshell.run_cell('%sosput a b')
        self.assertDictValue('a', 12345)
        self.assertDictValue('b', "12345")

    def testSoSGet(self):
        '''test %get'''
        self.ipshell.run_cell('%sos a = 12345')
        self.ipshell.run_cell('%sos b = "12345"')
        self.ipshell.run_cell('%sosget a b')
        self.assertEqual(self.ipshell.user_ns['a'], 12345)
        self.assertEqual(self.ipshell.user_ns['b'], "12345")

    def testSigil(self):
        '''test change sigil'''
        self.ipshell.run_cell('''%%sos
%set_options sigil='[ ]'
''')
        self.ipshell.run_cell('%sos a = 100')
        self.ipshell.run_cell('%sos b = "${a}"')
        self.ipshell.run_cell('%sosget b')
        self.assertEqual(self.ipshell.user_ns['b'], '${a}')
        #
        self.ipshell.run_cell('''%%sos
%set_options sigil=None
''')
        self.ipshell.run_cell('%sos b = "${a}"')
        self.ipshell.run_cell('%sosget b')
        self.assertEqual(self.ipshell.user_ns['b'], '${a}')
        #
        self.ipshell.run_cell('''%%sos
%set_options sigil='${ }'
''')

if __name__ == '__main__':
    unittest.main()
