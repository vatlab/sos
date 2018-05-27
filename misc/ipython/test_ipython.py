#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

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

if __name__ == '__main__':
    unittest.main()
