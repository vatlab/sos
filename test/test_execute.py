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

class TestRun(unittest.TestCase):
    def testSignature(self):
        '''Test recognizing the format of SoS script'''
        script = SoS_Script(r"""
[*_0]
output: 'temp/a.txt', 'temp/b.txt'

run('''echo "a.txt" > 'temp/a.txt' ''')
run('''echo "b.txt" > 'temp/b.txt' ''')

[1: output_alias=oa]
dest = ['temp/c.txt', 'temp/d.txt']
input: group_by='single', labels='dest'
output: dest

run(''' cp ${input} ${_dest} ''')
""")
        wf = script.workflow('default:0')
        wf.run()
        # not the default value of 1.0
        self.assertTrue(os.path.isfile('temp/a.txt'))
        self.assertTrue(os.path.isfile('temp/b.txt'))
        self.assertTrue(open('temp/a.txt').read(), 'a.txt')
        self.assertTrue(open('temp/b.txt').read(), 'b.txt')
        #
        wf = script.workflow('default')
        wf.run()
        # not the default value of 1.0
        self.assertTrue(os.path.isfile('temp/c.txt'))
        self.assertTrue(os.path.isfile('temp/d.txt'))
        self.assertTrue(open('temp/c.txt').read(), 'a.txt')
        self.assertTrue(open('temp/d.txt').read(), 'b.txt')
        self.assertEqual(wf.locals['oa'], ['temp/c.txt', 'temp/d.txt'])
        

    def testForEach(self):
        '''Test for_each option of input'''
        script = SoS_Script(r"""
[0]
files = ['a.txt', 'b.txt']
names = ['a', 'b', 'c']
c = ['1', '2']
counter = 0
all_names = ''
all_loop = ''

input: 'a.pdf', files, group_by='single', labels='names', for_each='c'

all_names += _names[0] + " "
all_loop += _c + " "

counter = counter + 1
""")
        wf = script.workflow('default')
        wf.run()
        self.assertEqual(wf.locals['counter'], 6)
        self.assertEqual(wf.locals['all_names'], "a b c a b c ")
        self.assertEqual(wf.locals['all_loop'], "1 1 1 2 2 2 ")

    def testAlias(self):
        script = SoS_Script(r"""
[0: input_alias=ia, output_alias=oa]
files = ['a.txt', 'b.txt']
names = ['a', 'b', 'c']
c = ['1', '2']
counter = "0"

input: 'a.pdf', files, group_by='single', labels='names', for_each='c'

counter = str(int(counter) + 1)
""")
        wf = script.workflow('default')
        wf.run()
        self.assertEqual(wf.locals['ia'], ["a.pdf", 'a.txt', 'b.txt'])

if __name__ == '__main__':
    unittest.main()
