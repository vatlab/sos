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

class TestRemove(unittest.TestCase):
    def setUp(self):
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        os.chdir('temp')
        with open('test.sos', 'w') as script:
            script.write('''
[0]
output:  't_f1'
sh:
    touch t_f1

[1]
output:  't_d1/t_f2'
sh:
    touch t_d1/t_f2
    touch t_d1/ut_f4

[2]
output:  't_d2/t_d3/t_f3'
sh:
    touch t_d2/t_d3/t_f3

''')
        subprocess.call('sos run test', shell=True)
        # create some other files and directory
        for d in ('ut_d1', 'ut_d2', 'ut_d2/ut_d3'):
            os.mkdir(d)
        for f in ('ut_f1', 'ut_d1/ut_f2', 'ut_d2/ut_d3/ut_f3'):
            with open(f, 'w') as tf:
                tf.write(f)

    def assertExists(self, fdlist):
        for fd in fdlist:
            self.assertTrue(os.path.exists(fd), '{} does not exist'.format(fd))

    def assertNonExists(self, fdlist):
        for fd in fdlist:
            self.assertFalse(os.path.exists(fd), '{} still exists'.format(fd))

    def testSetup(self):
        self.assertExists(['ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_f1', 'ut_d1/ut_f2', 'ut_d2/ut_d3/ut_f3'])
        self.assertExists(['t_f1', 't_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2'])
        # this is the tricky part, directory containing untracked file should remain
        self.assertExists(['t_d1', 't_d1/ut_f4'])

    def testRemoveAllTracked(self):
        '''test list files'''
        subprocess.call('sos remove . -t -y', shell=True)
        self.assertExists(['ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_f1', 'ut_d1/ut_f2', 'ut_d2/ut_d3/ut_f3'])
        self.assertNonExists(['t_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2'])
        # this is the tricky part, directory containing untracked file should remain
        self.assertExists(['t_d1', 't_f1',  't_d1/ut_f4'])
        
    def testRemoveSpecificTracked(self):
        # note the t_f1, which is under current directory and has to be remove specifically.
        subprocess.call('sos remove t_f1 ut_f1 t_d2 ut_d2 -t -y', shell=True)
        self.assertExists(['ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_f1', 'ut_d1/ut_f2', 'ut_d2/ut_d3/ut_f3', 't_d1/t_f2', 't_d1', 't_d1/ut_f4'])
        self.assertNonExists(['t_f1', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2'])

    def testRemoveAllUntracked(self):
        '''test list files'''
        subprocess.call('sos remove . -u -y', shell=True)
        self.assertNonExists(['ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_d1/ut_f2', 't_d1/ut_f4', 'ut_d2/ut_d3/ut_f3'])
        self.assertExists(['t_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2', 't_d1', 't_f1'])
        # this is the tricky part, files under the current directory are not removed
        self.assertExists(['ut_f1'])

    def testRemoveSpecificUntracked(self):
        # note the t_f1, which is under current directory and has to be remove specifically.
        subprocess.call('sos remove t_f1 ut_f1 ut_d1/ut_f2 t_d1 -u -y', shell=True)
        self.assertNonExists(['ut_f1', 'ut_d1/ut_f2', 't_d1/ut_f4'])
        self.assertExists(['t_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2', 't_d1', 't_f1'])
        self.assertExists(['ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_d2/ut_d3/ut_f3'])

    def testRemoveAll(self):
        subprocess.call('sos remove ut_d1 t_d1 ut_d2/ut_d3 -y', shell=True)
        self.assertNonExists(['ut_d1', 't_d1', 'ut_d2/ut_d3'])
        self.assertExists(['t_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2', 't_f1'])

    def tearDown(self):
        os.chdir('..')
        shutil.rmtree('temp')

if __name__ == '__main__':
    unittest.main()

