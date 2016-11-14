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

class TestPack(unittest.TestCase):
    def setUp(self):
        if os.path.isdir('temp'):
            shutil.rmtree('temp')
        os.mkdir('temp')
        os.chdir('temp')
        with open('test.sos', 'w') as script:
            script.write('''
%from included include *
parameter: name='t_f1'
[0]
output:  name
sh:
    dd if=/dev/urandom of=${name} count=10000

[1: sigil='[ ]']
output:  't_d1/t_f2'
task:
sh:
    dd if=/dev/urandom of=[output] count=50000
    dd if=/dev/urandom of=t_d1/ut_f4 count=500

[2]
output:  't_d2/t_d3/t_f3'
sh:
    dd if=/dev/urandom of=${output} count=6000

''')
        with open('included.sos', 'w') as script:
            script.write('''
# does nothing
a = 1
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

    def testDryrun(self):
        '''Test dryrun mode'''
        self.assertEqual(subprocess.call('sos pack -o b.sar -i t_d1/ut_f4 --dryrun', shell=True), 0)
        self.assertFalse(os.path.isfile('b.sar'))

    def testPackUnpack(self):
        '''Test pack command'''
        self.assertEqual(subprocess.call('sos pack -o a.sar', shell=True), 0)
        # extra file
        self.assertEqual(subprocess.call('sos pack -o b.sar -i t_d1/ut_f4', shell=True), 0)
        # extra directory
        self.assertEqual(subprocess.call('sos pack -o b.sar -i t_d1 -y', shell=True), 0)
        # unpack
        self.assertEqual(subprocess.call('sos unpack a.sar', shell=True), 0)
        # unpack to a different directory
        self.assertEqual(subprocess.call('sos unpack a.sar -d tmp', shell=True), 0)
        # list content
        self.assertEqual(subprocess.call('sos unpack a.sar -l', shell=True), 0)

    def testUnpackScript(self):
        '''Test -s option of unpack'''
        self.assertEqual(subprocess.call('sos pack -o a.sar', shell=True), 0)
        os.remove('test.sos')
        os.remove('included.sos')
        # unpack
        self.assertEqual(subprocess.call('sos unpack a.sar', shell=True), 0)
        self.assertFalse(os.path.isfile('test.sos'))
        self.assertFalse(os.path.isfile('included.sos'))
        # unpack to a different directory
        self.assertEqual(subprocess.call('sos unpack a.sar -s -y', shell=True), 0)
        self.assertTrue(os.path.isfile('test.sos'))
        self.assertTrue(os.path.isfile('included.sos'))

    def testUnpackSelected(self):
        # unpack selected file
        self.assertEqual(subprocess.call('sos pack -o a.sar -i t_d1/ut_f4', shell=True), 0)
        shutil.rmtree('.sos')
        shutil.rmtree('t_d1')
        shutil.rmtree('t_d2')
        self.assertEqual(subprocess.call('sos unpack a.sar ut_f4', shell=True), 0)
        self.assertTrue(os.path.isfile('t_d1/ut_f4'))
        self.assertFalse(os.path.exists('t_d2'))

    def tearDown(self):
        os.chdir('..')
        shutil.rmtree('temp')

if __name__ == '__main__':
    unittest.main()

