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

#from sos.sos_script import SoS_Script
from sos.utils import env
from sos.target import FileTarget
#from sos.rq.sos_executor import RQ_Executor

import unittest
import subprocess

class TestRQQueue(unittest.TestCase):
    def setUp(self):
        env.reset()
        self.temp_files = []
        subprocess.call('redis-server &', shell=True)
        subprocess.call('rq worker &', shell=True)

    def tearDown(self):
        #
        for f in self.temp_files:
            FileTarget(f).remove('both')
        subprocess.call('redis-cli shutdown', shell=True)
        # how to kill rq worker?

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

#     def testSharedVar(self):
#         '''Test shared var with rq queue'''
#         script = SoS_Script(
# '''
# [work_1: shared = {'data': 'output'}]
# input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
# output: expand_pattern('{_name}.out')
# task: concurrent = True
# run:
#   touch ${_output}
# 
# [work_2]
# input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}', paired_with = ['data']
# output: expand_pattern('{_name}.out2')
# task: concurrent = True
# run:
#   touch ${_data} ${_output}
# 
# [default]
# sos_run("work:1+work:2")
# ''')
#         self.touch(['1.txt', '2.txt'])
#         subprocess.call('sos remove . -t -y', shell=True)
#         wf = script.workflow()
#         RQ_Executor(wf).run()
#         for f in ['1.out', '1.out2', '2.out', '2.out2']:
#             self.assertTrue(FileTarget(f).exists('target'))
#             FileTarget(f).remove('both')


if __name__ == '__main__':
    unittest.main()
