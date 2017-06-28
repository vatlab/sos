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
import sys
import unittest
from ipykernel.tests.utils import assemble_output, execute
from sos.jupyter.test_utils import sos_kernel

class TestPreview(unittest.TestCase):
    def setUp(self):
        self.resource_dir = os.path.abspath(os.path.split(__file__)[0])

    @unittest.skipIf(sys.platform == 'win32', 'pysam is not installed under windows')
    def testMagicPreview(self):
        with sos_kernel() as kc:
            # preview bam file
            iopub = kc.iopub_channel
            execute(kc=kc, code='''
%preview -n {}/sim_reads_aligned.bam
'''.format(self.resource_dir))
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            self.assertTrue('PG' in stdout)

            # preview sam file
            execute(kc=kc, code='''
%preview -n {}/sim_reads_aligned.sam
'''.format(self.resource_dir))
            stdout, stderr = assemble_output(iopub)
            self.assertEqual(stderr, '')
            self.assertTrue('PG' in stdout)


if __name__ == '__main__':
    unittest.main()
