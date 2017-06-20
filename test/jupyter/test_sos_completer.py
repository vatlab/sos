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

#
# NOTE: for some namespace reason, this test can only be tested using
# nose.
#
# % nosetests test_kernel.py
#
#
import os
import unittest
from ipykernel.tests.utils import execute, wait_for_idle
from sos.jupyter.test_utils import sos_kernel, flush_channels

def get_reply(kc, text):
    flush_channels()
    kc.complete(text, len(text))
    reply = kc.get_shell_msg(timeout=2)
    return reply['content']

class TestKernel(unittest.TestCase):
    def testCompleter(self):
        with sos_kernel() as kc:
            # match magics
            self.assertTrue('%get ' in get_reply(kc, '%g')['matches'])
            self.assertTrue('%get ' in get_reply(kc, '%')['matches'])
            self.assertTrue('%with ' in get_reply(kc, '%w')['matches'])
            # path complete
            for m in get_reply(kc, '!ls ')['matches'][:5]:
                self.assertTrue(os.path.exists(m))
            #
            wait_for_idle(kc)
            # variable complete
            execute(kc=kc, code='alpha=5')
            wait_for_idle(kc)
            self.assertTrue('alpha' in get_reply(kc, 'al')['matches'])
            self.assertTrue('all(' in get_reply(kc, 'al')['matches'])

if __name__ == '__main__':
    unittest.main()
