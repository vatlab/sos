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
from ipykernel.tests.utils import execute, wait_for_idle
from sos.jupyter.test_utils import sos_kernel, flush_channels

def get_completions(kc, text):
    flush_channels()
    kc.complete(text, len(text))
    reply = kc.get_shell_msg(timeout=2)
    return reply['content']

class TestKernel(unittest.TestCase):
    def testCompleter(self):
        with sos_kernel() as kc:
            # match magics
            self.assertTrue('%get ' in get_completions(kc, '%g')['matches'])
            self.assertTrue('%get ' in get_completions(kc, '%')['matches'])
            self.assertTrue('%with ' in get_completions(kc, '%w')['matches'])
            # path complete
            self.assertGreater(len(get_completions(kc, '!ls ')['matches']), 0)
            self.assertEqual(len(get_completions(kc, '!ls SOMETHING')['matches']), 0)
            #
            wait_for_idle(kc)
            # variable complete
            execute(kc=kc, code='alpha=5')
            wait_for_idle(kc)
            execute(kc=kc, code='%use Python3')
            wait_for_idle(kc)
            self.assertTrue('alpha' in get_completions(kc, 'al')['matches'])
            self.assertTrue('all(' in get_completions(kc, 'al')['matches'])
            # for no match
            self.assertEqual(len(get_completions(kc, 'alphabetatheta')['matches']), 0)
            # get with all variables in
            self.assertTrue('alpha' in get_completions(kc, '%get ')['matches'])
            self.assertTrue('alpha' in get_completions(kc, '%get al')['matches'])
            # with use and restart has kernel name
            self.assertTrue('Python3' in get_completions(kc, '%with ')['matches'])
            self.assertTrue('Python3' in get_completions(kc, '%use ')['matches'])
            self.assertTrue('Python3' in get_completions(kc, '%shutdown ')['matches'])
            self.assertTrue('Python3' in get_completions(kc, '%shutdown ')['matches'])
            self.assertTrue('Python3' in get_completions(kc, '%use Py')['matches'])
            #
            self.assertEqual(len(get_completions(kc, '%use SOME')['matches']), 0)
            #
            wait_for_idle(kc)
            execute(kc=kc, code='%use SoS')
            wait_for_idle(kc)

if __name__ == '__main__':
    unittest.main()
