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

def inspect(kc, name, pos=0):
    flush_channels()
    kc.inspect(name, pos)
    reply = kc.get_shell_msg(timeout=2)
    return reply['content']

class TestKernel(unittest.TestCase):
    def testCompleter(self):
        with sos_kernel() as kc:
            # match magics
            ins_print = inspect(kc, 'print')['data']['text/plain']
            self.assertTrue('built-in function' in ins_print,
                    'Returned: {}'.format(ins_print))
            wait_for_idle(kc)
            #
            execute(kc=kc, code='alpha=5')
            wait_for_idle(kc)
            execute(kc=kc, code='%use Python3')
            wait_for_idle(kc)
            
            ins_alpha = inspect(kc, 'alpha')['data']['text/plain']
            self.assertTrue('5' in ins_alpha, 'Returned: {}'.format(ins_alpha))
            wait_for_idle(kc)
            ins_get = inspect(kc, '%get', 2)['data']['text/plain']
            self.assertTrue('usage: %get' in ins_get, 'Returned: {}'.format(ins_get))
            wait_for_idle(kc)
            execute(kc=kc, code='%use SoS')
            wait_for_idle(kc)

if __name__ == '__main__':
    unittest.main()
