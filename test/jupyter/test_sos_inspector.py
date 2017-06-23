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

def is_complete(kc, code):
    flush_channels()
    kc.is_complete(code)
    reply = kc.get_shell_msg(timeout=2)
    return reply['content']

class TestSoSCompleter(unittest.TestCase):
    def testCompleter(self):
        with sos_kernel() as kc:
            # match magics
            ins_print = inspect(kc, 'print')['data']['text/plain']
            self.assertTrue('built-in function' in ins_print,
                    'Returned: {}'.format(ins_print))
            wait_for_idle(kc)
            #
            # keywords
            ins_depends = inspect(kc, 'depends:')['data']['text/plain']
            self.assertTrue('dependent targets' in ins_depends,
                    'Returned: {}'.format(ins_depends))
            wait_for_idle(kc)
            #
            execute(kc=kc, code='alpha=5')
            wait_for_idle(kc)
            execute(kc=kc, code='%use Python3')
            wait_for_idle(kc)
            # action
            ins_run = inspect(kc, 'run:')['data']['text/plain']
            self.assertTrue('sos.actions' in ins_run,
                    'Returned: {}'.format(ins_run))
            wait_for_idle(kc)
            #
            ins_alpha = inspect(kc, 'alpha')['data']['text/plain']
            self.assertTrue('5' in ins_alpha, 'Returned: {}'.format(ins_alpha))
            wait_for_idle(kc)
            for magic in ('get', 'run', 'set', 'sosrun', 'rerun', 'toc'):
                ins_magic = inspect(kc, '%' + magic, 2)['data']['text/plain']
                self.assertTrue('usage: %' + magic in ins_magic, 'Returned: {}'.format(ins_magic))
            wait_for_idle(kc)
            execute(kc=kc, code='%use SoS')
            wait_for_idle(kc)

    def testIsComplete(self):
        with sos_kernel() as kc:
            # match magics
            status = is_complete(kc, 'prin')
            self.assertEqual(status['status'], 'incomplete')
            #
            status = is_complete(kc, 'a=1')
            self.assertEqual(status['status'], 'incomplete')
            #
            status = is_complete(kc, '')
            self.assertEqual(status['status'], 'complete')
            #
            status = is_complete(kc, 'input:\n a=1,')
            self.assertEqual(status['status'], 'incomplete')
            #
            status = is_complete(kc, '%dict -r')
            self.assertEqual(status['status'], 'complete')
            wait_for_idle(kc)

if __name__ == '__main__':
    unittest.main()
