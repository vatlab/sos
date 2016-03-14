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

class TestActions(unittest.TestCase):
    def testCheckCommand(self):
        '''Test action check_command'''
        script = SoS_Script(r"""
[0]
check_command('cat')
""")
        wf = script.workflow()
        # should be ok
        wf.run()
        #
        script = SoS_Script(r"""
[0]
check_command('catmouse')
""")
        env.run_mode = 'dryrun'
        wf = script.workflow()
        # should fail in dryrun mode
        self.assertRaises(RuntimeError, wf.run)
        #
        env.run_mode = 'run'
        wf = script.workflow()
        # should fail also in run mode
        self.assertRaises(RuntimeError, wf.run)

    def testFailIf(self):
        '''Test action fail if'''
        script = SoS_Script(r"""
[0]
input: 'a.txt'
fail_if(len(_step.input) == 1)
""")
        env.run_mode = 'dryrun'
        wf = script.workflow()
        # should fail in dryrun mode
        self.assertRaises(RuntimeError, wf.run)
        script = SoS_Script(r"""
[0]
input: 'a.txt', 'b.txt'
fail_if(len(_step.input) == 1)
""")
        env.run_mode = 'dryrun'
        wf = script.workflow()
        # should be ok
        wf.run()

    def testWarnIf(self):
        '''Test action fail if'''
        script = SoS_Script(r"""
[0]
input: 'a.txt'
warn_if(len(_step.input) == 1, 'Expect to see a warning message')
""")
        env.run_mode = 'dryrun'
        wf = script.workflow()
        # should see a warning message.
        wf.run()
        #self.assertRaises(RuntimeError, wf.run)
        script = SoS_Script(r"""
[0]
input: 'a.txt', 'b.txt'
warn_if(len(_step.input) == 1)
""")
        env.run_mode = 'dryrun'
        wf = script.workflow()
        # should be silent
        wf.run()

    def testCheckOutput(self):
        '''Test action check_output'''
        script = SoS_Script(r"""
[0]
check_output('cat test_actions.py', 'abcde' + 'fgh')
""")
        wf = script.workflow()
        # should raise an error
        self.assertRaises(RuntimeError, wf.run)
        #
        script = SoS_Script(r"""
check_output('cat test_actions.py', 'testCheckOutput')
""")
        wf = script.workflow()
        wf.run()

if __name__ == '__main__':
    unittest.main()
