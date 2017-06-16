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
from ipykernel.tests.utils import assemble_output, execute, wait_for_idle
from sos.jupyter.test.utils import sos_kernel, get_result, get_display_data
from sos.target import FileTarget

import nose.tools as nt

TIMEOUT = 60
def long_execute(code='', kc=None, **kwargs):
    """wrapper for doing common steps for validating an execution request"""
    from ipykernel.tests.test_message_spec import validate_message
    if kc is None:
        kc = KC
    msg_id = kc.execute(code=code, **kwargs)
    reply = kc.get_shell_msg(timeout=TIMEOUT)
    validate_message(reply, 'execute_reply', msg_id)
    busy = kc.get_iopub_msg(timeout=TIMEOUT)
    validate_message(busy, 'status', msg_id)
    nt.assert_equal(busy['content']['execution_state'], 'busy')

    if not kwargs.get('silent'):
        execute_input = kc.get_iopub_msg(timeout=TIMEOUT)
        validate_message(execute_input, 'execute_input', msg_id)
        nt.assert_equal(execute_input['content']['code'], code)

    return msg_id, reply['content']

class TestKernel(unittest.TestCase):
    #
    # Beacuse these tests would be called from sos/test, we
    # should switch to this directory so that some location
    # dependent tests could run successfully
    #
    def setUp(self):
        self.olddir = os.getcwd()
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

    def tearDown(self):
        os.chdir(self.olddir)

    def testForceTask(self):
        '''Test the execution of tasks with -s force'''
        with sos_kernel() as kc:
            iopub = kc.iopub_channel
            # the cell will actually be executed several times
            # with automatic-reexecution
            code = """
%set -v1 -s force
[10]
input: for_each={'i': range(1)}
task:
run:
   echo this is "${i}"
   sleep ${i}

[20]
input: for_each={'i': range(2)}
task:
run:
   echo this aa is "${i}"
   sleep ${i}

"""
            # these should be automatically rerun by the frontend
            long_execute(kc=kc, code=code)
            wait_for_idle(kc)

if __name__ == '__main__':
    unittest.main()
