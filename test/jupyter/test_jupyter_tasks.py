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
import unittest
import subprocess
from ipykernel.tests.utils import wait_for_idle, execute
from sos.jupyter.test_utils import sos_kernel, KC, get_display_data

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

class TestJupyterTasks(unittest.TestCase):
    #
    # Beacuse these tests would be called from sos/test, we
    # should switch to this directory so that some location
    # dependent tests could run successfully
    #
    def setUp(self):
        subprocess.call(['sos', 'purge'])

    def tearDown(self):
        pass

    def testForceTask(self):
        '''Test the execution of tasks with -s force'''
        with sos_kernel() as kc:
            # the cell will actually be executed several times
            # with automatic-reexecution
            code = """
%set -v1
%run -s force
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

    def testPendingTask(self):
        '''Test the execution of tasks with -s force'''
        with sos_kernel() as kc:
            # the cell will actually be executed several times
            # with automatic-reexecution
            code = """
%run -s force -W
[10]
input: for_each={'i': range(2)}
task:
run:
   echo this is jupyter pending test "${i}"
   sleep  ${10+i}

"""
            # these should be automatically rerun by the frontend
            execute(kc=kc, code=code)
            wait_for_idle(kc)
            # check for task?
            execute(kc=kc, code='%tasks')
            res = get_display_data(kc.iopub_channel)
            # get IDs
            # table_localhost_ac755352394584f797cebddf2c0b8ca7"
            tid = res.split('table_localhost_')[-1].split('"')[0]
            # now we have the tid, we can check task info
            execute(kc=kc, code='%taskinfo ' + tid)
            res = get_display_data(kc.iopub_channel)
            self.assertTrue(tid in res)
            # there should be two tasks
            lines = subprocess.check_output(['sos', 'status']).decode().splitlines()
            self.assertGreaterEqual(len(lines), 2)

if __name__ == '__main__':
    unittest.main()
