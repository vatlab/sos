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

###############################################################################
#                                                                             #
# To test support for a language:                                             #
#                                                                             #
# 1. Import needed functions from ipykernel.tests.utils and                   #
#    sos.jupyter.test.utils.                                                  #
# 2. Define a test case with functions to get and put variables               #
#                                                                             #
# The following is a placeholder class for a new unittest. Please refer to    #
# tests for other languages for details.                                      #
#                                                                             #
###############################################################################

# import unittest
# from ipykernel.tests.utils import execute, wait_for_idle
# from sos.jupyter.test.utils import sos_kernel, get_result
#
# class TestKernel(unittest.TestCase):
#     def setUp(self):
#         self.olddir = os.getcwd()
#         if os.path.dirname(__file__):
#             os.chdir(os.path.dirname(__file__))
#
#     def tearDown(self):
#         os.chdir(self.olddir)
#
#     def testGetData(self):
#         with sos_kernel() as kc:
#             iopub = kc.iopub_channel
#             msg_id, content = execute(kc=kc, code='')
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code="%use LAN")
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code="%get VAR")
#             wait_for_idle(kc)
#
#     def testPutData(self):
#         with sos_kernel() as kc:
#             iopub = kc.iopub_channel
#             msg_id, content = execute(kc=kc, code="%use LAN")
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code='')
#             wait_for_idle(kc)
#             msg_id, content = execute(kc=kc, code="%put VAR")
#             wait_for_idle(kc)
#
#
# if __name__ == '__main__':
#     unittest.main()
