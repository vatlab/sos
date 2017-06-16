#!/usr/bin/env python
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
# To test additional actions for a language:                                  #
#                                                                             #
# 1. Import needed functions from sos.sos_script, sos.utils and               #
#    sos.sos_executor                                                         #
# 2. Define a test case with functions to test each new action                #
#                                                                             #
# The following is a placeholder class for a new unittest. Please refer to    #
# tests for other actions for details.                                        #
#                                                                             #
###############################################################################

#
# import unittest
#
# from sos.sos_script import SoS_Script
# from sos.utils import env
# from sos.sos_executor import Base_Executor, ExecuteError
# from sos.target import FileTarget
#
# class TestActions(unittest.TestCase):
#     def setUp(self):
#         env.reset()
#
#     def testACTION(self):
#         script = SoS_Script(r'''
# [0]
# ACTION:
# ''')
#         wf = script.workflow()
#         Base_Executor(wf).run()
#
# if __name__ == '__main__':
#     unittest.main()
