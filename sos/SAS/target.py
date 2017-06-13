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
# To define new targets for your particular language                          #
#                                                                             #
# 1. Import BaseTarget and other needed functions from sos.target             #
# 2. Define a class that is derived from BaseTarget. The class should         #
#    define function exists(mode=any|target|signature) to test if the         #
#    target itself and/or its signature exists, function name to return       #
#    its name, and a signature (hash string) for the target.                  #
# 3. Modify setup.py and add the target to the [sos_targets] section          #
#    of entry_points.                                                         #
#                                                                             #
# The following is a placeholder class for a new target. Please refer to      #
# implementation of existing targets and the 'Extending_SoS' section of       #
# the SoS documentation (vatlab.github.io/SOS/) for details.                  #
#                                                                             #
###############################################################################


# from sos.target import BaseTarget, textMD5
#
# class My_Target(BaseTarget):
#     def __init__(self, name):
#         super(Py_Module, self).__init__()
#         self._name = name
#
#     def exists(self, mode='any'):
#         return True
#
#    def name(self):
#        return self._name
#
#    def signature(self):
#        return textMD5(self._module)
