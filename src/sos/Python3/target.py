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
from sos.target import BaseTarget, textMD5
import importlib

class Py_Module(BaseTarget):
    '''A target for a Python module.'''

    LIB_STATUS_CACHE = {}

    def __init__(self, module):
        super(Py_Module, self).__init__()
        self._module = module

    def _install(self, name):
        '''Check existence of Python module and install it using command
        pip install if necessary.'''
        spam_spec = importlib.util.find_spec(name)
        if spam_spec is not None:
            return True
        # try to install it?
        import subprocess
        ret = subprocess.call(['pip', 'install', self._module])
        return ret == 0

    def exists(self, mode='any'):
        if self._module in self.LIB_STATUS_CACHE:
            return self.LIB_STATUS_CACHE[self._module]
        else:
            ret = self._install(self._module)
            self.LIB_STATUS_CACHE[self._module] = ret
            return ret

    def name(self):
        return self._module

    def signature(self, mode='any'):
        # we are supposed to get signature of the module, but we cannot
        return textMD5('Python module ' + self._module)
