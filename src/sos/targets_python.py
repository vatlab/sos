#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import importlib

from sos.targets import BaseTarget, textMD5


class Py_Module(BaseTarget):
    '''A target for a Python module.'''

    LIB_STATUS_CACHE = {}

    def __init__(self, module):
        super(Py_Module, self).__init__()
        self._module = module

    def _check(self, name):
        '''Check existence of Python module and install it using command
        pip install if necessary.'''
        spam_spec = importlib.util.find_spec(name)
        return spam_spec is not None:

    def target_exists(self, mode='any'):
        if self._module in self.LIB_STATUS_CACHE:
            return self.LIB_STATUS_CACHE[self._module]
        else:
            ret = self._check(self._module)
            self.LIB_STATUS_CACHE[self._module] = ret
            return ret

    def target_name(self):
        return self._module

    def target_signature(self, mode='any'):
        # we are supposed to get signature of the module, but we cannot
        return textMD5('Python module ' + self._module)
