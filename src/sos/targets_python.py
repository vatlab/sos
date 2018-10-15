#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import importlib
import pkg_resources

from .targets import BaseTarget, textMD5
from .utils import env


class Py_Module(BaseTarget):
    '''A target for a Python module.'''

    LIB_STATUS_CACHE = {}

    def __init__(self, module, version=None, autoinstall=False):
        super(Py_Module, self).__init__()
        self._module = module
        self._version = version
        self._autoinstall = autoinstall

    def _install(self, name, autoinstall):
        '''Check existence of Python module and install it using command
        pip install if necessary.'''
        spam_spec = importlib.util.find_spec(name)
        if spam_spec is not None:
            if self._version:
                if hasattr(spam_spec, '__version__'):
                    ver = spam_spec.__version__
                else:
                    try:
                        ver = pkg_resources.get_distribution(name).version
                    except Exception as e:
                        env.logger.debug(f'Failed to get version of {name}: {e}')
                        return True
                if pkg_resources.parse_version(ver) >= pkg_resources.parse_version(self._version):
                    return True
                else:
                    env.logger.error(f'Version {ver} of installed {name} does not match specified version {self._version}')
                    return False
            return True

        if not autoinstall:
            return False
        # try to install it?
        import subprocess
        cmd = ['pip', 'install', self._module if self._autoinstall is True else self._autoinstall]
        env.logger.info(f'Installing python module {name} with command {" ".join(cmd)}')
        ret = subprocess.call(cmd)
        # try to check version
        return ret == 0 and self._install(name, False)

    def target_exists(self, mode='any'):
        if self._module in self.LIB_STATUS_CACHE:
            return self.LIB_STATUS_CACHE[(self._module, self._autoinstall)]
        else:
            ret = self._install(self._module, self._autoinstall)
            self.LIB_STATUS_CACHE[(self._module, self._autoinstall)] = ret
            return ret

    def target_name(self):
        return self._module

    def target_signature(self, mode='any'):
        # we are supposed to get signature of the module, but we cannot
        return textMD5('Python module ' + self._module)
