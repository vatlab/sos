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


import os
import unittest

from sostestcase import SoSTestCase

from pysos.utils import *

class TestUtils(SoSTestCase):
    def testLogger(self):
        'Test logging level'
        for verbosity in ['0', '1', '2', '3']:
            env.verbosity = verbosity
            env.logger.trace('Verbosity {}:trace message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.debug('Verbosity {}:debug message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.info('Verbosity {}:info message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.warning('Verbosity {}:warning message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.error('Verbosity {}:error message with ``empahsized text`` in between'.format(env.verbosity))
        # log
        if os.path.isfile('test.log'):
            os.remove('test.log')
        env.logfile = 'test.log'
        for verbosity in ['0', '1', '2', '3']:
            env.verbosity = verbosity
            env.logger.trace('Verbosity {}:trace message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.debug('Verbosity {}:debug message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.info('Verbosity {}:info message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.warning('Verbosity {}:warning message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.error('Verbosity {}:error message with ``empahsized text`` in between'.format(env.verbosity))
        # log file should not have any color codes
        with open('test.log') as logfile:
            line_count = 0
            for line in logfile:
                line_count += 1
                self.assertFalse('\033[' in line)
            # 4 lines for all logging level (logging level of logfile is fixed to DEBUG)
            self.assertEqual(line_count, 16)

if __name__ == '__main__':
    unittest.main()
