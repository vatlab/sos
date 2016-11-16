#!/usr/bin/env python3
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
import re
import unittest

def importTests():
    tests = unittest.TestSuite()
    for file in os.listdir('.'):
        match = re.match("^(test_(.*))\\.py$", file)
        if match:
            m = match.group(1)
            print("Adding test cases in %s" % m)
            module = __import__(m)
            tests.addTest(unittest.defaultTestLoader.loadTestsFromModule( module ))
    return tests

if __name__ == '__main__':
    test_runner = unittest.TextTestRunner(verbosity=2)
    try:
        import nose
        # we use nose for testing because the ipython tests have some namespace
        # conflict with unittest.
        nose.core.run(importTests())
    except ImportError:
        test_runner.run(importTests())
