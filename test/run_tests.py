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
import sys
import re
import unittest

def importTests():
    tests = unittest.TestSuite()
    for root, dirs, files in os.walk('..'):
        files[:] = [x for x in files if re.match("^(test_(.*))\\.py$", x)]
        for file in files:
            match = re.match("^(test_(.*))\\.py$", file)
            m = match.group(1)
            print("Adding test cases in {}/{}".format(root, file))
            sys.path.insert(0, root)
            module = __import__(m)
            tests.addTest(unittest.defaultTestLoader.loadTestsFromModule( module ))
        dirs[:] = [x for x in dirs if not x.startswith('.') and x not in ('dist', 'build', 'development')]
    return tests

if __name__ == '__main__':
    test_runner = unittest.TextTestRunner(verbosity=2)
    try:
        import nose
        # we use nose for testing because the ipython tests have some namespace
        # conflict with unittest.
        sys.exit(0 if nose.core.run(importTests()) else 1)
    except ImportError:
        sys.exit(0 if test_runner.run(importTests()) else 1)
