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
    def testFileFormat(self):
        '''Test recognizing the format of SoS script'''
        pass

    def testSections(self):
        '''Test section definitions'''
        pass

    def testGlovalVariables(self):
        '''Test definition of variables'''
        pass

    def testParameters(self):
        '''Test parameters section'''
        pass

    def testSectionVariables(self):
        '''Test variables in sections'''
        pass

    def testSectionDirectives(self):
        '''Test directives of sections'''
        pass

    def testSectionActions(self):
        '''Test actions of sections'''
        pass

if __name__ == '__main__':
    unittest.main()
