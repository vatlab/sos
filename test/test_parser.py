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

# passing string as unicode to python 2 version of SoS
# to ensure compatibility
from __future__ import unicode_literals

import os
import unittest

from sostestcase import SoS_TestCase

from pysos import *


class TestParser(unittest.TestCase): #SoSTestCase):
    def testFileFormat(self):
        '''Test recognizing the format of SoS script'''
        # file format must be 'fileformat=SOSx.x'
        script = SoS_Script_Parser()
        self.assertRaises(ParsingError, script.parse,
            '#fileformat=SS2')
        #
        script = SoS_Script_Parser()
        script.read('scripts/section1.sos')
        # not the default value of 1.0
        self.assertEqual(script.format_version, '1.1')

    def testSections(self):
        '''Test section definitions'''
        script = SoS_Script_Parser()
        script.read('scripts/section1.sos')
        # not the default value of 1.0
        for name in ('parameters', 'section_1', 'section_2', 'section_3, section_4'):
            self.assertTrue('parameters' in [x[0] for x in script.sections.keys()])

    def testGlobalVariables(self):
        '''Test definition of variables'''
        script = SoS_Script_Parser()
        script.read('scripts/section1.sos')
        # not the default value of 1.0

    def testParameters(self):
        '''Test parameters section'''
        pass

    def testSectionVariables(self):
        '''Test variables in sections'''
        pass

    def testSectionDirectives(self):
        '''Test directives of sections'''
        script = SoS_Script_Parser()
        self.assertRaises(ParsingError, script.parse,
            '''input: 'filename' ''')

    def testSectionActions(self):
        '''Test actions of sections'''
        pass

if __name__ == '__main__':
    unittest.main()
