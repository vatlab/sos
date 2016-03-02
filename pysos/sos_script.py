#!/usr/bin/env python
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
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

import re
from collections import OrderedDict

class SoS_Step:
    #
    # A single sos step
    #
    def __init__(self, action):
        self.action = action


class SoS_Workflow:
    #
    # A SoS workflow with multiple steps
    #
    def __init__(self):
        pass

class SoS_Script_Parser:
    # Regular expressions for parsing section headers and options
    _SECT_TMPL = r'''
        \[                                 # [
        (?P<header>[^]]+)                  # anything but ], can be accessed by name 'header'
        \]                                 # ]
        '''
    _OPT_TMPL = r'''
        (?P<option>.*?)                    # non-greedy
        \s*(?P<vi>=)\s*                    # any number of space/tab,
                                           # followed by any of the
                                           # allowed delimiters,
                                           # followed by any space/tab
        (?P<value>.*)$                     # everything up to eol
        '''
    SECTCRE = re.compile(_SECT_TMPL, re.VERBOSE)
    OPTCRE = re.compile(_OPT_TMPL, re.VERBOSE)

    def __init__(self):
        """Parse a sectioned SoS script file.

        Each section in a SoS script contains a header in square brackets ('[]'). The
        header contains a comma separated section name, followed by comma seperated
        key=value options. Section name and options should be separated by a colon (':').
        
        Each section contains, in any order, either comments, a directive (name : values),
        a expression (key = value), or an action (func(...)).

        Values can span multiple lines, as long as they are indented deeper
        than the first line of the value. The action can also span multiple lines
        until it reaches a blank line. Newlines in triple quotes ('''  ''' and """ """)
        are part of the string though. Because an action is a valid python function
        call, it is actually parsed by a Python tokenizer.
        """
        pass

    def read(self, filename):
        with open(filename) as fp:
            for lineno, line in enumerate(fp, start=1):
 

class SoS_Script:
    #
    # A SoS script with multiple scripts
    #
    def __init__(self, script_file):
        pass

    def _parse(self, script_file):
        pass
        
