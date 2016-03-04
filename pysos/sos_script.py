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

from .utils import env
import pprint

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
    _SECTION_HEADER_TMPL = r'''
        \[                                 # [
        (?P<section_name>[\d\w_,\s]+)      # digit, alphabet, _ and ,
        (:\s*                           # optional Options
        (?P<section_option>[^]]*)           # section options 
        )?
        \]                                 # ]
        '''

    _FORMAT_LINE_TMPL = r'''                   
        ^                                  # from first column
        \#fileformat=SOS                    # starts with #fileformat=SOS
        (?P<format_version>[\d\.]+)        # any number and .
        \s*$                               # till end of line
        '''

    _DIRECTIVE_TMPL = r'''
        ^                                   # from start of line
        (?P<directive_name>input|output|depends)  # can be input, output or depends
        \s*:\s*                             # followed by :
        (?P<directive_value>.*)           # and values
        '''

    _ASSIGNMENT_TMPL = r'''
        ^                                   # from start of line
        (?P<var_name>[\w_][\d\w_]+)         # variable name
        \s*=\s*                             # assignment
        (?P<var_value>.*)                   # variable content
        '''

    SECTION_HEADER = re.compile(_SECTION_HEADER_TMPL, re.VERBOSE)
    FORMAT_LINE = re.compile(_FORMAT_LINE_TMPL, re.VERBOSE)
    DIRECTIVE = re.compile(_DIRECTIVE_TMPL, re.VERBOSE)
    ASSIGNMENT = re.compile(_ASSIGNMENT_TMPL, re.VERBOSE)

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
            self._read(fp)

    def _read(self, fp):
        # a very crude parser
        self.sections = {}
        cursect = None
        in_action = False
        self.format_version = '1.0'

        for lineno, line in enumerate(fp, start=1):
            comment_block = 1
            if line.startswith('#'):
                if comment_block == 1:
                    # look for format information
                    mo = self.FORMAT_LINE.match(line)
                    if mo:
                        self.format_version = mo.group('format_version')
            # a blank line
            elif not line.strip():
                comment_block += 1
                if not in_action:
                    curitem = None
            # a continuation of previous item?
            elif line[0].isspace() and cursect is not None and self.sections[cursect]:
                value = line.strip()
                if value:
                    self.sections[cursect][-1][-1].append(value)
            else:
                # section header?
                mo = self.SECTION_HEADER.match(line)
                if mo:
                    section_name = mo.group('section_name').strip()
                    section_option = mo.group('section_option')
                    cursect = (section_name, section_option)
                    self.sections[cursect] = []
                    continue
                # assignment?
                mo = self.ASSIGNMENT.match(line)
                if mo:
                    if cursect is None:
                        cursect = '__global__'
                        self.sections[cursect] = []
                    var_name = mo.group('var_name')
                    var_value = mo.group('var_value')
                    if not self.sections[cursect] or self.sections[cursect][-1][0] == '=':
                        self.sections[cursect].append(
                            ['=', var_name, [var_value]])
                    elif self.sections[cursect][-1][0] != '!':
                        self.sections[cursect].append(
                            ['!', ['{} = {}\n'.format(var_name, var_value)]])
                    else:
                        self.sections[cursect][-1][-1].append(
                            '{} = {}\n'.format(var_name, var_value))

                    continue
                # directive?
                mo = self.DIRECTIVE.match(line)
                if mo:
                    if cursect is None:
                        raise RuntimeError('directives are not allowed in global section')
                    if self.sections[cursect]:
                        if self.sections[cursect][-1][0] == '!':
                            raise RuntimeError('Cannot define a directive after action')
                    directive_name = mo.group('directive_name')
                    directive_value = mo.group('directive_value')
                    self.sections[cursect].append(
                        [':', directive_name, [directive_value]])
                    continue
                # all others?
                if not cursect or cursect == '__global__':
                    raise RuntimeError('Line not recognized: {}: cursect {}'.format(line, cursect))
                #
                if not self.sections[cursect]:
                    raise RuntimeError('Line not recognized in section {}: {}'.format(cursect, line))
                #
                # It should be an action
                if not self.sections[cursect] or self.sections[cursect][-1][0] != '!':
                    self.sections[cursect].append( ['!', [line]])
                else:
                    self.sections[cursect][-1][-1].append(line)
        #
        # now, let us merge the multi-line strings and check if the syntax is correct
        #
        for header, content in self.sections.items():
            for item in content:
                # expression
                if item[0] == '=':
                    value = '\n'.join(item[2]).strip()
                    try:
                        compile(value, filename='<string>', mode='eval')
                        item[2] = value
                    except Exception as e:
                        env.logger.error('Invalid assignment of "{}" to variable "{}": {}'
                            .format(value, item[1], e))
                if item[0] == ':':
                    value = '({})'.format('\n'.join(item[2]).strip())
                    try:
                        compile(value, filename='<string>', mode='eval')
                        item[2] = value
                    except Exception as e:
                        env.logger.error('Invalid directive {} with value "{}": {}'
                            .format(item[1], ' '.join(item[2]).strip(), e))
                if item[0] == '!':
                    value = ''.join(item[1]).strip()
                    try:
                        compile(value, filename='<string>', mode='exec')
                        item[1] = value
                    except Exception as e:
                        env.logger.error('Invalid step action with value "{}": {}'
                            .format(value, e))
        #pp = pprint.PrettyPrinter()
        #pp.pprint(self.sections)
        #
 

class SoS_Script:
    #
    # A SoS script with multiple scripts
    #
    def __init__(self, script_file, args):
        script = SoS_ScriptParser()
        script.read(script_file)
        # 
        # this will update values in the default section with 
        # values read from command line
        self.parseCommandLineArgs(script, args)
        #
        # parse header and get a series of workflows
        self.parseWorkflows(script)

    def parseCommandLineArgs(self, script, args):
        #
        # look for parameters section
        parameters = []
        if 'parameters' in [x[0] for x in script.sections.keys()]:
            parameter_section = [y for x,y in script.sections.items() if x[0] == 'parameters'][0]
            for items in parameter_section:
                try:
                    parameters.append([items[1], eval(items[2])])
                    # FIXME, check type
                except Exception as e:
                    raise RuntimeError('Incorrect initial value for parameter {}'.format(items[2]))
        else:
            return
        parser = argparse.ArgumentParser()
        for var, defvalue in parameters:
            parser.add_argument('--{}'.format(var), 
                nargs='1' if isinstance(defvalue, str) else '*', default=defvalue)
        #
        args = vars(parser.parse_args(args))
        #
        # now change the value with passed values
        for items in parameter_section:
            if items[1] in args:
                items[2] = args[items[1]]
            else:
                items[2] = eval(items[2])
            
    def parseWorkflows(self, script):
        headers = [x if isinstance(x, str) else x[0] for x in script.keys()]
        # FIXME: process wild card characters
        # FIXME: separate headers shared by multiple steps
        # FIXME: parse header names

        # we temporarily say the script defines a single default 
        # pipeline and we do not worry about its order now
        #workflow = [y for y 
        self.workflows = {
                'default': 
                    script
            }
