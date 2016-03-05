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
# Python 2.7 should also have this module
from io import StringIO

from .utils import env
import pprint

# exception classes
class Error(Exception):
    '''Base class for SoS_ScriptParser exceptions.'''

    def _get_message(self):
        '''Getter for 'message'; needed only to override deprecation in
        BaseException.'''
        return self.__message

    def _set_message(self, value):
        '''Setter for 'message'; needed only to override deprecation in
        BaseException.'''
        self.__message = value

    # BaseException.message has been deprecated since Python 2.6.  To prevent
    # DeprecationWarning from popping up over this pre-existing attribute, use
    # a new property that takes lookup precedence.
    message = property(_get_message, _set_message)

    def __init__(self, msg=''):
        self.message = msg
        Exception.__init__(self, msg)

    def __repr__(self):
        return self.message

    __str__ = __repr__

class DuplicateSectionError(Error):
    """Raised when a section is multiply-created."""

    def __init__(self, section):
        Error.__init__(self, "Section %r already exists" % section)
        self.section = section
        self.args = (section, )

class InterpolationError(Error):
    """Base class for interpolation-related exceptions."""

    def __init__(self, option, section, msg):
        Error.__init__(self, msg)
        self.option = option
        self.section = section
        self.args = (option, section, msg)

class ParsingError(Error):
    """Raised when a configuration file does not follow legal syntax."""

    def __init__(self, filename):
        Error.__init__(self, 'File contains parsing errors: %s' % filename)
        self.filename = filename
        self.errors = []
        self.args = (filename, )

    def append(self, lineno, line, msg):
        self.errors.append((lineno, line))
        self.message += '\n\t[line %2d]: %s\n%s' % (lineno, line, msg)

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
    _DIRECTIVES = ['input', 'output', 'depends']
    _PARAMETERS_SECTION = 'parameters'

    # Regular expressions for parsing section headers and options
    _SECTION_HEADER_TMPL = r'''
        \[                                 # [
        (?P<section_name>[\d\w_,\s]+)      # digit, alphabet, _ and ,
        (:\s*                              # :
        (?P<section_option>[^]]*)          # section options 
        )?                                 # optional 
        \]                                 # ]
        '''

    _FORMAT_LINE_TMPL = r'''
        ^                                  # from first column
        \#fileformat\s*=\s*                # starts with #fileformat=SOS
        (?P<format_name>.*)                # format name
        \s*$                               # till end of line
        '''

    _FORMAT_VERSION_TMPL = r'''
        ^                                  # from first column
        (?P<format_name>[a-zA-Z]+)         # format name
        (?P<format_version>[\d\.]+)        # any number and .
        \s*$                               # till end of line
        '''

    _DIRECTIVE_TMPL = r'''
        ^                                  # from start of line
        (?P<directive_name>{})             # can be input, output or depends
        \s*:\s*                            # followed by :
        (?P<directive_value>.*)            # and values
        '''.format('|'.join(_DIRECTIVES))

    _ASSIGNMENT_TMPL = r'''
        ^                                   # from start of line
        (?P<var_name>[\w_][\d\w_]*)         # variable name
        \s*=\s*                             # assignment
        (?P<var_value>.*)                   # variable content
        '''

    SECTION_HEADER = re.compile(_SECTION_HEADER_TMPL, re.VERBOSE)
    FORMAT_LINE = re.compile(_FORMAT_LINE_TMPL, re.VERBOSE)
    FORMAT_VERSION = re.compile(_FORMAT_VERSION_TMPL, re.VERBOSE)
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

    def parse(self, content):
        '''Parse specified content as string for parsing specified text directly.'''
        with StringIO(content) as fp:
            self._read(fp, '<string>')

    def read(self, filename):
        '''Read a SoS script and parse it '''
        with open(filename) as fp:
            self._read(fp, filename)

    def _isValidExpr(self, expr):
        try:
            compile('(' + ''.join(expr) + ')', filename='<string>', mode='eval')
            return True
        except Exception as e:
            return False

    def _isValidStmt(self, stmt):
        try:
            compile(''.join(stmt), filename='<string>', mode='exec')
            return True
        except Exception as e:
            return False

    def _read(self, fp, fpname):
        self.sections = {}
        self.format_version = '1.0'
        self.workflow_descriptions = []
        #
        comment_block = 1
        cursect = None
        last_expression = []
        last_statement = []
        #
        # this ParsingError is a container for all parsing errors. It will be
        # raised after parsing if there is at least one parsing error.
        parsing_errors = ParsingError(fpname)
        for lineno, line in enumerate(fp, start=1):
            #
            # comments in SoS scripts are mostly informative
            if line.startswith('#'):
                # Comment blocks before any section
                if cursect is None:
                    if comment_block == 1:
                        # look for format information
                        mo = self.FORMAT_LINE.match(line)
                        if mo:
                            format_name = mo.group('format_name')
                            if not format_name.upper().startswith('SOS'):
                                parsing_errors.append(lineno, line,
                                    'Unrecognized file format name {}. Expecting SOS.'.format(format_name))
                            mo = self.FORMAT_VERSION.match(format_name)
                            if mo:
                                self.format_version = mo.group('format_version')
                            else:
                                parsing_errors.append(lineno, line,
                                    'Unrecognized file format version in {}.'.format(format_name))
                    elif comment_block > 1:
                        # anything before the first section can be pipeline
                        # description.
                        self.workflow_descriptions[-1].append(line)
                else:
                    if cursect[0] == self._PARAMETERS_SECTION:
                        # in the parameter section, the comments are description
                        # of parameters
                        #
                        # multi-line comment
                        if self.sections[cursect] and self.sections[cursect][-1][0] == '#':
                            self.sections[cursect][-1][-1] += ' ' + line.strip()
                        else:
                            # first line of the comment
                            self.sections[cursect].append(['#', line.strip()])
                    elif comment_block == 1:
                        # in a regular section, we only record the first comment block
                        if self.sections[cursect] and self.sections[cursect][-1][0] == '#':
                            self.sections[cursect][-1][-1] += ' ' + line.strip()
                        else:
                            # first line of the comment
                            self.sections[cursect].append(['#', line.strip()])
                continue
            elif not line.strip():
                # a blank line start a new comment block if we are still
                # in the front of the script
                if cursect is None:
                    comment_block += 1
                    self.workflow_descriptions.append([])
                else:
                    if self.sections[cursect] and self.sections[cursect][-1][0] == '#':
                        comment_block += 1
                continue
            #
            # a continuation of previous item?
            if line[0].isspace() and cursect is not None and self.sections[cursect]:
                value = line.strip()
                if value:
                    self.sections[cursect][-1][-1].append(value)
                    if last_expression:
                        last_expression.append(value)
                    if last_statement:
                        last_statement.append(value)
                continue
            #
            # a new line (start from first column)
            # 
            # section header?
            mo = self.SECTION_HEADER.match(line)
            if mo:
                # check previous expression before a new assignment
                if last_expression and not self._isValidExpr(last_expression):
                    parsing_errors.append(lineno -1 , ''.join(last_expression), 'Invalid expression')
                if last_statement and not self._isValidStmt(last_statement):
                    parsing_errors.append(lineno -1 , ''.join(last_statement), 'Invalid statement')
                # start a new section
                section_name = mo.group('section_name').strip()
                section_option = mo.group('section_option')
                cursect = (section_name, section_option)
                self.sections[cursect] = []
                last_statement = []
                continue
            #
            # assignment?
            mo = self.ASSIGNMENT.match(line)
            if mo:
                if cursect is None:
                    cursect = ('__global__', None)
                    self.sections[cursect] = []
                # check previous expression before a new assignment
                if last_expression and not self._isValidExpr(last_expression):
                    parsing_errors.append(lineno -1 , ''.join(last_expression), 'Invalid expression')
                #
                var_name = mo.group('var_name')
                var_value = mo.group('var_value')
                # if first line of the section, or following another assignment
                # this is assignment
                if not self.sections[cursect] or self.sections[cursect][-1][0] == '=':
                    self.sections[cursect].append(['=', var_name, [var_value]])
                    last_expression = [var_value]
                # 
                # if following a directive, this must be start of an action
                elif self.sections[cursect][-1][0] == ':':
                    self.sections[cursect].append(
                        ['!', ['{} = {}\n'.format(var_name, var_value)]])
                    last_expression = []
                    last_statement = ['{} = {}\n'.format(var_name, var_value)]
                else:
                    #
                    # otherwise it is an continuation of the existing action
                    self.sections[cursect][-1][-1].append(
                        '{} = {}\n'.format(var_name, var_value))
                    last_expression = []
                    last_statement.append('{} = {}\n'.format(var_name, var_value))
                continue
            #
            # directive?
            mo = self.DIRECTIVE.match(line)
            if mo:
                # check previous expression before a new directive
                if last_expression and not self._isValidExpr(last_expression):
                    parsing_errors.append(lineno -1 , ''.join(last_expression), 'Invalid expression')
                #
                directive_name = mo.group('directive_name')
                directive_value = mo.group('directive_value')
                if cursect is None:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed out side of a SoS step'.format(directive_name))
                    continue
                if cursect[0] == self._PARAMETERS_SECTION:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed in {} section'.format(directive_name, self._PARAMETERS_SECTION))
                    continue
                if self.sections[cursect] and self.sections[cursect][-1][0] == '!':
                    parsing_errors.append(lineno, line, 'Directive {} should be be defined before step action'.format(directive_name))
                    continue
                self.sections[cursect].append([':', directive_name, [directive_value]])
                last_expression = [directive_value]
                continue
            # 
            # is it a continuation of uncompleted assignment or directive?
            if last_expression and not self._isValidExpr(last_expression):
                last_expression.append(line)
                self.sections[cursect][-1][-1].append(line)
                continue
            if last_statement and not self._isValidStmt(last_statement):
                last_statement.append(line)
                self.sections[cursect][-1][-1].append(line)
                continue                
            #
            # all others?
            if not cursect or cursect == ('__global__', None):
                parsing_errors.append(lineno, line, 'Only variable assignment is allowed before section definitions.')
                continue
            #
            # It should be an action
            if cursect is None:
                parsing_errors.append(lineno, line, 'Action statement is not allowed in global section')
                continue
            elif cursect[0] == self._PARAMETERS_SECTION:
                parsing_errors.append(lineno, line, 'Action statement is not allowed in {} section'.format(self._PARAMETERS_SECTION))
                continue

            if not self.sections[cursect] or self.sections[cursect][-1][0] != '!':
                self.sections[cursect].append( ['!', [line]])
                last_statement = [line]
            else:
                self.sections[cursect][-1][-1].append(line)
                last_statement.append(line)
            last_expression = []
        #
        # check the last expression before a new directive
        if last_expression and not self._isValidExpr(last_expression):
            parsing_errors.append(lineno -1 , ''.join(last_expression), 'Invalid expression')
        # check the last statement before a new directive
        if last_statement and not self._isValidStmt(last_statement):
            parsing_errors.append(lineno -1 , ''.join(last_statement), 'Invalid statement')
        # if there is any parsing error, raise an exception
        if parsing_errors.errors:
            raise parsing_errors
        #
        # now, let us merge the multi-line strings and check if the syntax is correct
        #
        for header, content in self.sections.items():
            for item in content:
                # expression
                if item[0] == '=':
                    item[2] = '\n'.join(item[2]).strip()
                if item[0] == ':':
                    item[2] = '({})'.format('\n'.join(item[2]).strip())
                if item[0] == '!':
                    item[1] = ''.join(item[1]).strip()
        #
        # pp = pprint.PrettyPrinter()
        # pp.pprint(self.sections)
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
