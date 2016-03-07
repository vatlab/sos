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

import os
import re
import copy
import argparse
from collections import OrderedDict, defaultdict
# Python 2.7 should also have this module
from io import StringIO

from .utils import env, Error
import pprint


class ArgumentError(Error):
    """Raised when an invalid argument is passed."""
    def __init__(self, msg):
        Error.__init__(self, msg)
        self.args = (msg, )

class DuplicateSectionError(Error):
    """Raised when a section is multiply-created."""
    def __init__(self, section):
        Error.__init__(self, "Section %r already exists" % section)
        self.section = section
        self.args = (section, )

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
    def __init__(self, names=[], options=[], is_global=False, is_parameters=False):
        # A step will not have a name and index until it is copied to separate workflows
        self.name = None
        self.index = None
        # it initially hold multiple names with/without wildcard characters
        self.names = names
        self.options = options
        self.comment = ''
        self.parameters = []
        self.assignments = []
        self.directives = []
        self.statements = []
        # is it global section?
        self.is_global = is_global
        # is it the parameters section?
        self.is_parameters = is_parameters
        # indicate the type of input of the last line
        self.last_line = None
    
    def empty(self):
        '''If there is no content (comment does not count)'''
        return self.last_line is None

    def extend(self, line):
        if self.last_line == ':':
            self.add_directive(None, line)
        elif self.last_line == '=':
            self.add_assignment(None, line)
        else:
            self.add_statement(line)

    def add_comment(self, line):
        '''Add comment line'''
        self.comment += ' ' + line.lstrip('#').strip()

    def add_assignment(self, key, value):
        '''Assignments are items with '=' type '''
        if key is None:
            # continuation of multi-line assignment
            if self.is_parameters:
                self.parameters[-1][1] += value
            else:
                self.assignments[-1][1] += value
        else:
            # new assignment
            if self.is_parameters:
                # in assignment section, comments belong to their following
                # parameter definition
                self.parameters.append([key, value, self.comment])
                self.comment = ''
            else:
                self.assignments.append([key, value])
            self.last_line = '='

    def add_directive(self, key, value):
        '''Assignments are items with ':' type '''
        if key is None:
            # continuation of multi-line directive
            self.directives[-1][1] += value
        else:
            # new directive
            self.directives.append([key, value])
            self.last_line = ':'

    def add_statement(self, line):
        '''Assignments are items with ':' type '''
        # there can be only one statement block
        self.statements.append(line)
        self.last_line = '!'

    def __repr__(self):
        result = ''
        if self.is_global:
            result += '## global definitions ##\n'
        elif self.is_parameters:
            result += '[parameters]\n'
        else:
            result += '[{}:{}]'.format(','.join('{}_{}'.format(x,y) if y else x for x,y in self.names),
                ','.join('{}={}'.format(x,y) for x,y in self.options))
        result += self.comment + '\n'
        for key, value, comment in self.parameters:
            # FIXME: proper line wrap
            result += '# {}\n'.format(comment)
            result += '{} = {}\n'.format(key, value)
        for key, value in self.assignments:
            result += '{} = {}\n'.format(key, value)
        for key, value in self.directives:
            result += '{} = {}\n'.format(key, value)
        for line in self.statements:
            result += line
        result += '\n'
        return result


class SoS_Workflow:
    #
    # A SoS workflow with multiple steps
    #
    def __init__(self, workflow_name, sections, description):
        '''create a workflow from its name and a list of SoS_Sections (using name matching)'''
        self.description = description
        self.sections = []
        self.global_section = None
        self.parameters_section = None
        self.auxillary_sections = []
        #
        for section in sections:
            if section.is_global:
                self.global_section = copy.deepcopy(section)
                continue
            elif section.is_parameters:
                self.parameters_section = copy.deepcopy(section)
                continue
            for name, index in section.names:
                if index is None:
                    self.auxillary_sections.append(copy.deepcopy(section))
                elif '*' in name:
                    # check if name match. There should be no special character in section name
                    # so no worry about invalid regular expression
                    pattern = re.compile(name.replace('*', '.*'))
                    if pattern.match(workflow_name):
                        self.sections.append(copy.deepcopy(section))
                        self.sections[-1].name = workflow_name
                        self.sections[-1].index = int(index)
                elif name == workflow_name:
                    self.sections.append(copy.deepcopy(section))
                    self.sections[-1].name = workflow_name
                    self.sections[-1].index = int(index)
        #
        # sort sections by index
        self.sections.sort(key=lambda x: x.index)

    def __repr__(self):
        result = '__WORKFLOW__\n'
        # FIXME: proper line wrap
        result += '# ' + self.description
        if self.global_section:
            result += repr(self.global_section)
        if self.parameters_section:
            result += repr(self.parameters_section)
        for sect in self.sections:
            result += repr(sect)
        return result 


class ExprStack:
    '''This is a helper class that keeps track of partial expressions,
    value of directives, and step actions. We use a stack because these
    values can expand to several lines so we have to keep reading until
    a valid expression is obtained.'''
    def __init__(self):
        # expression type can be DIRECTIVE, EXPRESSION, and STATEMENTS
        # 
        # EXPRESSION should be a valid Python expression
        # STATEMENTS should be valid Python statements
        # DIRECTIVE should be values plus optional keyword arguments, similar
        #     to parameters to function calls.
        self.category = None
        self.values = []

    def clear(self):
        self.category = None
        self.values = []

    def set(self, expr, category):
        # For safety, we require checking and clearing stack before setting a new one.
        # This ensures all parsing error be checked.
        if self.values:
            raise ValueError('Please manually clear expression stack before setting a new one.')
        self.values = [expr]
        self.category = category

    def push(self, value):
        if self.category is None:
            raise RuntimeError('Value cannot be added to ExprStack without being initialized')
        self.values.append(value)

    def isValid(self):
        if not self.values:
            return True
        try:
            if self.category == 'expression':
                compile(''.join(self.values), filename='<string>', mode='eval')
            elif self.category == 'directive':
                # we add func() because the expression can be multi-line and
                # can have keyword-argument like options
                #
                # However, python considers
                #
                #     func('value', )
                #
                # a valid syntax but we do want , to continue to the next line
                if self.values[-1].strip().endswith(','):
                    return False
                compile('func(' + ''.join(self.values) + ')', filename='<string>', mode='eval')
            elif self.category == 'statements':
                compile(''.join(self.values), filename='<string>', mode='exec')
            else:
                sys.exit('Invalid category {}'.format(self.category))
            return True
        except Exception as e:
            return False


class SoS_Script:
    _DIRECTIVES = ['input', 'output', 'depends']
    _SECTION_OPTIONS = ['input_alias', 'output_alias', 'nonconcurrent', 
        'skip', 'blocking', 'sigil', 'target']
    _PARAMETERS_SECTION_NAME = 'parameters'

    # Regular expressions for parsing section headers and options
    _SECTION_HEADER_TMPL = r'''
        \[                                 # [
        (?P<section_name>[\d\w_,*\s]+)     # digit, alphabet, _ and ,
        (:\s*                              # :
        (?P<section_option>[^]]*)          # section options 
        )?                                 # optional 
        \]                                 # ]
        '''

    _SECTION_NAME_TMPL = '''
        ^\s*                               # start
        (?P<name>                          # optional name
        [a-zA-Z*]                          # alphabet or '*'
        ([\w\d_*]*?                        # followed by alpha numeric or '*'
        [a-zA-Z\d*])?                      # but last character cannot be _
        )?                                 # name is optional
        (?(name)                           # if there is name
        (_(?P<index>\d+))?                 #   optional _index
        |(?P<default_index>\d+))           # no name, then index
        \s*$                                  
        '''
    
    _SECTION_OPTION_TMPL = '''
        ^\s*                               # start
        (?P<name>{})                       # one of the option names
        (\s*=\s*                           # =
        (?P<value>.+)                      # value
        )?                                 # value is optional
        \s*$                                  
        '''.format('|'.join(_SECTION_OPTIONS))

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
    SECTION_NAME = re.compile(_SECTION_NAME_TMPL, re.VERBOSE)
    SECTION_OPTION = re.compile(_SECTION_OPTION_TMPL, re.VERBOSE)
    FORMAT_LINE = re.compile(_FORMAT_LINE_TMPL, re.VERBOSE)
    FORMAT_VERSION = re.compile(_FORMAT_VERSION_TMPL, re.VERBOSE)
    DIRECTIVE = re.compile(_DIRECTIVE_TMPL, re.VERBOSE)
    ASSIGNMENT = re.compile(_ASSIGNMENT_TMPL, re.VERBOSE)

    def __init__(self, content, args=[]):
        '''Parse a sectioned SoS script file. Please refer to the SoS manual 
        for detailed specification of this format.
        
        Parameter `content` can be either a filename or a content of a 
        SoS script in unicode.
        '''
        if os.path.isfile(content):
            with open(content) as fp:
                self._read(fp, content)
        else:
            with StringIO(content) as fp:
                self._read(fp, '<string>')
        # 
        # this will update values in the default section with 
        # values read from command line
        self._parse_cla(args)
        #
        self._create_workflows()
        
    def _read(self, fp, fpname):
        self.sections = []
        self.format_version = '1.0'
        self.workflow_descriptions = []
        #
        comment_block = 1
        # cursect always point to the last section
        cursect = None
        last_expression = []
        last_statement = []
        #
        # this ParsingError is a container for all parsing errors. It will be
        # raised after parsing if there is at least one parsing error.
        parsing_errors = ParsingError(fpname)
        stck = ExprStack()
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
                        self.workflow_descriptions[-1] += line.lstrip('#')
                else:
                    if cursect.is_parameters:
                        # in the parameter section, the comments are description
                        # of parameters and are all significant
                        cursect.add_comment(line)
                    elif comment_block == 1 and cursect.empty():
                        # in a regular section, we only record the first comment block
                        cursect.add_comment(line)
                continue
            elif not line.strip():
                # a blank line start a new comment block if we are still
                # in the front of the script
                if cursect is None:
                    comment_block += 1
                    self.workflow_descriptions.append('')
                elif cursect.comment:
                    comment_block += 1
                continue
            #
            # a continuation of previous item?
            if line[0].isspace() and cursect is not None and not cursect.empty():
                if line.strip():
                    cursect.extend(line)
                    stck.push(line)
                continue
            # 
            # is it a continuation of uncompleted assignment or directive?
            if not stck.isValid():
                stck.push(line)
                cursect.extend(line)
                continue
            #
            # a new line (start from first column)
            # 
            # section header?
            mo = self.SECTION_HEADER.match(line)
            if mo:
                # check previous expression before a new assignment
                if not stck.isValid():
                    parsing_errors.append(lineno -1 , ''.join(stck.values), 'Invalid ' + stck.category)
                stck.clear()
                # start a new section
                section_name = mo.group('section_name').strip()
                section_option = mo.group('section_option')
                step_names = []
                step_options = []
                for name in section_name.split(','):
                    mo = self.SECTION_NAME.match(name)
                    if mo:
                        n, i, di = mo.group('name', 'index', 'default_index')
                        if n:
                            step_names.append((n, i))
                        if di:
                            step_names.append(('default', di))
                    else:
                        parsing_errors.append(lineno - 1, line, 'Invalid section name')
                if section_option is not None:
                    for option in section_option.split(','):
                        mo = self.SECTION_OPTION.match(option)
                        if mo:
                            step_options.append(mo.group('name', 'value'))
                        else:
                            parsing_errors.append(lineno - 1, line, 'Invalid section option')
                self.sections.append(SoS_Step(step_names, step_options, is_parameters= step_names and step_names[0][0] == self._PARAMETERS_SECTION_NAME))
                cursect = self.sections[-1]
                continue
            #
            # assignment?
            mo = self.ASSIGNMENT.match(line)
            if mo:
                if cursect is None:
                    self.sections.append(SoS_Step(is_global=True))
                    cursect = self.sections[-1]
                # check previous expression before a new assignment
                if not stck.isValid():
                    parsing_errors.append(lineno -1 , ''.join(stck.values), 'Invalid ' + stck.category)
                stck.clear()
                #
                var_name = mo.group('var_name')
                var_value = mo.group('var_value')
                # if first line of the section, or following another assignment
                # this is assignment
                if cursect.empty() or cursect.last_line == '=':
                    cursect.add_assignment(var_name, var_value)
                    stck.set(var_value, 'expression')
                # 
                # if following a directive, this must be start of an action
                elif cursect.last_line == ':':
                    cursect.add_statement('{} = {}\n'.format(var_name, var_value))
                    stck.set('{} = {}\n'.format(var_name, var_value), 'statements')
                else:
                    # otherwise it is an continuation of the existing action
                    cursect.extend('{} = {}\n'.format(var_name, var_value))
                    stck.set('{} = {}\n'.format(var_name, var_value), 'statements')
                continue
            #
            # directive?
            mo = self.DIRECTIVE.match(line)
            if mo:
                # check previous expression before a new directive
                if not stck.isValid():
                    parsing_errors.append(lineno -1 , ''.join(stck.values), 'Invalid ' + stck.category)
                stck.clear()
                #
                directive_name = mo.group('directive_name')
                directive_value = mo.group('directive_value')
                if cursect is None:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed out side of a SoS step'.format(directive_name))
                    continue
                if cursect.is_parameters:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed in {} section'.format(directive_name, self._PARAMETERS_SECTION_NAME))
                    continue
                if not cursect.empty() and cursect.last_line == '!':
                    parsing_errors.append(lineno, line, 'Directive {} should be be defined before step action'.format(directive_name))
                    continue
                cursect.add_directive(directive_name, directive_value)
                stck.set(directive_value, 'directive')
                continue
            #
            # all others?
            if not cursect or cursect.is_global:
                parsing_errors.append(lineno, line, 'Only variable assignment is allowed before section definitions.')
                continue
            #
            # It should be an action
            if cursect is None:
                parsing_errors.append(lineno, line, 'Action statement is not allowed in global section')
                continue
            elif cursect.is_parameters:
                parsing_errors.append(lineno, line, 'Action statement is not allowed in {} section'.format(self._PARAMETERS_SECTION_NAME))
                continue
            #
            if cursect.empty() or cursect.last_line != '!':
                # new statement
                cursect.add_statement(line)
                stck.clear()
                stck.set(line, 'statements')
            else:
                # existing one
                cursect.extend(line)
                stck.push(line)
        #
        # check the last expression before a new directive
        if not stck.isValid():
            parsing_errors.append(lineno -1 , ''.join(stck.values), 'Invalid ' + stck.category)
        #
        # if there is any parsing error, raise an exception
        if parsing_errors.errors:
            raise parsing_errors

    def _parse_error(self, msg):
        raise ArgumentError(msg)

    def _parse_cla(self, args):
        '''Parse command line arguments and set values to parameters section'''
        if not args:
            return
        # look for parameters section
        sections = [x for x in self.sections if x.is_parameters]
        if len(sections) == 0:
            return
        elif len(sections) > 1:
            raise DuplicateSectionError(self._PARAMETERS_SECTION_NAME)
        #
        parser = argparse.ArgumentParser()
        for key, defvalue, _ in sections[0].parameters:
            try:
                # FIXME: proper evaluation
                defvalue = eval(defvalue)
            except Exception as e:
                raise RuntimeError('Incorrect initial value {} for parameter {}: {}'.format(defvalue, key, e))
            parser.add_argument('--{}'.format(key), 
                nargs='?' if isinstance(defvalue, basestring) else '*', 
                default=defvalue)
        #
        parser.error = self._parse_error
        #
        args = vars(parser.parse_args(args))
        print(args)
        # now change the value with passed values
        for idx in range(len(sections[0].parameters)):
            if sections[0].parameters[idx][0] in args:
                # FIXME: proper passing
                sections[0].parameters[idx][1] = repr(args[sections[0].parameters[idx][0]])

    def _create_workflows(self):
        #
        self.script_description = ''
        # now, let check what workflows have been defined
        section_steps = sum([x.names for x in self.sections], [])
        # (name, None) is auxiliary steps
        workflow_names = set([x[0] for x in section_steps if x[1] is not None and '*' not in x[0]])
        #
        # separate descriptions by workflow names
        cur_description = None
        workflow_description = defaultdict(str)
        for block in self.workflow_descriptions:
            for name in workflow_names:
                if block.lstrip().startswith(name):
                    cur_description = name
                    break
            if cur_description is None:
                self.script_description += block + '\n'
            else:
                workflow_description[cur_description] += block + '\n'
        # create workflows
        self.workflows = {name: SoS_Workflow(name, self.sections, workflow_description[name])
            for name in workflow_names}
    
    def __repr__(self):
        result = 'SOS Script (version {}\n'.format(self.format_version)
        result += 'workflows:\n    ' + '\n    '.join(self.workflows.keys())
        return result

    #
    # for testing purposes
    # 
    def parameter(self, name):
        sections = [x for x in self.sections if x.is_parameters]
        if len(sections) == 0:
            return None
        for key, value, _ in sections[0].parameters:
            if key == name:
                return eval(value)
        return None
        
