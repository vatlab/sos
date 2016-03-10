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
import sys
import re
import copy
import argparse
import textwrap
from collections import OrderedDict, defaultdict, Sequence
from itertools import tee, combinations
try:
    # python 2.7
    import izip
except ImportError:
    # python 3.x
    izip = zip

# Python 2.7 should also have this module
from io import StringIO

from .utils import env, Error, _WorkflowDict, SoS_eval, SoS_exec, RuntimeInfo, dehtml, getTermWidth
from .actions import *

class ArgumentError(Error):
    """Raised when an invalid argument is passed."""
    def __init__(self, msg):
        Error.__init__(self, msg)
        self.args = (msg, )

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
        # comment at the end of a section that could be a workflow description
        self.back_comment = ''
        self.parameters = []
        self.assignments = []
        self.directives = []
        self.statements = []
        # is it global section?
        self.is_global = is_global
        # is it the parameters section?
        self.is_parameters = is_parameters
        # indicate the type of input of the last line
        self.category = None
        self.values = []

    def empty(self):
        '''If there is no content (comment does not count)'''
        return self.category is None

    def extend(self, line):
        if self.category == 'directive':
            self.add_directive(None, line)
        elif self.category == 'expression':
            self.add_assignment(None, line)
        else:
            self.add_statement(line)
        self.back_comment = ''

    def add_comment(self, line):
        '''Add comment line'''
        # in parameter section, comments will always be kept
        if self.is_parameters or self.empty():
            self.comment += line.lstrip('#').lstrip()
        else:
            self.back_comment += line.lstrip('#').lstrip()

    def add_assignment(self, key, value):
        '''Assignments are items with '=' type '''
        if key is None:
            # continuation of multi-line assignment
            if self.is_parameters:
                self.parameters[-1][1] += value
            else:
                self.assignments[-1][1] += value
            self.values.append(value)
        else:
            # new assignment
            if self.is_parameters:
                # in assignment section, comments belong to their following
                # parameter definition
                self.parameters.append([key, value, self.comment])
                self.comment = ''
            else:
                self.assignments.append([key, value])
            self.category = 'expression'
            self.values = [value]
        self.back_comment = ''

    def add_directive(self, key, value):
        '''Assignments are items with ':' type '''
        if key is None:
            # continuation of multi-line directive
            self.directives[-1][1] += value
            self.values.append(value)
        else:
            # new directive
            self.directives.append([key, value])
            self.category = 'directive'
            self.values = [value]
        self.back_comment = ''

    def add_statement(self, line):
        '''Assignments are items with ':' type '''
        # there can be only one statement block
        self.statements.append(line)
        if self.category != 'statements':
            self.values = [line]
        else:
            self.values.append(line)
        self.category = 'statements'
        self.back_comment = ''

    def _get_groups(self, ifiles, group_by):
        # for this special case, the step is skipped
        if group_by == 'single':
            return [[x] for x in ifiles]
        elif group_by == 'all':
            return [ifiles]
        elif group_by == 'paired':
            if len(ifiles) // 2 != 0:
                raise ValueError('Paired group_by has to have even number of input files')
            return list(zip(ifiles[:len(ifiles)//2], ifiles[len(ifiles)//2:]))
        elif group_by == 'pairwise':
            f1, f2 = tee(ifiles)
            next(f2, None)
            return [list(x) for x in izip(f1, f2)]
        elif group_by == 'combinations':
            return [list(x) for x in combinations(ifiles, 2)]

    def _get_labels(self, labels):
        if labels is None or not labels:
            labels = []
        elif isinstance(labels, str):
            labels = [labels]
        elif isinstance(labels, list):
            labels = labels
        else:
            raise ValueError('Unacceptable value for parameter labels: {}'.format(labels))
        #
        ifiles = self.locals['step_input']
        set_vars = [{} for x in self._groups]
        for wv in labels:
            values = self.locals[wv]
            if not isinstance(values, list):
                raise ValueError('with_var variable {} is not a list ("{}")'.format(wv, values))
            if len(values) != len(ifiles):
                raise ValueError('Length of variable {} (length {}) should match the number of input files (length {}).'
                    .format(wv, len(values), len(ifiles)))
            file_map = {x:y for x,y in zip(ifiles, values)}
            for idx,val in enumerate(values):
                set_vars[idx]['_' + wv] = [file_map[x] for x in self._groups[idx]]
        return set_vars

    def _get_for_each(self, for_each):
        if for_each is None or not for_each:
            for_each = []
        elif isinstance(for_each, str):
            for_each = [for_each]
        elif isinstance(for_each, list):
            for_each = for_each
        else:
            raise ValueError('Unacceptable value for parameter for_each: {}'.format(for_each))
        #
        for fe_all in for_each:
            loop_size = None
            for fe in fe_all.split(','):
                values = self.locals[fe]
                if not isinstance(values, list):
                    raise ValueError('for_each variable {} is not a list ("{}")'.format(fe, values))
                if loop_size is None:
                    loop_size = len(values)
                elif loop_size != len(values):
                    raise ValueError('Length of variable {} (length {}) should match the length of variable {} (length {}).'
                        .format(fe, len(values), fe_all.split(',')[0], loop_size))
            # expand
            self._groups = self._groups * loop_size
            tmp = []
            for vidx in range(loop_size):
                for idx in range(len(self._vars)):
                    for fe in fe_all.split(','):
                        self._vars[idx]['_' + fe] = self.locals[fe][vidx]
                tmp.extend(copy.deepcopy(self._vars))
            self._vars = tmp

    # directive functions
    def _directive_input(self, *args, **kwargs):
        # first *args are filenames
        if args:
            ifiles = []
            for arg in args:
                if isinstance(arg, basestring):
                    ifiles.append(arg)
                elif isinstance(arg, list):
                    if not all(isinstance(x, basestring) for x in arg):
                        raise RuntimeError('Invalid input file: {}'.format(arg))
                    ifiles.extend(arg)
            self.locals['step_input'] = ifiles
        else:
            ifiles = self.locals['step_input']
        #
        # handle group_by
        if 'group_by' in kwargs:
            self._groups = self._get_groups(ifiles, kwargs['group_by'])
        else:
            self._groups = [ifiles]
        # handle labels
        if 'labels' in kwargs:
            self._vars = self._get_labels(kwargs['labels'])
        else:
            self._vars = [{} for x in self._groups]
        # handle for_each
        if 'for_each' in kwargs:
            self._get_for_each(kwargs['for_each'])


    def _directive_depends(self, *args, **kwargs):
        dfiles = []
        for arg in args:
            if isinstance(arg, basestring):
                dfiles.append(arg)
            elif isinstance(arg, list):
                if not all(isinstance(x, basestring) for x in arg):
                    raise RuntimeError('Invalid input file: {}'.format(arg))
                dfiles.extend(arg)
        self.locals['step_depends'] = dfiles

    def _directive_output(self, *args, **kwargs):
        ofiles = []
        for arg in args:
            if isinstance(arg, basestring):
                ofiles.append(arg)
            elif isinstance(arg, list):
                if not all(isinstance(x, basestring) for x in arg):
                    raise RuntimeError('Invalid input file: {}'.format(arg))
                ofiles.extend(arg)
        self.locals['step_output'] = ofiles

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
                raise RuntimeError('Unrecognized expression type {}'.format(self.category))
            return True
        except Exception as e:
            return False

    def run(self, globals, locals):
        if isinstance(self.index, int):
            locals['workflow_index'] = str(self.index)
        #
        # assignment
        #
        if self.is_parameters:
            for key, value, _ in self.parameters:
                try:
                    locals[key] = SoS_eval(value, globals, locals)
                except Exception as e:
                    raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
        else:
            for key, value in self.assignments:
                try:
                    locals[key] = SoS_eval(value, globals, locals)
                except Exception as e:
                    raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
        #
        # directives
        #
        globals['self'] = self
        self.locals = locals
        self.globals = globals
        locals['step_output'] = []
        locals['step_depends'] = []
        self._groups = [[]]
        self._vars = [[]]
        for key, value in self.directives:
            try:
                ret = SoS_eval('self._directive_{}({})'.format(key, value), globals, locals)
            except Exception as e:
                raise RuntimeError('Failed to process directive {}: {}'.format(key, e))
        # input alias
        if 'input_alias' in self.options:
            locals[self.options['input_alias']] = locals['step_input']
        if 'output_alias' in self.options:
            locals[self.options['output_alias']] = locals['step_output']

        if locals['step_output']:
            signature = RuntimeInfo(locals['step_output'])
            for ofile in locals['step_output']:
                if not os.path.isdir(os.path.split(ofile)[0]):
                    os.makedirs(os.path.split(ofile)[0])
            if signature.validate(''.join(self.statements),
                locals['step_input'], locals['step_output'], locals['step_depends']):
                # everything matches
                env.logger.info('Reusing existing output files {}'.format(', '.join(locals['step_output'])))
                return
        else:
            signature = None
        #
        for g, v in zip(self._groups, self._vars):
            locals.update(v)
            locals['input'] = g
            # action
            try:
                env.logger.debug('Running \n{}'.format(''.join(self.statements)))
                SoS_exec(''.join(self.statements), globals, locals)
            except Exception as e:
                raise RuntimeError('Failed to execute statement\n\n{}\n\n{}'.format(''.join(self.statements), e))
        for ofile in locals['step_output']:
            if not os.path.isfile(ofile):
                raise RuntimeError('Output file {} does not exist after completion of action'.format(ofile))
        if signature:
            signature.write(''.join(self.statements),
                locals['step_input'], locals['step_output'], locals['step_depends'])

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
    def __init__(self, workflow_name, allowed_steps, sections, description):
        '''create a workflow from its name and a list of SoS_Sections (using name matching)'''
        self.name = workflow_name
        self.description = description
        self.sections = []
        self.global_section = None
        self.parameters_section = None
        self.auxillary_sections = []
        #
        for section in sections:
            if section.is_global:
                # section global is shared by all workflows
                self.global_section = section
                continue
            elif section.is_parameters:
                # section parameters is shared by all workflows
                self.parameters_section = section
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
        #
        if allowed_steps:
            all_steps = {x.index:False for x in self.sections}
            #
            for item in allowed_steps.split(','):
                # remove space
                item = ''.join([x for x in item if x != ' '])
                if item.isdigit():
                    # pipeline:100
                    all_steps[int(item)] = True
                elif '-' in item and item.count('-') == 1:
                    l, u = item.split('-')
                    if (l and not l.isdigit()) or (u and not u.isdigit()) or \
                        (l and u and int(l) > int(u)):
                        raise ValueError('Invalid pipeline step item {}'.format(item))
                    # pipeline:-100, pipeline:100+ or pipeline:10-100
                    if not l:
                        l = min(all_steps.keys())
                    if not u:
                        u = max(all_steps.keys())
                    #
                    for key in all_steps.keys():
                        if key >= int(l) and key <= int(u):
                            all_steps[key] = True
                else:
                    raise ValueError('Invalid pipeline step item {}'.format(item))
            # keep only selected steps
            self.sections = [x for x in self.sections if all_steps[x.index]]

    def prepareVars(self):
        '''Prepare global variables '''
        self.globals = globals()
        self.locals = _WorkflowDict()
        # initial values
        try:
            self.locals['home'] = os.environ['HOME']
        except:
            self.locals['home'] = '.'
        #
        self.locals['workflow_name'] = self.name
        self.locals['workdir'] = os.path.abspath('.')
        #
        if self.global_section:
            self.global_section.run(self.globals, self.locals)

    def run(self):
        '''Very preliminary run function
        '''
        self.prepareVars()
        # There is no input initially
        self.locals['step_input'] = []
        #
        if self.parameters_section:
            self.parameters_section.run(self.globals, self.locals)
        for section in self.sections:
            section.run(self.globals, self.locals)
            if 'step_output' in self.locals:
                # passing step output to step_input of next step
                self.locals['step_input'] = self.locals['step_output']
                self.locals.pop('step_output')

    def show(self, parameters=True):
        textWidth = max(60, getTermWidth())
        paragraphs = dehtml(self.description).split('\n\n')
        print('\n' + '\n'.join(textwrap.wrap('{}:  {}'
            .format(self.name, paragraphs[0]), width=textWidth)))
        for paragraph in paragraphs[1:]:
            print('\n'.join(textwrap.wrap(paragraph, width=textWidth)))
        for step in self.sections:
            # hide a step if there is no comment
            text = '{:<22}'.format('  {}_{}:'.format(step.name, step.index)) + step.comment
            print('\n'.join(textwrap.wrap(text, width=textWidth, subsequent_indent=' '*22)))
        if parameters and self.parameters_section:
            print('\nParameters:')
            for k,v,c in self.parameters_section.parameters:
                #
                text = '  ' + k + \
                    (' '*(22-len(k)-2) if len(k)<20 else ' ') + \
                    (c + ' ' if c else '') + \
                    ('(default: {})'.format(v) if v else '')
                print('\n'.join(textwrap.wrap(text, subsequent_indent=' '*22,
                    width=textWidth)))


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
        [a-zA-Z\d*])??                     # but last character cannot be _
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
        # workflows in this script
        section_steps = sum([x.names for x in self.sections], [])
        # (name, None) is auxiliary steps
        self.workflows = list(set([x[0] for x in section_steps if x[1] is not None and '*' not in x[0]]))
        if not self.workflows:
            self.workflows = ['default']
        # this will update values in the default section with
        # values read from command line
        self._parse_args(args)

    def _read(self, fp, fpname):
        self.sections = []
        self.format_version = '1.0'
        self.description = ''
        self.workflow_descriptions = []
        #
        comment_block = 1
        # cursect always point to the last section
        cursect = None
        last_expression = []
        last_statement = []
        all_step_names = []
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
                        self.workflow_descriptions[-1] += line.lstrip('#').lstrip()
                else:
                    # in the parameter section, the comments are description
                    # of parameters and are all significant
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
                continue
            #
            # is it a continuation of uncompleted assignment or directive?
            if cursect and not cursect.isValid():
                cursect.extend(line)
                continue
            #
            # a new line (start from first column)
            #
            # section header?
            mo = self.SECTION_HEADER.match(line)
            if mo:
                # check previous expression before a new assignment
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(lineno -1 , ''.join(cursect.values), 'Invalid ' + cursect.category)
                    cursect.category = None
                    cursect.values = []
                # start a new section
                section_name = mo.group('section_name').strip()
                section_option = mo.group('section_option')
                step_names = []
                step_options = {}
                for name in section_name.split(','):
                    mo = self.SECTION_NAME.match(name)
                    if mo:
                        n, i, di = mo.group('name', 'index', 'default_index')
                        if n:
                            if i is None and '*' in n:
                                parsing_errors.append(lineno - 1, line, 'Auxillary section name cannot contain wildcard character (*).')
                            step_names.append((n, i))
                        if di:
                            step_names.append(('default', di))
                    else:
                        parsing_errors.append(lineno - 1, line, 'Invalid section name')
                if section_option is not None:
                    for option in section_option.split(','):
                        mo = self.SECTION_OPTION.match(option)
                        if mo:
                            opt_name, opt_value = mo.group('name', 'value')
                            if opt_name in step_options:
                                parsing_errors.append(lineno - 1, line, 'Duplicate options')
                            step_options[opt_name] = opt_value
                        else:
                            parsing_errors.append(lineno - 1, line, 'Invalid section option')
                for name in step_names:
                    prev_workflows = [x[0] for x in all_step_names if '*' not in x[0]]
                    for prev_name in all_step_names:
                        # auxillary step
                        if name[1] is None and prev_name[1] is None and name[0] != prev_name[0]:
                            continue
                        # index not euqal (one of them can be None)
                        if name[1] != prev_name[1]:
                            continue
                        # index equal and one of them have wild card character
                        if '*' in name[0]:
                            names = [x for x in prev_workflows if re.match(name[0].replace('*', '.*'), x)]
                        else:
                            names = [name[0]]
                        if '*' in prev_name:
                            prev_names = [x for x in prev_workflows if re.match(prev_name[0].replace('*', '.*'), x)]
                        else:
                            prev_names = [prev_name[0]]
                        if len(set(prev_names) & set(names)):
                            parsing_errors.append(lineno - 1, line, 'Duplicate section names')
                all_step_names.extend(step_names)
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
                if not cursect.isValid():
                    parsing_errors.append(lineno -1 , ''.join(cursect.values), 'Invalid ' + cursect.category)
                    continue
                cursect.values = []
                #
                var_name = mo.group('var_name')
                if var_name in self._DIRECTIVES:
                    parsing_errors.append(lineno - 1, line, 'directive name cannot be used as variables')
                    continue
                # newline should be kept for multi-line assignment
                var_value = mo.group('var_value') + '\n'
                # if first line of the section, or following another assignment
                # this is assignment
                if cursect.empty() or cursect.category == 'expression':
                    cursect.add_assignment(var_name, var_value)
                #
                # if following a directive, this must be start of an action
                elif cursect.category == 'directive':
                    cursect.add_statement('{} = {}\n'.format(var_name, var_value))
                else:
                    # otherwise it is an continuation of the existing action
                    cursect.extend('{} = {}\n'.format(var_name, var_value))
                continue
            #
            # directive?
            mo = self.DIRECTIVE.match(line)
            if mo:
                # check previous expression before a new directive
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(lineno -1 , ''.join(cursect.values), 'Invalid ' + cursect.category)
                    cursect.values = []
                #
                directive_name = mo.group('directive_name')
                # newline should be kept in case of multi-line directive
                directive_value = mo.group('directive_value') + '\n'
                if cursect is None:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed out side of a SoS step'.format(directive_name))
                    continue
                if cursect.is_parameters:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed in {} section'.format(directive_name, self._PARAMETERS_SECTION_NAME))
                    continue
                if not cursect.empty() and cursect.category == 'statements':
                    parsing_errors.append(lineno, line, 'Directive {} should be be defined before step action'.format(directive_name))
                    continue
                cursect.add_directive(directive_name, directive_value)
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
            if cursect.empty() or cursect.category != 'statements':
                # new statement
                cursect.add_statement(line)
            else:
                # existing one
                cursect.extend(line)
        #
        # check the last expression before a new directive
        if cursect and not cursect.isValid():
            parsing_errors.append(lineno -1 , ''.join(cursect.values), 'Invalid ' + cursect.category)
        #
        # if there is any parsing error, raise an exception
        if parsing_errors.errors:
            raise parsing_errors

    def _parse_error(self, msg):
        '''This function will replace error() function in argparse module so that SoS
        can hijack errors raised from it.'''
        raise ArgumentError(msg)

    def _parse_args(self, args):
        '''Parse command line arguments and set values to parameters section'''
        # first, we need to look for the global section and evaluate it because
        # the parameter section might use variables defined in it.
        wf = self.workflow(self.workflows[0] + ':0')
        wf.prepareVars()
        if not wf.parameters_section:
            return
        #
        if not args:
            return
        #
        parser = argparse.ArgumentParser()
        for key, defvalue, _ in wf.parameters_section.parameters:
            try:
                # FIXME: proper evaluation
                defvalue = SoS_eval(defvalue, wf.globals, wf.locals)
            except Exception as e:
                raise RuntimeError('Incorrect initial value {} for parameter {}: {}'.format(defvalue, key, e))
            parser.add_argument('--{}'.format(key), type=type(defvalue),
                nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, basestring) else '?',
                default=defvalue)
        #
        parser.error = self._parse_error
        #
        args = vars(parser.parse_args(args))
        # now change the value with passed values
        for idx in range(len(wf.parameters_section.parameters)):
            if wf.parameters_section.parameters[idx][0] in args:
                wf.parameters_section.parameters[idx][1] = repr(args[wf.parameters_section.parameters[idx][0]])

    def run(self, wf_name=None):
        wf = self.workflow(wf_name)
        wf.run()

    def workflow(self, wf_name):
        '''Return a workflow with name:step specified in wf_name'''
        allowed_steps = None
        if not wf_name:
            wf_name = ''
        else:
            if ':' in wf_name:
                wf_name, allowed_steps = wf_name.split(':', 1)
        if not wf_name:
            if len(self.workflows) == 1:
                wf_name = list(self.workflows)[0]
            elif 'default' in self.workflows:
                wf_name = 'default'
            else:
                raise ValueError('Name of workflow should be specified because '
                    'the script defines more than one pipelines without a default one. '
                    'Available pipelines are: {}.'.format(', '.join(self.workflows)))
        elif wf_name not in self.workflows:
            raise ValueError('Workflow {} is undefined. Available workflows are: {}'.format(wf_name,
                ', '.join(self.workflows)))
        #
        cur_description = None
        description = ''
        for block in self.workflow_descriptions:
            lines = [x for x in block.split('\n') if x.strip()]
            if not lines:
                continue
            for name in self.workflows:
                if lines[0].strip() == name:
                    cur_description = name
                    break
            if cur_description == wf_name:
                description += '\n'.join(lines[1:] if lines[0].strip() == wf_name else lines) + '\n'
            elif cur_description is None:
                self.description += block + '\n'
        for section in self.sections:
            lines = [x for x in section.back_comment.split('\n') if x.strip()]
            if lines and lines[0].strip() == wf_name:
                description += '\n'.join(lines[1:]) + '\n'
        # create workflows
        return SoS_Workflow(wf_name, allowed_steps, self.sections, description)

    def __repr__(self):
        result = 'SOS Script (version {}\n'.format(self.format_version)
        result += 'workflows:\n    ' + '\n    '.join(self.workflows)
        return result

    def show(self):
        textWidth = max(60, getTermWidth())
        if self.description:
            # separate \n\n
            for paragraph in dehtml(self.description).split('\n\n'):
                print('\n'.join(textwrap.wrap(paragraph, width=textWidth)))
        #
        text = 'Available workflows: {}'.format(', '.join(sorted(self.workflows)))
        print('\n' + '\n'.join(textwrap.wrap(text, width=textWidth, subsequent_indent=' '*8)))
        #
        # parameters
        for wf_name in self.workflows:
            #print('\n{}:'.format(wf_name))
            wf = self.workflow(wf_name)
            wf.show(parameters=False)
        #
        if wf.parameters_section:
            print('\nParameters:')
            for k,v,c in wf.parameters_section.parameters:
                #
                text = '  ' + k + \
                    (' '*(22-len(k)-2) if len(k)<20 else ' ') + \
                    (c + ' ' if c else '') + \
                    ('(default: {})'.format(v) if v else '')
                print('\n'.join(textwrap.wrap(text, subsequent_indent=' '*22,
                    width=textWidth)))
    #
    # for testing purposes
    #
    def parameter(self, name):
        wf = self.workflow(self.workflows[0] + ':0')
        if not wf.parameters_section:
            return
        wf.prepareVars()
        for key, value, _ in wf.parameters_section.parameters:
            if key == name:
                return SoS_eval(value, wf.globals, wf.locals)
        return None

def sos_show(args, argv):
    try:
        script = SoS_Script(args.script, argv)
        if args.workflow:
            workflow = script.workflow(args.workflow)
            workflow.show()
        else:
            script.show()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)

def sos_run(args, argv):
    try:
        script = SoS_Script(args.script, argv)
        workflow = script.workflow(args.workflow)
        workflow.run()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
