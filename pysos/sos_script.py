#!/usr/bin/env python3
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
import fnmatch
import textwrap

from io import StringIO
from collections import defaultdict
from collections.abc import Sequence

from .utils import env, Error, dehtml, getTermWidth, locate_script, text_repr
from .sos_eval import Undetermined
from .sos_syntax import SOS_FORMAT_LINE, SOS_FORMAT_VERSION, SOS_SECTION_HEADER, \
    SOS_SECTION_NAME, SOS_SECTION_OPTION, SOS_PARAMETERS_SECTION_NAME, \
    SOS_DIRECTIVE, SOS_DIRECTIVES, SOS_ASSIGNMENT, SOS_SUBWORKFLOW, SOS_REPORT_PREFIX

__all__ = ['SoS_Script']

class ParsingError(Error):
    '''Raised when a configuration file does not follow legal syntax.'''
    def __init__(self, filename):
        Error.__init__(self, 'File contains parsing errors: %s' % filename)
        self.filename = filename
        self.errors = []
        self.args = (filename, )

    def append(self, lineno, line, msg):
        self.errors.append((lineno, line))
        self.message += '\n\t[line %2d]: %s\n%s' % (lineno, line, msg)

class SoS_Step:
    '''Parser of a SoS step. This class accepts strings sent by the parser, determine
    their types and add them to appropriate sections (directive, assignment, statement,
    scripts etc) '''
    def __init__(self, context=None, names=[], options={}, is_global=False, is_parameters=False):
        '''A sos step '''
        self.context = context
        # A step will not have a name and index until it is copied to separate workflows
        self.name = None
        self.index = None
        # it initially hold multiple names with/without wildcard characters
        self.names = names
        self.options = options
        self.comment = ''
        # comment at the end of a section that could be a workflow description
        self.back_comment = ''
        # parameters for parameters section
        self.parameters = []
        # everything before step process
        self.statements = []
        # step processes
        self.global_def = ''
        self.task = ''
        # is it global section? This is a temporary indicator because the global section
        # will be inserted to each step of the workflow.
        self.is_global = is_global
        # is it the parameters section?
        self.is_parameters = is_parameters
        # indicate the type of input of the last line
        self.values = []
        self.lineno = None
        #
        self.runtime_options = {}
        if 'sigil' in self.options:
            self.sigil = self.options['sigil']
        else:
            self.sigil = '${ }'
        #
        # string mode to collect all strings as part of an action
        self._action = None
        self._action_options = ''
        self._script = ''

    def category(self):
        '''Determine the category of existing statement'''
        if self.statements:
            if self.statements[-1][0] == '=':
                return 'expression'
            elif self.statements[-1][0] == '%':
                return 'report'
            elif self.statements[-1][0] == ':':
                # a hack. ... to avoid calling isValid recursively
                def validDirective():
                    if not self.values:
                        return True
                    if self.values[-1].strip().endswith(','):
                        return False
                    try:
                        compile('func(' + ''.join(self.values) + ')', filename='<string>', mode='eval')
                    except:
                        return False
                    return True
                if validDirective() and self._action is not None:
                    return 'script'
                else:
                    return 'directive'
            else:
                return 'statements'
        elif self.parameters:
            return 'expression'
        else:
            return None

    def isValid(self):
        '''Determine if the statement, expression or directive is valid. Otherwise
        the parser will continue until a valid multi-line expression or statement
        can be found.'''
        if not self.values:
            return True
        try:
            if self.category() == 'expression':
                compile(''.join(self.values), filename='<string>', mode='eval')
            elif self.category() == 'directive':
                # we add func() because the expression can be multi-line and
                # can have keyword-argument like options
                #
                # However, python considers
                #
                #     func('value', )
                #
                # a valid syntax but we do want , to continue to the next line
                if self.values[-1].strip().endswith(','):
                    self.error_msg = 'Trailing ,'
                    return False
                compile('func(' + ''.join(self.values) + ')', filename='<string>', mode='eval')
            elif self.category() == 'statements':
                compile(''.join(self.values), filename='<string>', mode='exec')
            elif self.category() == 'script':
                #
                # We are talking about this script here
                #
                # [0]
                # input: 'filename',  'filename2', opt=value==1
                # python3:
                #
                # with open('something') as e:
                #   e.write("""
                # [section]
                # """)
                #
                # The script is obviously valid but [section] takes priority so python3 cannot
                # get the complete script and fail. We can potentially fix this problem but
                # we will see another bug, namely, "${}" expand to something that is not
                # string in the script. That is to say
                #
                # [0]
                # input: 'filename',  'filename2', opt=value==1
                # python3:
                #
                # with open('something') as e:
                #   a = ${input}
                #   e.write("""
                # [section]
                # """)
                #
                # would fail because a=${input} is not a valid statement.
                #
                #if self._action in ['python3', 'task']:
                #    # we only know how to parse python script, but that is good enough
                #    try:
                #        compile(textwrap.dedent(self._script), filename='<string>', mode='exec')
                #        return True
                #    except Exception as e:
                #        #self.error_msg = repr(e)
                #        return False
                #else:
                return True
            elif self.category() == 'report':
                return True
            else:
                raise RuntimeError('Unrecognized expression type {}'.format(self.category()))
            return True
        except Exception as e:
            self.error_msg = repr(e)
            return False

    def empty(self):
        '''If there is no content (comment does not count)'''
        return self.category() is None

    def extend(self, line):
        '''Extend the current directive, expression or script'''
        if self.category() == 'directive':
            self.add_directive(None, line)
        elif self.category() == 'expression':
            self.add_assignment(None, line)
        elif self.category() == 'script':
            self._script += line
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

    def add_report(self, line):
        self.wrap_script()
        if self.statements and self.statements[-1][0] == '%':
            self.statements[-1][-1] += line
        else:
            self.statements.append(['%', line])

    def add_assignment(self, key, value, lineno=None):
        '''Assignments are items with '=' type '''
        if key is None:
            # continuation of multi-line assignment
            if self.is_parameters:
                self.parameters[-1][1] += value
            else:
                self.statements[-1][-1] += value
            self.values.append(value)
        else:
            # new assignment
            if self.is_parameters:
                # in assignment section, comments belong to their following
                # parameter definition
                self.parameters.append([key, value, self.comment])
                self.comment = ''
            else:
                self.statements.append(['=', key, value])
            self.values = [value]
        self.back_comment = ''
        if lineno:
            self.lineno = lineno

    def add_directive(self, key, value, lineno=None):
        '''Assignments are items with ':' type '''
        if key is None:
            # continuation of multi-line directive
            self.statements[-1][-1] += value
            self.values.append(value)
            if self._action is not None:
                self._action_options += value
        else:
            # new directive
            self.statements.append([':', key, value])
            self.values = [value]
        self.back_comment = ''
        if lineno:
            self.lineno = lineno

    def add_script(self, key, value, lineno=None):
        '''script starts with key: value'''
        # we need a fake directive here because the : directive can be multi-line and
        # we need to borrow the syntax checking of directives here.
        self.statements.append([':', '__script__', ''])
        self.values = [value]
        self._action = key
        self._action_options = value
        self.back_comment = ''
        if lineno:
            self.lineno = lineno

    def add_statement(self, line, lineno=None):
        '''Assignments are items with ':' type '''
        # there can be only one statement block
        if self.category() != 'statements':
            self.values = [line]
        else:
            self.values.append(line)
        if self.statements and self.statements[-1][0] == '!':
            self.statements[-1][-1] += line
        else:
            self.statements.append(['!', line])
        self.back_comment = ''
        if lineno:
            self.lineno = lineno

    def wrap_script(self):
        '''convert action: script to task: action(script)'''
        if self._action is None:
            return
        # _action options can contain both runtime option and action options
        opt = self._action_options.strip()
        if self.statements[-1][0] != ':' or self.statements[-1][1] != '__script__':
            raise RuntimeError('Failed to parse script')
        self.statements[-1] = ['!', '{}({}{})\n'.format(self._action, text_repr(textwrap.dedent(self._script)), (', ' + opt) if opt else '')]
        self._action = None
        self._action_options = None
        self._script = ''

    def finalize(self):
        ''' split statement and task by last directive '''
        self.wrap_script()
        if not self.statements:
            self.task = ''
            return
        # convert all ! statement to report action
        for idx, statement in enumerate(self.statements):
            if statement[0] == '%':
                # report should be converted to action
                self.statements[idx] = ['!', 'report({})\n'.format(text_repr(statement[1].
                    replace(SOS_REPORT_PREFIX + '\n', '\n').replace(SOS_REPORT_PREFIX + ' ', '')))]
        #
        task_directive = [idx for idx, statement in enumerate(self.statements) if statement[0] == ':' and statement[1] == 'task']
        if not task_directive:
            self.task = ''
            return
        start_task = task_directive[0] + 1
        # convert statement to task
        self.task = ''
        for statement in self.statements[start_task:]:
            if statement[0] == '=':
                self.task += '{} = {}'.format(statement[1], statement[2])
            elif statement[0] == ':':
                if statement[1] in ('input', 'output', 'depends'):
                    raise ValueError('Step task should be defined as the last item in a SoS step')
                elif statement[1] == 'task':
                    raise ValueError('Only one task is allowed for a step')
                # ignore ...
                self.task += '\n'
            else:
                self.task += statement[1]
        # remove task from self.statement
        if task_directive:
            self.statements = self.statements[:start_task]

    def show(self):
        '''Output for command sos show'''
        textWidth = max(60, getTermWidth())
        if self.is_parameters:
            print('Accepted parameters:')
            for k,v,c in self.parameters:
                # paragraphs = dehtml(c).split('\n\n')
                # FIXME: print paragraphs one by one...
                text = '{:<16}'.format(k + ':') + (c + ' ' if c else '') + \
                     ('(default: {})'.format(v) if v else '')
                print('\n'.join(
                     textwrap.wrap(text,
                     initial_indent=' '*2,
                     subsequent_indent=' '*18,
                     width=textWidth)
                     ))
        else:
            text = '  {:<20}'.format('Step {}_{}:'.format(self.name, self.index)) + self.comment
            print('\n'.join(
                textwrap.wrap(text,
                    width=textWidth,
                    initial_indent='',
                    subsequent_indent=' '*22)))


class SoS_Workflow:
    '''A SoS workflow with multiple steps. It is created from multiple sections of a SoS script
    and consists of multiple SoS_Step.'''
    def __init__(self, workflow_name, allowed_steps, sections, description):
        '''create a workflow from its name and a list of SoS_Sections (using name matching)'''
        self.name = workflow_name
        self.description = description
        self.sections = []
        self.auxillary_sections = []
        #
        for section in sections:
            if section.is_parameters:
                self.sections.append(copy.deepcopy(section))
                self.sections[-1].name = workflow_name
                # for ordering purpose, this section is always after global
                self.sections[-1].index = -1
                continue
            for name, index in section.names:
                if 'target' in section.options:
                    self.auxillary_sections.append(section)
                elif fnmatch.fnmatch(workflow_name, name):
                    self.sections.append(copy.deepcopy(section))
                    self.sections[-1].name = workflow_name
                    self.sections[-1].index = 0 if index is None else int(index)
        #
        # sort sections by index
        self.sections.sort(key=lambda x: x.index)
        #
        # disable some disallowed steps
        if allowed_steps:
            all_steps = {x.index:False for x in self.sections if x.index >= 0}
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
            # keep only selected steps (and the global and parameters section)
            self.sections = [x for x in self.sections if x.index < 0 or all_steps[x.index]]
        #
        env.logger.debug('Workflow {} created with {} sections: {}'
            .format(workflow_name, len(self.sections),
            ', '.join('{}_{}'.format(section.name,
                    'global' if section.index == -2 else ('parameters' if section.index == -1 else section.index))
            for section in self.sections)))

    def extend(self, workflow):
        '''Append another workflow to existing one to created a combined workflow'''
        # all sections are simply appended ...
        self.sections.extend(workflow.sections)

    def show(self, parameters=True):
        textWidth = max(60, getTermWidth())
        paragraphs = dehtml(self.description).split('\n\n')
        print('\n'.join(
            textwrap.wrap('{} {}:  {}'.format(
                'Workflow',
                self.name, paragraphs[0]),
                width=textWidth)
            ))
        for paragraph in paragraphs[1:]:
            print('\n'.join(
            textwrap.wrap(paragraph, width=textWidth)
            ))
        for section in [x for x in self.sections if not x.is_parameters]:
            section.show()
        # parameters display at last
        if parameters:
            print('\n')
            for section in [x for x in self.sections if x.is_parameters]:
                section.show()

class SoS_ScriptContent:
    '''A small class to record the script information to be used by nested
    workflow.'''
    def __init__(self, content='', filename=None):
        self.content = content
        self.filename = filename

class SoS_Script:
    def __init__(self, content='', filename=None, transcript=None):
        '''Parse a sectioned SoS script file. Please refer to the SoS manual
        for detailed specification of this format.

        Parameter `content` can be either a filename or a content of a
        SoS script in unicode, which is convenient for passing scripts for
        testing purposes.

        Parameter `filename` should be used if the content should be read
        from a file.
        '''
        if filename:
            content, self.sos_script = locate_script(filename, start='.')
            self.content = SoS_ScriptContent(content, self.sos_script)
        else:
            self.sos_script = '<string>'
            self.content = SoS_ScriptContent(content, None)
        # save a parsed version of the script for displaying purpose only
        self.transcript = transcript
        self.global_def = ''
        # open the file
        if content:
            with StringIO(content) as fp:
                self._read(fp)
        else:
            with open(self.sos_script) as fp:
                self._read(fp)
        #
        # workflows in this script, from sections that are not skipped.
        section_steps = sum([x.names for x in self.sections if \
            not (x.is_global or x.is_parameters) and \
            not ('skip' in x.options and (x.options['skip'] is None or x.options['skip'])) and \
            not ('target' in x.options)], [])
        # (name, None) is auxiliary steps
        self.workflows = list(set([x[0] for x in section_steps if '*' not in x[0]]))
        if not self.workflows:
            self.workflows = ['default']
        #
        # now we need to record the workflows to the global and parameters section so
        # that we know if which has to be included when a subworkflow is used.
        for section in self.sections:
            if section.is_global or section.is_parameters:
                section.names = self.workflows
        #
        # get script descriptions
        cur_description = None
        self.description = ''
        self.workflow_descriptions = defaultdict(str)
        for block in self.descriptions:
            lines = [x for x in block.split('\n') if x.strip()]
            if not lines:
                continue
            for name in self.workflows:
                if lines[0].strip() == name:
                    cur_description = name
                    break
            if cur_description:
                self.workflow_descriptions[cur_description] += '\n'.join(lines[1:] if lines[0].strip() == cur_description else lines) + '\n'
            else:
                self.description += block + '\n'
        for section in self.sections:
            lines = [x for x in section.back_comment.split('\n') if x.strip()]
            for name in self.workflows:
                if lines and lines[0].strip() == name:
                    self.workflow_descriptions[name] += '\n'.join(lines[1:]) + '\n'

    def _read(self, fp):
        self.sections = []
        self.format_version = '1.0'
        self.descriptions = []
        self.gloal_def = ''
        #
        comment_block = 1
        # cursect always point to the last section
        cursect = None
        all_step_names = []
        #
        # this ParsingError is a container for all parsing errors. It will be
        # raised after parsing if there is at least one parsing error.
        parsing_errors = ParsingError(self.sos_script)
        for lineno, line in enumerate(fp, start=1):
            #
            # comments in SoS scripts are mostly informative
            if line.startswith('#'):
                # Comment blocks before any section
                if cursect is None:
                    if comment_block == 1:
                        # look for format information
                        mo = SOS_FORMAT_LINE.match(line)
                        if mo:
                            format_name = mo.group('format_name')
                            if not format_name.upper().startswith('SOS'):
                                parsing_errors.append(lineno, line,
                                    'Unrecognized file format name {}. Expecting SOS.'.format(format_name))
                            mo = SOS_FORMAT_VERSION.match(format_name)
                            if mo:
                                self.format_version = mo.group('format_version')
                            else:
                                parsing_errors.append(lineno, line,
                                    'Unrecognized file format version in {}.'.format(format_name))
                    elif comment_block > 1:
                        # anything before the first section can be pipeline
                        # description.
                        self.descriptions[-1] += line.lstrip('#').lstrip()
                    if self.transcript:
                        self.transcript.write('COMMENT\t{}\t{}'.format(lineno, line))
                else:
                    # in the parameter section, the comments are description
                    # of parameters and are all significant
                    if cursect.is_parameters or cursect.empty():
                        cursect.add_comment(line)
                        if self.transcript:
                            self.transcript.write('COMMENT\t{}\t{}'.format(lineno, line))
                    elif cursect.category() == 'script':
                        cursect.extend(line)
                        if self.transcript:
                            self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
                    elif cursect.category() in ('statement', 'expression') and cursect.isValid():
                        # this can be comment or back comment
                        cursect.add_comment(line)
                        if self.transcript:
                            self.transcript.write('COMMENT\t{}\t{}'.format(lineno, line))
                    else:
                        # ignored.
                        if self.transcript:
                            self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
                continue
            elif not line.strip():
                # a blank line start a new comment block if we are still
                # in the front of the script
                if cursect is None:
                    comment_block += 1
                    self.descriptions.append('')
                else:
                    if cursect.category() in ('statements', 'script'):
                        cursect.extend(line)
                    elif cursect.comment:
                        comment_block += 1
                if self.transcript:
                    self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
                continue
            elif line.startswith(SOS_REPORT_PREFIX):
                if cursect is None:
                    # global section can have reports, but the reports will be
                    # written repeatedly, which is not good.
                    self.sections.append(SoS_Step(is_global=True))
                    cursect = self.sections[-1]
                elif cursect.is_parameters:
                    parsing_errors.append(lineno, line, 'Report line is not allowed in parameters section')
                    if self.transcript:
                        self.transcript.write('ERROR\t{}\t{}'.format(lineno, line))
                    continue
                elif not cursect.isValid():
                    parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid {}: {}'.format(cursect.category(), cursect.error_msg))
                    if self.transcript:
                        self.transcript.write('ERROR\t{}\t{}'.format(lineno, line))
                    continue
                #elif cursect.category() == 'script':
                #    cursect.extend(line)
                #    if self.transcript:
                #        self.transcript.write('SCRIPT\t{}\t{}'.format(lineno, line))
                #    continue
                #
                if len(line) > 1 and (line[1] != ' ' and line[1] != '\n'):
                    parsing_errors.append(lineno, line, 'Invalid report line: {} symbol should be followed by a space.'.format(SOS_REPORT_PREFIX))
                    if self.transcript:
                        self.transcript.write('ERROR\t{}\t{}'.format(lineno, line))
                    continue
                cursect.add_report(line)
                if self.transcript:
                    self.transcript.write('REPORT\t{}\t{}'.format(lineno, line))
                continue
            #
            # a continuation of previous item?
            if line[0].isspace() and cursect is not None and not cursect.empty():
                cursect.extend(line)
                if self.transcript:
                    self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
                continue
            #
            # is it a continuation of uncompleted assignment or directive?
            if cursect and not cursect.isValid():
                cursect.extend(line)
                if self.transcript:
                    self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
                continue
            #
            # a new line (start from first column)
            #
            # section header?
            mo = SOS_SECTION_HEADER.match(line)
            if mo:
                # check previous expression before a new assignment
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid {}: {}'.format(cursect.category(), cursect.error_msg))
                    cursect.values = []
                    cursect.finalize()
                # start a new section
                section_name = mo.group('section_name').strip()
                section_option = mo.group('section_option')
                step_names = []
                step_options = {}
                #
                for name in section_name.split(','):
                    mo = SOS_SECTION_NAME.match(name)
                    if mo:
                        n, i, di = mo.group('name', 'index', 'default_index')
                        if n:
                            if i is None and '*' in n:
                                parsing_errors.append(lineno, line, 'Unindexed section name cannot contain wildcard character (*).')
                            step_names.append((n, i))
                        if di:
                            step_names.append(('default', di))
                    else:
                        parsing_errors.append(lineno, line, 'Invalid section name')
                if section_option is not None:
                    # this does not work directly because list parameter can have ,
                    # without having to really evaluate all complex expressions, we
                    # have to try to put syntax correctly pieces together.
                    pieces = section_option.split(',')
                    idx = 0
                    while True:
                        try:
                            # test current group
                            compile(pieces[idx].strip(), filename = '<string>', mode='exec' if '=' in pieces[idx] else 'eval')
                            # if it is ok, go next
                            idx += 1
                            if idx == len(pieces):
                                break
                        except Exception as e:
                            # error happens merge the next piece
                            if idx < len(pieces) - 1:
                                pieces[idx] += ',' + pieces[idx + 1]
                                # error happens merge the next piece
                                pieces.pop(idx + 1)
                            else:
                                # if no next group, expand previously correct one
                                if idx == 0:
                                    parsing_errors.append(lineno, line, 'Invalid section option')
                                    break
                                # break myself again
                                pieces = pieces[: idx] + pieces[idx].split(',') + pieces[idx+1:]
                                # go back
                                idx -= 1
                                pieces[idx] += '\n' + pieces[idx + 1]
                                pieces.pop(idx+1)
                    #
                    for option in pieces:
                        mo = SOS_SECTION_OPTION.match(option)
                        if mo:
                            opt_name, opt_value = mo.group('name', 'value')
                            #
                            #opt_value should also be a valid Python expression (with quote etc)
                            #which is most likely a string.
                            if opt_value:
                                try:
                                    # now, the expression might depend on some globle varialbe
                                    # or even option so we might not be able to parse it at parsing time
                                    opt_value = eval(opt_value)
                                except Exception as e:
                                    if opt_name == 'sigil':
                                        parsing_errors.append(lineno, line, e)
                                    else:
                                        env.logger.debug('Step option {}={} (line {}) cannot be resolved during parsing.'.format(opt_name, opt_value, lineno))
                                        opt_value = Undetermined(opt_value)
                            if opt_name == 'sigil':
                                if opt_value.count(' ') != 1 or opt_value[0] in (' ', "'") or \
                                    opt_value[-1] in (' ', "'") or \
                                    opt_value.split(' ')[0] == opt_value.split(' ')[1]:
                                    parsing_errors.append(lineno, line, 'Incorrect sigil "{}"'.format(opt_value))
                            if opt_name in step_options:
                                parsing_errors.append(lineno, line, 'Duplicate options')
                            step_options[opt_name] = opt_value
                        else:
                            parsing_errors.append(lineno, line, 'Invalid section option')
                    env.logger.trace('Header parsed with names {} and options {}'
                        .format(step_names, step_options))
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
                            parsing_errors.append(lineno, line, 'Duplicate section names')
                all_step_names.extend(step_names)
                self.sections.append(SoS_Step(self.content, step_names, step_options, is_parameters= step_names and step_names[0][0] == SOS_PARAMETERS_SECTION_NAME))
                cursect = self.sections[-1]
                if self.transcript:
                    self.transcript.write('SECTION\t{}\t{}'.format(lineno, line))
                continue
            #
            # directive?
            mo = SOS_DIRECTIVE.match(line)
            if mo:
                # check previous expression before a new directive
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid {}: {}'.format(cursect.category(), cursect.error_msg))
                    cursect.values = []
                    # allow multiple process-style actions
                    cursect.wrap_script()
                #
                directive_name = mo.group('directive_name')
                # newline should be kept in case of multi-line directive
                directive_value = mo.group('directive_value') + '\n'
                if cursect and cursect.is_parameters:
                    parsing_errors.append(lineno, line, 'Directive {} is not allowed in {} section'.format(directive_name, SOS_PARAMETERS_SECTION_NAME))
                    continue
                # is it an action??
                if directive_name in SOS_DIRECTIVES:
                    if cursect is None:
                        parsing_errors.append(lineno, line, 'Directive {} is not allowed outside of a SoS step'.format(directive_name))
                        continue
                    cursect.add_directive(directive_name, directive_value, lineno)
                    if self.transcript:
                        self.transcript.write('DIRECTIVE\t{}\t{}'.format(lineno, line))
                #
                elif directive_name == 'process':
                    env.logger.warning('Keyword "process" is depredated and will be removed in a later release. Please use "task" instead.')
                    if cursect is None:
                        parsing_errors.append(lineno, line, 'Directive {} is not allowed outside of a SoS step'.format(directive_name))
                        continue
                    cursect.add_directive('task', directive_value, lineno)
                    if self.transcript:
                        self.transcript.write('DIRECTIVE\t{}\t{}'.format(lineno, line))
                else:
                    # should be in script mode, which is ok for global section
                    if cursect is None:
                        self.sections.append(SoS_Step(is_global=True))
                        cursect = self.sections[-1]
                    cursect.add_script(directive_name, directive_value, lineno)
                    if self.transcript:
                        self.transcript.write('SCRIPT_{}\t{}\t{}'.format(directive_name, lineno, line))
                continue
            # if section is string mode?
            if cursect and cursect.isValid() and cursect.category() == 'script':
                cursect.extend(line)
                if self.transcript:
                    self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
                continue
            #
            # assignment?
            mo = SOS_ASSIGNMENT.match(line)
            if mo:
                if cursect is None:
                    self.sections.append(SoS_Step(is_global=True))
                    cursect = self.sections[-1]
                # check previous expression before a new assignment
                if not cursect.isValid():
                    parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid {}: {}'.format(cursect.category(), cursect.error_msg))
                    if self.transcript:
                        self.transcript.write('ERROR\t{}\t{}'.format(lineno, line))
                    continue
                cursect.values = []
                #
                var_name = mo.group('var_name')
                if var_name in SOS_DIRECTIVES:
                    parsing_errors.append(lineno, line, 'directive name cannot be used as variables')
                    if self.transcript:
                        self.transcript.write('ERROR\t{}\t{}'.format(lineno, line))
                    continue
                # newline should be kept for multi-line assignment
                var_value = mo.group('var_value') + '\n'
                # if first line of the section, or following another assignment
                # this is assignment
                if cursect.empty() or cursect.category() == 'expression':
                    cursect.add_assignment(var_name, var_value, lineno)
                #
                # if following a directive, this can be start of an action or between directives
                elif cursect.category() == 'directive':
                    cursect.add_assignment(var_name, var_value, lineno)
                else:
                    # otherwise it is an continuation of the existing action
                    cursect.extend('{} = {}\n'.format(var_name, var_value))
                if self.transcript:
                    self.transcript.write('ASSIGNMENT\t{}\t{}'.format(lineno, line))
                continue
            #
            # all others?
            if not cursect:
                self.sections.append(SoS_Step(is_global=True))
                cursect = self.sections[-1]
                cursect.add_statement(line, lineno)
                if self.transcript:
                    self.transcript.write('STATEMENT\t{}\t{}'.format(lineno, line))
                continue
            elif cursect.is_parameters:
                parsing_errors.append(lineno, line, 'Action statement is not allowed in {} section'.format(SOS_PARAMETERS_SECTION_NAME))
                if self.transcript:
                    self.transcript.write('ERROR\t{}\t{}'.format(lineno, line))
                continue
            #
            if cursect.empty() or cursect.category() != 'statements':
                # new statement
                cursect.add_statement(line, lineno)
                if self.transcript:
                    self.transcript.write('STATEMENT\t{}\t{}'.format(lineno, line))
            else:
                # existing one
                cursect.extend(line)
                if self.transcript:
                    self.transcript.write('FOLLOW\t{}\t{}'.format(lineno, line))
        #
        # check the last expression before a new directive
        if cursect:
            if not cursect.isValid():
                parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid {}: {}'.format(cursect.category(), cursect.error_msg))
            else:
                cursect.finalize()
        #
        # if there is any parsing error, raise an exception
        if parsing_errors.errors:
            raise parsing_errors
        #
        # as the last step, let us insert the global section to all sections
        global_section = [(idx,x) for idx,x in enumerate(self.sections) if x.is_global]
        if global_section:
            for statement in global_section[0][1].statements:
                if statement[0] == '=':
                    self.global_def += '{} = {}\n'.format(statement[1], statement[2])
                    env.readonly_vars.add(statement[1])
                else:
                    self.global_def += statement[1]
            #
            for section in self.sections:
                section.global_def = self.global_def
            # remove the global section after inserting it to each step of the process
            self.sections.pop(global_section[0][0])

    def workflow(self, workflow_name=None, source={}):
        '''Return a workflow with name_step+name_step specified in wf_name
        This function might be called recursively because of nested
        workflow. Additional workflows can be specified in the source parameter.'''
        allowed_steps = None
        if not workflow_name:
            wf_name = ''
        else:
            # if consists of multiple workflows
            if '+' in workflow_name:
                wfs = []
                for wf in workflow_name.split('+'):
                    if not SOS_SUBWORKFLOW.match(wf):
                        raise ValueError('Incorrect workflow name {}'.format(workflow_name))
                    # if this is a combined workflow, extra_section might be specied.
                    wfs.append(self.workflow(wf, source=source))
                combined_wf = wfs[0]
                for wf in wfs[1:]:
                    combined_wf.extend(wf)
                combined_wf.name = workflow_name
                return combined_wf
            # if a single workflow
            # workflow_10-15 etc
            mo = SOS_SUBWORKFLOW.match(workflow_name)
            if not mo:
                raise ValueError('Incorrect workflow name {}'.format(workflow_name))
            wf_name, allowed_steps = mo.group('name', 'steps')
        # check source
        source_scripts = []
        if source:
            if isinstance(source, str):
                source = {'': [source]}
            elif isinstance(source, dict):
                source = {x: ([y] if isinstance(y, str) else y) for x,y in source.items()}
            elif isinstance(source, Sequence):
                source = {'': source}
            else:
                raise RuntimeError('Invalid value for option source {}'.format(source))
            #
            for key in source.keys():
                source_scripts = []
                for sos_file in source[key]:
                    try:
                        if self.sos_script and self.sos_script != '<string>':
                            content = locate_script(sos_file, start=os.path.split(self.sos_script)[0])
                        else:
                            content = locate_script(sos_file)
                    except Exception as e:
                        raise RuntimeError('Source file for nested workflow {} does not exist: {}'.format(sos_file, e))
                    source_scripts.append(SoS_Script(*content))
                #
                source[key] = source_scripts
        # get workflow name from source files
        extra_workflows = {name: list(set(sum([x.workflows for x in scripts], []))) for name,scripts in source.items()}
        extra_workflow_names = []
        for k,v in extra_workflows.items():
            extra_workflow_names.extend(['{}.{}'.format(k, x).lstrip('.') for x in v])
        extra_sections = {name: sum([x.sections for x in scripts], []) for name,scripts in source.items()}
        if extra_sections:
            env.logger.debug('Importing workflows {}'.format(', '.join(extra_workflow_names)))
        #
        if not wf_name:
            if len(self.workflows) == 1:
                wf_name = list(self.workflows)[0]
            elif 'default' in self.workflows:
                wf_name = 'default'
            else:
                raise ValueError('Name of workflow should be specified because '
                    'the script defines more than one pipelines without a default one. '
                    'Available pipelines are: {}.'.format(', '.join(self.workflows)))
        elif wf_name not in self.workflows + extra_workflow_names:
            raise ValueError('Workflow {} is undefined. Available workflows are: {}'.format(wf_name,
                ', '.join(self.workflows + extra_workflow_names)))
        # do not send extra parameters of ...
        sections = []
        # look for relevant sections in self.sections and extra sections from another script
        if '.' in wf_name:
            candidate_sections = extra_sections[wf_name.split('.')[0]]
        else:
            candidate_sections = self.sections + extra_sections.get('', [])
        for section in candidate_sections:
            # skip, skip=True, skip=1 etc are all allowed.
            if 'skip' in section.options and (section.options['skip'] is None or section.options['skip'] is True):
                continue
            if section.is_parameters:
                # include parameter only if they apply to wf_name
                if wf_name.split('.')[-1] in section.names:
                    sections.append(section)
                continue
            if 'target' in section.options:
                # section global is shared by all workflows
                sections.append(section)
                continue
            for name, index in section.names:
                # exact match or filename like match if name contains * etc
                if fnmatch.fnmatch(wf_name.split('.')[-1], name):
                    sections.append(section)
                    break
        return SoS_Workflow(wf_name.split('.')[-1], allowed_steps, sections, self.workflow_descriptions.get(wf_name, ''))

    def show(self):
        textWidth = max(60, getTermWidth())
        if self.description:
            # separate \n\n
            for paragraph in dehtml(self.description).split('\n\n'):
                print('\n'.join(textwrap.wrap(paragraph, width=textWidth)))
        #
        text = 'Available workflows: {}'.format(', '.join(sorted(self.workflows)))
        print('\n' + '\n'.join(textwrap.wrap(text, width=textWidth, subsequent_indent=' '*8)))
        for idx, workflow in enumerate(sorted(self.workflows)):
            wf = self.workflow(workflow)
            # does not define parameters
            wf.show(parameters=idx == len(self.workflows)-1)
            print('')

