#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
import os
import ast
import copy
import fnmatch
import textwrap
import shutil
import types

# used by structural directive
import sys
assert sys
import glob
assert glob

from io import StringIO
from tokenize import generate_tokens
from uuid import uuid4

from .utils import env, Error, locate_script, text_repr
from .eval import on_demand_options
from .targets import textMD5, file_target, sos_targets, path, paths
from .syntax import SOS_FORMAT_LINE, SOS_FORMAT_VERSION, SOS_SECTION_HEADER, \
    SOS_SECTION_NAME, SOS_SECTION_OPTION, SOS_DIRECTIVE, SOS_DIRECTIVES, \
    SOS_SUBWORKFLOW, SOS_INCLUDE, SOS_FROM_INCLUDE, SOS_AS, SOS_STRU, SOS_IF, \
    SOS_ELIF, SOS_ELSE, SOS_ENDIF, SOS_CELL, SOS_MAGIC, INDENTED

__all__ = ['SoS_Script']

class ParsingError(Error):
    '''Raised when a configuration file does not follow legal syntax.'''
    def __init__(self, filename):
        Error.__init__(self, f'File contains parsing errors: {filename if filename != "<string>" else ""}')
        self.filename = filename
        self.errors = []
        self.args = (filename, )

    def append(self, lineno, line, msg):
        if (lineno, line) in self.errors:
            return
        self.errors.append((lineno, line))
        self.message += f'\n\t[line {lineno:2d}]: {line}\n{msg}'

def get_type_hint(stmt):
    try:
        ns = {'file_target': file_target, 'sos_targets': sos_targets,
                'path': path, 'paths': paths}
        # let us grab the part before =
        exec(stmt.split('=', 1)[0], ns)
        # if it can compile, it can be typetrait, or something like
        #
        # python: input='a'
        #
        # where input is recognied
        #
        if '__annotations__' in ns and all(not isinstance(x, types.BuiltinFunctionType) and callable(x) for x in ns['__annotations__'].values()):
            return ns.popitem()[1]
        return None
    except:
        # if cannot compile, not type hint
        #
        # python: args='whatever'
        return None

def extract_option_from_arg_list(options, optname, default_value):
    if not options:
        return default_value, options
    try:
        args = list(ast.iter_fields(ast.parse(f"f({options})", mode='eval')))[0][1].keywords
        for idx,field in enumerate(args):
            if field.arg == optname:
                try:
                    value = eval(compile(ast.Expression(body=field.value), filename="<ast>", mode="eval"))
                    new_options = ','.join([x for x in options.split(',') if not x.strip().startswith(optname)])
                    return value, new_options.strip()
                except:
                    raise ValueError(f"A constant value is expected for option {optname}: {options} provided.")
        return default_value, options
    except SyntaxError as e:
        raise ValueError(f"Expect a list of keyword arguments: {options} provided")

def separate_options(options):
    pieces = options.split(',')
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
                    raise ValueError('Invalid section option')
                # break myself again
                pieces = pieces[: idx] + pieces[idx].split(',') + pieces[idx+1:]
                # go back
                idx -= 1
                pieces[idx] += '\n' + pieces[idx + 1]
                pieces.pop(idx+1)
    return pieces

def replace_sigil(text, sigil):
    if sigil == '{ }':
        return text
    if sigil is not None and (sigil.count(' ') != 1 or sigil[0] in (' ', "'") or \
        sigil[-1] in (' ', "'") or \
        sigil.split(' ')[0] == sigil.split(' ')[1]):
        raise ValueError(f'Incorrect sigil "{sigil}"')
    # then we need to replace left sigil as { and right sigil asn }
    l, r = sigil.split(' ')
    # now that we have the correct sigil
    # first, we need to replace all { as {{ and } as }}
    lp = re.compile(re.escape(l))
    rp = re.compile(re.escape(r))
    final_text = ''
    while True:
        pieces = lp.split(text, 1)
        if len(pieces) == 1:
            # cannot split
            final_text += text.replace('{', '{{').replace('}', '}}')
            break
        else:
            # find right sigil
            right_pieces = rp.split(pieces[1], 1)
            if len(right_pieces) == 1:
                raise ValueError(f'Missing right sigil in {text}')
            final_text += pieces[0].replace('{', '{{').replace('}', '}}') + \
                '{' + right_pieces[0].replace('{', '{{').replace('}', '}}') + '}'
            text = right_pieces[1]
    # finally, replace LSIGIL etc
    return final_text



class SoS_Step:
    '''Parser of a SoS step. This class accepts strings sent by the parser, determine
    their types and add them to appropriate sections (directive, statement,
    scripts etc) '''
    def __init__(self, context=None, names=None, options=None, is_global=False):
        '''A sos step '''
        self.context = context
        # A step will not have a name and index until it is copied to separate workflows
        self.name = None
        self.index = None
        self.alias = None
        # it initially hold multiple names with/without wildcard characters
        self.names = [] if names is None else names
        self.comment = ''
        self.comment_ended = False
        # everything before step process
        self.statements = []
        self.parameters = {}
        # step processes
        self.global_def = ''
        self.task = ''
        # is it global section? This is a temporary indicator because the global section
        # will be inserted to each step of the workflow.
        self.is_global = is_global
        # indicate the type of input of the last line
        self.values = []
        self.lineno = 0
        #
        self.runtime_options = {}
        self.options = on_demand_options(options)
        #
        # string mode to collect all strings as part of an action
        self._action = None
        self._action_options = ''
        self._script = ''

    def has_external_task(self):
        return self.task != ''

    def step_name(self, alias = False):
        if not self.name:
            n, i, a = self.names[0]
            if alias and a:
                return a
            else:
                return n + \
                    (f'_{i}' if isinstance(i, str) and i.isdigit() else '')
        else:
            if alias and self.alias:
                return self.alias
            else:
                return self.name + \
                    (f'_{self.index}' if isinstance(self.index, int) else '')

    def match(self, step_name):
        # if this step provides name...
        for name, index, _ in self.names:
            if step_name == name or step_name == f'{name}_{index}':
                return True
        return False

    def indented_script(self):
        ''' check self._script and see if it is indented '''
        # get all leading space, tab and newline
        leading = INDENTED.match(self._script)
        return leading is not None and leading.group(2)

    def category(self):
        '''Determine the category of existing statement'''
        if self.statements:
            if self.statements[-1][0] == ':':
                # a hack. ... to avoid calling isValid recursively
                def validDirective():
                    if not self.values:
                        return True
                    if self.values[-1].strip().endswith(','):
                        return False
                    try:
                        compile('func(' + ''.join(self.values) + ')', filename='<string>', mode='eval')
                    except Exception:
                        return False
                    return True
                if validDirective() and self._action is not None:
                    return 'script'
                else:
                    return 'directive'
            else:
                return 'statements'
        else:
            return None

    def isValid(self):
        '''Determine if the statement, expression or directive is valid. Otherwise
        the parser will continue until a valid multi-line expression or statement
        can be found.'''
        if not self.values:
            return True
        try:
            if self.category() == 'directive':
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
                # to allow type trait, we will have to test the expression as if in a function
                # definition, with something like "def func(a : str, b : list=[])"
                try:
                    compile('func(' + ''.join(self.values) + ')', filename='<string>', mode='eval')
                except:
                    compile('def func(' + ''.join(self.values) + '):\n  pass', filename='<string>', mode='exec')
            elif self.category() == 'statements':
                compile((''.join(self.values)), filename='<string>', mode='exec')
            elif self.category() == 'script':
                #
                # A valid script has an identation defined at the first line. That is to say
                #
                # line 1
                # line 2
                #
                # is allowed
                #
                #     line 1
                #     line 2
                #
                # line 3
                #
                # is not so the addition of line 3 would fail. However, the last line
                # will be tested before inserted so this function will always return True
                return True
            else:
                raise RuntimeError(f'Unrecognized expression type {self.category()}')
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
        elif self.category() == 'script':
            self._script += line
        else:
            self.add_statement(line)

    def add_comment(self, line):
        '''Add comment line'''
        if self.empty() and not self.comment_ended:
            self.comment += (' ' if self.comment else '') + line.lstrip('#').strip()

    def end_comment(self):
        self.comment_ended = True

    def add_directive(self, key, value, lineno=None):
        '''Assignments are items with ':' type '''
        if key is None:
            # continuation of multi-line directive
            self.statements[-1][2] += value
            self.values.append(value)
            if self._action is not None:
                self._action_options += value
        else:
            # new directive, the comment before it are used
            self.statements.append([':', key, value])
            self.values = [value]
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
        if lineno:
            self.lineno = lineno

    def add_statement(self, line, lineno=None):
        '''statements are regular python statements'''
        # there can be only one statement block
        if self.category() != 'statements':
            self.values = [line]
        else:
            self.values.append(line)
        if self.statements and self.statements[-1][0] == '!':
            self.statements[-1][-1] += line
        else:
            self.statements.append(['!', line])
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
        # under window, the lines will be ended with \r\n, which will cause
        # trouble with textwrap.dedent.
        self._script = '\n'.join(self._script.splitlines()) + '\n'
        # the script will be considered a f-string, but we will need to handle sigil option here
        self._script = text_repr(textwrap.dedent(self._script))
        # let us look for 'expand=""' in options
        prefix = ''
        if 'expand' in opt:
            sigil, opt = extract_option_from_arg_list(opt, 'expand', None)
            if sigil is None or sigil is False:
                pass
            elif sigil is True:
                prefix = 'f'
            else:
                prefix = 'f'
                self._script = replace_sigil(self._script, sigil)
        self.statements[-1] = ['!', f'{self._action}({prefix}{self._script}{(", " + opt) if opt else ""})\n']
        self.values = []
        self._action = None
        self._action_options = None
        self._script = ''

    def get_tokens(self):
        '''Get tokens after input statement'''
        def _get_tokens(statement):
            return [x[1] for x in generate_tokens(StringIO(statement).readline)]

        tokens = []
        for statement in self.statements:
            if statement[0] == ':':
                # we only keep statements after input
                if statement[1] == 'input':
                    tokens = []
                continue
            tokens.extend(_get_tokens(statement[1]))

        if self.task:
            tokens.extend(_get_tokens(self.task))

        return ' '.join(tokens)

    def finalize(self):
        ''' split statement and task by last directive '''
        self.wrap_script()
        if not self.statements:
            self.task = ''
            return
        # handle parameter
        for idx, statement in enumerate(self.statements):
            if statement[0] == ':' and statement[1] == 'parameter':
                if '=' not in statement[2]:
                    if ':' in statement[2]:
                        if not get_type_hint(statement[2]):
                            raise ValueError(f'Invalid type trait in parameter specification {statement[2]}')
                        name, value = statement[2].split(':')
                    else:
                        name = statement[2]
                        value = 'str'
                else:
                    name, value = statement[2].split('=', 1)
                    # ignore type trait if a default value is specified
                    name = name.split(':')[0]
                name = name.strip()
                if name.startswith('_'):
                    raise ValueError(f'Invalid parameter name {name}: names with leading underscore is not allowed.')
                if name in SOS_DIRECTIVES:
                    raise ValueError(f'Invalid parameter name {name}: {name} is a SoS keyword')
                if not value.strip():
                    raise ValueError(f'{self.step_name()}: Invalid parameter definition: {statement[2]}')
                self.statements[idx] = ['!',
                                        f'if "sos_handle_parameter_" in globals():\n    {name} = sos_handle_parameter_({name.strip()!r}, {value})\n', statement[2].strip()]
                self.parameters[name] = value
        # handle tasks
        task_directive = [idx for idx, statement in enumerate(self.statements) if statement[0] == ':' and statement[1] == 'task']
        if not task_directive:
            self.task = ''
        else:
            start_task = task_directive[0] + 1
            # convert statement to task
            self.task = ''
            for statement in self.statements[start_task:]:
                if statement[0] == ':':
                    if statement[1] in ('input', 'output', 'depends'):
                        raise ValueError(f'{self.step_name()}: Step task should be defined as the last item in a SoS step')
                    elif statement[1] == 'task':
                        raise ValueError(f'{self.step_name()}: Only one task is allowed for a step')
                    elif statement[1] == 'parameter':
                        raise ValueError(f'{self.step_name()}: Parameters should be defined before step task')
                    # ignore ...
                    self.task += '\n'
                else:
                    self.task += statement[1]
                self.statements = self.statements[:start_task]
        # merge multiple statments at the end
        if self.statements[-1][0] == '!' and len(self.statements) > 1:
            starting = len(self.statements) - 1
            for idx in range(starting-1,-1,-1):
                if self.statements[idx][0] == '!':
                    starting = idx
                else:
                    break
            # merge
            for idx in range(starting + 1, len(self.statements)):
                self.statements[starting][1] += self.statements[idx][1]
            # remove the rest of the statements
            self.statements = self.statements[:starting + 1]
        #
        # auto provides #859
        if not any(opt in self.options for opt in ('provides', 'shared')) and \
                len([x for x in self.statements if x[0] == ':' and x[1] == 'output']) == 1:
            # let us check if output is a "plain output"
            output_stmt = [x for x in self.statements if x[0] == ':' and x[1] == 'output'][0][2]
            try:
                plain_output = eval(output_stmt)
                if isinstance(plain_output, str) or \
                    (isinstance(plain_output, (list, tuple, set)) and \
                            all(isinstance(x, str) for x in plain_output)):
                    self.options['autoprovides'] = output_stmt
            except:
                # if otuput has options and rely on anything, it cannot be treated as
                # auto output
                pass


    def show(self):
        '''Output for command sos show'''
        textWidth = max(60, shutil.get_terminal_size((80, 20)).columns)
        text = f'  {self.step_name() + ":":<20} ' + self.comment
        print('\n'.join(
            textwrap.wrap(text,
                width=textWidth,
                initial_indent='',
                subsequent_indent=' '*22)))
        for statement in self.statements:
            if statement[0] == '!' and statement[1].startswith('if "sos_handle_parameter_" '):
                print(f'    Parameter: {statement[2]}')


class SoS_Workflow:
    '''A SoS workflow with multiple steps. It is created from multiple sections of a SoS script
    and consists of multiple SoS_Step.'''
    def __init__(self, content, workflow_name, allowed_steps, sections, global_def):
        '''create a workflow from its name and a list of SoS_Sections (using name matching)'''
        self.content = content
        self.name = workflow_name
        self.sections = []
        self.auxiliary_sections = []
        self.global_def = global_def
        #
        for section in sections:
            for name, index, alias in section.names:
                #if 'provides' in section.options or 'shared' in section.options:
                self.auxiliary_sections.append(section)
                self.auxiliary_sections[-1].name = section.names[0][0]
                self.auxiliary_sections[-1].index = None if index is None else int(index)
                self.auxiliary_sections[-1].alias = alias
                self.auxiliary_sections[-1].uuid = uuid4()
                # an auxiliary step can also serve as a regular step
                # as long as it matches workflow name.
                if fnmatch.fnmatch(workflow_name, name):
                    self.sections.append(copy.deepcopy(section))
                    self.sections[-1].name = workflow_name
                    self.sections[-1].index = 0 if index is None else int(index)
                    self.sections[-1].alias = alias
                    self.sections[-1].uuid = uuid4()
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
                        raise ValueError(f'Invalid pipeline step item {item}')
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
                    raise ValueError(f'Invalid pipeline step item {item}')
            # keep only selected steps (and the global section)
            self.sections = [x for x in self.sections if x.index < 0 or all_steps[x.index]]
        #
        env.logger.debug('Workflow {} created with {} sections: {}'
            .format(workflow_name, len(self.sections),
            ', '.join(f'{section.name}_{"global" if section.index == -2 else section.index}'
                      for section in self.sections)))

    def section_by_id(self, uuid):
        for section in self.sections:
            if uuid == section.uuid:
                return section
        for section in self.auxiliary_sections:
            if uuid == section.uuid:
                return section
        raise RuntimeError(f'Failed to find section with uuid {uuid}')

    def extend(self, workflow):
        '''Append another workflow to existing one to created a combined workflow'''
        # all sections are simply appended ...
        self.sections.extend(workflow.sections)

    def has_external_task(self):
        return any(x.has_external_task() for x in self.sections) or \
                any(x.has_external_task() for x in self.auxiliary_sections)

    def parameters(self):
        # collect parameters defined by `parameter:` of steps
        par = {}
        for x in self.sections + self.auxiliary_sections:
            par.update(x.parameters)
        return par

class SoS_ScriptContent:
    '''A small class to record the script information to be used by nested
    workflow.'''
    def __init__(self, content='', filename=None):
        self.content = content
        self.filename = filename
        self.included = []
        self.md5 = self.calc_md5()

    def calc_md5(self):
        if self.content:
            cnt = self.content
        else:
            with open(self.filename) as script:
                cnt = script.read()
        # additional files
        for script in self.included:
            cnt += script[0]
        #
        return textMD5(cnt)

    def add(self, content='', filename=None):
        if filename and not content:
            with open(filename) as script:
                content = script.read()
        if filename not in [x[1] for x in self.included]:
            self.included.append((content, filename))
            self.md5 = self.calc_md5()

    def __repr__(self):
        return f'{self.md5}: filename: {self.filename}, content: {self.content}'


    def __eq__(self, other):
        return self.md5 == other.md5

    def __ne__(self, other):
        return self.md5 != other.md5

class SoS_Script:
    def __init__(self, content='', filename=None, transcript=None):
        '''Parse a sectioned SoS script file. Please refer to the SoS manual
        for detailed specification of this format.

        Parameter `content` can be either a filename or a content of a
        SoS script in unicode, which is convenient for passing scripts for
        testing purposes.

        Parameter `filename` should be used if the content should be read
        from a file.

        A SoS_Script generated from a .sos file contains the following items

        sections:
            A list of SoS_Step objects that present all sections read.
            Among other things, section.names
        '''
        if not content:
            try:
                content, self.sos_script = locate_script(filename, start='.')
            except Exception:
                # try to add .sos extension?
                if not filename.endswith('.sos'):
                    try:
                        content, self.sos_script = locate_script(filename + '.sos', start='.')
                    except Exception:
                        if not filename.endswith('.ipynb'):
                            try:
                                content, self.sos_script = locate_script(filename + '.ipynb', start='.')
                            except Exception as e:
                                env.logger.debug(e)
                                env.logger.error(f'Failed to locate {filename}, {filename}.sos, or {filename}.ipynb')
                                sys.exit(1)
                        else:
                            raise
                else:
                    raise
            # Is this script in sos or ipynb format?
            ext = os.path.splitext(self.sos_script)[-1]
            if ext == '.ipynb':
                # convert ipynb to sos
                from sos_notebook.converter import notebook_to_script
                with StringIO() as script:
                    notebook_to_script(self.sos_script, script)
                    content = script.getvalue()
                self.sos_script = '<string>'
                self.content = SoS_ScriptContent(content, None)
            else:
                self.content = SoS_ScriptContent(content, self.sos_script)
        else:
            self.sos_script = '<string>'
            self.content = SoS_ScriptContent(content, None)
        # save a parsed version of the script for displaying purpose only
        self.transcript = transcript
        self.global_def = ''

        self.description = []

        # open the file
        if content:
            with StringIO(content) as fp:
                self._read(fp)
        else:
            with open(self.sos_script) as fp:
                self._read(fp)
        #
        # workflows in this script, from sections that are not skipped.
        all_section_steps = sum([x.names for x in self.sections], [])
        forward_section_steps = sum([x.names for x in self.sections if \
            not any(opt in x.options for opt in ('provides', 'shared'))], [])
        # (name, None) is auxiliary steps
        self.workflows = list(set([x[0] for x in all_section_steps if '*' not in x[0]]))
        forward_workflows = list(set([x[0] for x in forward_section_steps if '*' not in x[0]]))
        if not forward_workflows:
            self.workflows.append('default')
            self.default_workflow = 'default'
        elif len(forward_workflows) == 1:
            self.default_workflow = forward_workflows[0]
        else:
            self.default_workflow = None
        #
        # now we need to record the workflows to the global section so
        # that we know if which has to be included when a subworkflow is used.
        for section in self.sections:
            if section.is_global:
                section.names = self.workflows

    def _find_include_file(self, sos_file):
        # we could almost use SoS_script directly but we need to be able to start searching
        # from path of the master file.
        try:
            if self.sos_script and self.sos_script != '<string>':
                start_path = os.path.split(self.sos_script)[0]
            else:
                start_path = ''
            try:
                content, script_file = locate_script(sos_file + '.sos', start=start_path)
            except Exception:
                content, script_file = locate_script(sos_file + '.ipynb', start=start_path)
                # convert ipynb to sos
                from sos_notebook.converter import notebook_to_script
                with StringIO() as script:
                    notebook_to_script(script_file, script)
                    content = script.getvalue()
        except Exception:
            raise RuntimeError(
                f'Source file for nested workflow {sos_file} with extension .sos or .ipynb does not exist')

        return content, script_file

    def _include_namespace(self, sos_file, alias):
        content, script_file = self._find_include_file(sos_file)
        self.content.add(content, script_file)
        script = SoS_Script(content, script_file)
        if not alias:
            alias = sos_file
        # section names are changed from A to sos_file.A
        for section in script.sections:
            for idx in range(len(section.names)):
                section.names[idx] = [alias + '.' + section.names[idx][0]] + list(section.names[idx][1:])
        # The global definition of sos_file should be accessible as
        # sos_file.name
        self.sections.extend(script.sections)
        self.global_def += f"{alias} = sos_namespace_({repr(script.global_def)})\n"

    def _include_content(self, sos_file, name_map):
        content, script_file = self._find_include_file(sos_file)
        #
        self.content.add(content, script_file)
        script = SoS_Script(content, script_file)
        if not name_map:
            self.sections.extend(script.sections)
            self.global_def += script.global_def
        else:
            # if name_map, we only include selected names from the script
            for section in script.sections:
                for name, _, _ in section.names:
                    if any(fnmatch.fnmatch(x, name) for x in name_map):
                        # match ...
                        self.sections.append(section)
            # global_def is more complicated
            self.global_def += f"__{sos_file} = sos_namespace_({repr(script.global_def)})\n"
            #
            self.global_def += f'''
for __n, __v in {repr(name_map)}.items():
    if hasattr(__{sos_file}, __n):
        globals()[__v if __v else __n] = getattr(__{sos_file}, __n)
'''
            # env.logger.trace(self.global_def)

    def _read(self, fp):
        self.sections = []
        self.format_version = '1.0'
        self.gloal_def = ''
        #
        comment_block = 1
        # cursect always point to the last section
        cursect = None
        all_step_names = []
        #
        condition_ignore = False
        condition_met = None
        # this ParsingError is a container for all parsing errors. It will be
        # raised after parsing if there is at least one parsing error.
        parsing_errors = ParsingError(self.sos_script)
        for lineno, line in enumerate(fp, start=1):
            #
            # for structural lines
            if SOS_MAGIC.match(line):
                # ignore cell directive in batch mode
                if self.transcript:
                    self.transcript.write(f'COMMENT\t{lineno}\t{line}')

                continue

            if SOS_STRU.match(line):
                # ignore cell directive in batch mode
                if self.transcript:
                    self.transcript.write(f'COMMENT\t{lineno}\t{line}')

                if SOS_INCLUDE.match(line) or SOS_FROM_INCLUDE.match(line):
                    if cursect is not None:
                        parsing_errors.append(lineno, line, 'include magic can only be defined before any other statemetns.')

                    # handle import
                    mo = SOS_INCLUDE.match(line)
                    if mo:
                        sos_files = [x.strip() for x in mo.group('sos_files').split(',')]
                        for sos_file in sos_files:
                            ma = SOS_AS.match(sos_file)
                            self._include_namespace(ma.group('name'), alias=ma.group('alias'))
                        continue
                    mo = SOS_FROM_INCLUDE.match(line)
                    if mo:
                        sos_file = mo.group('sos_file')
                        name_map = {}
                        if mo.group('names') != '*':
                            for wf in mo.group('names').split(','):
                                ma = SOS_AS.match(wf.strip())
                                if not ma:
                                    parsing_errors.append(lineno, line, f'unacceptable include name "{wf}"')
                                else:
                                    name_map[ma.group('name')] = ma.group('alias')
                        self._include_content(sos_file, name_map=name_map)
                        continue

                mo = SOS_CELL.match(line)
                if mo:
                    continue

                mo = SOS_IF.match(line)
                if mo:
                    cond = mo.group('condition')
                    try:
                        cond_value = eval(cond)
                    except Exception as e:
                        parsing_errors.append(lineno, line, f'Invalid expression {cond}: {e}')
                        continue
                    #
                    if cond_value:
                        condition_met = True
                        condition_ignore = False
                    else:
                        condition_met = False
                        condition_ignore = True
                    continue

                mo = SOS_ELIF.match(line)
                if mo:
                    if condition_met is None:
                        parsing_errors.append(lineno, line, f'%elif not following %if: {line}')
                        continue

                    if condition_met:
                        condition_ignore = True
                        continue

                    cond = mo.group('condition')
                    try:
                        cond_value = eval(cond)
                    except Exception as e:
                        parsing_errors.append(lineno, line, f'Invalid expression {cond}: {e}')

                    if cond_value:
                        condition_met = True
                        condition_ignore = False
                    continue

                mo = SOS_ELSE.match(line)
                if mo:
                    if condition_met is None:
                        parsing_errors.append(lineno, line, f'%else not following %if: {line}')
                        continue

                    condition_ignore = condition_met
                    continue

                mo = SOS_ENDIF.match(line)
                if mo:
                    condition_met = None
                    condition_ignore = False
                    continue

                else:
                    parsing_errors.append(lineno, line, f'Unrecognized SoS magic statement: {line}')
                    continue

            if condition_ignore:
                # ignore cell directive in batch mode
                if self.transcript:
                    self.transcript.write(f'COMMENT\t{lineno}\t{line}')
                continue

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
                                                      f'Unrecognized file format name {format_name}. Expecting SOS.')
                            mo = SOS_FORMAT_VERSION.match(format_name)
                            if mo:
                                self.format_version = mo.group('format_version')
                            else:
                                parsing_errors.append(lineno, line,
                                                      f'Unrecognized file format version in {format_name}.')
                    elif comment_block > 1:
                        # anything before the first section can be pipeline
                        # description.
                        if cursect is None:
                            self.description.append(line)
                    if self.transcript:
                        self.transcript.write(f'COMMENT\t{lineno}\t{line}')
                else:
                    # this is description of the section
                    if cursect.empty():
                        cursect.add_comment(line)
                        if self.transcript:
                            self.transcript.write(f'COMMENT\t{lineno}\t{line}')
                    # this is comment in scripts (and perhaps not even comment)
                    elif cursect.category() == 'script':
                        if cursect.indented_script():
                            # if the script is indented and encounters a comment
                            # from first column, switch to comment mode
                            try:
                                cursect.wrap_script()
                            except Exception as e:
                                parsing_errors.append(lineno, line, e)
                            cursect.add_comment(line)
                            if self.transcript:
                                self.transcript.write(f'COMMENT\t{lineno}\t{line}')
                        else:
                            cursect.extend(line)
                            if self.transcript:
                                self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                    # this can be comment or back comment
                    elif cursect.category() == 'statements':
                        if cursect.isValid():
                             # this can be comment or back comment
                            cursect.add_comment(line)
                            if self.transcript:
                                self.transcript.write(f'COMMENT\t{lineno}\t{line}')
                        else:
                            cursect.extend(line)
                            if self.transcript:
                                self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                    else:
                        # ignored.
                        if self.transcript:
                            self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                continue
            elif not line.strip():
                # a blank line start a new comment block if we are still
                # in the front of the script
                if cursect is None:
                    comment_block += 1
                else:
                    if cursect.category() in ('statements', 'script'):
                        cursect.extend(line)
                    elif cursect.comment:
                        cursect.end_comment()
                if self.transcript:
                    self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                continue

            #
            # a continuation of previous item?
            if line[0].isspace() and cursect is not None and not cursect.empty():
                cursect.extend(line)
                if self.transcript:
                    self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                continue
            #
            # is it a continuation of uncompleted directive?
            if cursect and not cursect.isValid():
                cursect.extend(line)
                if self.transcript:
                    self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                continue
            #
            # a new line (start from first column)
            #
            # section header?
            mo = SOS_SECTION_HEADER.match(line)
            if mo:
                # check previous expression before a new section
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]),
                                              f'Invalid {cursect.category()}: {cursect.error_msg}')
                    cursect.values = []
                    try:
                        cursect.finalize()
                    except Exception as e:
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), e)
                # start a new section
                section_name = mo.group('section_name').strip()
                section_option = mo.group('section_option')
                step_names = []
                step_options = {}
                #
                for name in section_name.split(','):
                    mo = SOS_SECTION_NAME.match(name)
                    if mo:
                        n, i, di, al = mo.group('name', 'index', 'default_index', 'alias')
                        if n == 'global' and i is not None:
                            parsing_errors.append(lineno, line, 'Invalid global section definition')
                        if n:
                            if i is None and '*' in n:
                                parsing_errors.append(lineno, line, 'Unindexed section name cannot contain wildcard character (*).')
                            step_names.append([n, i, al])
                        if di:
                            step_names.append(['default', di, al])
                    else:
                        parsing_errors.append(lineno, line, 'Invalid section name')
                if 'global' in [x[0] for x in step_names] and len(step_names) > 1:
                    parsing_errors.append(lineno, line, 'Global section cannot be shared with another step')
                if section_option is not None:
                    # this does not work directly because list parameter can have ,
                    # without having to really evaluate all complex expressions, we
                    # have to try to put syntax correctly pieces together.
                    try:
                        pieces = separate_options(section_option)
                        for option in pieces:
                            mo = SOS_SECTION_OPTION.match(option)
                            if mo:
                                opt_name, opt_value = mo.group('name', 'value')
                                if opt_name in step_options:
                                    parsing_errors.append(lineno, line, 'Duplicate options')
                                step_options[opt_name] = opt_value
                            else:
                                parsing_errors.append(lineno, line, 'Invalid section option')
                    except Exception as e:
                        parsing_errors.append(lineno, line, e)
                    env.logger.trace('Header parsed with names {} and options {}'
                        .format(step_names, step_options))
                for name in step_names:
                    prev_workflows = [x[0] for x in all_step_names if '*' not in x[0]]
                    for prev_name in all_step_names:
                        # auxiliary step
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
                        if len(set(prev_names) & set(names)) and 'global' not in names:
                            parsing_errors.append(lineno, line, f'Duplicate section name {names}')
                all_step_names.extend(step_names)
                if 'global' in [x[0] for x in step_names]:
                    if step_options:
                        parsing_errors.append(lineno, line, 'Global section does not accept any option')
                    self.sections.append(SoS_Step(is_global=True))
                else:
                    self.sections.append(SoS_Step(self.content, step_names, step_options))
                cursect = self.sections[-1]
                if self.transcript:
                    self.transcript.write(f'SECTION\t{lineno}\t{line}')
                continue
            #
            # directive?
            mo = SOS_DIRECTIVE.match(line)
            if mo and not get_type_hint(line):
                # check previous expression before a new directive
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]),
                                              f'Invalid {cursect.category()}: {cursect.error_msg}')
                    cursect.values = []
                    # allow multiple process-style actions
                    try:
                        cursect.wrap_script()
                    except Exception as e:
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), e)
                else:
                    self.sections.append(SoS_Step(is_global=True))
                    cursect = self.sections[-1]
                #
                directive_name = mo.group('directive_name')
                # newline should be kept in case of multi-line directive
                directive_value = mo.group('directive_value') + '\n'
                # is it an action??
                if directive_name in SOS_DIRECTIVES and directive_name != 'parameter':
                    if cursect is None:
                        parsing_errors.append(lineno, line,
                                              f'Directive {directive_name} is not allowed outside of a SoS step')
                        continue
                    cursect.add_directive(directive_name, directive_value, lineno)
                    if self.transcript:
                        self.transcript.write(f'DIRECTIVE\t{lineno}\t{line}')
                else:
                    if directive_name == 'parameter':
                        cursect.add_directive(directive_name, directive_value, lineno)
                        if self.transcript:
                            self.transcript.write(f'DIRECTIVE\t{lineno}\t{line}')
                    else:
                        # let us check if this is an acture action, or a type hint

                        cursect.add_script(directive_name, directive_value, lineno)
                        if self.transcript:
                            self.transcript.write(f'SCRIPT_{directive_name}\t{lineno}\t{line}')
                continue
            # if section is in script mode?
            if cursect and cursect.isValid() and cursect.category() == 'script':
                # if the script is indented and the line is not, the script
                # is ended.
                if not line[0].isspace() and cursect.indented_script():
                    cursect.wrap_script()
                    try:
                        cursect.wrap_script()
                    except Exception as e:
                        parsing_errors.append(lineno, line, e)
                else:
                    cursect.extend(line)
                    if self.transcript:
                        self.transcript.write(f'FOLLOW\t{lineno}\t{line}')
                    continue

            # all others?
            if not cursect:
                self.sections.append(SoS_Step(is_global=True))
                cursect = self.sections[-1]
                cursect.add_statement(line, lineno)
                if self.transcript:
                    self.transcript.write(f'STATEMENT\t{lineno}\t{line}')
                continue
            #
            if cursect.empty() or cursect.category() != 'statements':
                # new statement
                cursect.add_statement(line, lineno)
                if self.transcript:
                    self.transcript.write(f'STATEMENT\t{lineno}\t{line}')
            else:
                # existing one
                cursect.extend(line)
                if self.transcript:
                    # it is possible that we are switching from SCRIPT mode to regular
                    # statement, in which case the new mode should be STATEMENT, not FOLLOW
                    self.transcript.write(
                        f'{"STATEMENT" if cursect.category() == "statements" else "FOLLOW"}\t{lineno}\t{line}')
        #
        # check the last expression before a new directive
        if cursect:
            if not cursect.isValid():
                parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]),
                                      f'Invalid {cursect.category()}: {cursect.error_msg}')
            else:
                try:
                    cursect.finalize()
                except Exception as e:
                    parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), e)

        # non-matching %if ...
        if condition_met is not None:
            parsing_errors.append(lineno, '', 'Non-matching %if and %endif')

        #
        # if there is any parsing error, raise an exception
        #
        # if there is no section in the script, we create a default section with global
        # definition being the content.
        global_parameters = {}
        if not [x for x in self.sections if not x.is_global]:
            self.sections.append(SoS_Step(self.content, [('default', None, None)]))
            for section in [x for x in self.sections if x.is_global]:
                if self.sections[-1].task != '':
                    parsing_errors.append(cursect.lineno, 'Invalid section',
                                'Cannot define multiple default sections with a task in between.')
                self.sections[-1].statements.extend(section.statements)
                self.sections[-1].task = section.task
                self.global_def = ''
                global_parameters.update(section.parameters)
            # The sections should have been finalized so there is no need to finalize
            # again. In particular, finalizing a section would reset existing task #833
            #self.sections[-1].finalize()
        else:
            # as the last step, let us insert the global section to all sections
            for idx,sec in [(idx,x) for idx,x in enumerate(self.sections) if x.is_global]:
                for statement in sec.statements:
                    if statement[0] == ':':
                        parsing_errors.append(cursect.lineno, f'{statement[1]}:{statement[2]}',
                                'Global section cannot contain sos input, ouput, and task statements')
                    else:
                        self.global_def += statement[1]
                global_parameters.update(sec.parameters)
        # remove the global section after inserting it to each step of the process
        self.sections = [x for x in self.sections if not x.is_global]
        #
        if parsing_errors.errors:
            raise parsing_errors
        #
        for section in self.sections:
            # for nested / included sections, we need to keep their own global definition
            if '.' not in section.names[0][0]:
                section.global_def = self.global_def
                section.parameters.update(global_parameters)
            #
            section.tokens = section.get_tokens()
            section.md5 = textMD5(self.content.md5 + section.tokens)

    def workflow(self, workflow_name=None, use_default=True):
        '''Return a workflow with name_step+name_step specified in wf_name
        This function might be called recursively because of nested
        workflow.'''
        if workflow_name is None and not use_default:
            return SoS_Workflow(self.content, '', '',
                [section for section in self.sections if any(x in section.options for x in ('provides', 'shared', 'autoprovides'))], self.global_def)
        allowed_steps = None
        if not workflow_name:
            wf_name = ''
        else:
            # if consists of multiple workflows
            if '+' in workflow_name:
                wfs = []
                for wf in workflow_name.split('+'):
                    if not SOS_SUBWORKFLOW.match(wf):
                        raise ValueError(f'Incorrect workflow name {workflow_name}')
                    # if this is a combined workflow, extra_section might be specied.
                    wfs.append(self.workflow(wf))
                combined_wf = wfs[0]
                for wf in wfs[1:]:
                    combined_wf.extend(wf)
                combined_wf.name = workflow_name
                return combined_wf
            # if a single workflow
            # workflow_10:15 etc
            mo = SOS_SUBWORKFLOW.match(workflow_name)
            if not mo:
                raise ValueError(f'Incorrect workflow name {workflow_name}')
            wf_name, allowed_steps = mo.group('name', 'steps')
        # check source
        if not wf_name:
            if len(self.workflows) == 1:
                wf_name = list(self.workflows)[0]
            elif 'default' in self.workflows:
                wf_name = 'default'
            elif self.default_workflow:
                wf_name = self.default_workflow
            else:
                raise ValueError('Name of workflow should be specified because '
                    'the script defines more than one pipelines without a default one. '
                    'Available pipelines are: {}.'.format(', '.join(self.workflows)))
        elif wf_name not in self.workflows:
            raise ValueError(f'Workflow {wf_name} is undefined. Available workflows are: {", ".join(self.workflows)}')
        # do not send extra parameters of ...
        #sections = []
        #for section in self.sections:
        #    # skip, skip=True, skip=1 etc are all allowed.
        #    if 'provides' in section.options or 'shared' in section.options:
        #        # section global is shared by all workflows
        #        sections.append(section)
        #        continue
        #    for name, index, _ in section.names:
        #        # exact match or filename like match if name contains * etc
        #        if fnmatch.fnmatch(wf_name, name):
        #            sections.append(section)
        #            break
        return SoS_Workflow(self.content, wf_name, allowed_steps, self.sections, self.global_def)

    def print_help(self):
        '''print a help message from the script'''
        description = [x.lstrip('# ').strip() for x in self.description]
        try:
            # this does not work, but hopefully the author of mdv can make
            # the trivial change to make mdv compatible with python 3
            from mdv import main
            print(main(md=textwrap.dedent('\n'.join(description))))
        except Exception:
            print(textwrap.dedent('\n'.join(description)))
        #
        print('\nAvailable workflows')
        print('  ' + '\n  '.join(self.workflows))
        #
        print('\nSections')
        for section in self.sections:
            section.show()
