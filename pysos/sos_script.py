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
import sys
import re
import copy
import atexit
import glob
import fnmatch
import keyword
import argparse
import textwrap
#
# billiard is a Celery extension to multiprocessing
# 
# When we execute a step with subprocesses in parallel, the multiprocessing.pool
# module creates processes in daemon mode, which disallows the creation of new
# pools from these processes, which means we have to execute subworkflow
# separately from the step process. This forces us to get result from
# step.process and send to step.subworkflow, which can be really tedious
#
# The multiprocessing extension in Celery overcomes this problem and
# allows use to create new processes from a pooled process.
import billiard as mp
from billiard.pool import ApplyResult

from io import StringIO
from collections import OrderedDict, defaultdict
from collections.abc import Sequence, Iterable 
from itertools import tee, combinations

from . import __version__
from .utils import env, Error, WorkflowDict, SoS_eval, SoS_exec, RuntimeInfo, \
    dehtml, getTermWidth, interpolate, shortRepr, extract_pattern, expand_pattern, \
    print_traceback, pickleable, ProgressBar, frozendict

__all__ = ['SoS_Script']

#
# Note: we import actions here because we execute some actions in the main SoS 
# namespace. If we are moving all steps to their separate processes, this should
# be removed.
from .actions import *

# 
# global definitions of SoS syntax
# 
_INPUT_OPTIONS = ['group_by', 'skip', 'filetype', 'paired_with', 'for_each', 'pattern', 'dynamic']
_OUTPUT_OPTIONS = ['pattern', 'dynamic']
_DEPENDS_OPTIONS = ['pattern', 'dynamic']
_RUNTIME_OPTIONS = ['workdir', 'concurrent', 'docker', 'docker_options']


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


class StepInfo(object):
    '''A simple class to hold input, output, and index of step. Its attribute can
    only be set using an interface, and cannot be assigned. This is to make sure
    such information is not changed freely by malicious scripts'''
    def __init__(self):
        pass

    def set(self, key, value):
        object.__setattr__(self, key, value)

    def __setattr__(self, key, value):
        raise RuntimeError('Changing of step info {} is prohibited.'.format(key))

    def __repr__(self):
        return '{' + ', '.join('{}: {!r}'.format(x,y) for x,y in self.__dict__.items()) + '}'

def execute_step_process(step_process, global_process, subworkflow, sos_dict, sigil, signature, workdir):
    '''A function that has a local dictionary (from SoS env.sos_dict),
    a global process, and a step process. The processes are executed 
    in a separate process, independent of SoS. This makes it possible
    to execute the processes in background, on cluster, or submit as
    Celery tasks.'''
    env.register_process(os.getpid(), 'spawned_job with {} {}'
        .format(', '.join(sos_dict['_input']), ', '.join(sos_dict['_output'])))
    try:
        os.chdir(workdir)
        # switch context to the new dict and switch back once the with
        # statement ends (or if an exception is raised)
        env.sos_dict.quick_update(sos_dict)
        SoS_exec('import os, sys, glob')
        SoS_exec('from pysos import *')
        # global
        if global_process:
            SoS_exec(global_process, sigil)
        # step process
        SoS_exec(step_process, sigil)
        # subworkflow
        if subworkflow:
            subworkflow.run(nested=True)
        if signature:
            signature.write()
    except KeyboardInterrupt:
        raise RuntimeError('KeyboardInterrupt from {}'.format(os.getpid()))
    env.deregister_process(os.getpid())
    return 0

#
# Functions to handle directives
#

def handle_input_group_by(ifiles, group_by):
    '''Handle input option group_by'''
    if group_by == 'single':
        return [[x] for x in ifiles]
    elif group_by == 'all':
        # default option
        return [ifiles]
    elif group_by == 'pairs':
        if len(ifiles) % 2 != 0:
            raise ValueError('Paired group_by has to have even number of input files: {} provided'
                .format(len(ifiles)))
        return list(list(x) for x in zip(ifiles[:len(ifiles)//2], ifiles[len(ifiles)//2:]))
    elif group_by == 'pairwise':
        f1, f2 = tee(ifiles)
        next(f2, None)
        return [list(x) for x in zip(f1, f2)]
    elif group_by == 'combinations':
        return [list(x) for x in combinations(ifiles, 2)]


def handle_input_paired_with(paired_with, _groups, _vars):
    '''Handle input option paired_with'''
    if paired_with is None or not paired_with:
        paired_with = []
    elif isinstance(paired_with, str):
        paired_with = [paired_with]
    elif isinstance(paired_with, Iterable):
        paired_with = paired_with
    else:
        raise ValueError('Unacceptable value for parameter paired_with: {}'.format(paired_with))
    #
    ifiles = env.sos_dict['__step_input__']
    for wv in paired_with:
        if '.' in wv:
            if wv.split('.')[0] not in env.sos_dict:
                raise ValueError('Variable {} does not exist.'.format(wv))
            values = getattr(env.sos_dict[wv.split('.')[0]], wv.split('.', 1)[-1])
        else:
            if wv not in env.sos_dict:
                raise ValueError('Variable {} does not exist.'.format(wv))
            values = env.sos_dict[wv]
        if isinstance(values, str) or not isinstance(values, Iterable):
            raise ValueError('with_var variable {} is not a sequence ("{}")'.format(wv, values))
        if len(values) != len(ifiles):
            raise ValueError('Length of variable {} (length {}) should match the number of input files (length {}).'
                .format(wv, len(values), len(ifiles)))
        file_map = {x:y for x,y in zip(ifiles, values)}
        for idx, grp in enumerate(_groups):
            _vars[idx]['_' + wv.split('.')[0]] = [file_map[x] for x in grp]


def handle_input_pattern(pattern, ifiles, _groups, _vars):
    '''Handle input option pattern'''
    if pattern is None or not pattern:
        patterns = []
    elif isinstance(pattern, str):
        patterns = [pattern]
    elif isinstance(pattern, Iterable):
        patterns = pattern
    else:
        raise ValueError('Unacceptable value for parameter pattern: {}'.format(pattern))
    #
    for pattern in patterns:
        res = extract_pattern(pattern, ifiles)
        # now, assign the variables to env
        for k, v in res.items():
            if k in ('input', 'output', 'depends') or k.startswith('_'):
                raise RuntimeError('Pattern defined variable {} is not allowed'.format(k))
            env.sos_dict[k] = v
        # also make k, v pair with _input
        handle_input_paired_with(res.keys(), _groups, _vars)
    

def handle_input_for_each(for_each, _groups, _vars):
    if for_each is None or not for_each:
        for_each = []
    elif isinstance(for_each, str):
        for_each = [for_each]
    elif isinstance(for_each, Sequence):
        for_each = for_each
    else:
        raise ValueError('Unacceptable value for parameter for_each: {}'.format(for_each))
    #
    for fe_all in for_each:
        loop_size = None
        for fe in fe_all.split(','):
            values = env.sos_dict[fe]
            if not isinstance(values, Sequence):
                raise ValueError('for_each variable {} is not a sequence ("{}")'.format(fe, values))
            if loop_size is None:
                loop_size = len(values)
            elif loop_size != len(values):
                raise ValueError('Length of variable {} (length {}) should match the length of variable {} (length {}).'
                    .format(fe, len(values), fe_all.split(',')[0], loop_size))
        # expand
        _tmp_groups = copy.deepcopy(_groups)
        _groups.clear()
        for i in range(loop_size):
            _groups.extend(_tmp_groups)
        #
        _tmp_vars = copy.deepcopy(_vars)
        _vars.clear()
        for vidx in range(loop_size):
            for idx in range(len(_tmp_vars)):
                for fe in fe_all.split(','):
                    if '.' in fe:
                        if fe.split('.')[0] not in env.sos_dict:
                            raise ValueError('Variable {} does not exist.'.format(fe))
                        _tmp_vars[idx]['_' + fe.split('.')[0]] = getattr(env.sos_dict[fe.split('.')[0]], fe.split('.', 1)[-1])[vidx]
                    else:
                        if fe not in env.sos_dict:
                            raise ValueError('Variable {} does not exist.'.format(fe))
                        _tmp_vars[idx]['_' + fe] = env.sos_dict[fe][vidx]
            _vars.extend(copy.deepcopy(_tmp_vars))

def __null_func__(*args, **kwargs):
    '''This is a utility function for the parser'''
    return args, kwargs

# directive input
def directive_input(*args, **kwargs):
    # first *args are filenames
    for k in kwargs.keys():
        if k not in _INPUT_OPTIONS:
            raise RuntimeError('Unrecognized input option {}'.format(k))
    #
    if args:
        ifiles = []
        for arg in args:
            if isinstance(arg, str):
                ifiles.append(arg)
            elif isinstance(arg, Iterable):
                # in case arg is a Generator, check its type will exhaust it
                arg = list(arg)
                if not all(isinstance(x, str) for x in arg):
                    raise RuntimeError('Invalid input file: {}'.format(arg))
                ifiles.extend(arg)
            else:
                raise ValueError('Unrecognizable input type {}'.format(arg))
        env.sos_dict.set('__step_input__', ifiles)
    else:
        ifiles = env.sos_dict['__step_input__']
    # expand files with wildcard characters and check if files exist
    tmp = []
    for ifile in ifiles:
        if os.path.isfile(os.path.expanduser(ifile)):
            tmp.append(ifile)
        elif env.run_mode == 'run':
            # in this mode file must exist
            expanded = glob.glob(os.path.expanduser(ifile))
            if not expanded:
                raise RuntimeError('{} not exist.'.format(ifile))
            tmp.extend(expanded)
        elif env.run_mode == 'dryrun':
            # FIXME: this should be the 'dynamic' mode
            expanded = glob.glob(os.path.expanduser(ifile))
            if expanded:
                tmp.extend(expanded)
            else:
                tmp.append(ifile)
    #
    ifiles = tmp
    #
    if 'filetype' in kwargs:
        if isinstance(kwargs['filetype'], str):
            ifiles = fnmatch.filter(ifiles, kwargs['filetype'])
        elif isinstance(kwargs['filetype'], Iterable):
            ifiles = [x for x in ifiles if any(fnmatch.fnmatch(x, y) for y in kwargs['filetype'])]
        elif callable(kwargs['filetype']):
            ifiles = [x for x in ifiles if kwargs['filetype'](x)]
    #             
    # handle group_by
    if 'group_by' in kwargs:
        _groups = handle_input_group_by(ifiles, kwargs['group_by'])
    else:
        _groups = [ifiles]
    #
    _vars = [{} for x in _groups]
    # handle paired_with
    if 'paired_with' in kwargs:
        handle_input_paired_with(kwargs['paired_with'], _groups, _vars)
    # handle pattern
    if 'pattern' in kwargs:
        handle_input_pattern(kwargs['pattern'], ifiles, _groups, _vars)
    # handle for_each
    if 'for_each' in kwargs:
        handle_input_for_each(kwargs['for_each'], _groups, _vars)
    #
    if 'skip' in kwargs:
        if callable(kwargs['skip']):
            _tmp_groups = []
            _tmp_vars = []
            for g, v in zip(_groups, _vars):
                try:
                    res = kwargs['skip'](g, **v)
                except Exception as e:
                    raise RuntimeError('Failed to apply skip function to group {}: {}'.format(', '.join(g), e))
                if res:
                   env.logger.info('Input group {} is skipped.'.format(', '.join(g)))
                else:
                    _tmp_groups.append(g)
                    _tmp_vars.append(v)
            _groups = _tmp_groups
            _vars = _tmp_vars
        else:
            _groups = []
            _vars = []
    return _groups, _vars


def directive_depends(*args, **kwargs):
    '''handle directive depends'''
    for k in kwargs.keys():
        if k not in _DEPENDS_OPTIONS:
            raise RuntimeError('Unrecognized depends option {}'.format(k))
    # first *args are filenames
    dfiles = []
    for arg in args:
        if isinstance(arg, str):
            dfiles.append(arg)
        elif isinstance(arg, Iterable):
            arg = list(arg)
            if not all(isinstance(x, str) for x in arg):
                raise RuntimeError('Invalid dependent file: {}'.format(arg))
            dfiles.extend(arg)
        else:
            raise ValueError('Unrecognizable output type {}'.format(arg))

    env.sos_dict.set('_depends', dfiles)

def handle_output_pattern(pattern, ofiles):
    # 
    if pattern is None or not pattern:
        patterns = []
    elif isinstance(pattern, str):
        patterns = [pattern]
    elif isinstance(pattern, Iterable):
        patterns = pattern
    else:
        raise ValueError('Unacceptable value for parameter pattern: {}'.format(pattern))
    #
    for pattern in patterns:
        ofiles.extend(expand_pattern(pattern))
        
def directive_output(*args, **kwargs):
    for k in kwargs.keys():
        if k not in _OUTPUT_OPTIONS:
            raise RuntimeError('Unrecognized output option {}'.format(k))
    ofiles = []
    for arg in args:
        if isinstance(arg, str):
            ofiles.append(arg)
        elif isinstance(arg, Iterable):
            arg = list(arg)
            if not all(isinstance(x, str) for x in arg):
                raise RuntimeError('Invalid output file: {}'.format(arg))
            ofiles.extend(arg)
        else:
            raise ValueError('Unrecognizable output type {}'.format(arg))
    #
    if 'pattern' in kwargs:
        handle_output_pattern(kwargs['pattern'], ofiles)
    for ofile in ofiles:
        parent_dir = os.path.split(os.path.expanduser(ofile))[0]
        if parent_dir and not os.path.isdir(parent_dir):
            os.makedirs(parent_dir)
    env.sos_dict.set('_output', ofiles)

def directive_process(**kwargs):
    for k in kwargs.keys():
        if k not in _RUNTIME_OPTIONS:
            raise RuntimeError('Unrecognized runtime option {}'.format(k))
    #
    if '_runtime' not in env.sos_dict:
        env.sos_dict.set('_runtime', kwargs)

#
# handling of parameters step
#



class SoS_Step:
    #
    # A single sos step
    #
    def __init__(self, names=[], options={}, is_global=False, is_parameters=False):
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
        self.global_process = ''
        self.process = ''
        # subworkflow of a step
        self.subworkflow = None
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
        self._script = ''

    def category(self):
        if self.statements:
            if self.statements[-1][0] == '=':
                return 'expression'
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

    #
    # Parsing input
    #
    def empty(self):
        '''If there is no content (comment does not count)'''
        return self.category() is None

    def extend(self, line):
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

    def add_directive(self, key, value, lineno=None, action=None):
        '''Assignments are items with ':' type '''
        if key is None:
            # continuation of multi-line directive
            self.statements[-1][-1] += value
            self.values.append(value)
        else:
            # new directive
            self.statements.append([':', key, value])
            self.values = [value]
            if action is not None:
                self._action = action
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
        '''convert action: script to process: action(script)'''
        if self._action is None:
            return
        self.statements.append(['!', '{}({!r})'.format(self._action, textwrap.dedent(self._script))])
        self._action = None
        self._script = ''

    def finalize(self):
        ''' split statement and process by last directive '''
        self.wrap_script()
        if not self.statements:
            self.process = ''
            return
        process_directive = [idx for idx, statement in enumerate(self.statements) if statement[0] == ':' and statement[1] == 'process']
        if not process_directive:
            self.process = ''
            return
        start_process = process_directive[0] + 1
        # convert statement to process
        self.process = ''
        for statement in self.statements[start_process:]:
            if statement[0] == '=':
                self.process += '{} = {}'.format(statement[1], statement[2])
            elif statement[0] == ':':
                if statement[1] in ('input', 'output', 'depends'):
                    raise ValueError('Step process should be defined as the last item in a SoS step')
                elif statement[2].strip():
                    raise ValueError('Runtime options are not allowed for second or more actions in a SoS step')
                # ignore ...
                self.process += '\n'
            else:
                self.process += statement[1]
        # remove process from self.statement
        if process_directive:
            self.statements = self.statements[:start_process]

    def isValid(self):
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
                    return False
                compile('func(' + ''.join(self.values) + ')', filename='<string>', mode='eval')
            elif self.category() == 'statements':
                compile(''.join(self.values), filename='<string>', mode='exec')
            elif self.category() == 'script':
                return True
            else:
                raise RuntimeError('Unrecognized expression type {}'.format(self.category()))
            return True
        except Exception as e:
            return False

    def _parse_error(self, msg):
        '''This function will replace error() function in argparse module so that SoS
        can hijack errors raised from it.'''
        raise ArgumentError(msg)

    def parse_args(self, args, check_unused=False, cmd_name=''):
        '''Parse command line arguments and set values to parameters section'''
        env.logger.info('Execute ``{}_parameters``'.format(self.name))
        env.sos_dict.set('step_name', '{}_parameters'.format(self.name))
        if self.global_process:
            try:
                SoS_exec(self.global_process)
            except Exception as e:
                if env.verbosity > 2:
                    print_traceback()
                raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(self.global_process, e))
     
        def str2bool(v):
            if v.lower() in ('yes', 'true', 't', '1'):
                return True
            elif v.lower() in ('no', 'false', 'f', '0'):
                return False
            else:
                raise ArgumentError('Invalid value for bool argument "{}" (only yes,no,true,false,t,f,0,1 are allowed)'.format(v))
        #
        parser = argparse.ArgumentParser(prog='sos-runner {}'.format(cmd_name))
        parser.register('type', 'bool', str2bool)
        arguments = {}
        for key, defvalue, comment in self.parameters:
            try:
                defvalue = SoS_eval(defvalue, self.sigil)
                arguments[key] = defvalue
            except Exception as e:
                raise RuntimeError('Incorrect default value {} for parameter {}: {}'.format(defvalue, key, e))
            if isinstance(defvalue, type):
                if defvalue == bool:
                    parser.add_argument('--{}'.format(key), type='bool', help=comment, required=True, nargs='?') 
                else:
                    # if only a type is specified, it is a required document of required type
                    parser.add_argument('--{}'.format(key), type=str if hasattr(defvalue, '__iter__') else defvalue,
                        help=comment, required=True, nargs='+' if hasattr(defvalue, '__iter__') else '?')
            else:
                if isinstance(defvalue, bool):
                    parser.add_argument('--{}'.format(key), type='bool', help=comment,
                        nargs='?', default=defvalue)
                else:
                    if isinstance(defvalue, str):
                        deftype = str
                    elif isinstance(defvalue, Sequence):
                        if len(defvalue) > 0:
                            deftype = type(defvalue[0])
                        else:
                            deftype = str
                    else:
                        deftype = type(defvalue)
                    parser.add_argument('--{}'.format(key), type=deftype, help=comment,
                        nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, str) else '?',
                        default=defvalue)
        #
        parser.error = self._parse_error
        # 
        # because of the complexity of having combined and nested workflows, we cannot know how
        # many parameters section a workflow has and therfore have to assume that the unknown parameters
        # are for other sections.
        if check_unused:
            parsed = parser.parse_args(args)
        else:
            parsed, unknown = parser.parse_known_args(args)
            if unknown:
                env.logger.warning('Unparsed arguments [{}] that might be processed by another combined or nested workflow'
                    .format(' '.join(unknown)))
        #
        arguments.update(vars(parsed))
        # now change the value with passed values
        for k, v in arguments.items():
            env.sos_dict[k] = v
            # protect variables from being further modified
            env.readonly_vars.add(k)

    #
    # Execution
    #
    def step_signature(self):
        '''return everything that might affect the execution of the step 
        namely, global process, step definition etc to create a unique 
        signature that might will be changed with the change of SoS script.'''
        result = self.global_process
        for statement in self.statements:
            if statement[0] in (':', '='):
                result += '{}: {}\n'.format(statement[1], statement[2])
            else:
                result += statement[1] + '\n'
        result += self.process
        return re.sub(r'\s+', ' ', result)

    def run_with_queue(self, queue):
        '''Execute the step in a separate process and return the results through a 
        queue '''
        try:
            res = self.run()
            # shared variables will be sent from subprocess to the master
            for key in env.shared_vars:
                if key in env.sos_dict and key not in res:
                    res[key] = env.sos_dict[key]
            queue.put(res)
        except Exception as e:
            queue.put(e)

    def run(self):
        '''Execute a single step and return results '''
        # only results will be sent back to the master process
        result = {'__step_output__': []}
        # handle these two sections differently
        env.logger.info('Execute ``{}_{}``: {}'.format(self.name, self.index, self.comment.strip()))
        env.sos_dict.set('step_name', '{}_{}'.format(self.name, self.index))
        # 
        # the following is a quick hack to allow directive_input function etc to access 
        # the workflow dictionary
        #
        # these are temporary variables that should be removed if exist
        for var in ('input', 'output', 'depends', '_input', '_depends', '_output'):
            env.sos_dict.pop(var, '')
        #
        # default input groups and vars, might be reset by directive input
        if '__step_input__' in env.sos_dict:
            self._groups = [env.sos_dict['__step_input__']]
        else:
            self._groups = [[]]
        self._vars = [{}]
        #
        input_idx = [idx for idx,x in enumerate(self.statements) if x[0] == ':' and x[1] == 'input']
        if not input_idx:
            input_idx = None
        elif len(input_idx) == 1:
            input_idx = input_idx[0]
        else:
            raise RuntimeError('Only one step input is allowed')
        #
        # execute global process
        if self.global_process:
            try:
                SoS_exec(self.global_process)
            except Exception as e:
                if env.verbosity > 2:
                    print_traceback()
                raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(self.global_process, e))
        #
        public_vars = set()
        if input_idx is not None:
            # execute before input stuff
            for statement in self.statements[:input_idx]:
                if statement[0] == '=':
                    key, value = statement[1:]
                    public_vars.add(key)
                    try:
                        env.sos_dict[key] = SoS_eval(value, self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
                elif statement[0] == ':':
                    raise RuntimeError('Step input should be specified before others')
                else:
                    try:
                        SoS_exec(statement[1], self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to process statement {}: {}'.format(statement[1], e))
            # input
            key, value = self.statements[input_idx][1:]
            try:
                args, kwargs = SoS_eval('__null_func__({})'.format(value), self.sigil)
                self._groups, self._vars = directive_input(*args, **kwargs)
            except Exception as e:
                raise RuntimeError('Failed to process step {} : {} ({})'.format(key, value.strip(), e))
            input_idx += 1
        else:
            # assuming everything starts from 0 is after input
            input_idx = 0
        
        env.sos_dict.set('input', list(OrderedDict.fromkeys(sum(self._groups, []))))
        step_info = StepInfo()
        step_info.set('step_name', env.sos_dict['step_name'])
        step_info.set('input', env.sos_dict['input'])
        if 'alias' in self.options:
            # the step might be skipped 
            for var in public_vars:
                # there is a slight possibility that var is deleted
                if var in env.sos_dict and pickleable(env.sos_dict[var]):
                    step_info.set(var, env.sos_dict[var])
            result[self.options['alias']] = copy.deepcopy(step_info)
        if not self._groups:
            env.logger.info('Step {} is skipped'.format(self.index))
            return result
        # post input
        # output and depends can be processed many times
        step_sig = self.step_signature()
        env.logger.info('input:   ``{}``'.format(shortRepr(env.sos_dict['input'])))
        self._outputs = []
        self._depends = []
        for idx, (g, v) in enumerate(zip(self._groups, self._vars)):
            # other variables
            env.sos_dict.update(v)
            env.sos_dict.set('_input', g)
            env.sos_dict.set('_index', idx)
            for key in ('_output', '_depends'):
                if key in env.sos_dict:
                    env.sos_dict.pop(key)
            for statement in self.statements[input_idx:]:
                if statement[0] == '=':
                    key, value = statement[1:]
                    public_vars.add(key)
                    try:
                        env.sos_dict[key] = SoS_eval(value, self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
                elif statement[0] == ':':
                    key, value = statement[1:]
                    # output, depends, and process can be processed multiple times
                    try:
                        args, kwargs = SoS_eval('__null_func__({})'.format(value), self.sigil)
                        eval('directive_' + key)(*args, **kwargs)
                    except Exception as e:
                        raise RuntimeError('Failed to process step {}: {} ({})'.format(key, value.strip(), e))
                else:
                    old_run_mode = env.run_mode
                    if '_output' in env.sos_dict:
                        signature = RuntimeInfo(step_sig, env.sos_dict['_input'], env.sos_dict['_output'],
                            env.sos_dict.get('_depends', []))
                        if env.sig_mode == 'default' and signature.validate():
                            env.logger.info('Execute statement in dryrun mode and reuse existing output files {}'.format(', '.join(env.sos_dict['_output'])))
                            env.run_mode = 'dryrun'
                    try:
                        SoS_exec(statement[1], self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to process statement {}: {}'.format(statement[1], e))
                    env.run_mode = old_run_mode

            # collect _output and _depends
            if '_output' in env.sos_dict:
                self._outputs.append(env.sos_dict['_output'])
            else:
                self._outputs.append([])
            if '_depends' in env.sos_dict:
                self._depends.append(env.sos_dict['_depends'])
            else:
                self._depends.append([])
        #
        # if no output directive, assuming no output for each step
        if not self._outputs:
            self._outputs = [[] for x in self._groups]
        # if no depends directive, assuming no dependent files for each step
        if not self._depends:
            self._depends = [[] for x in self._groups]
        # we need to reduce output files in case they have been processed multiple times.
        env.sos_dict.set('output', list(OrderedDict.fromkeys(sum(self._outputs, []))))
        env.sos_dict.set('depends', list(OrderedDict.fromkeys(sum(self._depends, []))))
        env.logger.info('output:  ``{}``'.format(shortRepr(env.sos_dict['output'])))
        if env.sos_dict['depends']:
            env.logger.info('depends: ``{}``'.format(shortRepr(env.sos_dict['depends'])))
        result['__step_output__'] = env.sos_dict['output']
        if 'alias' in self.options:
            step_info.set('output', env.sos_dict['output'])
            step_info.set('depends', env.sos_dict['depends'])
            for var in public_vars:
                if var in env.sos_dict and pickleable(env.sos_dict[var]):
                    step_info.set(var, env.sos_dict[var])
            result[self.options['alias']] = copy.deepcopy(step_info)
        #
        # if the signature matches, the whole step is ignored, including subworkflows
        if env.sos_dict['output']:
            signature = RuntimeInfo(step_sig, 
                env.sos_dict['input'], env.sos_dict['output'], env.sos_dict['depends'])
            if env.run_mode == 'run':
                if env.sig_mode == 'default':
                    if signature.validate():
                        # everything matches
                        env.logger.info('Reusing existing output files {}'.format(', '.join(env.sos_dict['output'])))
                        return result
                elif env.sig_mode == 'assert':
                    if not signature.validate():
                        raise RuntimeError('Signature mismatch.')
                # elif env.sig_mode == 'ignore'
        else:
            signature = None
        #
        proc_results = []
        signatures = []
        if '_runtime' in env.sos_dict:
            self.runtime_options = env.sos_dict['_runtime']
        concurrent = env.max_jobs > 1 and len(self._groups) > 1 and 'concurrent' in self.runtime_options and self.runtime_options['concurrent']
        if concurrent:
            pool = mp.Pool(min(env.max_jobs, len(self._groups)))
        for idx, (g, v, o, d) in enumerate(zip(self._groups, self._vars, self._outputs, self._depends)):
            env.sos_dict.update(v)
            env.sos_dict.set('_input', g)
            env.sos_dict.set('_output', o)
            env.sos_dict.set('_depends', d)
            env.sos_dict.set('_index', idx)
            env.logger.debug('_idx: ``{}``'.format(idx))
            env.logger.debug('_input: ``{}``'.format(shortRepr(env.sos_dict['_input'])))
            env.logger.debug('_output: ``{}``'.format(shortRepr(env.sos_dict['_output'])))
            #
            # action
            # If the users specifies output files for each loop (using ${input} etc, we
            # can try to see if we can create partial signature. This would help if the
            # step is interrupted in the middle.
            partial_signature = None
            if env.sos_dict['_output'] and env.sos_dict['_output'] != env.sos_dict['output'] and env.run_mode == 'run':
                partial_signature = RuntimeInfo(step_sig, env.sos_dict['_input'], env.sos_dict['_output'], 
                    env.sos_dict['_depends'])
                if env.sig_mode == 'default':
                    if partial_signature.validate():
                        # everything matches
                        env.logger.info('Reusing existing output files {}'.format(', '.join(env.sos_dict['_output'])))
                        continue
                elif env.sig_mode == 'assert':
                    if not partial_signature.validate():
                        raise RuntimeError('Signature mismatch for input {} and output {}'.format(
                            ', '.join(env.sos_dict['_input']), ', '.join(env.sos_dict['_output'])))
            # now, if output file has already been generated using non-process statement
            # so that no process need to be run, or if the output need help from a subworkflow
            # we create signature from outside.
            if (not self.process and not self.subworkflow):
                if partial_signature is not None:
                    partial_signature.write()
                continue
            #
            try:
                if concurrent:
                    proc_results.append(pool.apply_async(
                        execute_step_process,   # function
                        (self.process,          # process
                        self.global_process,    # global process
                        self.subworkflow,       # subworkflow
                        env.sos_dict.clone_pickleable(),
                        self.sigil,
                        # if subworkflow contribute to outcome, we can not save signature here
                        partial_signature,
                        self.runtime_options['workdir'] if 'workdir' in self.runtime_options else os.getcwd())))
                else:
                    # execute in existing process
                    proc_results.append(
                        execute_step_process(   # function
                        self.process,           # process
                        '',                     # local execusion, no need to re-run global
                        self.subworkflow,       # subworkflow
                        # do not clone dict
                        env.sos_dict,
                        self.sigil,
                        partial_signature,
                        self.runtime_options['workdir'] if 'workdir' in self.runtime_options else os.getcwd()))
            except Exception as e:
                # FIXME: cannot catch exception from subprocesses
                if env.verbosity > 2:
                    print_traceback()
                raise RuntimeError('Failed to execute process\n"{}"\n{}'.format(self.process, e))
        # check results? This is only meaningful for pool
        if concurrent:
            try:
                proc_results = [res.get() if isinstance(res, ApplyResult) else res for res in proc_results]
            except KeyboardInterrupt:
                # if keyboard interrupt
                pool.terminate()
                pool.join()
                raise RuntimeError('KeyboardInterrupt fro m {} (master)'.format(os.getpid()))
            except Exception as e:
                # if keyboard interrupt etc
                env.logger.error('Caught {}'.format(e))
                pool.terminate()
                pool.join()
                raise
        if not all(x==0 for x in proc_results):
            raise RuntimeError('Step process returns non-zero value')
        if env.run_mode == 'run':
            for ofile in env.sos_dict['output']:
                if not os.path.isfile(os.path.expanduser(ofile)): 
                    raise RuntimeError('Output file {} does not exist after completion of action'.format(ofile))
        if signature and env.run_mode == 'run':
            signature.write()
        if concurrent:
            # finally, write results back to the master process
            pool.close()
            pool.join()
        return result

    def show(self, indent):
        '''Output for command sos show'''
        textWidth = max(60, getTermWidth())
        if self.is_parameters:
            print(indent + '  Parameters:')
            for k,v,c in self.parameters:
                 text = '{:<18}'.format(k) + (c + ' ' if c else '') + \
                     ('(default: {})'.format(v) if v else '')
                 print('\n'.join(
                     textwrap.wrap(text, 
                     initial_indent=indent + ' '*4,
                     subsequent_indent=indent + ' '*22,
                     width=textWidth)
                     ))
        else:
            # hide a step if there is no comment
            text = '  {:<20}'.format('Step {}_{}:'.format(self.name, self.index)) + self.comment
            print('\n'.join(
                textwrap.wrap(text, 
                    width=textWidth, 
                    initial_indent=indent,
                    subsequent_indent=indent + ' '*22)))
            if self.subworkflow:
                self.subworkflow.show(indent+ '   ' + chr(124), nested=True)
          

class SoS_Workflow:
    #
    # A SoS workflow with multiple steps
    #
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
                # parameters and global section will not have subworkflow
                self.sections[-1].subworkflow = None
                continue
            for name, index, subworkflow in section.names:
                if 'target' in section.options:
                    self.auxillary_sections.append(section)
                elif fnmatch.fnmatch(workflow_name, name):
                    self.sections.append(copy.deepcopy(section))
                    self.sections[-1].name = workflow_name
                    self.sections[-1].index = 0 if index is None else int(index)
                    self.sections[-1].subworkflow = subworkflow
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
        for section in self.sections:
            if section.subworkflow is not None and isinstance(section.subworkflow, str):
                raise RuntimeError('Subworkflow {} not executable most likely because of recursice expansion.'.format(section.subworkflow))
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

    def run(self, args=[], nested=False, cmd_name='', config_file=None):
        '''Execute a workflow with specified command line args. If sub is True, this 
        workflow is a nested workflow and be treated slightly differently.
        '''
        if nested:
            # if this is a subworkflow, we use _input as step input the workflow
            env.sos_dict.set('__step_input__', env.sos_dict['_input'])
        else:
            # Because this workflow might belong to a combined workflow, we do not clear
            # locals before the execution of workflow.
            # Need to choose what to inject to globals
            env.sos_dict = WorkflowDict()
            #
            # inject a few things
            env.sos_dict.set('__null_func__', __null_func__)
            for k, v in globals().items():
                # if this is an action inject
                if callable(v) and hasattr(v, 'run_mode'):
                    env.sos_dict.set(k, v)
                elif k in ('os', 'globl', 'sys', 'SoS_Action', 'interpolate'):
                    env.sos_dict.set(k, v)
            # initial values
            env.sos_dict.set('SOS_VERSION', __version__)
            cfg = {}
            if config_file is not None:
                if not os.path.isfile(config_file):
                    raise RuntimeError('Config file {} not found'.format(config_file))
                try:
                    import yaml
                    with open(config_file) as config:
                        cfg = yaml.safe_load(config)
                except Exception as e:
                    raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format?'.format(config_file))
            #
            env.sos_dict.set('CONFIG', frozendict(cfg))
            py_version = sys.version_info
            env.sos_dict.set('__step_input__', [])
        #
        # process step of the pipelinp
        #
        num_parameters_sections = len([x for x in self.sections if x.is_parameters or x.subworkflow])
        if num_parameters_sections == 0 and args:
            raise ArgumentError('Unused parameter {}'.format(' '.join(args)))
        #
        # the steps can be executed in the pool (Not implemented)
        # if nested = true, start a new progress bar
        prog = ProgressBar(self.name, len(self.sections))
        for idx, section in enumerate(self.sections):
            # global section will not change _step etc
            if section.is_parameters:
                # if there is only one parameters section and no nested workflow, check unused section
                section.parse_args(args, num_parameters_sections == 1, cmd_name=cmd_name)
                prog.progress(1)
                continue
            # 
            # execute section with specified input
            # 1. for first step of workflow, _step.input=[]
            # 2. for subworkflow, _step.input = _input
            # 3. for second to later step, _step.input = _step.output
            # each section can use a separate process
            queue = mp.Queue()
            proc = mp.Process(target=section.run_with_queue,
                args=(queue,))
            proc.start()
            proc.join()
            res = queue.get()
            # if the job is failed
            if isinstance(res, Exception):
                # error must have been displayed.
                raise RuntimeError(res)
            #res = section.run()
            for k, v in res.items():
                if k == '__step_output__':
                    env.sos_dict.set('__step_input__', v)
                env.sos_dict.set(k, v)
            prog.progress(1)
        prog.done()

    def show(self, indent = '', nested=False):
        textWidth = max(60, getTermWidth())
        paragraphs = dehtml(self.description).split('\n\n')
        print('\n'.join(
            textwrap.wrap('{} {}:  {}'.format(
                'Nested workflow' if nested else 'Workflow', 
                self.name, paragraphs[0]),
                initial_indent = indent,
                subsequent_indent = indent,
                width=textWidth)
            ))
        for paragraph in paragraphs[1:]:
            print('\n'.join(
            textwrap.wrap(paragraph, width=textWidth,
                initial_indent = indent,
                subsequent_indent = indent)
            ))
        for section in self.sections:
            section.show(indent)

class SoS_Script:
    _DIRECTIVES = ['input', 'output', 'depends', 'process']
    _SECTION_OPTIONS = ['alias', 'skip', 'sigil', 'target', 'source']
    _PARAMETERS_SECTION_NAME = 'parameters'

    # Regular expressions for parsing section headers and options
    _SECTION_HEADER_TMPL = r'''
        ^\[\s*                             # [
        (?P<section_name>[\d\w_,+=*\s-]+)  # digit, alphabet, _ and ,
        (:\s*                              # :
        (?P<section_option>.*)             # section options
        )?                                 # optional
        \]\s*$                             # ]
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
        (\s*=\s*
        (?P<subworkflow>                   # = subworkflow
        [\w\d_+\s-]+                       # subworkflow specification
        ))?                                # optional
        \s*$
        '''

    _SUBWORKFLOW_TMPL = '''
        ^\s*                               # leading space
        (?P<name>                          # name
        [a-zA-Z*]                          # cannot start with _ etc
        ([\w\d_]*?))                       # can have _ and digit
        (_(?P<steps>                       # index start from _
        [\d\s-]+))?                        # with - and digit
        \s*$                               # end
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
        (?P<directive_name>                # 
        (?!({})\s*:)                       # not a python keyword followed by : (can be input)
        ({}                                # name of directive
        |[a-zA-Z][\w\d_]*))                #    or action
        \s*:\s*                            # followed by :
        (?P<directive_value>.*)            # and values
        '''.format('|'.join(keyword.kwlist), '|'.join(_DIRECTIVES))

    _ASSIGNMENT_TMPL = r'''
        ^                                  # from start of line
        (?P<var_name>[\w_][\d\w_]*)        # variable name
        \s*=\s*                            # assignment
        (?P<var_value>.*)                  # variable content
        '''

    SECTION_HEADER = re.compile(_SECTION_HEADER_TMPL, re.VERBOSE)
    SECTION_NAME = re.compile(_SECTION_NAME_TMPL, re.VERBOSE)
    SUBWORKFLOW = re.compile(_SUBWORKFLOW_TMPL, re.VERBOSE)
    SECTION_OPTION = re.compile(_SECTION_OPTION_TMPL, re.VERBOSE)
    FORMAT_LINE = re.compile(_FORMAT_LINE_TMPL, re.VERBOSE)
    FORMAT_VERSION = re.compile(_FORMAT_VERSION_TMPL, re.VERBOSE)
    DIRECTIVE = re.compile(_DIRECTIVE_TMPL, re.VERBOSE)
    ASSIGNMENT = re.compile(_ASSIGNMENT_TMPL, re.VERBOSE)

    def __init__(self, content='', filename=None):
        '''Parse a sectioned SoS script file. Please refer to the SoS manual
        for detailed specification of this format.

        Parameter `content` can be either a filename or a content of a
        SoS script in unicode, which is convenient for passing scripts for
        testing purposes.

        Parameter `filename` should be used if the content should be read 
        from a file.
        '''
        if filename:
            self.sos_script = os.path.expanduser(filename)
            if not os.path.isfile(self.sos_script):
                raise ValueError('{} does not exist'.format(filename))
            with open(self.sos_script) as fp:
                self._read(fp)
        else:
            self.sos_script = '<string>'
            with StringIO(content) as fp:
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
        parsing_errors = ParsingError(self.sos_script)
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
                        self.descriptions[-1] += line.lstrip('#').lstrip()
                else:
                    # in the parameter section, the comments are description
                    # of parameters and are all significant
                    if cursect.isValid():
                        cursect.add_comment(line)
                    # comment add to script
                    else:
                        cursect.extend(line)
                continue
            elif not line.strip():
                # a blank line start a new comment block if we are still
                # in the front of the script
                if cursect is None:
                    comment_block += 1
                    self.descriptions.append('')
                else:
                    if cursect.category() == 'statements':
                        cursect.extend(line)
                    elif cursect.comment:
                        comment_block += 1
                continue
            #
            # a continuation of previous item?
            if line[0].isspace() and cursect is not None and not cursect.empty():
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
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid ' + cursect.category())
                    cursect.values = []
                    cursect.finalize()
                # start a new section
                section_name = mo.group('section_name').strip()
                section_option = mo.group('section_option')
                step_names = []
                step_options = {}
                for name in section_name.split(','):
                    mo = self.SECTION_NAME.match(name)
                    if mo:
                        n, i, di, s = mo.group('name', 'index', 'default_index', 'subworkflow')
                        if n:
                            if i is None and '*' in n:
                                parsing_errors.append(lineno, line, 'Unindexed section name cannot contain wildcard character (*).')
                            step_names.append((n, i, s))
                        if di:
                            step_names.append(('default', di, s))
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
                        mo = self.SECTION_OPTION.match(option)
                        if mo:
                            opt_name, opt_value = mo.group('name', 'value')
                            #
                            #opt_value should also be a valid Python expression (with quote etc)
                            #which is most likely a string.
                            if opt_value:
                                try:
                                    opt_value = eval(opt_value)
                                except Exception as e:
                                    parsing_errors.append(lineno, line, e)
                            if opt_name == 'sigil':
                                if opt_value.count(' ') != 1 or opt_value[0] in (' ', "'") or \
                                    opt_value[-1] in (' ', "'") or \
                                    opt_value.split(' ')[0] == opt_value.split(' ')[1]:
                                    parsing_errors.append(lineno, line, 'Incorrect sigil "{}"'.format(opt_value))
                            if opt_name == 'source':
                                opt_value = [opt_value] if isinstance(opt_value, str) else opt_value
                                for sos_file in opt_value:
                                    if not os.path.isfile(sos_file) and not os.path.isfile(os.path.join(os.path.split(self.sos_script)[0], sos_file)):
                                        env.logger.warning('Source file for nested workflow {} does not exist'.format(sos_file))
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
                self.sections.append(SoS_Step(step_names, step_options, is_parameters= step_names and step_names[0][0] == self._PARAMETERS_SECTION_NAME))
                cursect = self.sections[-1]
                continue
            #
            # directive?
            mo = self.DIRECTIVE.match(line)
            if mo:
                # check previous expression before a new directive
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid ' + cursect.category())
                    cursect.values = []
                    # allow multiple process-style actions
                    cursect.wrap_script()
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
                # is it an action??
                if directive_name in self._DIRECTIVES:
                    cursect.add_directive(directive_name, directive_value, lineno)
                else:
                    # should be in string mode ...
                    cursect.add_directive('process', directive_value, lineno, action=directive_name)
                continue
            # if section is string mode?
            if cursect and cursect.isValid() and cursect.category() == 'script':
                cursect.extend(line)
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
                    parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid ' + cursect.category())
                    continue
                cursect.values = []
                #
                var_name = mo.group('var_name')
                if var_name in self._DIRECTIVES:
                    parsing_errors.append(lineno, line, 'directive name cannot be used as variables')
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
                continue
            #
            # all others?
            if not cursect:
                self.sections.append(SoS_Step(is_global=True))
                cursect = self.sections[-1]
                cursect.add_statement(line, lineno)
                continue
            elif cursect.is_parameters:
                parsing_errors.append(lineno, line, 'Action statement is not allowed in {} section'.format(self._PARAMETERS_SECTION_NAME))
                continue
            #
            if cursect.empty() or cursect.category() != 'statements':
                # new statement
                cursect.add_statement(line, lineno)
            else:
                # existing one
                cursect.extend(line)

        #
        # check the last expression before a new directive
        if cursect:
            if not cursect.isValid():
                parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid ' + cursect.category())
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
            global_process = ''
            for statement in global_section[0][1].statements:
                if statement[0] == '=':
                    global_process += '{} = {}\n'.format(statement[1], statement[2])
                    env.readonly_vars.add(statement[1])
                else:
                    global_process += statement[1]
            #
            for section in self.sections:
                section.global_process = global_process
            # remove the global section after inserting it to each step of the process
            self.sections.pop(global_section[0][0])


    def workflow(self, workflow_name=None, extra_sections=[]):
        '''Return a workflow with name_step+name_step specified in wf_name
        This function might be called recursively because of nested
        workflow. extra_sections are sections read from other sos script
        when the source section option is encountered. '''
        allowed_steps = None
        if not workflow_name:
            wf_name = ''
        else:
            # if consists of multiple workflows
            if '+' in workflow_name:
                wfs = []
                for wf in workflow_name.split('+'):
                    if not self.SUBWORKFLOW.match(wf):
                        raise ValueError('Incorrect workflow name {}'.format(workflow_name))
                    # if this is a combined workflow, extra_section might be specied.
                    wfs.append(self.workflow(wf, extra_sections))
                combined_wf = wfs[0]
                for wf in wfs[1:]:
                    combined_wf.extend(wf)
                combined_wf.name = workflow_name
                return combined_wf
            # if a single workflow
            # workflow:15,-10 etc
            mo = self.SUBWORKFLOW.match(workflow_name)
            if not mo:
                raise ValueError('Incorrect workflow name {}'.format(workflow_name))
            wf_name, allowed_steps = mo.group('name', 'steps')
        # get workflow name from extra_sections
        extra_section_steps = sum([x.names for x in extra_sections if not \
            ('skip' in x.options and \
                (x.options['skip'] is None or x.options['skip'])) \
            and not ('target' in x.options)], [])
        extra_workflows = list(set([x[0] for x in extra_section_steps if '*' not in x[0]]))
        if extra_sections:
            env.logger.debug('Workflows {} imported'.format(', '.join(extra_workflows)))
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
        elif wf_name not in self.workflows + extra_workflows:
            raise ValueError('Workflow {} is undefined. Available workflows are: {}'.format(wf_name,
                ', '.join(self.workflows + extra_workflows)))
        # do not send extra parameters of ...
        sections = []
        # look for relevant sections in self.sections and extra sections from another script
        for section in self.sections + extra_sections:
            # skip, skip=True, skip=1 etc are all allowed.
            if 'skip' in section.options and (section.options['skip'] is None or section.options['skip']):
                continue
            if section.is_parameters:
                # include parameter only if they apply to wf_name
                if wf_name in section.names:
                    sections.append(section)
                continue
            if 'target' in section.options:
                # section global is shared by all workflows
                sections.append(section)
                continue
            for name, index, subworkflow in section.names:
                # exact match or filename like match if name contains * etc
                if fnmatch.fnmatch(wf_name, name):
                    # if there is an source option, ...
                    imported_sections = []
                    if 'source' in section.options:
                        for sos_file in section.options['source']:
                            if not os.path.isfile(sos_file):
                                sos_file = os.path.join(os.path.split(self.sos_script)[0], sos_file)
                            if not os.path.isfile(sos_file):
                                raise RuntimeError('Source file for nested workflow {} does not exist'.format(sos_file))
                            script = SoS_Script(filename=sos_file)
                            # this includes all global and parameters sections
                            imported_sections.extend(script.sections)
                    # expand subworkflow if needed
                    for i in range(len(section.names)):
                        if section.names[i][2] and isinstance(section.names[i][2], str):
                            # need to expand the workflow, but it is possible that it is nested ...
                            # so let us figure out what workflows we are expanding
                            expanded = []
                            for tmp in section.names[i][2].split('+'):
                                if ':' in tmp:
                                    tmp = tmp.split(':')[0]
                                if '_' in tmp and tmp.rsplit('_')[-1].isdigit():
                                    tmp = tmp.rsplit('_')[0]
                                expanded.append(tmp)
                            if wf_name in expanded:
                                env.logger.debug('NOT expanding {} because of potential loop'.format(section.names[i][2]))
                            else:
                                env.logger.debug('Expanding subworkflow {}'.format(section.names[i][2]))
                                # expand nested workflow
                                section.names[i] = (section.names[i][0], section.names[i][1], 
                                    self.workflow(section.names[i][2], extra_sections=extra_sections + imported_sections))
                    # now copy the step over
                    sections.append(section)
                    break
        return SoS_Workflow(wf_name, allowed_steps, sections, self.workflow_descriptions[wf_name])


    def show(self):
        textWidth = max(60, getTermWidth())
        if self.description:
            # separate \n\n
            for paragraph in dehtml(self.description).split('\n\n'):
                print('\n'.join(textwrap.wrap(paragraph, width=textWidth)))
        #
        text = 'Available workflows: {}'.format(', '.join(sorted(self.workflows)))
        print('\n' + '\n'.join(textwrap.wrap(text, width=textWidth, subsequent_indent=' '*8)))
        for workflow in sorted(self.workflows):
            wf = self.workflow(workflow)
            wf.show()
            print('')
#
# subcommmand show
#
def sos_show(args, workflow_args):
    try:
        script = SoS_Script(filename=args.script)
        if args.workflow:
            workflow = script.workflow(args.workflow)
            workflow.show()
        else:
            script.show()
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            print_traceback()
        env.logger.error(e)
        sys.exit(1)

#
# subcommand run
#
def sos_run(args, workflow_args):
    env.verbosity = args.verbosity
    env.max_jobs = args.__max_jobs__
    # kill all remainging processes when the master process is killed.
    atexit.register(env.cleanup)
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        if args.__dryrun__:
            env.run_mode = 'dryrun'
        if args.__rerun__:
            env.sig_mode = 'ignore'
        #
        workflow.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            print_traceback()
        env.logger.error(e)
        sys.exit(1)

#
# subcommand dryrun
#
def sos_dryrun(args, workflow_args):
    args.__max_jobs__ = 1
    args.__dryrun__ = True
    args.__rerun__ = False
    args.__config__ = None
    sos_run(args, workflow_args)
