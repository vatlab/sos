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
import glob
import fnmatch
import argparse
import textwrap
import traceback
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

from . import __version__
from .actions import *
from .utils import env, Error, WorkflowDict, SoS_eval, SoS_exec, RuntimeInfo, \
    dehtml, getTermWidth, interpolate, shortRepr, glob_wildcards, apply_wildcards

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

class RuntimeEnvironment:
    """Context manager for changing the current working directory"""
    def __init__(self, runtime={}):
        self.runtime = runtime
        if 'workdir' in runtime:
            self.workdir = os.path.expanduser(runtime['workdir'])
        else:
            self.workdir = None

    def __enter__(self):
        if self.workdir:
            self.saved_dir = os.getcwd()
            os.chdir(self.workdir)

    def __exit__(self, etype, value, traceback):
        if self.workdir:
            os.chdir(self.saved_dir)

class StepInfo:
    '''A simple class to hold input, output, and index of step'''
    def __init__(self):
        pass

class SoS_Step:
    #
    # A single sos step
    #
    _INPUT_OPTIONS = ['group_by', 'skip', 'filetype', 'paired_with', 'for_each', 'pattern', 'dynamic']
    _OUTPUT_OPTIONS = ['pattern', 'dynamic']
    _DEPENDS_OPTIONS = ['pattern', 'dynamic']
    _RUNTIME_OPTIONS = ['workdir']
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
        # step process
        self.process = ''
        # subworkflow of a step
        self.subworkflow = None
        # is it global section?
        self.is_global = is_global
        # is it the parameters section?
        self.is_parameters = is_parameters
        # indicate the type of input of the last line
        self.values = []
        self.lineno = None
        #
        self.runtime = RuntimeEnvironment()
        if 'sigil' in self.options:
            self.sigil = self.options['sigil']
        else:
            self.sigil = '${ }'

    def category(self):
        if self.statements:
            if self.statements[-1][0] == '=':
                return 'expression'
            elif self.statements[-1][0] == ':':
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

    def add_directive(self, key, value, lineno=None):
        '''Assignments are items with ':' type '''
        if key is None:
            # continuation of multi-line directive
            self.statements[-1][-1] += value
            self.values.append(value)
        else:
            # new directive
            self.statements.append([':', key, value])
            self.values = [value]
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
        self.statements.append(['!', line])
        self.back_comment = ''
        if lineno:
            self.lineno = lineno

    def finalize(self):
        ''' split statement and process by last directive '''
        if not self.statements:
            self.process = ''
            return
        directives = [idx for idx, statement in enumerate(self.statements) if statement[0] == ':']
        if not directives:
            start_process = 0
        else:
            start_process = directives[-1] + 1
        # convert statement to process
        self.process = ''
        for statement in self.statements[start_process:]:
            if statement[0] == '=':
                self.process += '{} = {}'.format(statement[1], statement[2])
            else:
                self.process += statement[1]
        # remove process from self.statement
        if not directives:
            self.statements = []
        else:
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
            else:
                raise RuntimeError('Unrecognized expression type {}'.format(self.category()))
            return True
        except Exception as e:
            return False


    #
    # Execution: input options
    #
    def _handle_input_group_by(self, ifiles, group_by):
        # for this special case, the step is skipped
        if group_by == 'single':
            return [[x] for x in ifiles]
        elif group_by == 'all':
            # default option
            return [ifiles]
        elif group_by == 'pairs':
            if len(ifiles) % 2 != 0:
                raise ValueError('Paired group_by has to have even number of input files: {} provided'
                    .format(len(ifiles)))
            env.logger.error( list(zip(ifiles[:len(ifiles)//2], ifiles[len(ifiles)//2:])))
            return list(list(x) for x in zip(ifiles[:len(ifiles)//2], ifiles[len(ifiles)//2:]))
        elif group_by == 'pairwise':
            f1, f2 = tee(ifiles)
            next(f2, None)
            return [list(x) for x in izip(f1, f2)]
        elif group_by == 'combinations':
            return [list(x) for x in combinations(ifiles, 2)]

    def _handle_input_pattern(self, pattern):
        ifiles = env.locals['_step'].input
        # 
        if pattern is None or not pattern:
            patterns = []
        elif isinstance(pattern, basestring):
            patterns = [pattern]
        elif isinstance(pattern, list):
            patterns = pattern
        else:
            raise ValueError('Unacceptable value for parameter pattern: {}'.format(pattern))
        #
        for pattern in patterns:
            res = glob_wildcards(pattern, [])
            for ifile in ifiles:
                matched = glob_wildcards(pattern, [ifile])
                for key in matched.keys():
                    if not matched[key]:
                        env.logger.warning('Filename {} does not match pattern {}. None returned.'.format(ifile, pattern))
                        res[key].append(None)
                    else:
                        res[key].extend(matched[key])
            # now, assign the variables to env
            for k, v in res.items():
                env.locals[k] = v
            # also make k, v pair with _input
            self._handle_input_paired_with(res.keys())
        
    def _handle_input_paired_with(self, paired_with):
        if paired_with is None or not paired_with:
            paired_with = []
        elif isinstance(paired_with, basestring):
            paired_with = [paired_with]
        elif isinstance(paired_with, list):
            paired_with = paired_with
        else:
            raise ValueError('Unacceptable value for parameter paired_with: {}'.format(paired_with))
        #
        ifiles = env.locals['_step'].input
        if not hasattr(self, '_vars'):
            self._vars = [{} for x in self._groups]
        for wv in paired_with:
            values = env.locals[wv]
            if isinstance(values, basestring) or not isinstance(values, Sequence):
                raise ValueError('with_var variable {} is not a sequence ("{}")'.format(wv, values))
            if len(values) != len(ifiles):
                raise ValueError('Length of variable {} (length {}) should match the number of input files (length {}).'
                    .format(wv, len(values), len(ifiles)))
            file_map = {x:y for x,y in zip(ifiles, values)}
            for idx, grp in enumerate(self._groups):
                self._vars[idx]['_' + wv] = [file_map[x] for x in grp]

    def _handle_input_for_each(self, for_each):
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
                values = env.locals[fe]
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
                        self._vars[idx]['_' + fe] = env.locals[fe][vidx]
                tmp.extend(copy.deepcopy(self._vars))
            self._vars = tmp

    #
    # Execution; directives
    #
    def _directive_input(self, *args, **kwargs):
        # first *args are filenames
        for k in kwargs.keys():
            if k not in self._INPUT_OPTIONS:
                raise RuntimeError('Unrecognized input option {}'.format(k))
        #
        if args:
            ifiles = []
            for arg in args:
                if isinstance(arg, basestring):
                    ifiles.append(arg)
                elif isinstance(arg, list):
                    if not all(isinstance(x, basestring) for x in arg):
                        raise RuntimeError('Invalid input file: {}'.format(arg))
                    ifiles.extend(arg)
            env.locals['_step'].input = ifiles
        else:
            ifiles = env.locals['_step'].input
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
        if 'skip' in kwargs and kwargs['skip']:
            self._groups = []
            self._vars = []
            return
        if 'filetype' in kwargs:
            if isinstance(kwargs['filetype'], basestring):
                ifiles = fnmatch.filter(ifiles, kwargs['filetype'])
            elif isinstance(kwargs['filetype'], (list, tuple)):
                ifiles = [x for x in ifiles if any(fnmatch.fnmatch(x, y) for y in kwargs['filetype'])]
            elif callable(kwargs['filetype']):
                ifiles = [x for x in ifiles if kwargs['filetype'](x)]
        #             
        # handle group_by
        if 'group_by' in kwargs:
            self._groups = self._handle_input_group_by(ifiles, kwargs['group_by'])
        else:
            self._groups = [ifiles]
        #
        self._vars = [{} for x in self._groups]
        # handle paired_with
        if 'paired_with' in kwargs:
            self._handle_input_paired_with(kwargs['paired_with'])
        # handle pattern
        if 'pattern' in kwargs:
            self._handle_input_pattern(kwargs['pattern'])
        # handle for_each
        if 'for_each' in kwargs:
            self._handle_input_for_each(kwargs['for_each'])


    def _directive_depends(self, *args, **kwargs):
        '''handle directive depends'''
        for k in kwargs.keys():
            if k not in self._DEPENDS_OPTIONS:
                raise RuntimeError('Unrecognized depends option {}'.format(k))
        # first *args are filenames
        dfiles = []
        for arg in args:
            if isinstance(arg, basestring):
                dfiles.append(arg)
            elif isinstance(arg, list):
                if not all(isinstance(x, basestring) for x in arg):
                    raise RuntimeError('Invalid dependent file: {}'.format(arg))
                dfiles.extend(arg)
        env.locals.set('_depends', dfiles)

    def _handle_output_pattern(self, pattern, ofiles):
        # 
        if pattern is None or not pattern:
            patterns = []
        elif isinstance(pattern, basestring):
            patterns = [pattern]
        elif isinstance(pattern, list):
            patterns = pattern
        else:
            raise ValueError('Unacceptable value for parameter pattern: {}'.format(pattern))
        #
        for pattern in patterns:
            sz = None
            res = glob_wildcards(pattern, [])
            sz = None
            wildcard = [{}]
            for key in res.keys():
                if key not in env.locals:
                    raise ValueError('Undefined variable {} in pattern {}'.format(key, pattern))
                if not isinstance(env.locals[key], basestring) and isinstance(env.locals[key], Sequence):
                    if sz is None:
                        sz = len(env.locals[key])
                        wildcard = [{} for x in range(sz)]
                    elif sz != len(env.locals[key]):
                        raise ValueError('Variables in output pattern should have the same length (other={}, len({})={})'
                            .format(sz, key, len(env.locals[key])))
                    for idx, value in enumerate(env.locals[key]):
                        wildcard[idx][key] = value
                else:
                    for v in wildcard:
                        v[key] = env.locals[key]
            #
            for card in wildcard:
                ofiles.append(apply_wildcards(pattern, card, fill_missing=False,
                   fail_dynamic=False, dynamic_fill=None, keep_dynamic=False))
            
    def _directive_output(self, *args, **kwargs):
        for k in kwargs.keys():
            if k not in self._OUTPUT_OPTIONS:
                raise RuntimeError('Unrecognized output option {}'.format(k))
        ofiles = []
        for arg in args:
            if isinstance(arg, basestring):
                ofiles.append(arg)
            elif isinstance(arg, list):
                if not all(isinstance(x, basestring) for x in arg):
                    raise RuntimeError('Invalid output file: {}'.format(arg))
                ofiles.extend(arg)
        #
        if 'pattern' in kwargs:
            self._handle_output_pattern(kwargs['pattern'], ofiles)
        env.locals.set('_output', ofiles)

    def _directive_runtime(self, **kwargs):
        for k in kwargs.keys():
            if k not in self._RUNTIME_OPTIONS:
                raise RuntimeError('Unrecognized runtime option {}'.format(k))
        #
        self.runtime = RuntimeEnvironment(kwargs)

    #
    # handling of parameters step
    #
    def _parse_error(self, msg):
        '''This function will replace error() function in argparse module so that SoS
        can hijack errors raised from it.'''
        raise ArgumentError(msg)

    def parse_args(self, args, check_unused=False, cmd_name=''):
        '''Parse command line arguments and set values to parameters section'''
        env.logger.info('Execute ``{}_parameters``'.format(self.name))
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
                    parser.add_argument('--{}'.format(key), type=type(defvalue), help=comment,
                        nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, basestring) else '?',
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
            env.locals[k] = v


    def run_global(self):
        '''Execute the global section. The biggest difference is that there is no directive
        and no changing or _step ...'''
        env.logger.info('Execute ``{}_global``'.format(self.name))
        try:
            SoS_exec(self.process)
        except Exception as e:
            if env.verbosity > 2:
                print_traceback()
            raise RuntimeError('Failed to execute process\n{}\n{}'.format(self.process, e))

    def execute_single(self):
        '''Execute a single step 
        '''
        # If the users specifies output files for each loop (using ${input} etc, we
        # can try to see if we can create partial signature. This would help if the
        # step is interrupted in the middle.
        partial_signature = None
        if env.locals['_output'] and env.locals['_output'] != env.locals['_step'].output and env.run_mode == 'run':
            partial_signature = RuntimeInfo(self.process, env.locals['_input'], env.locals['_output'], 
                env.locals['_depends'])
            if partial_signature.validate():
                # everything matches
                env.logger.info('Reusing existing output files {}'.format(', '.join(env.locals['_output'])))
                return
        with self.runtime:
            try:
                # in case of nested workflow, these names might be changed and need to be reset
                env.locals['_step'].name = '{}_{}'.format(self.name, self.index)
                env.locals['_step'].index = self.index
                #
                SoS_exec(self.process, self.sigil)
                #
                # if there is subworkflow, the subworkflow is considered part of the process
                # and will be executed mutliple times if necessary. The partial signature 
                # stuff works for these though
                if self.subworkflow:
                    self.subworkflow.run(nested=True)
            except Exception as e:
                if env.verbosity > 2:
                    print_traceback()
                raise RuntimeError('Failed to execute process\n{}\n{}'.format(self.process, e))
        if partial_signature:
            partial_signature.write()
        
    #
    # Execution
    #
    def run(self):
        # handle these two sections differently
        if isinstance(self.index, int):
            env.locals.set('_workflow_index', self.index)
        #
        env.logger.info('Execute ``{}_{}``: {}'.format(self.name, self.index, self.comment))
        #
        # 
        # the following is a quick hack to allow _directive_input function etc to access 
        # the workflow dictionary
        env.globals['self'] = self
        env.locals['_step'].name = '{}_{}'.format(self.name, self.index)
        env.locals['_step'].index = self.index
        env.locals['_step'].output = []
        env.locals['_step'].depends = []
        #
        # these are temporary variables that should be removed if exist
        for var in ('_input', '_depends', '_output'):
            env.locals.pop(var, '')
        #
        # default input groups and vars
        if hasattr(env.locals['_step'], 'input'):
            self._groups = [env.locals['_step'].input]
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
        # if there is no input, easily
        self._outputs = []
        self._depends = []
        if input_idx is None:
            # no input? _input is _step.input
            env.locals.set('_input', env.locals['_step'].input)
            for statement in self.statements:
                if statement[0] == '=':
                    key, value = statement[1:]
                    try:
                        env.locals[key] = SoS_eval(value, self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
                elif statement[0] == ':':
                    key, value = statement[1:]
                    # input and runtime is processed only once
                    try:
                        SoS_eval('self._directive_{}({})'.format(key, value), self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to process step {}: {}'.format(key, e))
            # these assignments could affect _output and _depends and affect
            # saved _outputs and _depends list
            if '_output' in env.locals:
                self._outputs.append(env.locals['_output'])
            else:
                self._outputs.append([])
            if '_depends' in env.locals:
                self._depends.append(env.locals['_depends'])
            else:
                self._depends.append([])
        else:
            # if there is input directive
            for statement in self.statements[:input_idx]:
                if statement[0] == '=':
                    key, value = statement[1:]
                    try:
                        env.locals[key] = SoS_eval(value, self.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
                elif statement[0] == ':':
                    raise RuntimeError('Step input should be specified before others')
            # input
            key, value = self.statements[input_idx][1:]
            try:
                SoS_eval('self._directive_{}({})'.format(key, value), self.sigil)
            except Exception as e:
                raise RuntimeError('Failed to process directive {}: {}'.format(key, e))
            # post input
            # output and depends can be processed many times
            for g, v in zip(self._groups, self._vars):
                # other variables
                env.locals.update(v)
                env.locals.set('_input', g)
                for statement in self.statements[input_idx+1:]:
                    if statement[0] == '=':
                        key, value = statement[1:]
                        try:
                            env.locals[key] = SoS_eval(value, self.sigil)
                        except Exception as e:
                            raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
                    elif statement[0] == ':':
                        key, value = statement[1:]
                        # input and runtime is processed only once
                        try:
                            SoS_eval('self._directive_{}({})'.format(key, value), self.sigil)
                        except Exception as e:
                            raise RuntimeError('Failed to process directive {}: {}'.format(key, e))
                # collect _output and _depends
                if '_output' in env.locals:
                    self._outputs.append(env.locals['_output'])
                    env.locals.pop('_output')
                else:
                    self._outputs.append([])
                if '_depends' in env.locals:
                    self._depends.append(env.locals['_depends'])
                    env.locals.pop('_depends')
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
        env.locals['_step'].output = list(OrderedDict.fromkeys(sum(self._outputs, [])))
        env.locals['_step'].depends = list(OrderedDict.fromkeys(sum(self._depends, [])))
        env.logger.info('_step.input: ``{}``'.format(shortRepr(env.locals['_step'].input)))
        env.logger.info('_step.output: ``{}``'.format(shortRepr(env.locals['_step'].output)))
        env.logger.info('_step.depends: ``{}``'.format(shortRepr(env.locals['_step'].depends)))
        #
        # input alias
        # if the step is ignored
        if not self._groups:
            env.logger.info('Step {} is skipped'.format(self.index))
            return
        if 'alias' in self.options:
            env.locals[self.options['alias']] = copy.deepcopy(env.locals['_step'])
        #
        if env.locals['_step'].output:
            signature = RuntimeInfo(self.process, 
                env.locals['_step'].input, env.locals['_step'].output, env.locals['_step'].depends)
            if env.run_mode == 'run':
                for ofile in env.locals['_step'].output:
                    parent_dir = os.path.split(ofile)[0]
                    if parent_dir and not os.path.isdir(parent_dir):
                        os.makedirs(parent_dir)
                if signature.validate():
                    # everything matches
                    env.logger.info('Reusing existing output files {}'.format(', '.join(env.locals['_step'].output)))
                    return
        else:
            signature = None
        #
        for idx, (g, v, o, d) in enumerate(zip(self._groups, self._vars, self._outputs, self._depends)):
            env.locals.update(v)
            env.locals['_input'] = g
            env.locals['_output'] = o
            env.locals['_depends'] = d
            env.locals['_index'] = idx
            env.logger.info('_idx: ``{}``'.format(idx))
            env.logger.info('_input: ``{}``'.format(shortRepr(env.locals['_input'])))
            env.logger.info('_output: ``{}``'.format(shortRepr(env.locals['_output'])))
            #
            # action
            self.execute_single()
        if env.run_mode == 'run':
            for ofile in env.locals['_step'].output:
                if not os.path.isfile(os.path.expanduser(ofile)): 
                    raise RuntimeError('Output file {} does not exist after completion of action'.format(ofile))
        if signature and env.run_mode == 'run':
            signature.write()

    def __repr__(self):
        result = ''
        if self.is_global:
            result += '## global definitions ##\n'
        elif self.is_parameters:
            result += '[parameters]\n'
        else:
            result += '[{}:{}]'.format(','.join('{}_{}'.format(x,y) if y else x for x,y,z in self.names),
                ','.join('{}={}'.format(x,y) for x,y in self.options.items()))
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
        if self.subworkflow:
            result += str(self.subworkflow)
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
        self.auxillary_sections = []
        #
        for section in sections:
            if section.is_global or section.is_parameters:
                self.sections.append(copy.deepcopy(section))
                self.sections[-1].name = workflow_name
                # for ordering purpose, this section is always after global
                self.sections[-1].index = -2 if section.is_global else -1
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
            if section.subworkflow is not None and isinstance(section.subworkflow, basestring):
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

    def run(self, args=[], nested=False, cmd_name=''):
        '''Execute a workflow with specified command line args. If sub is True, this 
        workflow is a nested workflow and be treated slightly differently.
        '''
        if nested:
            # if this is a subworkflow, we use _input as step input the workflow
            env.locals['_step'].input = env.locals['_input']
        else:
            # other wise we need to prepare running environment.
            env.globals = globals()
            # Because this workflow might belong to a combined workflow, we do not clear
            # locals before the execution of workflow.
            if not hasattr(env, 'locals'):
                env.locals = WorkflowDict()
            # initial values
            env.locals.set('SOS_VERSION', __version__)
            py_version = sys.version_info
            env.locals.set('SOS_PYTHON_VERSION', (py_version.major, py_version.minor, py_version.micro))

            #
            env.locals.set('_step', StepInfo())
            # initially there is no input, output, or depends
            env.locals['_step'].input = []
            env.locals['_step'].name = self.name
            env.locals['_step'].output = []
            env.locals['_step'].depends = []
        #
        # process step of the pipelinp
        #
        start = True
        num_parameters_sections = len([x for x in self.sections if x.is_parameters or x.subworkflow])
        if num_parameters_sections == 0 and args:
            raise ArgumentError('Unused parameter {}'.format(' '.join(args)))
        for idx, section in enumerate(self.sections):
            # global section will not change _step etc
            if section.is_global:
                section.run_global()
                continue
            elif section.is_parameters:
                # if there is only one parameters section and no nested workflow, check unused section
                section.parse_args(args, num_parameters_sections == 1, cmd_name=cmd_name)
                continue
            # 
            # execute section with specified input
            # 1. for first step of workflow, _step.input=[]
            # 2. for subworkflow, _step.input = _input
            # 3. for second to later step, _step.input = _step.output
            if not start:
                # set the output to the input of the next step
                if hasattr(env.locals['_step'], 'output'):
                    # passing step output to _step.input of next step
                    env.locals['_step'].input = env.locals['_step'].output
                    # step_output and depends are temporary
                else:
                    env.locals['_step'].input = []
                #
                env.locals['_step'].output = []
                env.locals['_step'].depends = []
            else:
                # use initial values
                start = False
            #
            section.run()

    def show(self, indent = 0):
        textWidth = max(60, getTermWidth())
        paragraphs = dehtml(self.description).split('\n\n')
        print('\n')
        print('\n'.join(
            textwrap.wrap('{}:  {}'.format(self.name, paragraphs[0]),
                initial_indent = ' '*indent,
                subsequent_indent = ' '*indent,
                width=textWidth)
            ))
        for paragraph in paragraphs[1:]:
            print('\n'.join(
            textwrap.wrap(paragraph, width=textWidth,
                initial_indent = ' '*indent,
                subsequent_indent = ' '*indent)
            ))
        for step in self.sections:
            if step.is_parameters:
                print('\n{}Parameters:'.format(' '*indent))
                for k,v,c in step.parameters:
                    #
                    text = '{:<20}'.format(k) + (c + ' ' if c else '') + \
                        ('(default: {})'.format(v) if v else '')
                    print('\n'.join(
                        textwrap.wrap(text, 
                        initial_indent=' '*(2 + indent),
                        subsequent_indent=' '*(22 + indent),
                        width=textWidth)
                        ))
            elif step.index >= 0:
                # hide a step if there is no comment
                text = '{:<20}'.format('{}_{}:'.format(step.name, step.index)) + step.comment
                print('\n'.join(
                    textwrap.wrap(text, 
                        width=textWidth, 
                        initial_indent=' '*indent,
                        subsequent_indent=' '*(22 + indent))))
                if step.subworkflow:
                    step.subworkflow.show(indent+4)


class SoS_Script:
    _DIRECTIVES = ['input', 'output', 'depends', 'runtime']
    _SECTION_OPTIONS = ['alias', 'nonconcurrent',
        'skip', 'blocking', 'sigil', 'target', 'source']
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
        (?P<directive_name>{})             # name of directive
        \s*:\s*                            # followed by :
        (?P<directive_value>.*)            # and values
        '''.format('|'.join(_DIRECTIVES))

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
            if not os.path.isfile(os.path.expanduser(filename)):
                raise ValueError('{} does not exist'.format(filename))
            with open(filename) as fp:
                self._read(fp, filename)
        else:
            if os.path.isfile(content):
                with open(content) as fp:
                    self._read(fp, content)
            else:
                with StringIO(content) as fp:
                    self._read(fp, '<string>')
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


    def _read(self, fp, fpname):
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
                        self.descriptions[-1] += line.lstrip('#').lstrip()
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
                                opt_value = [opt_value] if isinstance(opt_value, basestring) else opt_value
                                for sos_file in opt_value:
                                    if not os.path.isfile(sos_file):
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
            # directive?
            mo = self.DIRECTIVE.match(line)
            if mo:
                # check previous expression before a new directive
                if cursect:
                    if not cursect.isValid():
                        parsing_errors.append(cursect.lineno, ''.join(cursect.values[:5]), 'Invalid ' + cursect.category())
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
                if not cursect.empty() and cursect.category() == 'statements':
                    parsing_errors.append(lineno, line, 'Directive {} should be be defined before step action'.format(directive_name))
                    continue
                cursect.add_directive(directive_name, directive_value, lineno)
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
            if section.is_global or section.is_parameters:
                # include global or parameter only if they apply to wf_name
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
                                raise RuntimeError('Source file for nested workflow {} does not exist'.format(sos_file))
                            script = SoS_Script(filename=sos_file)
                            # this includes all global and parameters sections
                            imported_sections.extend(script.sections)
                    # expand subworkflow if needed
                    for i in range(len(section.names)):
                        if section.names[i][2] and isinstance(section.names[i][2], basestring):
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
        print('\nUse command "sos show script workflow_name" to display details of specific workflow.')

#
# subcommmand show
#
def print_traceback():
    exc_type, exc_value, exc_traceback = sys.exc_info()
    #print "*** print_tb:"
    traceback.print_tb(exc_traceback, limit=1, file=sys.stderr)
    #print "*** print_exception:"
    traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=5, file=sys.stderr)
    #print "*** print_exc:"
    #traceback.print_exc()
    #print "*** format_exc, first and last line:"
    #formatted_lines = traceback.format_exc().splitlines()
    #print formatted_lines[0]
    #print formatted_lines[-1]
    #print "*** format_exception:"
    #print repr(traceback.format_exception(exc_type, exc_value,
    #                                      exc_traceback))
    #print "*** extract_tb:"
    #print repr(traceback.extract_tb(exc_traceback))
    #print "*** format_tb:"
    #print repr(traceback.format_tb(exc_traceback))
    #print "*** tb_lineno:", exc_traceback.tb_lineno

def sos_show(args):
    # options such as -v might be put at the end and we need to 
    # extract them from args.options
    mini_parser = argparse.ArgumentParser()
    mini_parser.add_argument('-v', '--verbosity', type=int, choices=range(5))
    remaining_args, remainder = mini_parser.parse_known_args(args.options)
    for k, v in remaining_args.__dict__.items():
        setattr(args, k, v)
    args.options = remainder
    env.verbosity = args.verbosity
    try:
        script = SoS_Script(filename=args.script)
        # do not parse args
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
def sos_run(args):
    # options such as -d might be put at the end and we need to 
    # extract them from args.options
    mini_parser = argparse.ArgumentParser()
    mini_parser.add_argument('-v', '--verbosity', type=int, choices=range(5))
    mini_parser.add_argument('-d', action='store_true', dest='__dryrun__')
    remaining_args, remainder = mini_parser.parse_known_args(args.options)
    for k, v in remaining_args.__dict__.items():
        setattr(args, k, v)
    args.options = remainder
    env.verbosity = args.verbosity
    env.max_jobs = args.__max_jobs__
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        if args.__dryrun__:
            env.run_mode = 'dryrun'
        workflow.run(args.options, cmd_name='{} {}'.format(args.script, args.workflow))
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            print_traceback()
        env.logger.error(e)
        sys.exit(1)

#
# subcommand dryrun
#
def sos_dryrun(args):
    args.options.append('-d')
    sos_run(args)
