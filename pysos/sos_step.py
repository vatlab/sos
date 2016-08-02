#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
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
import os
import sys
import re
import copy
import glob
import fnmatch
import multiprocessing as mp
from multiprocessing.pool import AsyncResult

from collections.abc import Sequence, Iterable
from collections import OrderedDict
from itertools import tee, combinations

from .utils import env, Error, short_repr, get_traceback, pickleable, transcribe
from .pattern import extract_pattern
from .sos_eval import  SoS_eval, SoS_exec, Undetermined
from .signature import  RuntimeInfo, textMD5
from .sos_syntax import SOS_INPUT_OPTIONS, SOS_DEPENDS_OPTIONS, SOS_OUTPUT_OPTIONS, \
    SOS_RUNTIME_OPTIONS

__all__ = []


class ExecuteError(Error):
    '''Raised when there are errors in inspect mode. Such errors are not raised
    immediately, but will be collected and raised at the end'''
    def __init__(self, workflow):
        Error.__init__(self, 'SoS workflow contains errors: %s' % workflow)
        self.workflow = workflow
        self.errors = []
        self.traces = []
        self.args = (workflow, )

    def append(self, line, error):
        lines = [x for x in line.split('\n') if x.strip()]
        if not lines:
            short_line = '<empty>'
        else:
            short_line = lines[0][:40] if len(lines[0]) > 40 else lines[0]
        self.errors.append(short_line)
        self.traces.append(get_traceback())
        if isinstance(error, Exception):
            self.message += '\n[%s] %s:\n\t%s' % (short_line, error.__class__.__name__, error)
        else:
            self.message += '\n[%s]:\n\t%s' % (short_line, error)


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

def execute_task(task, global_def, sos_dict, sigil, signature, workdir):
    '''A function that execute specified task within a local dictionary 
    (from SoS env.sos_dict). This function should be self-contained in that
    it can be handled by a task manager, be executed locally in a separate
    process or remotely on a different machine.'''
    env.register_process(os.getpid(), 'spawned_job with {} {}'
        .format(sos_dict['_input'], sos_dict['_output']))
    try:
        os.chdir(os.path.expanduser(workdir))
        env.sos_dict.quick_update(sos_dict)
        SoS_exec('import os, sys, glob')
        SoS_exec('from pysos import *')
        # re-execute global definition because some of the definitions in the 
        # global section might not be pickaleable (e.g. functions) and cannot
        # be passed to this separate process.
        if global_def:
            SoS_exec(global_def, sigil)
        # step process
        SoS_exec(task, sigil)
        os.chdir(env.exec_dir)
    except KeyboardInterrupt:
        raise RuntimeError('KeyboardInterrupt from {}'.format(os.getpid()))
    env.deregister_process(os.getpid())
    env.logger.trace('Execution completed with output ``{}``'.format(short_repr(env.sos_dict['_output'])))
    return {'succ': 0, 'output': env.sos_dict['_output']}


class Base_Step_Executor:
    # This base class defines how steps are executed. The derived classes will reimplement
    # some function to behave differently in different modes.
    #
    def __init__(self, step):
        self.step = step
        self.step_signature = self.get_step_signature()
        self.step_id = textMD5(self.step_signature)

    #
    # The following functions should be redefined in an executor
    # because it may behave differently in different modes.
    #
    def expand_input_files(self, *args, **kwargs):
        '''Process input files (perhaps a pattern) to determine input
        files.

        ret: 
            Return a file list or Undertermined.
        '''
        raise RuntimeError('Undefined virtual function.')

    # Nested functions to handle different parameters of input directive
    @staticmethod
    def handle_group_by(ifiles, group_by):
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
        elif isinstance(group_by, int) or group_by.isdigit():
            group_by = int(group_by)
            if len(ifiles) % group_by != 0:
                raise ValueError('Size of group_by block has to be divisible by the number of input files: {} provided'
                    .format(len(ifiles)))
            group_by = max(1, group_by)
            return [ifiles[i:i + group_by] for i in range(0, len(ifiles), group_by)]
        else:
            raise ValueError('Unsupported group_by option ``{}``!'.format(group_by))
    
    @staticmethod
    def handle_paired_with(paired_with, ifiles, _groups, _vars):
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

    @staticmethod
    def handle_extract_pattern(pattern, ifiles, _groups, _vars):
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
            Base_Step_Executor.handle_input_paired_with(res.keys(), ifiles, _groups, _vars)

    @staticmethod
    def handle_for_each(for_each, _groups, _vars):
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

    # directive input
    def process_input_args(self, ifiles, **kwargs):
        '''This function handles directive input and all its parameters.
        It
            determines and set __step_input__
            determines and set pattern variables if needed
        returns
            _groups
            _vars
        which are groups of _input and related _vars
        '''
        if isinstance(ifiles, Undetermined):
            return [Undetermined], [{}]

        for k in kwargs.keys():
            if k not in SOS_INPUT_OPTIONS:
                raise RuntimeError('Unrecognized input option {}'.format(k))
        #
        if 'filetype' in kwargs:
            if isinstance(kwargs['filetype'], str):
                ifiles = fnmatch.filter(ifiles, kwargs['filetype'])
            elif isinstance(kwargs['filetype'], Iterable):
                ifiles = [x for x in ifiles if any(fnmatch.fnmatch(x, y) for y in kwargs['filetype'])]
            elif callable(kwargs['filetype']):
                ifiles = [x for x in ifiles if kwargs['filetype'](x)]
        #
        # input file is the filtered files
        env.sos_dict.set('__step_input__', ifiles)
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
            handle_input_paired_with(kwargs['paired_with'], ifiles,  _groups, _vars)
        # handle pattern
        if 'pattern' in kwargs:
            handle_input_extract_pattern(kwargs['pattern'], ifiles, _groups, _vars)
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
            if k not in SOS_DEPENDS_OPTIONS:
                raise RuntimeError('Unrecognized depends option {}'.format(k))
        # first *args are filenames
        tmp = []
        for arg in args:
            if isinstance(arg, str):
                tmp.append(os.path.expanduser(arg))
            elif isinstance(arg, Iterable):
                arg = list(arg)
                if not all(isinstance(x, str) for x in arg):
                    raise RuntimeError('Invalid dependent file: {}'.format(arg))
                tmp.extend(os.path.expanduser(x) for x in arg)
            else:
                raise ValueError('Unrecognizable dependent type {}'.format(arg))
        #
        # expand wild card variables
        dfiles = []
        for dfile in tmp:
            if not os.path.isfile(dfile):
                expanded = sorted(glob.glob(os.path.expanduser(dfile)))
                if not expanded:
                    if env.run_mode in ['interactive', 'run']:
                        raise RuntimeError('{} not exist.'.format(dfile))
                    else:
                        dfiles.append(dfile)
                else:
                    dfiles.extend(expanded)
            else:
                dfiles.append(dfile)
        env.sos_dict.set('_depends', dfiles)


    def directive_output(*args, **kwargs):
        for k in kwargs.keys():
            if k not in SOS_OUTPUT_OPTIONS:
                raise RuntimeError('Unrecognized output option {}'.format(k))
        tmp = []
        for arg in args:
            if isinstance(arg, str):
                tmp.append(os.path.expanduser(arg))
            elif isinstance(arg, Iterable):
                arg = list(arg)
                if not all(isinstance(x, str) for x in arg):
                    raise RuntimeError('Invalid output file: {}'.format(arg))
                tmp.extend(arg)
            else:
                raise ValueError('Unrecognizable output type {}'.format(arg))
        #
        # expand wild card variables
        ofiles = []
        for ofile in tmp:
            if '*' in ofile or '?' in ofile or ('[' in ofile and ']' in ofile):
                expanded = sorted(glob.glob(os.path.expanduser(ofile)))
                if not expanded:
                    env.logger.warning('{} does not expand to any valid file.'.format(ofile))
                ofiles.extend(expanded)
            else:
                ofiles.append(ofile)
        #
        for ofile in ofiles:
            parent_dir = os.path.split(os.path.expanduser(ofile))[0]
            if parent_dir and not os.path.isdir(parent_dir):
                os.makedirs(parent_dir)
        env.sos_dict.set('_output', ofiles)

    def directive_task(**kwargs):
        for k,v in kwargs.items():
            if k not in SOS_RUNTIME_OPTIONS:
                raise RuntimeError('Unrecognized runtime option {}={}'.format(k, v))
        env.sos_dict.set('_runtime', kwargs)



    def get_step_signature(self):
        '''returns a signature of the step. Change of the step content will
        lead to the invalidation of the signature, which will then cause the
        re-execution of the step for any result from the step. '''
        #
        # TBD: should we also include global definitions? These can affect all 
        # steps.
        result = ''
        for statement in self.step.statements:
            if statement[0] in (':', '='):
                result += '{}: {}\n'.format(statement[1], statement[2])
            else:
                result += statement[1] + '\n'
        result += self.step.task
        return re.sub(r'\s+', ' ', result)

    def log(self, stage=None, msg=None):
        raise RuntimeError('Please redefine the log function in derived step executor.')

    def assign(self, key, value):
        try:
            env.sos_dict[key] = SoS_eval(value, self.step.sigil)
        except Exception as e:
            raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))

    def execute(self, stmt):
        try:
            SoS_exec(stmt, self.step.sigil)
        except Exception as e:
            raise RuntimeError('Failed to process statement {}: {}'.format(short_repr(stmt), e))

    def run(self):
        '''Execute a single step and return results. The result for batch mode is the
        input, output etc returned as alias, and for interactive mode is the return value
        of the last expression. '''
        #
        self.log('start')
        # 
        # prepare environments, namely variables that can be used by the step
        #
        # * step_name:  name of the step, can be used by step process to determine 
        #               actions dynamically.
        env.sos_dict.set('step_name', '{}_{}'.format(self.step.name, self.step.index))

        # * input:      input files, which should be __step_output__ if it is defined, or
        #               None otherwise. 
        # * _input:     first batch of input, which should be input if no input statement is used
        # * output:     None at first, can be redefined by output statement
        # * _output:    None at first, can be redefined by output statement
        # * depends:    None at first, can be redefined by depends statement
        # * _depends:   None at first, can be redefined by depends statement
        #
        if '__step_output__' not in env.sos_dict:
            env.sos_dict.set('input', None)
        else:
            if env.sos_dict['__step_output__'] is not None and not isinstance(env.sos_dict['__step_output__'], list):
                raise RuntimeError('__step_output__ can only be None or a list of files.')
            env.sos_dict.set('input', copy.deepcopy(env.sos_dict['__step_output__']))

        env.sos_dict.set('_input', copy.deepcopy(env.sos_dict['input']))
        env.sos_dict.set('output', None)
        env.sos_dict.set('_output', None)
        env.sos_dict.set('depends', None)
        env.sos_dict.set('_depends', None)

        # look for input statement.
        input_statement_idx = [idx for idx,x in enumerate(self.step.statements) if x[0] == ':' and x[1] == 'input']
        if not input_statement_idx:
            input_statement_idx = None
        elif len(input_statement_idx) == 1:
            input_statement_idx = input_statement_idx[0]
        else:
            raise RuntimeError('More than one step input are specified in step {}_{}'.format(self.step.name, self.step.index))

        # if there is an input statement, execute the statements before it, and then the input statement
        if input_statement_idx is not None:
            # execute before input stuff
            for statement in self.step.statements[:input_statement_idx]:
                if statement[0] == '=':
                    self.assign(statement[1], statement[2])
                elif statement[0] == ':':
                    raise RuntimeError('Step input should be specified before others')
                else:
                    self.execute(statement[1])
            # input statement
            self.log('input statement', stmt)
            try:
                args, kwargs = SoS_eval('__null_func__({})'.format(stmt), self.step.sigil)
                input_files = self.expand_input_files(*args, **kwargs)
                self._groups, self._vars = self.process_input_args(input_files, **kwargs)
            except Exception as e:
                if '__execute_errors__' in env.sos_dict and env.sos_dict['__execute_errors__'].errors:
                    raise env.sos_dict['__execute_errors__']
                else:
                    raise RuntimeError('Failed to process step {} : {} ({})'.format(key, value.strip(), e))

            # we cannot do anything
            if not self._groups:
                env.sos_dict.set('input', None)
                env.sos_dict.set('output', None)
                env.sos_dict.set('depends', None)
                return self.collectResult([])
            elif None in self._groups:
                if not all(x is None for x in self._groups):
                    raise RuntimeError('Either none or all of input groups can be unknown.')
                env.sos_dict.set('input', None)
            else:
                env.sos_dict.set('input', env.sos_dict['__step_input__'])
            input_statement_idx += 1
        else:
            # default case
            self._groups = [[]]
            self._vars = [{}]
            # assuming everything starts from 0 is after input
            input_statement_idx = 0

        self.log('input')
        
        concurrent = env.max_jobs > 1 and len(self._groups) > 1 and 'concurrent' in self.runtime_options and self.runtime_options['concurrent']
        if concurrent:
            pool = mp.Pool(min(env.max_jobs, len(self._groups)))
        # run steps after input statement, which will be run multiple times for each input 
        # group.
        env.sos_dict.set('__num_groups__', len(self._groups))
        for idx, (g, v) in enumerate(zip(self._groups, self._vars)):
            # other variables
            env.sos_dict.update(v)
            env.sos_dict.set('_input', g)
            env.sos_dict.set('_index', idx)
            for key in ('_output', '_depends'):
                if key in env.sos_dict:
                    env.sos_dict.pop(key)
            for statement in self.step.statements[input_statement_idx:]:
                if statement[0] == '=':
                    key, value = statement[1:]
                    public_vars.add(key)
                    try:
                        env.sos_dict[key] = SoS_eval(value, self.step.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to assign {} to variable {}: {}'.format(value, key, e))
                elif statement[0] == ':':
                    key, value, _ = statement[1:]
                    # output, depends, and process can be processed multiple times
                    try:
                        args, kwargs = SoS_eval('__null_func__({})'.format(value), self.step.sigil)
                        # dynamic output or dependent files
                        if 'dynamic' in kwargs:
                            env.logger.trace('Handling dynamic {}'.format(key))
                            if key not in ('output', 'depends'):
                                raise RuntimeError('dynamic option is only allowed for step input, output or depends')
                            if not isinstance(kwargs['dynamic'], bool):
                                raise RuntimeError('Option dynamic can only be True or False. {} provided'.format(kwargs['dynamic']))
                            if kwargs['dynamic']:
                                if env.run_mode == 'run' and key == 'depends':
                                    # depends need to be resolved now at run mode
                                    eval('directive_' + key)(*args, **kwargs)
                                # in other cases, namely non-run mode and output in run mode
                                else:
                                    if len(kwargs) > 1:
                                        raise RuntimeError('dynamic {} does not accept other options'.format(key))
                                    env.sos_dict.set('_' + key, [Undetermined(value)])
                        else:
                            eval('directive_' + key)(*args, **kwargs)
                    except Exception as e:
                        raise RuntimeError('Failed to process step {}: {} ({})'.format(key, value.strip(), e))
                    #
                    # we need to reduce output files in case they have been processed multiple times.
                    #
                    # The following temporarily set 'output' as an accumulated version of _output
                    # It will be reset after all the loops are done
                    #
                    if key == 'output':
                        if '_output' in env.sos_dict and env.sos_dict['_output'] and not isinstance(env.sos_dict['_output'][0], Undetermined):
                            if 'output' not in env.sos_dict:
                                env.sos_dict.set('output', copy.deepcopy(env.sos_dict['_output']))
                            elif not self._outputs or env.sos_dict['_output'] != self._outputs[-1]:
                                env.sos_dict['output'].extend(env.sos_dict['_output'])
                    elif key == 'depends':
                        if '_depends' in env.sos_dict:
                            if 'depends' not in env.sos_dict:
                                env.sos_dict.set('depends', copy.deepcopy(env.sos_dict['_depends']))
                            elif not self._depends or env.sos_dict['_depends'] != self._depends[-1]:
                                env.sos_dict['depends'].extend(env.sos_dict['_depends'])
                else:
                    # in prepare mode, check signature and see if all results exist
                    if env.run_mode in ('prepare', 'run') and '_output' in env.sos_dict and env.sos_dict['_output'] is not None and env.sos_dict['_input'] is not None:
                        signature = RuntimeInfo(step_sig, env.sos_dict['_input'], env.sos_dict['_output'], env.sos_dict.get('_depends', []), index=idx)
                        if env.sig_mode == 'default':
                            res = signature.validate()
                            if res:
                                env.sos_dict.set('_output', res['output'])
                                env.logger.debug('_output: {}'.format(res['output']))
                                env.logger.debug('Reuse existing output files ``{}``'.format(short_repr(env.sos_dict['_output'])))
                    #
                    try:
                        SoS_exec(statement[1], self.step.sigil)
                    except Exception as e:
                        raise RuntimeError('Failed to process statement {}: {}'.format(short_repr(statement[1]), e))
            #
            if '_output' in env.sos_dict:
                self._outputs.append(env.sos_dict['_output'])
            else:
                self._outputs.append(None)
            if '_depends' in env.sos_dict:
                self._depends.append(env.sos_dict['_depends'])
            else:
                self._depends.append([])
            #
            dict_stack.append(env.sos_dict.clone_pickleable())

            if 'active' in self.runtime_options:
                if isinstance(self.runtime_options['active'], int):
                    if self.runtime_options['active'] >= 0 and env.sos_dict['_index'] != self.runtime_options['active']:
                        continue
                    if self.runtime_options['active'] < 0 and env.sos_dict['_index'] != self.runtime_options['active'] + env.sos_dict['__num_groups__']:
                        continue
                elif isinstance(self.runtime_options['active'], Sequence):
                    allowed_index = list([x if x >= 0 else env.sos_dict['__num_groups__'] + x for x in self.runtime_options['active']])
                    if env.sos_dict['_index'] not in allowed_index:
                        continue
                elif isinstance(self.runtime_options['active'], slice):
                    allowed_index = list(range(env.sos_dict['__num_groups__']))[self.runtime_options['active']]
                    if env.sos_dict['_index'] not in allowed_index:
                        continue
                else:
                    raise RuntimeError('Unacceptable value for option active: {}'.format(self.runtime_options['active']))
            #
            # If the users specifies output files for each loop (using ${input} etc, we
            # can try to see if we can create partial signature. This would help if the
            # step is interrupted in the middle.
            partial_signature = None
            if env.sos_dict['_output'] is not None and env.sos_dict['_output'] != env.sos_dict['output'] and env.run_mode == 'run':
                partial_signature = RuntimeInfo(step_sig, env.sos_dict['_input'], env.sos_dict['_output'], env.sos_dict['_depends'], index=idx)
                if env.sig_mode == 'default':
                    if partial_signature.validate():
                        # everything matches
                        env.logger.info('Reusing existing output files {}'.format(', '.join(env.sos_dict['_output'])))
                        continue
                elif env.sig_mode == 'assert':
                    if not partial_signature.validate():
                        raise RuntimeError('Signature mismatch for input {} and output {}'.format(
                            ', '.join(env.sos_dict['_input']), ', '.join(env.sos_dict['_output'])))
                elif env.sig_mode == 'construct':
                    try:
                        partial_signature.write()
                        env.logger.debug('Construct signature from existing output files {}'.format(short_repr(env.sos_dict['_output'])))
                        continue
                    except Exception as e:
                        env.logger.debug('Failed to reconstruct signature. {}'.format(e))
            # now, if output file has already been generated using non-process statement
            # so that no process need to be run, we create signature from outside.
            if not self.step.task:
                # if no process, we should be able to figure out undetermined output now
                if env.sos_dict['_output'] and isinstance(env.sos_dict['_output'][0], Undetermined) and env.run_mode == 'run':
                    value = env.sos_dict['_output'][0].expr
                    env.logger.trace('Processing output: {}'.format(value))
                    args, kwargs = SoS_eval('__null_func__({})'.format(value), self.step.sigil)
                    # now we should have _output
                    directive_output(*args)
                    env.logger.trace('Reset _output to {}'.format(env.sos_dict['_output']))
                    self._outputs[idx] = env.sos_dict['_output']
                if partial_signature is not None:
                    partial_signature.set(env.sos_dict['_output'], 'output')
                    partial_signature.write()
                continue
            #
            env.logger.trace('Executing step process')
            try:
                if concurrent:
                    proc_results.append(pool.apply_async(
                        execute_step_process,   # function
                        (self.step.task,          # process
                        self.step.global_def,    # global process
                        env.sos_dict.clone_pickleable(),
                        self.step.sigil,
                        partial_signature,
                        self.runtime_options['workdir'] if 'workdir' in self.runtime_options else os.getcwd())))
                else:
                    # execute in existing process
                    proc_results.append(
                        execute_step_process(   # function
                        self.step.task,           # process
                        '',                     # local execusion, no need to re-run global
                        # do not clone dict
                        env.sos_dict,
                        self.step.sigil,
                        partial_signature,
                        self.runtime_options['workdir'] if 'workdir' in self.runtime_options else os.getcwd()))
            except Exception as e:
                # FIXME: cannot catch exception from subprocesses
                if env.verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise RuntimeError('Failed to execute process\n"{}"\n{}'.format(short_repr(self.step.task), e))
        # check results? This is only meaningful for pool
        if concurrent:
            try:
                proc_results = [res.get() if isinstance(res, AsyncResult) else res for res in proc_results]
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
        if proc_results:
            if not all(x['succ']==0 for x in proc_results):
                raise RuntimeError('Step process returns non-zero value')
            for idx, res in enumerate(proc_results):
                if self._outputs[idx] and isinstance(self._outputs[idx][0], Undetermined):
                    env.logger.trace('Setting _output[{}] from proc output {}'.format(idx, short_repr(res['output'])))
                    self._outputs[idx] = res['output']
        env.logger.trace('Checking output files {}'.format(short_repr(env.sos_dict['output'])))
        if env.run_mode == 'run' and env.sos_dict['output'] is not None:
            if env.sos_dict['output'] and isinstance(env.sos_dict['output'][0], Undetermined):
                # at this point self._outputs should be expanded already.
                env.sos_dict.set('output', list(OrderedDict.fromkeys(sum(self._outputs, []))))
                env.logger.info('output:  ``{}``'.format(short_repr(env.sos_dict['output'], noneAsNA=True)))
            for ofile in env.sos_dict['output']:
                if not os.path.isfile(os.path.expanduser(ofile)):
                    raise RuntimeError('Output file {} does not exist after completion of action'.format(ofile))
        if signature and env.run_mode == 'run':
            signature.set(env.sos_dict['output'], 'output')
            signature.write()
        if concurrent:
            # finally, write results back to the master process
            pool.close()
            pool.join()
        return self.collectResult(public_vars)



    def prepare_report(self):
        if '__step_report__' not in env.sos_dict:
            env.sos_dict.set('__step_report__', os.path.join('.sos', 'report', '{}_{}.md'.format(self.step.name, self.step.index)))
        if os.path.isfile(env.sos_dict['__step_report__']):
            # truncate the file
            with open(env.sos_dict['__step_report__'], 'w'):
                pass



class Queued_Step_Executor(Base_Step_Executor):
    # this class execute the step in a separate process
    # and returns result using a queue
    def __init__(self, step, queue, DAG):
        Base_Step_Executor.__init__(self, step)
        self.queue = queue
        self.DAG = DAG

    def run(self):
        try:
            res = self.run(DAG)
            # shared variables will be sent from subprocess to the master
            for key in env.shared_vars:
                if key in env.sos_dict and key not in res:
                    res[key] = env.sos_dict[key]
            res['__dag__'] = DAG
            # in run mode, these variables must be valid (not Undetermined)
            if env.run_mode == 'run':
                for key in ('__step_input__', '__step_output__', '__step_depends__'):
                    if res[key] is not None and not isinstance(res[key], list):
                        raise RuntimeError('Step input, output or depends has to be None or a list of filenames')
                    if res[key] is not None:
                        for v in res[key]:
                            if isinstance(v, Undetermined):
                                raise RuntimeError('Step input, output, or depends cannot be undetermined in run mode')
            queue.put(res)
        except Exception as e:
            queue.put(e)

    def collectResult(self, public_vars):
        # only results will be sent back to the master process
        #
        # __step_input__:    input of this step
        # __steo_output__:   output of this step
        # __step_depends__:  dependent files of this step
        result = {
            '__step_input__': env.sos_dict['input'],
            '__step_output__': env.sos_dict['output'],
            '__step_depends__': env.sos_dict['depends'],
        }
        if '__execute_errors__' in env.sos_dict:
            result['__execute_errors__'] = env.sos_dict['__execute_errors__']
        if 'alias' in self.step.options:
            step_info = StepInfo()
            step_info.set('step_name', env.sos_dict['step_name'])
            step_info.set('input', env.sos_dict['input'])
            step_info.set('output', env.sos_dict['output'])
            step_info.set('depends', env.sos_dict['depends'])
            # the step might be skipped
            for var in public_vars:
                # there is a slight possibility that var is deleted
                if var in env.sos_dict and pickleable(env.sos_dict[var]):
                    step_info.set(var, env.sos_dict[var])
            #
            if isinstance(self.step.options['alias'], Undetermined):
                # it is time to evalulate this expression now
                self.step.options['alias'] = self.step.options['alias'].value(self.step.sigil)
            transcribe('{} = {}'.format(self.step.options['alias'], step_info))
            result[self.step.options['alias']] = copy.deepcopy(step_info)
        return result




class Inspect_Step_Executor(Base_Step_Executor):
    def __init__(self, step, queue):
        Base_Step_Executor.__init__(self, step)
        self.queue = queue

    def expand_input_files(self, *args, **kwargs):
        if 'dynamic' in kwargs:
            env.sos_dict.set('__step_input__', None)
            return Undetermined()
 
        if ifiles and isinstance(ifiles[0], Undetermined):
            return Undetermined()

        if args:
            # if args are specified, it overrides the default __step_input__
            ifiles = []
            for arg in args:
                if isinstance(arg, str):
                    ifiles.append(os.path.expanduser(arg))
                elif isinstance(arg, Iterable):
                    # in case arg is a Generator, check its type will exhaust it
                    arg = list(arg)
                    if not all(isinstance(x, str) for x in arg):
                        raise RuntimeError('Invalid input file: {}'.format(arg))
                    raise ValueError('Unrecognizable input type {}'.format(arg))
        else:
            # Otherwise, ifiles are from the system passed '__step_output__'.
            if '__step_output__' in env.sos_dict and env.sos_dict['__step_output__'] is not None:
                ifiles = env.sos_dict['__step_output__']
            else:
                ifiles = []

        # expand files with wildcard characters and check if files exist
        tmp = []
        for ifile in ifiles:
            if os.path.isfile(os.path.expanduser(ifile)):
                tmp.append(ifile)
            else:
                expanded = sorted(glob.glob(os.path.expanduser(ifile)))
                if not expanded:
                    tmp.append(ifile)
                else:
                    tmp.extend(expanded)
        return tmp

    def log(self, stage, msg=None):
        if stage == 'start':
            env.logger.trace('Inspecting ``{}_{}``: {}'.format(self.step.name, self.step.index, self.step.comment.strip()))
        

class Prepare_Step_Executor(Base_Step_Executor):
    def __init__(self, step):
        env.run_mode = 'prepare'
        Base_Step_Executor.__init__(self, step)

    def expand_input_files(self, *args, **kwargs):
        # this function is identical to expand_input_file for Inspect_Step_Executor
        if 'dynamic' in kwargs:
            env.sos_dict.set('__step_input__', None)
            return Undetermined()
 
        if args:
            # if args are specified, it overrides the default __step_input__
            ifiles = []
            for arg in args:
                if isinstance(arg, str):
                    ifiles.append(os.path.expanduser(arg))
                elif isinstance(arg, Iterable):
                    # in case arg is a Generator, check its type will exhaust it
                    arg = list(arg)
                    if not all(isinstance(x, str) for x in arg):
                        raise RuntimeError('Invalid input file: {}'.format(arg))
                    raise ValueError('Unrecognizable input type {}'.format(arg))
        else:
            # Otherwise, ifiles are from the system passed '__step_output__'.
            if '__step_output__' in env.sos_dict and env.sos_dict['__step_output__'] is not None:
                ifiles = env.sos_dict['__step_output__']
            else:
                ifiles = []

        if ifiles and isinstance(ifiles[0], Undetermined):
            return Undetermined()

        # expand files with wildcard characters and check if files exist
        tmp = []
        for ifile in ifiles:
            if os.path.isfile(os.path.expanduser(ifile)):
                tmp.append(ifile)
            else:
                expanded = sorted(glob.glob(os.path.expanduser(ifile)))
                if not expanded:
                    tmp.append(ifile)
                else:
                    tmp.extend(expanded)
        return tmp

    def log(self, stage=0):
        if stage == 'start':
            env.logger.trace('Preparing ``{}_{}``: {}'.format(self.step.name, self.step.index, self.step.comment.strip()))

class Run_Step_Executor(Base_Step_Executor):
    def __init__(self, step):
        env.run_mode = 'run'
        Base_Step_Executor.__init__(self, step)

    def log(self, stage=None, msg=None):
        if stage == 'start':
            env.logger.info('Executing ``{}_{}``: {}'.format(self.step.name, self.step.index, self.step.comment.strip()))
        elif stage == 'input statement':
            env.logger.trace('Handling input statement {}'.format(msg))
        elif stage == 'input':
            env.logger.info('input:   ``{}``'.format(short_repr(env.sos_dict['input'], noneAsNA=True)))

    def assign(self, key, value):
        Base_Step_Executor.assign(self, key, value)
        public_keys[key] = env.sos_dict[key]
        transcribe('{} = {}'.format(key, env.sos_dict[key]))

    def expand_input_files(self, *args, **kwargs):
        # We ignore 'dynamic' option in run mode
        # we are using step_output as our input
        if args:
            # if args are specified, it overrides the default __step_input__
            ifiles = []
            for arg in args:
                if isinstance(arg, str):
                    ifiles.append(os.path.expanduser(arg))
                elif isinstance(arg, Iterable):
                    # in case arg is a Generator, check its type will exhaust it
                    arg = list(arg)
                    if not all(isinstance(x, str) for x in arg):
                        raise RuntimeError('Invalid input file: {}'.format(arg))
                    raise ValueError('Unrecognizable input type {}'.format(arg))
        else:
            # Otherwise, ifiles are from the system passed '__step_output__'.
            if '__step_output__' in env.sos_dict and env.sos_dict['__step_output__'] is not None:
                ifiles = env.sos_dict['__step_output__']
            else:
                ifiles = [] 
        #
        if ifiles and isinstance(ifiles[0], Undetermined):
            raise RuntimeError('Run mode does not accept undetermined input')
        #
        # if files are determined,
        #
        # expand files with wildcard characters and check if files exist
        tmp = []
        for ifile in ifiles:
            if os.path.isfile(os.path.expanduser(ifile)):
                tmp.append(ifile)
            else:
                # in this mode file must exist
                expanded = sorted(glob.glob(os.path.expanduser(ifile)))
                if not expanded:
                    raise RuntimeError('{} not exist'.format(ifile))
                tmp.extend(expanded)
        return tmp

class Interactive_Step_Executor(Base_Step_Executor):
    def __init__(self, step):
        env.run_mode = 'interactive'
        Base_Step_Executor.__init__(self, step)
    
    def expand_input_files(self, *args, **kwargs):
        # We ignore 'dynamic' option in run mode
        # we are using step_output as our input
        if args:
            # if args are specified, it overrides the default __step_input__
            ifiles = []
            for arg in args:
                if isinstance(arg, str):
                    ifiles.append(os.path.expanduser(arg))
                elif isinstance(arg, Iterable):
                    # in case arg is a Generator, check its type will exhaust it
                    arg = list(arg)
                    if not all(isinstance(x, str) for x in arg):
                        raise RuntimeError('Invalid input file: {}'.format(arg))
                    raise ValueError('Unrecognizable input type {}'.format(arg))
        else:
            # Otherwise, ifiles are from the system passed '__step_output__'.
            if '__step_output__' in env.sos_dict and env.sos_dict['__step_output__'] is not None:
                ifiles = env.sos_dict['__step_output__']
            else:
                ifiles = [] 
        #
        if ifiles and isinstance(ifiles[0], Undetermined):
            raise RuntimeError('Run mode does not accept undetermined input')
        #
        # if files are determined,
        #
        # expand files with wildcard characters and check if files exist
        tmp = []
        for ifile in ifiles:
            if os.path.isfile(os.path.expanduser(ifile)):
                tmp.append(ifile)
            else:
                # in this mode file must exist
                expanded = sorted(glob.glob(os.path.expanduser(ifile)))
                if not expanded:
                    raise RuntimeError('{} not exist'.format(ifile))
                tmp.extend(expanded)
        return tmp

    def log(self, stage=None):
        return
