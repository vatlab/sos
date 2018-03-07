#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
import traceback
import contextlib
import psutil
from io import StringIO
from multiprocessing import Pool

from collections.abc import Sequence, Iterable, Mapping
from itertools import tee, combinations

from .utils import env, StopInputGroup, TerminateExecution, short_repr, stable_repr,\
    get_traceback, expand_size, format_HHMMSS, SlotManager
from .pattern import extract_pattern
from .eval import SoS_eval, SoS_exec, Undetermined, stmtHash
from .targets import BaseTarget, file_target, dynamic, remote, RuntimeInfo, UnknownTarget, RemovedTarget, UnavailableLock, sos_targets, path, paths
from .syntax import SOS_INPUT_OPTIONS, SOS_DEPENDS_OPTIONS, SOS_OUTPUT_OPTIONS, \
    SOS_RUNTIME_OPTIONS, SOS_TAG
from .tasks import TaskParams, MasterTaskParams

__all__ = []

class PendingTasks(Exception):
    def __init__(self, tasks, *args, **kwargs):
        super(PendingTasks, self).__init__(*args, **kwargs)
        self.tasks = tasks

def analyze_section(section, default_input=None):
    '''Analyze a section for how it uses input and output, what variables
    it uses, and input, output, etc.'''
    from .workflow_executor import __null_func__
    from ._version import __version__
    from .eval import accessed_vars

    # these are the information we need to build a DAG, by default
    # input and output and undetermined, and there are no variables.
    #
    # step input and output can be true "Undetermined", namely unspecified,
    # can be dynamic and has to be determined at run time, or undetermined
    # at this stage because something cannot be determined now.
    step_input = Undetermined()
    step_output = Undetermined()
    step_depends = []
    environ_vars = set()
    signature_vars = set()
    changed_vars = set()
    #
    # 1. execute global definition to get a basic environment
    #
    # FIXME: this could be made much more efficient
    if 'provides' in section.options:
        if '__default_output__' in env.sos_dict:
            step_output = env.sos_dict['__default_output__']
    else:
        #env.sos_dict = WorkflowDict()
        env.sos_dict.set('__null_func__', __null_func__)
        # initial values
        env.sos_dict.set('SOS_VERSION', __version__)
        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)

    #
    # Here we need to get "contant" values from the global section
    # Because parameters are considered variable, they has to be
    # removed. We achieve this by removing function sos_handle_parameter_
    # from the SoS_dict namespace
    #
    if section.global_def:
        try:
            SoS_exec('del sos_handle_parameter_\n' + section.global_def)
        except RuntimeError as e:
            if env.verbosity > 2:
                sys.stderr.write(get_traceback())
            raise RuntimeError(f'Failed to execute statements\n"{section.global_def}"\n{e}')
        finally:
            SoS_exec('from sos.runtime import sos_handle_parameter_', None)

    #
    # 2. look for input statement
    if 'shared' in section.options:
        svars = section.options['shared']
        if isinstance(svars, str):
            changed_vars.add(svars)
            svars = {svars: svars}
        elif isinstance(svars, Sequence):
            for item in svars:
                if isinstance(item, str):
                    changed_vars.add(item)
                elif isinstance(item, Mapping):
                    changed_vars |= set(item.keys())
                else:
                    raise ValueError(
                        f'Option shared should be a string, a mapping of expression, or list of string or mappings. {svars} provided')
        elif isinstance(svars, Mapping):
            changed_vars |= set(svars.keys())
        else:
            raise ValueError(
                f'Option shared should be a string, a mapping of expression, or list of string or mappings. {svars} provided')


    # look for input statement.
    input_statement_idx = [idx for idx,x in enumerate(section.statements) if x[0] == ':' and x[1] == 'input']
    if not input_statement_idx:
        input_statement_idx = None
    elif len(input_statement_idx) == 1:
        input_statement_idx = input_statement_idx[0]
    else:
        raise RuntimeError(f'More than one step input are specified in step {section.name}_{section.index}')

    # if there is an input statement, analyze the statements before it, and then the input statement
    if input_statement_idx is not None:
        # execute before input stuff
        for statement in section.statements[:input_statement_idx]:
            if statement[0] == ':':
                if statement[1] == 'depends':
                    environ_vars |= accessed_vars(statement[2])
                    key, value = statement[1:]
                    try:
                        args, kwargs = SoS_eval(f'__null_func__({value})')
                        if any(isinstance(x, (dynamic, remote)) for x in args):
                            step_depends = Undetermined()
                        else:
                            step_depends = _expand_file_list(True, *args)
                    except Exception as e:
                        env.logger.debug(f"Args {value} cannot be determined: {e}")
                else:
                    raise RuntimeError(f'Step input should be specified before {statement[1]}')
            else:
                environ_vars |= accessed_vars(statement[1])
        #
        # input statement
        stmt = section.statements[input_statement_idx][2]
        try:
            environ_vars |= accessed_vars(stmt)
            args, kwargs = SoS_eval(f'__null_func__({stmt})')

            if not args:
                if default_input is None:
                    step_input = Undetermined()
                else:
                    step_input = default_input
            elif not any(isinstance(x, (dynamic, remote)) for x in args):
                step_input = _expand_file_list(True, *args)
            env.sos_dict.set('input', Undetermined() if isinstance(step_input, Undetermined) else sos_targets(step_input))

            if 'paired_with' in kwargs:
                pw = kwargs['paired_with']
                if pw is None or not pw:
                    pass
                elif isinstance(pw, str):
                    environ_vars.add(pw)
                elif isinstance(pw, Iterable):
                    environ_vars |= set(pw)
                elif isinstance(pw, Iterable):
                    # value supplied, no environ var
                    environ_vars |= set()
                else:
                    raise ValueError(f'Unacceptable value for parameter paired_with: {pw}')
            if 'group_with' in kwargs:
                pw = kwargs['group_with']
                if pw is None or not pw:
                    pass
                elif isinstance(pw, str):
                    environ_vars.add(pw)
                elif isinstance(pw, Iterable):
                    environ_vars |= set(pw)
                elif isinstance(pw, Iterable):
                    # value supplied, no environ var
                    environ_vars |= set()
                else:
                    raise ValueError(f'Unacceptable value for parameter group_with: {pw}')
            if 'for_each' in kwargs:
                fe = kwargs['for_each']
                if fe is None or not fe:
                    pass
                elif isinstance(fe, str):
                    environ_vars |= set([x.strip() for x in fe.split(',')])
                elif isinstance(fe, Sequence):
                    for fei in fe:
                        environ_vars |= set([x.strip() for x in fei.split(',')])
                else:
                    raise ValueError(f'Unacceptable value for parameter fe: {fe}')
        except Exception as e:
            # if anything is not evalutable, keep Undetermined
            env.logger.debug(f'Input of step {section.name}_{section.index} is set to Undertermined: {e}')
            # expression ...
            step_input = Undetermined(stmt)
        input_statement_idx += 1
    else:
        # assuming everything starts from 0 is after input
        input_statement_idx = 0

    # other variables
    for statement in section.statements[input_statement_idx:]:
        # if input is undertermined, we can only process output:
        if statement[0] == '=':
            signature_vars |= accessed_vars('='.join(statement[1:3]))
        elif statement[0] == ':':
            key, value = statement[1:]
            #if key == 'depends':
            environ_vars |= accessed_vars(value)
            # output, depends, and process can be processed multiple times
            try:
                args, kwargs = SoS_eval(f'__null_func__({value})')
                if not any(isinstance(x, (dynamic, remote)) for x in args):
                    if key == 'output':
                        step_output = _expand_file_list(True, *args)
                    elif key == 'depends':
                        step_depends = _expand_file_list(True, *args)
            except Exception as e:
                env.logger.debug(f"Args {value} cannot be determined: {e}")
        else: # statement
            signature_vars |= accessed_vars(statement[1])
    # finally, tasks..
    if section.task:
        signature_vars |= accessed_vars(section.task)
    if '__default_output__' in env.sos_dict and not isinstance(step_output, Undetermined):
        for out in env.sos_dict['__default_output__']:
            if out not in step_output:
                raise ValueError(f'Defined output fail to produce expected output: {step_output} generated, {env.sos_dict["__default_output__"]} expected.')
 
    return {
        'step_name': f'{section.name}_{section.index}' if isinstance(section.index, int) else section.name,
        'step_input': step_input if isinstance(step_input, Undetermined) else sos_targets(step_input),
        'step_output': step_output if isinstance(step_output, Undetermined) else sos_targets(step_output),
        'step_depends': step_depends if isinstance(step_depends, Undetermined) else sos_targets(step_depends),
        # variables starting with __ are internals...
        'environ_vars': {x for x in environ_vars if not x.startswith('__')},
        'signature_vars': {x for x in signature_vars if not x.startswith('__')},
        'changed_vars': changed_vars
        }


class TaskManager:
    # manage tasks created by the step
    def __init__(self, trunk_size, trunk_workers):
        super(TaskManager, self).__init__()
        self.trunk_size = trunk_size
        self.trunk_workers = trunk_workers
        self._submitted_tasks = []
        self._unsubmitted_tasks = []
        # derived from _unsubmitted_tasks
        self._all_ids = []
        self._all_output = []
        #
        self._terminate = False

    def append(self, task_def):
        self._unsubmitted_tasks.append(task_def)
        if isinstance(task_def[2], Sequence):
            self._all_output.extend(task_def[2])
        self._all_ids.append(task_def[0])

    def has_task(self, task_id):
        return task_id in self._all_ids

    def has_output(self, output):
        if not isinstance(output, Sequence) or not self._unsubmitted_tasks:
            return False
        return any(x in self._all_output for x in output)

    def get_job(self, all_tasks=False):
        # save tasks
        if not self._unsubmitted_tasks:
            return None
        # single tasks
        if self.trunk_size == 1 or all_tasks:
            to_be_submitted = self._unsubmitted_tasks
            self._unsubmitted_tasks = []
        else:
            # save complete blocks
            num_tasks = len(self._unsubmitted_tasks) // self.trunk_size * self.trunk_size
            to_be_submitted = self._unsubmitted_tasks[: num_tasks]
            self._unsubmitted_tasks = self._unsubmitted_tasks[ num_tasks:]

        # save tasks
        ids = []
        if self.trunk_size == 1 or (all_tasks and len(self._unsubmitted_tasks) == 1):
            for task_id, taskdef, _ in to_be_submitted:
                job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.def')
                taskdef.save(job_file)
                ids.append(task_id)
        else:
            master = None
            for task_id, taskdef, _ in to_be_submitted:
                if master is not None and master.num_tasks() == self.trunk_size:
                    job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', master.ID + '.def')
                    ids.append(master.ID)
                    master.save(job_file)
                    master = None
                if master is None:
                    master = MasterTaskParams(self.trunk_workers)
                master.push(task_id, taskdef)
            # the last piece
            if master is not None:
                job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', master.ID + '.def')
                master.save(job_file)
                ids.append(master.ID)

        if not ids:
            return None

        self._submitted_tasks.extend(ids)
        return ids

    def clear_submitted(self):
        self._submitted_tasks = []


# overwrite concurrent_execute defined in Base_Step_Executor because sos notebook
# can only handle stdout/stderr from the master process
#
@contextlib.contextmanager
def stdoutIO():
    oldout = sys.stdout
    olderr = sys.stderr
    stdout = StringIO()
    stderr = StringIO()
    sys.stdout = stdout
    sys.stderr = stderr
    yield stdout, stderr
    sys.stdout = oldout
    sys.stderr = olderr

def concurrent_execute(stmt, proc_vars={}, sig=None, capture_output=False):
    '''Execute statements in the passed dictionary'''
    env.sos_dict.quick_update(proc_vars)
    outmsg = ''
    errmsg = ''
    try:
        if sig:
            sig.lock()
        if capture_output:
            with stdoutIO() as (out, err):
                SoS_exec(stmt, return_result=False)
                outmsg = out.getvalue()
                errmsg = err.getvalue()
        else:
            SoS_exec(stmt, return_result=False)
        if sig:
            sig.write()
        res = {'ret_code': 0}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    except (StopInputGroup, TerminateExecution, UnknownTarget, RemovedTarget, UnavailableLock, PendingTasks) as e:
        res = {'ret_code': 1, 'exception': e }
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    except (KeyboardInterrupt, SystemExit) as e:
        # Note that KeyboardInterrupt is not an instance of Exception so this piece is needed for
        # the subprocesses to handle keyboard interrupt. We do not pass the exception
        # back to the master process because the master process would handle KeyboardInterrupt
        # as well and has no chance to handle the returned code.
        procs = psutil.Process().children(recursive=True)
        if procs:
            if env.verbosity > 2:
                env.logger.info(f'{os.getpid()} interrupted. Killing subprocesses {" ".join(str(x.pid) for x in procs)}')
            for p in procs:
                p.terminate()
            gone, alive = psutil.wait_procs(procs, timeout=3)
            if alive:
                for p in alive:
                    p.kill()
            gone, alive = psutil.wait_procs(procs, timeout=3)
            if alive:
                for p in alive:
                    env.logger.warning(f'Failed to kill subprocess {p.pid}')
        elif env.verbosity > 2:
            env.logger.info(f'{os.getpid()} interrupted. No subprocess.')
        raise e
    except Exception as e:
        error_class = e.__class__.__name__
        cl, exc, tb = sys.exc_info()
        msg = ''
        for st in reversed(traceback.extract_tb(tb)):
            if st.filename.startswith('script_'):
                code = stmtHash.script(st.filename)
                line_number = st.lineno
                code = '\n'.join([f'{"---->" if i+1 == line_number else "     "} {x.rstrip()}' for i, x in enumerate(code.splitlines())][max(line_number - 3, 0):line_number + 3])
                msg += f'''\
{st.filename} in {st.name}
{code}
'''
        detail = e.args[0] if e.args else ''
        res = {'ret_code': 1, 'exception': RuntimeError(f'''
---------------------------------------------------------------------------
{error_class:42}Traceback (most recent call last)
{msg}
{error_class}: {detail}''') if msg else RuntimeError(f'{error_class}: {detail}')}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    finally:
        # release the lock even if the process becomes zombie? #871
        if sig:
            sig.release(quiet=True)

class Base_Step_Executor:
    # This base class defines how steps are executed. The derived classes will reimplement
    # some function to behave differently in different modes.
    #
    def __init__(self, step):
        self.step = step
        self.task_manager = None

    def expand_input_files(self, value, *args):
        if self.run_mode == 'dryrun' and any(isinstance(x, dynamic) for x in args):
            return Undetermined(value)

        # if unspecified, use __step_output__ as input (default)
        # resolve dynamic input.
        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        if not args:
            return env.sos_dict['step_input']
        else:
            return _expand_file_list(False, *args)

    def expand_depends_files(self, *args, **kwargs):
        '''handle directive depends'''
        if self.run_mode == 'dryrun' and any(isinstance(x, dynamic) for x in args):
            for k in args:
                if isinstance(k, dynamic):
                    env.logger.warning(f'Dependent target {k} is dynamic')
            return Undetermined()

        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        return _expand_file_list(False, *args)

    def expand_output_files(self, value, *args):
        '''Process output files (perhaps a pattern) to determine input files.
        '''
        if any(isinstance(x, dynamic) for x in args):
            return Undetermined(value)
        else:
            return _expand_file_list(True, *args)

    def verify_input(self):
        if self.run_mode == 'dryrun':
            return
        # now, if we are actually going to run the script, we
        # need to check the input files actually exists, not just the signatures
        for target in (env.sos_dict['_input'] if isinstance(env.sos_dict['_input'], list) else []) + \
            (env.sos_dict['_depends'] if isinstance(env.sos_dict['_depends'], list) else []):
            # if the file does not exist (although the signature exists)
            # request generation of files
            if isinstance(target, str):
                if not file_target(target).target_exists('target'):
                    # remove the signature and regenerate the file
                    file_target(target).remove_sig()
                    raise RemovedTarget(target)
            elif not target.target_exists('target'):
                target.remove_sig()
                raise RemovedTarget(target)

    def verify_output(self):
        if self.run_mode == 'dryrun':
            return
        if env.sos_dict['step_output'] is None:
            return
        if isinstance(env.sos_dict['step_output'], Undetermined):
            raise RuntimeError('Output of a completed step cannot be undetermined.')
        for target in env.sos_dict['step_output']:
            if isinstance(target, str):
                if not file_target(target).target_exists('any'):
                    raise RuntimeError(
                        f'Output target {target} does not exist after the completion of step {env.sos_dict["step_name"]} (curdir={os.getcwd()})')
            elif not target.target_exists('any'):
                raise RuntimeError(
                    f'Output target {target} does not exist after the completion of step {env.sos_dict["step_name"]}')

    def step_signature(self, index):
        '''returns a signature of the step. Change of the step content will
        lead to the invalidation of the signature, which will then cause the
        re-execution of the step for any result from the step. '''
        #
        if env.config['sig_mode'] == 'ignore' or env.config['run_mode'] == 'dryrun':
            return None
        env_vars = []
        for var in sorted(env.sos_dict['__signature_vars__']):
            if var in env.sos_dict and isinstance(env.sos_dict[var], (str, bool, int, float, complex, bytes, list, tuple, set, dict)):
                env_vars.append(f'{var} = {stable_repr(env.sos_dict[var])}\n')

        # env.logger.warning(''.join(env_vars) + '\n' + self.step.tokens)
        return ''.join(env_vars) + '\n' + self.step.tokens

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
                raise ValueError(f'Paired group_by has to have even number of input files: {len(ifiles)} provided')
            return list(list(x) for x in zip(ifiles[:len(ifiles)//2], ifiles[len(ifiles)//2:]))
        elif group_by == 'pairwise':
            f1, f2 = tee(ifiles)
            next(f2, None)
            return [list(x) for x in zip(f1, f2)]
        elif group_by == 'combinations':
            return [list(x) for x in combinations(ifiles, 2)]
        elif isinstance(group_by, int) or group_by.isdigit():
            group_by = int(group_by)
            if len(ifiles) % group_by != 0 and len(ifiles) > group_by:
                env.logger.warning(
                    f'Number of samples ({len(ifiles)}) is not a multiple of group_by ({group_by}). The last group would have less files than the other groups.')
            if group_by < 1:
                raise ValueError('Value of paramter group_by should be a positive number.')
            return [ifiles[i:i + group_by] for i in range(0, len(ifiles), group_by)]
        else:
            raise ValueError(f'Unsupported group_by option ``{group_by}``!')

    @staticmethod
    def handle_paired_with(paired_with, ifiles, _groups, _vars):
        '''Handle input option paired_with'''
        if paired_with is None or not paired_with:
            var_name = []
            var_value = []
        elif isinstance(paired_with, str):
            var_name = ['_' + paired_with]
            if paired_with not in env.sos_dict:
                raise ValueError(f'Variable {paired_with} does not exist.')
            var_value = [env.sos_dict[paired_with]]
        elif isinstance(paired_with, dict):
            var_name = []
            var_value = []
            for k,v in paired_with.items():
                var_name.append(k)
                var_value.append(v)
        elif isinstance(paired_with, Iterable):
            try:
                var_name = ['_'+x for x in paired_with]
            except Exception:
                raise ValueError(f'Invalud value for option paired_with {paired_with}')
            var_value = []
            for vn in var_name:
                if vn[1:] not in env.sos_dict:
                    raise ValueError(f'Variable {vn[1:]} does not exist.')
                var_value.append(env.sos_dict[vn[1:]])
        else:
            raise ValueError(f'Unacceptable value for parameter paired_with: {paired_with}')
        #
        for vn, vv in zip(var_name, var_value):
            if isinstance(vv, str) or not isinstance(vv, Iterable):
                raise ValueError(f'paired_with variable {vn} is not a sequence ("{vv}")')
            if len(vv) != len(ifiles):
                raise ValueError(
                    f'Length of variable {vn} (length {len(vv)}) should match the number of input files (length {len(ifiles)}).')
            file_map = {x:y for x,y in zip(ifiles, vv)}
            for idx, grp in enumerate(_groups):
                mapped_vars = [file_map[x] for x in grp]
                # 862. we make the paired variable the same type so that if the input is a paths or sos_targets,
                # the returned value is of the same type
                _vars[idx][vn] = type(vv)(mapped_vars)

    @staticmethod
    def handle_group_with(group_with, ifiles, _groups, _vars):
        '''Handle input option group_with'''
        if group_with is None or not group_with:
            var_name = []
            var_value = []
        elif isinstance(group_with, str):
            var_name = ['_' + group_with]
            if group_with not in env.sos_dict:
                raise ValueError(f'Variable {group_with} does not exist.')
            var_value = [env.sos_dict[group_with]]
        elif isinstance(group_with, dict):
            var_name = []
            var_value = []
            for k,v in group_with.items():
                var_name.append(k)
                var_value.append(v)
        elif isinstance(group_with, Iterable):
            try:
                var_name = ['_'+x for x in group_with]
            except Exception:
                raise ValueError(f'Invalud value for option group_with {group_with}')
            var_value = []
            for vn in var_name:
                if vn[1:] not in env.sos_dict:
                    raise ValueError(f'Variable {vn[1:]} does not exist.')
                var_value.append(env.sos_dict[vn[1:]])
        else:
            raise ValueError(f'Unacceptable value for parameter group_with: {group_with}')
        #
        for vn, vv in zip(var_name, var_value):
            if isinstance(vv, str) or not isinstance(vv, Iterable):
                raise ValueError(f'group_with variable {vn} is not a sequence ("{vv}")')
            if len(vv) != len(_groups):
                raise ValueError(
                    f'Length of variable {vn} (length {len(vv)}) should match the number of input groups (length {len(_groups)}).')
            for idx, grp in enumerate(_groups):
                _vars[idx][vn] = vv[idx]

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
            raise ValueError(f'Unacceptable value for parameter pattern: {pattern}')
        #
        for pattern in patterns:
            res = extract_pattern(pattern, ifiles)
            # now, assign the variables to env
            for k, v in res.items():
                if k in ('step_input', 'step_output', 'step_depends') or k.startswith('_'):
                    raise RuntimeError(f'Pattern defined variable {k} is not allowed')
                env.sos_dict[k] = v
            # also make k, v pair with _input
            Base_Step_Executor.handle_paired_with(res.keys(), ifiles, _groups, _vars)

    @staticmethod
    def handle_for_each(for_each, _groups, _vars):
        if for_each is None or not for_each:
            for_each = []
        elif isinstance(for_each, (str, dict)):
            for_each = [for_each]
        elif isinstance(for_each, Sequence):
            for_each = for_each
        else:
            raise ValueError(f'Unacceptable value for parameter for_each: {for_each}')
        #
        for fe_all in for_each:
            if isinstance(fe_all, dict):
                # in the format of {'name': value}
                fe_iter_names = []
                fe_values = []
                for k, v in fe_all.items():
                    if ',' in k:
                        names = [x.strip() for x in k.split(',')]
                        if isinstance(v, Iterable):
                            v = list(v)
                        if any(len(_v) != len(names) for _v in v):
                            raise ValueError(
                                f'Unable to unpack object {short_repr(v)} for variables {k} (of length {len(names)})')
                        fe_iter_names.extend(names)
                        fe_values.extend(list(zip(*v)))
                    else:
                        fe_iter_names.append(k)
                        fe_values.append(v)
            else:
                if ',' in fe_all:
                    fe_var_names = [x.strip() for x in fe_all.split(',')]
                    fe_iter_names = ['_' + x for x in fe_var_names]
                else:
                    fe_var_names = [fe_all]
                    fe_iter_names = ['_' + fe_all]
                # check iterator variable name
                for name in fe_iter_names:
                    if '.' in name:
                        raise ValueError(f'Invalid iterator variable {name}')
                # check variables
                fe_values = []
                for name in fe_var_names:
                    if name.split('.')[0] not in env.sos_dict:
                        raise ValueError(f'Variable {name} does not exist.')
                    if '.' in name:
                        fe_values.append(getattr(env.sos_dict[name.split('.')[0]], name.split('.', 1)[-1]))
                    else:
                        fe_values.append(env.sos_dict[name])

            # get loop size
            loop_size = None
            for name, values in zip(fe_iter_names, fe_values):
                if not isinstance(values, Sequence):
                    try:
                        import pandas as pd
                        if not isinstance(values, (pd.DataFrame, pd.Series, pd.Index)):
                            raise ValueError(f'Unacceptable for_each data type {values.__class__.__name__}')
                    except Exception as e:
                        raise ValueError(f'Cannot iterate through variable {name}: {e}')
                if loop_size is None:
                    loop_size = len(values)
                elif loop_size != len(values):
                    raise ValueError(
                        f'Length of variable {name} (length {len(values)}) should match the length of other variables (length {loop_size}).')
            # expand
            _tmp_groups = copy.deepcopy(_groups)
            _groups.clear()
            for _ in range(loop_size):
                _groups.extend(_tmp_groups)
            #
            _tmp_vars = copy.deepcopy(_vars)
            _vars.clear()
            for vidx in range(loop_size):
                for idx, _ in enumerate(_tmp_vars):
                    for var_name, values in zip(fe_iter_names, fe_values):
                        if isinstance(values, Sequence):
                            _tmp_vars[idx][var_name] = values[vidx]
                        elif isinstance(values, (pd.DataFrame, pd.Series)):
                            _tmp_vars[idx][var_name] = values.iloc[vidx]
                        elif isinstance(values, pd.Index):
                            _tmp_vars[idx][var_name] = values[vidx]
                        else:
                            raise ValueError(f'Failed to iterate through for_each variable {short_repr(values)}')
                _vars.extend(copy.deepcopy(_tmp_vars))

    # directive input
    def process_input_args(self, ifiles, **kwargs):
        """This function handles directive input and all its parameters.
        It
            determines and set __step_input__
            determines and set pattern variables if needed
        returns
            _groups
            _vars
        which are groups of _input and related _vars
        """
        if isinstance(ifiles, Undetermined):
            env.sos_dict.set('step_input', Undetermined())
            env.sos_dict.set('_input', Undetermined())
            # temporarily set depends and output to Undetermined because we cannot
            # go far with such input
            env.sos_dict.set('step_output', Undetermined())
            return [Undetermined()], [{}]

        for k in kwargs.keys():
            if k not in SOS_INPUT_OPTIONS:
                raise RuntimeError(f'Unrecognized input option {k}')
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
        env.sos_dict.set('step_input', sos_targets(ifiles))
        env.sos_dict.set('_input', sos_targets(ifiles))
        #
        # handle group_by
        if 'group_by' in kwargs:
            _groups = Base_Step_Executor.handle_group_by(ifiles, kwargs['group_by'])
        else:
            _groups = [ifiles]
        #
        _vars = [{} for x in _groups]
        # handle paired_with
        if 'paired_with' in kwargs:
            Base_Step_Executor.handle_paired_with(kwargs['paired_with'], ifiles,  _groups, _vars)
        # handle pattern
        if 'pattern' in kwargs:
            Base_Step_Executor.handle_extract_pattern(kwargs['pattern'], ifiles, _groups, _vars)
        # handle group_with
        if 'group_with' in kwargs:
            Base_Step_Executor.handle_group_with(kwargs['group_with'], ifiles,  _groups, _vars)
        # handle for_each
        if 'for_each' in kwargs:
            Base_Step_Executor.handle_for_each(kwargs['for_each'], _groups, _vars)
        return _groups, _vars

    def process_depends_args(self, dfiles, **kwargs):
        for k in kwargs.keys():
            if k not in SOS_DEPENDS_OPTIONS:
                raise RuntimeError(f'Unrecognized depends option {k}')
        if isinstance(dfiles, Undetermined):
            raise ValueError(r"FIXME depends needs to handle undetermined")
        if isinstance(dfiles, sos_targets):
            raise ValueError(r"FIXME depends should not be targets type")

        env.sos_dict.set('_depends', sos_targets(dfiles))
        if env.sos_dict['step_depends'] is None:
            env.sos_dict.set('step_depends', sos_targets(dfiles))
        elif env.sos_dict['step_depends'].targets() != dfiles:
            env.sos_dict['step_depends'].extend(dfiles)

    def process_output_args(self, ofiles, **kwargs):
        if isinstance(ofiles, sos_targets):
            raise ValueError(r"FIXME output should not be targets type")

        for k in kwargs.keys():
            if k not in SOS_OUTPUT_OPTIONS:
                raise RuntimeError(f'Unrecognized output option {k}')
        # create directory
        if not isinstance(ofiles, Undetermined):
            for ofile in ofiles:
                if isinstance(ofile, str):
                    parent_dir = os.path.split(os.path.expanduser(ofile))[0]
                    if parent_dir and not os.path.isdir(parent_dir):
                        os.makedirs(parent_dir)

        if 'group_by' in kwargs:
            _ogroups = Base_Step_Executor.handle_group_by(ofiles, kwargs['group_by'])
            if len(_ogroups) != len(self._groups):
                raise RuntimeError(
                    f'Output option group_by produces {len(_ogroups)} output groups which is different from the number of input groups ({len(self._groups)}).')
            ofiles = _ogroups[env.sos_dict['_index']]

        if isinstance(ofiles, sos_targets):
            raise ValueError(r"FIXME output should not be targets type")

        # set variables
        env.sos_dict.set('_output', ofiles if isinstance(ofiles, Undetermined) else sos_targets(ofiles))
        #
        if isinstance(env.sos_dict['step_output'], (type(None), Undetermined)):
            env.sos_dict.set('step_output', copy.deepcopy(ofiles))
        elif not isinstance(env.sos_dict['step_output'], Undetermined) and \
            env.sos_dict['step_output'].targets() != ofiles:
            env.sos_dict['step_output'].extend(ofiles)

    def process_task_args(self, **kwargs):
        env.sos_dict.set('_runtime', {})
        for k,v in kwargs.items():
            if k not in SOS_RUNTIME_OPTIONS:
                raise RuntimeError(f'Unrecognized runtime option {k}={v}')
            # standardize walltime to an integer
            if k == 'walltime':
                v = format_HHMMSS(v)
            elif k == 'mem':
                v = expand_size(v)
            env.sos_dict['_runtime'][k] = v

    def reevaluate_output(self):
        # re-process the output statement to determine output files
        if isinstance(env.sos_dict['step_output'], Undetermined):
            args, _ = SoS_eval(f'__null_func__({env.sos_dict["output"].expr})')
            # handle dynamic args
            args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
            env.sos_dict.set('step_output', sos_targets(self.expand_output_files('', *args)))
        else:
            # this needs to be re-visited
            for target in env.sos_dict['step_output'].targets():
                if not isinstance(target, Undetermined):
                    continue
                args, _ = SoS_eval(f'__null_func__({target.expr})')
                # handle dynamic args
                args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
                env.sos_dict.set('step_output', sos_targets(self.expand_output_files('', *args)))


    def prepare_task(self):
        env.sos_dict['_runtime']['cur_dir'] = os.getcwd()
        # we need to record the verbosity and sigmode of task during creation because
        # they might be changed while the task is in the queue waiting to be
        # submitted (this happens when tasks are submitted from Jupyter)
        env.sos_dict['_runtime']['verbosity'] = env.verbosity
        env.sos_dict['_runtime']['sig_mode'] = env.config.get('sig_mode', 'default')
        env.sos_dict['_runtime']['run_mode'] = env.config.get('run_mode', 'run')
        env.sos_dict['_runtime']['home_dir'] = os.path.expanduser('~')
        if 'workdir' in env.sos_dict['_runtime'] and not os.path.isdir(os.path.expanduser(env.sos_dict['_runtime']['workdir'])):
            try:
                os.makedirs(os.path.expanduser(env.sos_dict['_runtime']['workdir']))
            except Exception:
                raise RuntimeError(f'Failed to create workdir {env.sos_dict["_runtime"]["workdir"]}')

        # NOTE: we do not explicitly include 'step_input', 'step_output',
        # 'step_depends' and 'CONFIG'
        # because they will be included by env.sos_dict['__signature_vars__'] if they are actually
        # used in the task. (issue #752)
        task_vars = env.sos_dict.clone_selected_vars(env.sos_dict['__signature_vars__'] \
                    | {'_input', '_output', '_depends', '_index', '__args__', 'step_name', '_runtime',
                    '__signature_vars__', '__step_context__'
                    })

        task_tags = [env.sos_dict.get('step_name', ''), os.path.basename(env.sos_dict.get('__workflow_sig__', '')).rsplit('.', 1)[0]]
        if 'tags' in env.sos_dict['_runtime']:
            if isinstance(env.sos_dict['_runtime']['tags'], str):
                tags = [env.sos_dict['_runtime']['tags']]
            elif isinstance(env.sos_dict['_runtime']['tags'], Sequence):
                tags = list(env.sos_dict['_runtime']['tags'])
            else:
                env.logger.warning(f'Unacceptable value for parameter tags: {env.sos_dict["_runtime"]["tags"]}')
            #
            for tag in tags:
                if not tag.strip():
                    continue
                if not SOS_TAG.match(tag):
                    new_tag = re.sub(r'[^\w_.-]', '', tag)
                    if new_tag:
                        env.logger.warning(f'Invalid tag "{tag}" is added as "{new_tag}"')
                        task_tags.append(new_tag)
                    else:
                        env.logger.warning(f'Invalid tag "{tag}" is ignored')
                else:
                    task_tags.append(tag)

        # save task to a file
        taskdef = TaskParams(
            name = '{} (index={})'.format(self.step.step_name(), env.sos_dict['_index']),
            global_def = self.step.global_def,
            task = self.step.task,          # task
            sos_dict = task_vars,
            tags = task_tags
        )
        # if no output (thus no signature)
        # temporarily create task signature to obtain sig_id
        task_id = RuntimeInfo(self.step.md5, self.step.task, task_vars['_input'].targets(),
            task_vars['_output'].targets(), task_vars['_depends'].targets(),
            task_vars['__signature_vars__'], task_vars).sig_id

        # workflow ID should be included but not part of the signature, this is why it is included
        # after task_id is created.
        if '__workflow_sig__' in env.sos_dict and env.sos_dict['__workflow_sig__']:
            task_vars['__workflow_sig__'] = env.sos_dict['__workflow_sig__']


        if self.task_manager is None:
            if 'trunk_size' in env.sos_dict['_runtime']:
                if not isinstance(env.sos_dict['_runtime']['trunk_size'], int):
                    raise ValueError(
                        f'An integer value is expected for runtime option trunk, {env.sos_dict["_runtime"]["trunk_size"]} provided')
                trunk_size = env.sos_dict['_runtime']['trunk_size']
            else:
                trunk_size = 1
            if 'trunk_workers' in env.sos_dict['_runtime']:
                if not isinstance(env.sos_dict['_runtime']['trunk_workers'], int):
                    raise ValueError(
                        f'An integer value is expected for runtime option trunk_workers, {env.sos_dict["_runtime"]["trunk_workers"]} provided')
                trunk_workers = env.sos_dict['_runtime']['trunk_workers']
            else:
                trunk_workers = 0

            #if 'queue' in env.sos_dict['_runtime'] and env.sos_dict['_runtime']['queue']:
            #    host = env.sos_dict['_runtime']['queue']
            #else:
            #    # otherwise, use workflow default
            #    host = '__default__'

            self.task_manager = TaskManager(trunk_size, trunk_workers)

        #618
        # it is possible that identical tasks are executed (with different underlying random numbers)
        # we should either give a warning or produce different ids...
        if self.task_manager.has_task(task_id):
            raise RuntimeError(f'Identical task generated from _index={env.sos_dict["_index"]}.')
        elif self.task_manager.has_output(task_vars['_output']):
            raise RuntimeError(
                f'Task produces output files {", ".join(task_vars["_output"])} that are output of other tasks.')
        # if no trunk_size, the job will be submitted immediately
        # otherwise tasks will be accumulated and submitted in batch
        self.task_manager.append((task_id, taskdef, task_vars['_output'].targets()))
        tasks = self.task_manager.get_job()
        if tasks:
            self.submit_tasks(tasks)
        return task_id

    def wait_for_results(self):
        if self.concurrent_input_group:
            self.proc_results = [x.get() for x in self.proc_results]
            SlotManager().release(self.worker_pool._processes - 1)
            self.worker_pool.close()
            self.worker_pool.join()
            self.worker_pool = None
            return

        if self.task_manager is None:
            return {}

        # submit the last batch of tasks
        tasks = self.task_manager.get_job(all_tasks=True)
        if tasks:
            self.submit_tasks(tasks)

        # waiting for results of specified IDs
        results = self.wait_for_tasks(self.task_manager._submitted_tasks)
        self.task_manager.clear_submitted()
        for idx, task in enumerate(self.proc_results):
            # if it is done
            if isinstance(task, dict):
                continue
            if task in results:
                self.proc_results[idx] = results[task]
            else:
                # can be a subtask
                for _, mres in results.items():
                    if 'subtasks' in mres and task in mres['subtasks']:
                        self.proc_results[idx] = mres['subtasks'][task]
        #
        # check if all have results?
        if any(isinstance(x, str) for x in self.proc_results):
            raise RuntimeError(
                f'Failed to get results for tasks {", ".join(x for x in self.proc_results if isinstance(x, str))}')
        #
        # now, if the task has shared variable, merge to sos_dict
        shared = {}
        for res in self.proc_results:
            #
            # shared looks like: {0: {'a': 100}} where the first 0 is _index
            # we need to convert it to {'a': {0: 100}}
            env.logger.debug(f'Collect shared result {res["shared"]}')
            for idx, sh in res['shared'].items():
                for k, v in sh.items():
                    if k in shared:
                        shared[k][idx] = v
                    else:
                        shared[k] = {idx: v}
        env.sos_dict.update(shared)


    def log(self, stage=None, msg=None):
        if stage == 'start':
            env.logger.info(
                f'{"Checking" if self.run_mode == "dryrun" else "Executing"} ``{self.step.step_name(True)}``: {self.step.comment.strip()}')
        elif stage == 'input statement':
            env.logger.trace(f'Handling input statement {msg}')
        elif stage == '_input':
            if env.sos_dict['_input'] is not None:
                env.logger.debug(f'_input: ``{short_repr(env.sos_dict["_input"])}``')
        elif stage == '_depends':
            if env.sos_dict['_depends'] is not None:
                env.logger.debug(f'_depends: ``{short_repr(env.sos_dict["_depends"])}``')
        elif stage == 'input':
            if env.sos_dict['step_input'] is not None:
                env.logger.info(f'input:   ``{short_repr(env.sos_dict["step_input"])}``')
        elif stage == 'output':
            if env.sos_dict['step_output'] is not None:
                env.logger.info(f'output:   ``{short_repr(env.sos_dict["step_output"])}``')

    def execute(self, stmt, sig=None):
        try:
            if sig is None:
                env.sos_dict.set('__step_sig__', None)
            else:
                env.sos_dict.set('__step_sig__', os.path.basename(sig.proc_info).split('.')[0])
            self.last_res = SoS_exec(stmt, return_result=self.run_mode == 'interactive')
        except (StopInputGroup, TerminateExecution, UnknownTarget, RemovedTarget, UnavailableLock, PendingTasks):
            raise
        except Exception as e:
            error_class = e.__class__.__name__
            cl, exc, tb = sys.exc_info()
            msg = ''
            for st in reversed(traceback.extract_tb(tb)):
                if st.filename.startswith('script_'):
                    code = stmtHash.script(st.filename)
                    line_number = st.lineno
                    code = '\n'.join([f'{"---->" if i+1 == line_number else "     "} {x.rstrip()}' for i, x in enumerate(code.splitlines())][max(line_number - 3, 0):line_number + 3])
                    msg += f'''\
{st.filename} in {st.name}
{code}
'''
            detail = e.args[0] if e.args else ''
            if msg:
                raise RuntimeError(f'''
---------------------------------------------------------------------------
{error_class:42}Traceback (most recent call last)
{msg}
{error_class}: {detail}''')
            else:
                raise RuntimeError(f'{error_class}: {detail}')
        finally:
            env.sos_dict.set('__step_sig__', None)


    def collect_result(self):
        # only results will be sent back to the master process
        #
        # __step_input__:    input of this step
        # __steo_output__:   output of this step
        # __step_depends__:  dependent files of this step
        result = {
            '__step_input__': env.sos_dict['step_input'],
            '__step_output__': env.sos_dict['step_output'],
            '__step_depends__': env.sos_dict['step_depends'],
            '__step_name__': env.sos_dict['step_name'],
        }
        result['__last_res__'] = self.last_res
        result['__changed_vars__'] = set()
        result['__shared__'] = {}
        if 'shared' in self.step.options:
            rvars = self.step.options['shared']
            if isinstance(rvars, str):
                result['__changed_vars__'].add(rvars)
                result['__shared__'][rvars] = copy.deepcopy(env.sos_dict[rvars])
            elif isinstance(rvars, Mapping):
                result['__changed_vars__'] |= rvars.keys()
                for var in rvars.keys():
                    result['__shared__'][var] = copy.deepcopy(env.sos_dict[var])
            elif isinstance(rvars, Sequence):
                for item in rvars:
                    if isinstance(item, str):
                        result['__changed_vars__'].add(item)
                        result['__shared__'][item] = copy.deepcopy(env.sos_dict[item])
                    elif isinstance(item, Mapping):
                        result['__changed_vars__'] |= item.keys()
                        for var in item.keys():
                            result['__shared__'][var] = copy.deepcopy(env.sos_dict[var])
                    else:
                        raise ValueError(
                            f'Option shared should be a string, a mapping of expression, or a list of string or mappings. {rvars} provided')
            else:
                raise ValueError(
                    f'Option shared should be a string, a mapping of expression, or a list of string or mappings. {rvars} provided')

        if hasattr(env, 'accessed_vars'):
            result['__environ_vars__'] = self.environ_vars
            result['__signature_vars__'] = env.accessed_vars
        return result

    def run(self):
        '''Execute a single step and return results. The result for batch mode is the
        input, output etc returned as alias, and for interactive mode is the return value
        of the last expression. '''
        # return value of the last executed statement
        self.last_res = None
        #
        self.log('start')
        #
        # prepare environments, namely variables that can be used by the step
        #
        # * step_name:  name of the step, can be used by step process to determine
        #               actions dynamically.
        env.sos_dict.set('step_name', self.step.step_name())
        env.sos_dict.set('step_id', self.step.md5)
        env.sos_dict.set('master_id', env.config.get('master_md5', ''))
        try:
            env.sos_dict.set('workflow_id', os.path.split(env.sos_dict['__workflow_sig__'])[-1].split('.')[0])
        except Exception as e:
            env.logger.debug(f'Failed to set workflow_id for step {self.step.step_name()}: {e}')
            env.sos_dict.set('workflow_id', env.config.get('master_md5', ''))
        # used by nested workflow
        env.sos_dict.set('__step_context__', self.step.context)

        env.sos_dict.set('_runtime', {})

        # * input:      input files, which should be __step_output__ if it is defined, or
        #               None otherwise.
        # * _input:     first batch of input, which should be input if no input statement is used
        # * output:     None at first, can be redefined by output statement
        # * _output:    None at first, can be redefined by output statement
        # * depends:    None at first, can be redefined by depends statement
        # * _depends:   None at first, can be redefined by depends statement
        #
        if '__step_output__' not in env.sos_dict:
            env.sos_dict.set('step_input', sos_targets())
        else:
            if env.sos_dict['__step_output__'] is not None:
                if isinstance(env.sos_dict['__step_output__'], (Undetermined, sos_targets)):
                    env.sos_dict.set('step_input', copy.deepcopy(env.sos_dict['__step_output__']))
                elif isinstance(env.sos_dict['__step_output__'], list):
                    #env.logger.warning(f"__step_output__ should not be list")
                    env.sos_dict.set('step_input', sos_targets(env.sos_dict['__step_output__']))
                else:
                    raise RuntimeError('__step_output__ can only be None, Undetermined, or a list of files.')
            else:
                env.sos_dict.set('step_input', sos_targets())

        # input can be Undetermined from undetermined output from last step
        env.sos_dict.set('_input', copy.deepcopy(env.sos_dict['step_input']))
        if '__default_output__' in env.sos_dict:
            #if not isinstance(env.sos_dict['__default_output__'], sos_targets):
            #    env.logger.warning("__default_output__ should be sos_targets")
            env.sos_dict.set('step_output', sos_targets(env.sos_dict['__default_output__']))
            env.sos_dict.set('_output', sos_targets(env.sos_dict['__default_output__']))
        else:
            env.sos_dict.set('step_output', sos_targets())
            env.sos_dict.set('_output', sos_targets())
        env.sos_dict.set('step_depends', sos_targets())
        env.sos_dict.set('_depends', sos_targets())
        # _index is needed for pre-input action's active option and for debug output of scripts
        env.sos_dict.set('_index', 0)

        # look for input statement.
        input_statement_idx = [idx for idx,x in enumerate(self.step.statements) if x[0] == ':' and x[1] == 'input']
        if not input_statement_idx:
            input_statement_idx = None
        elif len(input_statement_idx) == 1:
            input_statement_idx = input_statement_idx[0]
        else:
            raise ValueError(f'More than one step input are specified in step {self.step.step_name()}')

        # if there is an input statement, execute the statements before it, and then the input statement
        self.concurrent_input_group = False
        if input_statement_idx is not None:
            # execute before input stuff
            for statement in self.step.statements[:input_statement_idx]:
                if statement[0] == ':':
                    key, value = statement[1:]
                    if key != 'depends':
                        raise ValueError(f'Step input should be specified before {key}')
                    try:
                        args, kwargs = SoS_eval(f'__null_func__({value})')
                        dfiles = self.expand_depends_files(*args)
                        # dfiles can be Undetermined
                        self.process_depends_args(dfiles, **kwargs)
                    except (UnknownTarget, RemovedTarget, UnavailableLock):
                        raise
                    except Exception as e:
                        raise RuntimeError(f'Failed to process step {key}: {value.strip()} ({e})')
                else:
                    try:
                        self.execute(statement[1])
                    except StopInputGroup as e:
                        if e.message:
                            env.logger.warning(e)
                        return self.collect_result()
            # input statement
            stmt = self.step.statements[input_statement_idx][2]
            self.log('input statement', stmt)
            try:
                args, kwargs = SoS_eval(f"__null_func__({stmt})")
                # Files will be expanded differently with different running modes
                input_files = self.expand_input_files(stmt, *args)
                self._groups, self._vars = self.process_input_args(input_files, **kwargs)
                self.concurrent_input_group = 'concurrent' in kwargs and kwargs['concurrent'] and len(self._groups) > 1
            except (UnknownTarget, RemovedTarget, UnavailableLock):
                raise
            except Exception as e:
                raise ValueError(f'Failed to process input statement {stmt}: {e}')

            input_statement_idx += 1
        else:
            # default case
            self._groups = [env.sos_dict['step_input'].targets()]
            self._vars = [{}]
            # assuming everything starts from 0 is after input
            input_statement_idx = 0

        self.log('input')

        self.proc_results = []
        # run steps after input statement, which will be run multiple times for each input
        # group.
        env.sos_dict.set('__num_groups__', len(self._groups))

        # determine if a single index or the whole step should be skipped
        skip_index = False
        # signatures of each index, which can remain to be None if no output
        # is defined.
        signatures = [None for x in self._groups]
        self.output_groups = [[] for x in self._groups]

        if self.concurrent_input_group:
            if self.step.task:
                self.concurrent_input_group = False
                env.logger.debug('Input groups are executed sequentially because of existence of tasks')
            else:
                conc_stmts = [x for x in self.step.statements[input_statement_idx:] if x[0] != ':']
                if len(conc_stmts) > 1:
                    self.concurrent_input_group = False
                    env.logger.debug('Input groups are executed sequentially because of existence of directives between statements.')
                elif any('sos_run' in x[1] for x in self.step.statements[input_statement_idx:]):
                    self.concurrent_input_group = False
                    env.logger.debug('Input groups are executed sequentially because of existence of nested workflow.')
                else:
                    sm = SlotManager()
                    # because the master process pool will count one worker in (step)
                    gotten = sm.acquire(len(self._groups) - 1, env.config.get('max_procs', max(int(os.cpu_count() / 2), 1)))
                    if gotten == 0:
                        env.logger.debug(f'Input group executed sequencially due to -j constraint')
                        self.concurrent_input_group = False
                    else:
                        env.logger.debug(f'Using process pool with size {gotten+1}')
                        self.worker_pool = Pool(gotten + 1)

        try:
            for idx, (g, v) in enumerate(zip(self._groups, self._vars)):
                # other variables
                #
                env.sos_dict.update(v)
                env.sos_dict.set('_input', sos_targets(g))
                self.log('_input')
                env.sos_dict.set('_index', idx)
                #
                pre_statement = []
                if not any(st[0] == ':' and st[1] == 'output' for st in self.step.statements[input_statement_idx:]) and \
                    '__default_output__' in env.sos_dict:
                    pre_statement = [[':', 'output', '_output']]

                for statement in pre_statement + self.step.statements[input_statement_idx:]:
                    # if input is undertermined, we can only process output:
                    if isinstance(g, Undetermined) and statement[0] != ':':
                        return self.collect_result()
                    if statement[0] == ':':
                        key, value = statement[1:]
                        # output, depends, and process can be processed multiple times
                        try:
                            args, kwargs = SoS_eval(f'__null_func__({value})')
                            # dynamic output or dependent files
                            if key == 'output':
                                # if output is defined, its default value needs to be cleared
                                if idx == 0:
                                    env.sos_dict.set('step_output', sos_targets())
                                ofiles = self.expand_output_files(value, *args)
                                if not isinstance(g, (type(None), Undetermined)) and not isinstance(ofiles, (type(None), Undetermined)):
                                    if any(x in g for x in ofiles):
                                        raise RuntimeError(
                                            f'Overlapping input and output files: {", ".join(repr(x) for x in ofiles if x in g)}')
                                # set variable _output and output
                                self.process_output_args(ofiles, **kwargs)
                                self.output_groups[idx] = env.sos_dict['_output'].targets()

                                # ofiles can be Undetermined
                                sg = self.step_signature(idx)
                                if sg is not None and not isinstance(g, Undetermined):
                                    signatures[idx] = RuntimeInfo(self.step.md5, sg, env.sos_dict['_input'].targets(),
                                        env.sos_dict['_output'].targets(), env.sos_dict['_depends'].targets(),
                                        env.sos_dict['__signature_vars__'])
                                    signatures[idx].lock()
                                    if env.config['sig_mode'] == 'default':
                                        # if users use sos_run, the "scope" of the step goes beyong names in this step
                                        # so we cannot save signatures for it.
                                        if 'sos_run' in env.sos_dict['__signature_vars__']:
                                            skip_index = False
                                        else:
                                            matched = signatures[idx].validate()
                                            if isinstance(matched, dict):
                                                # in this case, an Undetermined output can get real output files
                                                # from a signature
                                                env.sos_dict.set('_input', sos_targets(matched['input']))
                                                env.sos_dict.set('_depends', sos_targets(matched['depends']))
                                                env.sos_dict.set('_output', sos_targets(matched['output']))
                                                env.sos_dict.update(matched['vars'])
                                                env.logger.info(
                                                    f'``{self.step.step_name(True)}`` (index={idx}) is ``ignored`` due to saved signature')
                                                skip_index = True
                                            else:
                                                env.logger.debug(f'Signature mismatch: {matched}')
                                    elif env.config['sig_mode'] == 'assert':
                                        matched = signatures[idx].validate()
                                        if isinstance(matched, str):
                                            raise RuntimeError(f'Signature mismatch: {matched}')
                                        else:
                                            env.sos_dict.set('_input', sos_targets(matched['input']))
                                            env.sos_dict.set('_depends', sos_targets(matched['depends']))
                                            env.sos_dict.set('_output', sos_targets(matched['output']))
                                            env.sos_dict.update(matched['vars'])
                                            env.logger.info(
                                                f'Step ``{self.step.step_name(True)}`` (index={idx}) is ``ignored`` with matching signature')
                                            skip_index = True
                                    elif env.config['sig_mode'] == 'build':
                                        # build signature require existence of files
                                        if 'sos_run' in env.sos_dict['__signature_vars__']:
                                            skip_index = True
                                        elif signatures[idx].write(rebuild=True):
                                            env.logger.info(
                                                f'Step ``{self.step.step_name(True)}`` (index={idx}) is ``ignored`` with signature constructed')
                                            skip_index = True
                                    elif env.config['sig_mode'] == 'force':
                                        skip_index = False
                                    else:
                                        raise RuntimeError(f'Unrecognized signature mode {env.config["sig_mode"]}')
                                if skip_index:
                                    break
                            elif key == 'depends':
                                try:
                                    dfiles = self.expand_depends_files(*args)
                                    # dfiles can be Undetermined
                                    self.process_depends_args(dfiles, **kwargs)
                                    self.log('_depends')
                                except Exception as e:
                                    env.logger.info(e)
                                    raise
                            elif key == 'task':
                                self.process_task_args(*args, **kwargs)
                            else:
                                raise RuntimeError(f'Unrecognized directive {key}')
                        except (UnknownTarget, RemovedTarget, UnavailableLock):
                            raise
                        except Exception as e:
                            # if input is Undertermined, it is possible that output cannot be processed
                            # due to that, and we just return
                            if isinstance(g, Undetermined):
                                return self.collect_result()
                            raise RuntimeError(f'Failed to process step {key}: {value.strip()} ({e})')
                    else:
                        try:
                            self.verify_input()
                            if self.concurrent_input_group:
                                proc_vars = env.sos_dict.clone_selected_vars(env.sos_dict['__signature_vars__'] \
                                    | {'_input', '_output', '_depends', '_index', '__args__', 'step_name', '_runtime',
                                    '__signature_vars__', '__step_context__'
                                    })

                                if signatures[idx] is None:
                                    proc_vars['__step_sig__'] = None
                                else:
                                    proc_vars['__step_sig__'] = os.path.basename(signatures[idx].proc_info).split('.')[0]
                                    # we need to release the signature otherwise there can be too many opened
                                    # signatures for concurrent jobs
                                    signatures[idx].release()

                                self.proc_results.append(
                                        self.worker_pool.apply_async(concurrent_execute, kwds=dict(stmt=statement[1],
                                            proc_vars=proc_vars, sig=signatures[idx],
                                            capture_output= self.run_mode == 'interactive')))
                                # signature will be written by the concurrent executor
                                signatures[idx] = None
                            else:
                                self.execute(statement[1], signatures[idx])
                        except StopInputGroup as e:
                            self.output_groups[idx] = []
                            if e.message:
                                env.logger.info(e)
                            skip_index = True
                            break

                # if this index is skipped, go directly to the next one
                if skip_index:
                    skip_index = False
                    if signatures[idx]:
                        signatures[idx].release()
                        signatures[idx] = None
                    continue

                # if concurrent input group, there is no task
                if self.concurrent_input_group:
                    continue
                # finally, tasks..
                # now the regular step process is done and we are going to the task part
                # we should be able to release the signature because external task has its own signatures

                if signatures[idx] is not None:
                    signatures[idx].release()

                if not self.step.task:
                    if signatures[idx] is not None:
                        if 'sos_run' not in env.sos_dict['__signature_vars__']:
                            signatures[idx].write()
                        signatures[idx] = None
                    continue

                # check if the task is active
                if 'active' in env.sos_dict['_runtime']:
                    active = env.sos_dict['_runtime']['active']
                    if active is True:
                        pass
                    elif active is False:
                        continue
                    elif isinstance(active, int):
                        if active >= 0 and env.sos_dict['_index'] != active:
                            continue
                        if active < 0 and env.sos_dict['_index'] != active + env.sos_dict['__num_groups__']:
                            continue
                    elif isinstance(active, Sequence):
                        allowed_index = list([x if x >= 0 else env.sos_dict['__num_groups__'] + x for x in active])
                        if env.sos_dict['_index'] not in allowed_index:
                            continue
                    elif isinstance(active, slice):
                        allowed_index = list(range(env.sos_dict['__num_groups__']))[active]
                        if env.sos_dict['_index'] not in allowed_index:
                            continue
                    else:
                        raise RuntimeError(f'Unacceptable value for option active: {active}')

                #
                self.log('task')
                try:
                    task = self.prepare_task()
                    self.proc_results.append(task)
                except Exception as e:
                    # FIXME: cannot catch exception from subprocesses
                    if env.verbosity > 2:
                        sys.stderr.write(get_traceback())
                    raise RuntimeError(f'Failed to execute process\n"{short_repr(self.step.task)}"\n{e}')
                #
                # if not concurrent, we have to wait for the completion of the task
                if 'concurrent' in env.sos_dict['_runtime'] and env.sos_dict['_runtime']['concurrent'] is False:
                    self.wait_for_results()
                #
                # endfor loop for each input group
                #
            # check results? This is only meaningful for pool
            self.wait_for_results()
            for idx,res in enumerate(self.proc_results):
                if signatures[idx] is not None:
                    if res['ret_code'] == 0:
                        signatures[idx].write()
                    signatures[idx] = None
            # check results
            for proc_result in [x for x in self.proc_results if x['ret_code'] == 0]:
                if 'stdout' in proc_result and proc_result['stdout']:
                    sys.stdout.write(proc_result['stdout'])
                if 'stderr' in proc_result and proc_result['stderr']:
                    sys.stderr.write(proc_result['stderr'])
            for proc_result in [x for x in self.proc_results if x['ret_code'] != 0]:
                if 'stdout' in proc_result and proc_result['stdout']:
                    sys.stdout.write(proc_result['stdout'])
                if 'stderr' in proc_result and proc_result['stderr']:
                    sys.stderr.write(proc_result['stderr'])
                if 'exception' in proc_result:
                    raise proc_result['exception']
            # if output is Undetermined, re-evalulate it
            #
            # NOTE: dynamic output is evaluated at last, so it sets output,
            # not _output. For the same reason, signatures can be wrong if it has
            # Undetermined output.
            if env.config['run_mode'] in ('run', 'interactive'):
                if isinstance(env.sos_dict['step_output'], Undetermined) or \
                env.sos_dict['step_output'].has_undetermined():
                    self.reevaluate_output()
                    # if output is no longer Undetermined, set it to output
                    # of each signature
                    for sig in signatures:
                        if sig is not None:
                            sig.set(env.sos_dict['step_output'], 'output')
                else:
                    # finalize output from output_groups because some output might be skipped
                    # this is the final version of the output but we do maintain output
                    # during the execution of step, for compatibility.
                    env.sos_dict.set('step_output', sos_targets(self.output_groups[0]))
                    for og in self.output_groups[1:]:
                        if og != env.sos_dict['step_output'].targets():
                            env.sos_dict['step_output'].extend(og)

            self.log('output')
            # variables defined by the shared option needs to be available to be verified
            if 'shared' in self.step.options:
                if isinstance(self.step.options['shared'], Mapping):
                    for var, val in self.step.options['shared'].items():
                        if var == val:
                            continue
                        try:
                            env.sos_dict.set(var, SoS_eval(val))
                        except Exception as e:
                            raise RuntimeError(f'Failed to evaluate shared variable {var} from expression {val}: {e}')
                # if there are dictionaries in the sequence, e.g.
                # shared=['A', 'B', {'C':'D"}]
                elif isinstance(self.step.options['shared'], Sequence):
                    for item in self.step.options['shared']:
                        if isinstance(item, Mapping):
                            for var, val in item.items():
                                if var == val:
                                    continue
                                try:
                                    env.sos_dict.set(var, SoS_eval(val))
                                except Exception as e:
                                    raise RuntimeError(
                                        f'Failed to evaluate shared variable {var} from expression {val}: {e}')
            #
            self.verify_output()
            return self.collect_result()
        except KeyboardInterrupt:
            # if the worker_pool is not properly shutdown (e.g. interrupted by KeyboardInterrupt #871)
            # we try to kill all subprocesses. Because it takes time to shutdown all processes, impatient
            # users might hit Ctrl-C again, interrupting the shutdown process again. In this case we
            # simply catch the KeyboardInterrupt exception and try again.
            #
            if self.concurrent_input_group and self.worker_pool:
                if env.verbosity > 2:
                    env.logger.info(f'{os.getpid()} terminating worker pool')
                while self.worker_pool:
                    try:
                        self.worker_pool.terminate()
                        SlotManager().release(self.worker_pool._processes - 1)
                        self.worker_pool = None
                    except KeyboardInterrupt:
                        continue
        finally:
            # release all signatures
            for sig in signatures:
                if sig is not None:
                    try:
                        sig.release()
                    except:
                        pass


def _expand_file_list(ignore_unknown, *args):
    ifiles = []

    for arg in args:
        if arg is None:
            continue
        elif isinstance(arg, BaseTarget):
            ifiles.append(arg)
        elif isinstance(arg, (str, path)):
            ifiles.append(os.path.expanduser(arg))
        elif isinstance(arg, sos_targets):
            ifiles.extend(arg.targets())
        elif isinstance(arg, paths):
            ifiles.extend(arg.paths())
        elif isinstance(arg, Iterable):
            # in case arg is a Generator, check its type will exhaust it
            arg = list(arg)
            if not all(isinstance(x, (str, path, BaseTarget)) for x in arg):
                raise RuntimeError(f'Invalid target: {arg}')
            ifiles.extend(arg)
        else:
            raise RuntimeError(f'Unrecognized file: {arg} of type {type(arg).__name__}')

    # expand files with wildcard characters and check if files exist
    tmp = []
    for ifile in ifiles:
        if isinstance(ifile, BaseTarget):
            if ignore_unknown or ifile.target_exists():
                tmp.append(ifile)
            else:
                raise UnknownTarget(ifile)
        elif file_target(ifile).target_exists('target'):
            tmp.append(ifile)
        elif file_target(ifile).target_exists('signature'):
            env.logger.debug(f'``{ifile}`` exists in signature form (actual target has been removed).')
            tmp.append(ifile)
        elif isinstance(ifile, sos_targets):
            raise ValueError("sos_targets should not appear here")
        else:
            expanded = sorted(glob.glob(os.path.expanduser(ifile)))
            # no matching file ... but this is not a problem at the
            # inspection stage.
            #
            # NOTE: if a DAG is given, the input file can be output from
            # a previous step..
            #
            if not expanded:
                if not ignore_unknown:
                    raise UnknownTarget(ifile)
                else:
                    tmp.append(ifile)
            else:
                tmp.extend(expanded)
    return tmp



class Step_Executor(Base_Step_Executor):
    '''Single process step executor'''
    def __init__(self, step, pipe, mode='run'):
        self.run_mode = mode
        env.config['run_mode'] = mode
        if hasattr(env, 'accessed_vars'):
            delattr(env, 'accessed_vars')
        super(Step_Executor, self).__init__(step)
        self.pipe = pipe
        # because step is executed in a separate SoS_Worker process, this
        # __pipe__ is available to all the actions that will be executed
        # in the step
        env.__pipe__ = pipe

    def submit_tasks(self, tasks):
        env.logger.debug(f'Send {tasks}')
        if 'queue' in env.sos_dict['_runtime'] and env.sos_dict['_runtime']['queue']:
            host = env.sos_dict['_runtime']['queue']
        else:
            # otherwise, use workflow default
            host = '__default__'
        self.pipe.send(f'tasks {host} {" ".join(tasks)}')

    def wait_for_tasks(self, tasks):
        if not tasks:
            return {}
        # wait till the executor responde
        results = {}
        while True:
            res = self.pipe.recv()
            if res is None:
                sys.exit(0)
            results.update(res)
            # all results have been obtained.
            if len(results) == len(tasks):
                break
        return results

    def run(self):
        try:
            res = Base_Step_Executor.run(self)
            if self.pipe is not None:
                env.logger.debug(f'Step {self.step.step_name()} sends result {short_repr(res)}')
                self.pipe.send(res)
            else:
                return res
        except Exception as e:
            if env.verbosity > 2:
                sys.stderr.write(get_traceback())
            if self.pipe is not None:
                env.logger.debug(f'Step {self.step.step_name()} sends exception {e}')
                self.pipe.send(e)
            else:
                raise e





