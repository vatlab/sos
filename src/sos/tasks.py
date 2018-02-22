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
import pickle
import time
import copy
import threading
import lzma
import random
from io import StringIO
from tokenize import generate_tokens
from collections.abc import Sequence, Mapping
import concurrent.futures
import traceback

from .utils import env, short_repr, sample_of_file, tail_of_file, linecount_of_file, \
    format_HHMMSS, expand_time, expand_size, StopInputGroup
from .eval import SoS_exec, SoS_eval, stmtHash, cfg_interpolate

from .targets import textMD5, RuntimeInfo, Undetermined, file_target, UnknownTarget, remote, sos_step, sos_targets
from .eval import interpolate
from .monitor import ProcessMonitor

from collections import OrderedDict
import subprocess


monitor_interval = 5
resource_monitor_interval = 60

class TaskParams(object):
    '''A parameter object that encaptulates parameters sending to
    task executors. This would makes the output of workers, especially
    in the web interface much cleaner (issue #259)'''
    def __init__(self, name, global_def, task, sos_dict, tags=[]):
        self.name = name
        self.global_def = global_def
        self.task = task
        self.sos_dict = sos_dict
        self.tags = sorted(list(set(tags)))

    def save(self, job_file):
        with open(job_file, 'wb') as jf:
            jf.write(f'SOSTASK1.2\n{" ".join(self.tags)}\n'.encode())
            # remove __builtins__ from sos_dict #835
            if 'CONFIG' in self.sos_dict and '__builtins__' in self.sos_dict['CONFIG']:
                self.sos_dict['CONFIG'].pop('__builtins__')
            try:
                jf.write(lzma.compress(pickle.dumps(self)))
            except Exception as e:
                env.logger.warning(e)
                raise

    def __repr__(self):
        return self.name

class MasterTaskParams(TaskParams):
    def __init__(self, num_workers=0):
        self.ID = 'M_0'
        self.name = self.ID
        self.global_def = ''
        self.task = ''
        self.sos_dict = {'_runtime': {}, '_input': sos_targets(), '_output': sos_targets(), '_depends': sos_targets(),
                'step_input': sos_targets(), 'step_output':sos_targets(),
                'step_depends': sos_targets(), 'step_name': '',
                '_index': 0}
        self.num_workers = num_workers
        self.tags = []
        # a collection of tasks that will be executed by the master task
        self.task_stack = []

    def num_tasks(self):
        return len(self.task_stack)

    def push(self, task_id, params):
        # update walltime, cores, and mem
        # right now we require all tasks to have same resource requirment, which is
        # quite natural because they are from the same step
        #
        # update input, output, and depends
        #
        # walltime
        if not self.task_stack:
            for key in ('walltime', 'max_walltime', 'cores', 'max_cores', 'mem', 'max_mem', 'map_vars',
                        'name', 'cur_dir', 'home_dir', 'verbosity', 'sig_mode', 'run_mode'):
                if key in params.sos_dict['_runtime'] and params.sos_dict['_runtime'][key] is not None:
                    self.sos_dict['_runtime'][key] = params.sos_dict['_runtime'][key]
            self.sos_dict['step_name'] = params.sos_dict['step_name']
            self.tags = params.tags
        else:
            for key in ('walltime', 'max_walltime', 'cores', 'max_cores', 'mem', 'max_mem',
                        'name', 'cur_dir', 'home_dir'):
                val0 = self.task_stack[0][1].sos_dict['_runtime'].get(key, None)
                val = params.sos_dict['_runtime'].get(key, None)
                if val0 != val:
                    raise ValueError(f'All tasks should have the same resource {key}')
                #
                nrow = len(self.task_stack) if self.num_workers <= 1 else ((len(self.task_stack) + 1) // self.num_workers + (0 if (len(self.task_stack) + 1) % self.num_workers == 0 else 1))
                if self.num_workers == 0:
                    ncol = 1
                elif nrow > 1:
                    ncol = self.num_workers
                else:
                    ncol = len(self.task_stack) + 1

                if val0 is None:
                    continue
                elif key == 'walltime':
                    # if define walltime
                    self.sos_dict['_runtime']['walltime'] = format_HHMMSS(nrow * expand_time(val0))
                elif key == 'mem':
                    # number of columns * mem for each + 100M for master
                    self.sos_dict['_runtime']['mem'] = ncol * expand_size(val0) + (expand_size('100M') if self.num_workers > 0 else 0)
                elif key == 'cores':
                    # number of columns * cores for each + 1 for the master
                    self.sos_dict['_runtime']['cores'] = ncol * val0 + (1 if self.num_workers > 0 else 0)
                elif key == 'name':
                    self.sos_dict['_runtime']['name'] = f'{val0}_{len(self.task_stack) + 1}'

            self.tags.extend(params.tags)
        #
        # input, output, preserved vars etc
        for key in ['_input', '_output', '_depends']:
            if key in params.sos_dict and isinstance(params.sos_dict[key], list):
                # do not extend duplicated input etc
                self.sos_dict[key].extend(list(set(params.sos_dict[key]) - set(self.sos_dict[key])))
        #
        self.task_stack.append([task_id, params])
        self.tags = sorted(list(set(self.tags)))
        #
        self.ID = f'M{len(self.task_stack)}_{self.task_stack[0][0]}'
        self.name = self.ID

def loadTask(filename):
    try:
        with open(filename, 'rb') as task:
            try:
                header = task.readline().decode()
                if header.startswith('SOSTASK1.1'):
                    # ignore the tags
                    task.readline()
                    return pickle.load(task)
                elif header.startswith('SOSTASK1.2'):
                    task.readline()
                    try:
                        return pickle.loads(lzma.decompress(task.read()))
                    except:
                        # at some point, the task files were compressed with zlib
                        import zlib
                        return pickle.loads(zlib.decompress(task.read()))
                else:
                    raise ValueError('Try old format')
            except:
                # old format
                task.seek(0)
                param = pickle.load(task)
                # old format does not have tags
                param.tags = []
                return param
    except ImportError as e:
        raise RuntimeError(
            f'Failed to load task {os.path.basename(filename)}, which is likely caused by incompatible python modules between local and remote hosts: {e}')


def addTags(filename, new_tags):
    with open(filename, 'rb') as task:
        header = task.readline()
        # read the tags
        tags = task.readline().decode().strip().split(' ')
        if isinstance(new_tags, str):
            if new_tags in tags:
                return
            else:
                tags.append(new_tags)
        elif isinstance(new_tags, Sequence):
            new_tags = [tag for tag in new_tags if tag not in tags]
            if new_tags:
                tags.extend(new_tags)
            else:
                return
        else:
            raise ValueError(f'Cannot add tags {new_tags} to task {filename}')
        body = task.read()
    with open(filename, 'wb') as task:
        task.write(header)
        task.write((' '.join(tags) + '\n').encode())
        task.write(body)


def taskDuration(task):
    filename = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', f'{task}.task')
    return os.path.getatime(filename) - os.path.getmtime(filename)


def taskTags(task):
    filename = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', f'{task}.task')
    atime = os.path.getatime(filename)
    try:
        with open(filename, 'rb') as task:
            try:
                header = task.readline().decode()
                if header.startswith('SOSTASK'):
                    return task.readline().decode().strip()
                else:
                    return ''
            except:
                return ''
    except Exception as e:
        env.logger.warning(f'Failed to get tags for task {task}: {e}')
        return []
    finally:
        os.utime(filename, (atime, os.path.getmtime(filename)))


def collect_task_result(task_id, sos_dict):
    shared = {}
    if 'shared' in env.sos_dict['_runtime']:
        svars = env.sos_dict['_runtime']['shared']
        if isinstance(svars, str):
            if vars not in env.sos_dict:
                raise ValueError(f'Unavailable shared variable {svars} after the completion of task {task_id}')
            shared[svars] = copy.deepcopy(env.sos_dict[svars])
        elif isinstance(svars, Mapping):
            for var, val in svars.items():
                if var != val:
                    env.sos_dict.set(var, SoS_eval(val))
                if var not in env.sos_dict:
                    raise ValueError(f'Unavailable shared variable {var} after the completion of task {task_id}')
                shared[var] = copy.deepcopy(env.sos_dict[var])
        elif isinstance(svars, Sequence):
            # if there are dictionaries in the sequence, e.g.
            # shared=['A', 'B', {'C':'D"}]
            for item in svars:
                if isinstance(item, str):
                    if item not in env.sos_dict:
                        raise ValueError(f'Unavailable shared variable {item} after the completion of task {task_id}')
                    shared[item] = copy.deepcopy(env.sos_dict[item])
                elif isinstance(item, Mapping):
                    for var, val in item.items():
                        if var != val:
                            env.sos_dict.set(var, SoS_eval(val))
                        if var not in env.sos_dict:
                            raise ValueError(
                                f'Unavailable shared variable {var} after the completion of task {task_id}')
                        shared[var] = copy.deepcopy(env.sos_dict[var])
                else:
                    raise ValueError(
                        f'Option shared should be a string, a mapping of expression, or a list of string or mappings. {svars} provided')
        else:
            raise ValueError(
                f'Option shared should be a string, a mapping of expression, or a list of string or mappings. {svars} provided')
        env.logger.debug(f'task {task_id} (index={env.sos_dict["_index"]}) return shared variable {shared}')
    # the difference between sos_dict and env.sos_dict is that sos_dict (the original version) can have remote() targets
    # which should not be reported.
    if env.sos_dict['_output'] is None:
        output = {}
    elif isinstance(env.sos_dict['_output'], Undetermined):
        from .workflow_executor import __null_func__
        from .targets import dynamic
        from .step_executor import _expand_file_list
        env.sos_dict.set('__null_func__', __null_func__)
        # re-process the output statement to determine output files
        args, _ = SoS_eval(f'__null_func__({env.sos_dict["_output"].expr})')
        # handle dynamic args
        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        output = {x:file_target(x).target_signature() for x in _expand_file_list(True, *args)}
    elif sos_dict['_output'] is None:
        output = {}
    else:
        output = {x:file_target(x).target_signature() for x in sos_dict['_output'] if isinstance(x, (str, file_target))}

    input = {} if env.sos_dict['_input'] is None or sos_dict['_input'] is None else {x:file_target(x).target_signature() for x in sos_dict['_input'] if isinstance(x, (str, file_target))}
    depends = {} if env.sos_dict['_depends'] is None or sos_dict['_depends'] is None else {x:file_target(x).target_signature() for x in sos_dict['_depends'] if isinstance(x, (str, file_target))}
    return {'ret_code': 0, 'task': task_id, 'input': input, 'output': output, 'depends': depends,
            'shared': {env.sos_dict['_index']: shared} }

def execute_task(task_id, verbosity=None, runmode='run', sigmode=None, monitor_interval=5,
    resource_monitor_interval=60):
    res = _execute_task(task_id, verbosity, runmode, sigmode, monitor_interval, resource_monitor_interval)
    # write result file
    res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.res')
    with open(res_file, 'wb') as res_file:
        pickle.dump(res, res_file)
    if res['ret_code'] != 0 and 'exception' in res:
        with open(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err'), 'a') as err:
            err.write(f'Task {task_id} exits with code {res["ret_code"]}')
    return res['ret_code']

def _execute_task(task_id, verbosity=None, runmode='run', sigmode=None, monitor_interval=5,
    resource_monitor_interval=60):
    '''A function that execute specified task within a local dictionary
    (from SoS env.sos_dict). This function should be self-contained in that
    it can be handled by a task manager, be executed locally in a separate
    process or remotely on a different machine.'''
    # start a monitoring file, which would be killed after the job
    # is done (killed etc)
    if isinstance(task_id, str):
        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
        params = loadTask(task_file)
        subtask = False
    else:
        # subtask
        subtask = True
        task_id, params = task_id
        env.logger.trace(f'Executing subtask {task_id}')

    if hasattr(params, 'task_stack'):
        # pulse thread
        m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
            resource_monitor_interval=resource_monitor_interval,
            max_walltime=params.sos_dict['_runtime'].get('max_walltime', None),
            max_mem=params.sos_dict['_runtime'].get('max_mem', None),
            max_procs=params.sos_dict['_runtime'].get('max_procs', None))
        m.start()

        master_out = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.out')
        master_err = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err')
        # if this is a master task, calling each sub task
        with open(master_out, 'wb') as out, open(master_err, 'wb') as err:
            def copy_out_and_err(result):
                tid = result['task']
                out.write(f'{tid}: {"completed" if result["ret_code"] == 0 else "failed"}\n'.encode())
                if 'output' in result:
                    out.write(f'output: {result["output"]}\n'.encode())
                sub_out = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', tid + '.out')
                if os.path.isfile(sub_out):
                    with open(sub_out, 'rb') as sout:
                        out.write(sout.read())

                sub_err = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', tid + '.err')
                err.write(f'{tid}: {"completed" if result["ret_code"] == 0 else "failed"}\n'.encode())
                if os.path.isfile(sub_err):
                    with open(sub_err, 'rb') as serr:
                        err.write(serr.read())

            if params.num_workers > 1:
                from multiprocessing.pool import Pool
                p = Pool(params.num_workers)
                results = []
                for t in params.task_stack:
                    results.append(p.apply_async(_execute_task, (t, verbosity, runmode,
                        sigmode, monitor_interval, resource_monitor_interval), callback=copy_out_and_err))
                for idx,r in enumerate(results):
                    results[idx] = r.get()
                p.close()
                p.join()
                # we wait for all results to be ready to return or raise
                # but we only raise exception for one of the subtasks
                for res in results:
                    if 'exception' in res:
                        failed = [x.get("task", "") for x in results if "exception" in x]
                        env.logger.error(f'{task_id} ``failed`` due to failure of subtask{"s" if len(failed) > 1 else ""} {", ".join(failed)}')
                        return {'ret_code': 1, 'exception': res['exception'], 'task': task_id }
            else:
                results = []
                for tid, tdef in params.task_stack:
                    res = _execute_task((tid, tdef), verbosity=verbosity, runmode=runmode,
                        sigmode=sigmode, monitor_interval=monitor_interval,
                        resource_monitor_interval=resource_monitor_interval)
                    copy_out_and_err(res)
                    results.append(res)
                for res in results:
                    if 'exception' in res:
                        failed = [x.get("task", "") for x in results if "exception" in x]
                        env.logger.error(f'{task_id} ``failed`` due to failure of subtask{"s" if len(failed) > 1 else ""} {", ".join(failed)}')
                        return {'ret_code': 1, 'exception': res['exception'], 'task': task_id}
        #
        # now we collect result
        all_res = {'ret_code': 0, 'output': {}, 'subtasks': {}, 'shared': {}}
        for tid, x in zip(params.task_stack, results):
            all_res['ret_code'] += x['ret_code']
            all_res['output'].update(x['output'])
            all_res['subtasks'][tid[0]] = x
            all_res['shared'].update(x['shared'])
        return all_res

    global_def, task, sos_dict = params.global_def, params.task, params.sos_dict

    # task output
    env.sos_dict.set('__std_out__', os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.out'))
    env.sos_dict.set('__std_err__', os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err'))
    env.logfile = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err')
    # clear the content of existing .out and .err file if exists, but do not create one if it does not exist
    if os.path.exists(env.sos_dict['__std_out__']):
        open(env.sos_dict['__std_out__'], 'w').close()
    if os.path.exists(env.sos_dict['__std_err__']):
        open(env.sos_dict['__std_err__'], 'w').close()

    if verbosity is not None:
        env.verbosity = verbosity
    try:
        # global def could fail due to execution on remote host...
        # we also execute global_def way before others and allows variables set by
        # global_def be overwritten by other passed variables
        #
        # note that we do not handle parameter in tasks because values should already be
        # in sos_task dictionary
        SoS_exec('''\
import os, sys, glob
from sos.runtime import *
CONFIG = {}
del sos_handle_parameter_
''' + global_def)
    except Exception as e:
        env.logger.trace(f'Failed to execute global definition {short_repr(global_def)}: {e}')

    if '_runtime' not in sos_dict:
        sos_dict['_runtime'] = {}

    # pulse thread
    m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
        resource_monitor_interval=resource_monitor_interval,
        max_walltime=sos_dict['_runtime'].get('max_walltime', None),
        max_mem=sos_dict['_runtime'].get('max_mem', None),
        max_procs=sos_dict['_runtime'].get('max_procs', None))

    m.start()
    if sigmode is not None:
        env.config['sig_mode'] = sigmode
    env.config['run_mode'] = runmode
    #
    if subtask:
        env.logger.debug(f'{task_id} ``started``')
    else:
        env.logger.info(f'{task_id} ``started``')

    env.sos_dict.quick_update(sos_dict)

    # if targets are defined as `remote`, they should be resolved during task execution
    def resolve_remote(x):
        if isinstance(x, remote):
            x = x.resolve()
            if isinstance(x, str):
                x = interpolate(x, env.sos_dict._dict)
        return x

    for key in ['step_input', '_input',  'step_output', '_output', 'step_depends', '_depends']:
        if key in sos_dict and isinstance(sos_dict[key], (list, sos_targets)):
            # resolve remote() target
            env.sos_dict.set(key, sos_targets(resolve_remote(x) for x in sos_dict[key] if not isinstance(x, sos_step)))

    skipped = False
    if env.config['sig_mode'] == 'ignore':
        sig = None
    else:
        tokens = [x[1] for x in generate_tokens(StringIO(task).readline)]
        # try to add #task so that the signature can be different from the step
        # if everything else is the same
        sig = RuntimeInfo(textMD5('#task\n' + ' '.join(tokens)), task,
            env.sos_dict['_input'].targets(), env.sos_dict['_output'].targets(),
            env.sos_dict['_depends'].targets(), env.sos_dict['__signature_vars__'])
        sig.lock()

        idx = env.sos_dict['_index']
        if env.config['sig_mode'] == 'default':
            matched = sig.validate()
            if isinstance(matched, dict):
                # in this case, an Undetermined output can get real output files
                # from a signature
                env.sos_dict.set('_input', sos_targets(matched['input']))
                env.sos_dict.set('_depends', sos_targets(matched['depends']))
                env.sos_dict.set('_output', sos_targets(matched['output']))
                env.sos_dict.update(matched['vars'])
                env.logger.info(
                    f'Task ``{env.sos_dict["step_name"]}`` (index={idx}) is ``ignored`` due to saved signature')
                skipped = True
        elif env.config['sig_mode'] == 'assert':
            matched = sig.validate()
            if isinstance(matched, str):
                raise RuntimeError(f'Signature mismatch: {matched}')
            else:
                env.sos_dict.set('_input', sos_targets(matched['input']))
                env.sos_dict.set('_depends', sos_targets(matched['depends']))
                env.sos_dict.set('_output', sos_targets(matched['output']))
                env.sos_dict.update(matched['vars'])
                env.logger.info(
                    f'Step ``{env.sos_dict["step_name"]}`` (index={idx}) is ``ignored`` with matching signature')
                skipped = True
        elif env.config['sig_mode'] == 'build':
            # build signature require existence of files
            if sig.write(rebuild=True):
                env.logger.info(
                    f'Task ``{env.sos_dict["step_name"]}`` (index={idx}) is ``ignored`` with signature constructed')
                skipped = True
            else:
                env.logger.info(
                    f'Task ``{env.sos_dict["step_name"]}`` (index={idx}) is ``executed`` with failed signature constructed')
        elif env.config['sig_mode'] == 'force':
            skipped = False
        else:
            raise RuntimeError(f'Unrecognized signature mode {env.config["sig_mode"]}')

    if skipped:
        env.logger.info(f'{task_id} ``skipped``')
        return collect_task_result(task_id, sos_dict)

    # if we are to really execute the task, touch the task file so that sos status shows correct
    # execution duration.
    if not subtask:
        os.utime(task_file, None)

    try:
        # go to 'cur_dir'
        if '_runtime' in sos_dict and 'cur_dir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['cur_dir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
                    os.chdir(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
                except Exception as e:
                    # sometimes it is not possible to go to a "cur_dir" because of
                    # file system differences, but this should be ok if a work_dir
                    # has been specified.
                    env.logger.debug(f'Failed to create cur_dir {sos_dict["_runtime"]["cur_dir"]}')
            else:
                os.chdir(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
        #
        orig_dir = os.getcwd()

        if runmode != 'dryrun':
            # we will need to check existence of targets because the task might
            # be executed on a remote host where the targets are not available.
            for target in (sos_dict['_input'] if isinstance(sos_dict['_input'], list) else []) + \
                (sos_dict['_depends'] if isinstance(sos_dict['_depends'], list) else []):
                # if the file does not exist (although the signature exists)
                # request generation of files
                if isinstance(target, str):
                    if not file_target(target).target_exists('target'):
                        # remove the signature and regenerate the file
                        file_target(target).remove_sig()
                        raise UnknownTarget(target)
                # the sos_step target should not be checked in tasks because tasks are
                # independently executable units.
                elif not isinstance(target, sos_step) and not target.target_exists('target'):
                    target.remove_sig()
                    raise UnknownTarget(target)

        # create directory. This usually has been done at the step level but the task can be executed
        # on a remote host where the directory does not yet exist.
        ofiles = env.sos_dict['_output']
        if not isinstance(ofiles, (type(None), Undetermined)):
            for ofile in ofiles:
                if isinstance(ofile, str):
                    parent_dir = os.path.split(os.path.expanduser(ofile))[0]
                    if parent_dir and not os.path.isdir(parent_dir):
                        try:
                            os.makedirs(parent_dir)
                        except Exception as e:
                            # this can fail but we do not really care because the task itself might
                            # create this directory, or if the directory has already been created by other tasks
                            env.logger.warning(f'Failed to create directory {parent_dir}: {e}')


        # go to user specified workdir
        if '_runtime' in sos_dict and 'workdir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['workdir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['workdir']))
                except Exception as e:
                    raise RuntimeError(f'Failed to create workdir {sos_dict["_runtime"]["workdir"]}')
            os.chdir(os.path.expanduser(sos_dict['_runtime']['workdir']))
        # set environ ...
        # we join PATH because the task might be executed on a different machine
        if '_runtime' in sos_dict:
            if 'env' in sos_dict['_runtime']:
                for key, value in sos_dict['_runtime']['env'].items():
                    if 'PATH' in key and key in os.environ:
                        new_path = OrderedDict()
                        for p in value.split(os.pathsep):
                            new_path[p] = 1
                        for p in value.split(os.environ[key]):
                            new_path[p] = 1
                        os.environ[key] = os.pathsep.join(new_path.keys())
                    else:
                        os.environ[key] = value
            if 'prepend_path' in sos_dict['_runtime']:
                if isinstance(sos_dict['_runtime']['prepend_path'], str):
                    os.environ['PATH'] = sos_dict['_runtime']['prepend_path'] + os.pathsep + os.environ['PATH']
                elif isinstance(env.sos_dict['_runtime']['prepend_path'], Sequence):
                    os.environ['PATH'] = os.pathsep.join(sos_dict['_runtime']['prepend_path']) + os.pathsep + os.environ['PATH']
                else:
                    raise ValueError(
                        f'Unacceptable input for option prepend_path: {sos_dict["_runtime"]["prepend_path"]}')


        # step process
        SoS_exec(task)

        if subtask:
            env.logger.debug(f'{task_id} ``completed``')
        else:
            env.logger.info(f'{task_id} ``completed``')

    except StopInputGroup as e:
        # task ignored with stop_if exception
        if e.message:
            env.logger.warning(f'{task_id} ``stopped``: {e.message}')
        return {'ret_code': 0, 'task': task_id, 'input': [],
            'output': [], 'depends': [], 'shared': {}}
    except KeyboardInterrupt:
        env.logger.error(f'{task_id} ``interrupted``')
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
            env.logger.debug(f'''
---------------------------------------------------------------------------
{error_class:42}Traceback (most recent call last)
{msg}
{error_class}: {detail}''')
            env.logger.debug(f'{error_class}: {detail}')

        env.logger.error(f'{task_id} ``failed``: {error_class} {detail}')
        return {'ret_code': 1, 'exception': e, 'task': task_id, 'shared': {}}
    finally:
        env.sos_dict.set('__step_sig__', None)
        os.chdir(orig_dir)
        if not subtask:
            # after the task is completed, we change the access time
            # but keep the modify time of the task file, which serves
            # as the "starting" time of the task.
            os.utime(task_file, (time.time(), os.path.getmtime(task_file)))

    if sig:
        sig.write()
        sig.release()

    # the final result should be relative to cur_dir, not workdir
    # because output is defined outside of task
    return collect_task_result(task_id, sos_dict)

def check_task(task):
    #
    # status of the job, please refer to https://github.com/vatlab/SOS/issues/529
    # for details.
    #
    task_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.task')
    if not os.path.isfile(task_file):
        return 'missing'
    pulse_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.pulse')
    if not os.path.isfile(pulse_file):
        pulse_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.status')

    def has_pulse():
        return os.path.isfile(pulse_file) and os.stat(pulse_file).st_mtime >= os.stat(task_file).st_mtime

    res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.res')

    def has_res():
        return os.path.isfile(res_file) and os.stat(res_file).st_mtime >= os.stat(task_file).st_mtime

    job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.sh')

    def has_job():
        return os.path.isfile(job_file) and os.stat(job_file).st_mtime >= os.stat(task_file).st_mtime

    if has_res():
        try:
            from .targets import file_target
            with open(res_file, 'rb') as result:
                res = pickle.load(result)
            if ('ret_code' in res and res['ret_code'] == 0) or ('succ' in res and res['succ'] == 0):
                for var in ('input', 'output', 'depends'):
                    if var not in res or not isinstance(res[var], dict):
                        continue
                    for x,y in res[var].items():
                        if not file_target(x).target_exists() or file_target(x).target_signature() != y:
                            env.logger.debug(f'{x} not found or signature mismatch')
                            return 'signature-mismatch'
                return 'completed'
            else:
                return 'failed'
        except Exception as e:
            # sometimes the resfile is changed while we are reading it
            # so we wait a bit and try again.
            env.logger.warning(e)
            time.sleep(.5)
            return check_task(task)
    #
    if has_pulse():
        # dead?
        # if the status file is readonly
        if not os.access(pulse_file, os.W_OK):
            return 'aborted'
        start_stamp = os.stat(pulse_file).st_mtime
        elapsed = time.time() - start_stamp
        if elapsed < 0:
            env.logger.warning(f'{pulse_file} is created in the future. Your system time might be problematic')
        # if the file is within 5 seconds
        if elapsed < monitor_interval:
            return 'running'
        elif elapsed > 2 * monitor_interval:
            if has_res():
                # result file appears during sos tatus run
                return check_task(task)
            else:
                return 'aborted'
        # otherwise, let us be patient ... perhaps there is some problem with the filesystem etc
        time.sleep(2 * monitor_interval)
        end_stamp = os.stat(pulse_file).st_mtime
        # the process is still alive
        if has_res():
            return check_task(task)
        elif start_stamp != end_stamp:
            return 'running'
        else:
            return 'aborted'
    # if there is no status file
    if has_job():
        return 'submitted'
    else:
        return 'pending'

def check_tasks(tasks, verbosity=1, html=False, start_time=False, age=None, tags=None, status=None):
    # verbose is ignored for now
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    if not tasks:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getmtime(x)) for x in tasks]
        if not all_tasks:
            return
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', f'{t}*.task'))
            matched = [(os.path.basename(x)[:-5], os.path.getmtime(x)) for x in matched]
            if not matched:
                all_tasks.append((t, None))
            else:
                all_tasks.extend(matched)
    if age is not None:
        age = expand_time(age, default_unit='d')
        if age > 0:
            all_tasks = [x for x in all_tasks if time.time() - x[1] >= age]
        else:
            all_tasks = [x for x in all_tasks if time.time() - x[1] <= -age]

    all_tasks = sorted(list(set(all_tasks)), key=lambda x: 0 if x[1] is None else x[1])

    if tags:
        all_tasks = [x for x in all_tasks if any(x in tags for x in taskTags(x[0]).split(' '))]

    if not all_tasks:
        env.logger.info('No matching tasks')
        return
    # at most 20 threads
    p = Pool(min(20, len(all_tasks)))
    obtained_status = p.map(check_task, [x[0] for x in all_tasks])
    if status:
        all_tasks = [x for x, s in zip(all_tasks, obtained_status) if s in status]
        obtained_status = [x for x in obtained_status if x in status]
    #
    # automatically remove non-running tasks that are more than 30 days old
    to_be_removed = []
    #
    if html:
        verbosity = -1
    if verbosity == 0:
        print('\n'.join(obtained_status))
    elif verbosity == 1:
        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30*24*60*60 and s != 'running':
                to_be_removed.append(t)
                continue
            print(f'{t}\t{s}')
    elif verbosity == 2:
        from .utils import PrettyRelativeTime
        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30*24*60*60 and s != 'running':
                to_be_removed.append(t)
                continue
            if start_time:
                if d is None:
                    print(f'{t}\t{taskTags(t)}\t{time.time()}\t{s}')
                else:
                    print(f'{t}\t{taskTags(t)}\t{d}\t{s}')
            else:
                if d is None:
                    print(f'{t}\t{taskTags(t)}\t{"":>15}\t{s}')
                elif s in ('pending', 'submitted', 'running'):
                    print(f'{t}\t{taskTags(t)}\t{PrettyRelativeTime(time.time() - d):>15}\t{s}')
                else:
                    # completed or failed
                    print(f'{t}\t{taskTags(t)}\t{PrettyRelativeTime(taskDuration(t)):>15}\t{s}')
    elif verbosity > 2:
        from .utils import PrettyRelativeTime
        import pprint
        from .monitor import summarizeExecution

        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30*24*60*60 and s != 'running':
                to_be_removed.append(t)
                continue
            print(f'{t}\t{s}\n')
            task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                continue
            if d is not None:
                print(f'Started {PrettyRelativeTime(time.time() - d)} ago')
            if s not in ('pending', 'submitted', 'running'):
                print(f'Duration {PrettyRelativeTime(taskDuration(t))}')
            params = loadTask(task_file)
            print('TASK:\n=====')
            print(params.task)
            print('TAGS:\n=====')
            print(' '.join(params.tags))
            print()
            if params.global_def:
                print('GLOBAL:\n=======')
                print(params.global_def)
                print()
            print('ENVIRONMENT:\n============')
            job_vars = params.sos_dict
            for k in sorted(job_vars.keys()):
                v = job_vars[k]
                print(f'{k:22}{short_repr(v) if verbosity == 3 else pprint.pformat(v)}')
            print()
            print('EXECUTION STATS:\n================')
            print(summarizeExecution(t, status=s))
            if verbosity == 4:
                # if there are other files such as job file, print them.
                files = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.*'))
                for f in sorted([x for x in files if os.path.splitext(x)[-1] not in ('.res',
                    '.task', '.pulse', '.status', '.def')]):
                    print(f'{os.path.basename(f)}:\n{"="*(len(os.path.basename(f))+1)}')
                    try:
                        with open(f) as fc:
                            print(fc.read())
                    except Exception:
                        print('Binary file')
    else:
        # HTML output
        from .utils import PrettyRelativeTime, isPrimitive
        from .monitor import summarizeExecution
        import pprint
        print('<table width="100%" class="resource_table">')
        def row(th=None, td=None):
            if td is None:
                print(f'<tr><th align="right" width="30%">{th}</th><td></td></tr>')
            elif th is None:
                print(f'<tr><td colspan="2" align="left"  width="30%">{td}</td></tr>')
            else:
                print(
                    f'<tr><th align="right"  width="30%">{th}</th><td align="left"><div class="one_liner">{td}</div></td></tr>')
        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30*24*60*60:
                to_be_removed.append(t)
                continue
            row('ID', t)
            row('Status', s)
            task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                print('</table>')
                continue
            if d is not None:
                row('Start', f'{PrettyRelativeTime(time.time() - d):>15} ago')
            if s not in ('pending', 'submitted', 'running'):
                row('Duration', f'{PrettyRelativeTime(taskDuration(t)):>15}')
            task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                continue
            params = loadTask(task_file)
            row('Task')
            row(td=f'<pre style="text-align:left">{params.task}</pre>')
            row('Tags')
            row(td=f'<pre style="text-align:left">{" ".join(params.tags)}</pre>')
            if params.global_def:
                row('Global')
                row(td=f'<pre style="text-align:left">{params.global_def}</pre>')
            #row('Environment')
            job_vars = params.sos_dict
            for k in sorted(job_vars.keys()):
                v = job_vars[k]
                if not k.startswith('__') and not k == 'CONFIG':
                    if k == '_runtime':
                        for _k, _v in v.items():
                            if isPrimitive(_v) and _v not in (None, '', [], (), {}):
                                row(_k, _v)
                    elif isPrimitive(v) and v not in (None, '', [], (), {}):
                        row(k, f'<pre style="text-align:left">{pprint.pformat(v)}</pre>')
            summary = summarizeExecution(t, status=s)
            if summary:
                #row('Execution')
                for line in summary.split('\n'):
                    fields = line.split(None, 1)
                    if fields[0] == 'task':
                        continue
                    row(fields[0], '' if fields[1] is None else fields[1])
                # this is a placeholder for the frontend to draw figure
                row(td=f'<div id="res_{t}"></div>')
            #
            files = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.*'))
            for f in sorted([x for x in files if os.path.splitext(x)[-1] not in ('.def', '.res', '.task', '.pulse', '.status')]):
                numLines = linecount_of_file(f)
                row(os.path.splitext(f)[-1], '(empty)' if numLines == 0 else f'{numLines} lines{"" if numLines < 200 else " (showing last 200)"}')
                try:
                    row(td=f'<small><pre style="text-align:left">{tail_of_file(f, 200, ansi2html=True)}</pre></small>')
                except Exception:
                    row(td='<small><pre style="text-align:left">ignored.</pre><small>')
            print('</table>')
            #
            # supplement run time information
            pulse_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.pulse')
            if not os.path.isfile(pulse_file):
                return
            else:
                # A sample of 400 point should be enough to show the change of resources
                lines = sample_of_file(pulse_file, 400).splitlines()
                if len(lines) <= 2:
                    return
            # read the pulse file and plot it
            #time   proc_cpu        proc_mem        children        children_cpu    children_mem
            try:
                etime = []
                cpu = []
                mem = []
                for line in lines:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.split()
                    etime.append(float(fields[0]))
                    cpu.append(float(fields[1]) + float(fields[4]))
                    mem.append(float(fields[2]) / 1e6 + float(fields[5]) / 1e6)
                if not etime:
                    return
            except Exception:
                return
            #
            print('''
<script>
    function loadFiles(files, fn) {
        if (!files.length) {
            files = [];
        }
        var head = document.head || document.getElementsByTagName('head')[0];

        function loadFile(index) {
            if (files.length > index) {
                if (files[index].endsWith('.css')) {
                    var fileref = document.createElement('link');
                    fileref.setAttribute("rel", "stylesheet");
                    fileref.setAttribute("type", "text/css");
                    fileref.setAttribute("href", files[index]);
                } else {
                    var fileref = document.createElement('script');
                    fileref.setAttribute("type", "text/javascript");
                    fileref.setAttribute("src", files[index]);
                }
                console.log('Load ' + files[index]);
                head.appendChild(fileref);
                index = index + 1;
                // Used to call a callback function
                fileref.onload = function() {
                    loadFile(index);
                }
            } else if (fn) {
                fn();
            }
        }
        loadFile(0);
    }

function plotResourcePlot_''' + t + '''() {
    // get the item
    // parent element is a table cell, needs enlarge
    document.getElementById("res_''' + t + '''").parentElement.setAttribute("height", "300px;");
    $("#res_''' + t + '''").css("height", "300px");
    $("#res_''' + t + '''").css("width", "100%");
    $("#res_''' + t + '''").css("min-height", "300px");

    var cpu = [''' + ','.join([f'[{x*1000},{y}]' for x, y in zip(etime, cpu)]) + '''];
    var mem = [''' + ','.join([f'[{x*1000},{y}]' for x, y in zip(etime, mem)]) + '''];

    $.plot('#res_''' + t + '''', [{
                data: cpu,
                label: "CPU (%)"
            },
            {
                data: mem,
                label: "mem (M)",
                yaxis: 2
            }
        ], {
            xaxes: [{
                mode: "time"
            }],
            yaxes: [{
                min: 0
            }, {
                position: "right",
                tickFormatter: function(v, axis) {
                    return v.toFixed(1) + 'M';
                }
            }],
            legend: {
                position: "nw"
            }
        });
}

var dt = 100;
// the frontend might be notified before the table is inserted as results.
function showResourceFigure_''' + t + '''() {
    if ( $("#res_''' + t + '''").length === 0) {
          dt = dt * 1.5; // slow-down checks for datatable as time goes on;
          setTimeout(showResourceFigure_''' + t + ''', dt);
          return;
    } else {
        $("#res_''' + t + '''").css('width', "100%").css('height', "300px");
        loadFiles(["http://www.flotcharts.org/flot/jquery.flot.js",
             "http://www.flotcharts.org/flot/jquery.flot.time.js"
            ], plotResourcePlot_''' + t + ''');
    }
}
showResourceFigure_''' + t + '''()
</script>
''')
    # remove jobs that are older than 1 month
    if to_be_removed:
        purge_tasks(to_be_removed, verbosity = 0)

def kill_tasks(tasks, tags=None):
    #
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    if not tasks:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [os.path.basename(x)[:-5] for x in tasks]
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', f'{t}*.task'))
            matched = [os.path.basename(x)[:-5] for x in matched]
            if not matched:
                env.logger.warning(f'{t} does not match any existing task')
            else:
                all_tasks.extend(matched)
    if tags:
        all_tasks = [x for x in all_tasks if any(x in tags for x in taskTags(x).split(' '))]

    if not all_tasks:
        env.logger.warning('No task to kill')
        return
    all_tasks = sorted(list(set(all_tasks)))
    p = Pool(len(all_tasks))
    killed = p.map(kill_task, all_tasks)
    for s, t in zip(killed, all_tasks):
        print(f'{t}\t{s}')

def kill_task(task):
    status = check_task(task)
    if status == 'pending':
        return 'cancelled'
    # remove job file as well
    job_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.sh')
    if os.path.isfile(job_file):
        try:
            os.remove(job_file)
        except Exception:
            pass
    if status != 'running':
        return status
    # job is running
    pulse_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.pulse')
    from stat import S_IREAD, S_IRGRP, S_IROTH
    os.chmod(pulse_file, S_IREAD|S_IRGRP|S_IROTH)
    return 'killed'


def purge_tasks(tasks, purge_all=False, age=None, status=None, tags=None, verbosity=2):
    # verbose is ignored for now
    import glob
    if tasks:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', f'{t}*.task'))
            matched = [(os.path.basename(x)[:-5], os.path.getmtime(x)) for x in matched]
            all_tasks.extend(matched)
    else:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getmtime(x)) for x in tasks]
    #
    if age is not None:
        age = expand_time(age, default_unit='d')
        if age > 0:
            all_tasks = [x for x in all_tasks if time.time() - x[1] >= age]
        else:
            all_tasks = [x for x in all_tasks if time.time() - x[1] <= -age]

    if status:
        # at most 20 threads
        from multiprocessing.pool import ThreadPool as Pool
        p = Pool(min(20, len(all_tasks)))
        task_status = p.map(check_task, [x[0] for x in all_tasks])
        all_tasks = [x for x,s in zip(all_tasks, task_status) if s in status]

    if tags:
        all_tasks = [x for x in all_tasks if any(x in tags for x in taskTags(x[0]).split(' '))]
    #
    # remoe all task files
    all_tasks = set([x[0] for x in all_tasks])
    if all_tasks:
        #
        # find all related files, including those in nested directories
        from collections import defaultdict
        to_be_removed = defaultdict(list)
        for dirname, _, filelist in os.walk(os.path.join(os.path.expanduser('~'), '.sos', 'tasks')):
            for f in filelist:
                ID = os.path.basename(f).split('.', 1)[0]
                if ID in all_tasks:
                    to_be_removed[ID].append(os.path.join(dirname, f))
        #
        for task in all_tasks:
            removed = True
            for f in to_be_removed[task]:
                try:
                    if verbosity > 3:
                        env.logger.trace(f'Remove {f}')
                    os.remove(f)
                except Exception as e:
                    removed = False
                    if verbosity > 0:
                        env.logger.warning(f'Failed to purge task {task[0]}')
            if removed and verbosity > 1:
                env.logger.info(f'Task ``{task}`` removed.')
    elif verbosity > 1:
        env.logger.info('No matching tasks')
    if purge_all:
        matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*'))
        count = 0
        for f in matched:
            if os.path.isdir(f):
                import shutil
                try:
                    shutil.rmtree(f)
                    count += 1
                except Exception as e:
                    if verbosity > 0:
                        env.logger.warning(f'Failed to remove {f}')
            else:
                try:
                    os.remove(f)
                    count += 1
                except Exception as e:
                    if verbosity > 0:
                        env.logger.warning(f'Failed to remove {e}')
        if count > 0 and verbosity > 1:
            env.logger.info(f'{count} other files and directories are removed.')
    return ''

class TaskEngine(threading.Thread):
    def __init__(self, agent):
        threading.Thread.__init__(self)
        self.daemon = True
        #
        # agent is the agent that provides function
        #
        #    run_command
        #
        # to submit command, which can be a direct process call, or a call
        # on the remote server.
        #
        self.agent = agent
        self.config = agent.config
        self.alias = self.config['alias']

        self.engine_ready = threading.Event()

        self.running_tasks = []
        self.pending_tasks = []
        self.submitting_tasks = {}
        self.canceled_tasks = []
        self.resuming_tasks = set()

        self.task_status = OrderedDict()
        self.task_date = {}
        self.last_checked = None
        if 'status_check_interval' not in self.config:
            self.status_check_interval = 10
        else:
            self.status_check_interval = self.config['status_check_interval']
        #
        if env.config['max_running_jobs'] is not None:
            # override from command line
            self.max_running_jobs = env.config['max_running_jobs']
        elif 'max_running_jobs' in self.config:
            # queue setting
            self.max_running_jobs = self.config['max_running_jobs']
        else:
            # default
            self.max_running_jobs = max(os.cpu_count() // 2, 1)
        #
        # multiple thread job submission does not work because the threads share the
        # same namespace with variables such as sos_dict. Variables changed by
        # one thread can be changed again by another thread, which makes it unsafe
        # to submit jobs this way. Multi-processing is possible but that can be done
        # later.
        #
        self._thread_workers = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        self._status_checker = None
        #
        if env.config['wait_for_task'] is not None:
            self.wait_for_task = env.config['wait_for_task']
        elif 'wait_for_task' in self.config:
            self.wait_for_task = self.config['wait_for_task']
        else:
            # default
            self.wait_for_task = True
        #
        # how to stack jobs ... currently backgroun process queue
        # allows stacking of up to 1000 tasks, but PBS queue does not
        # allow stacking.
        self.batch_size = 1

    def notify(self, msg):
        # GUI ...
        if hasattr(env, '__task_notifier__'):
            if not isinstance(msg, str):
                env.__task_notifier__(msg)
        elif isinstance(msg, str):
            env.logger.info(msg)
        # text mode does not provide detailed message change information

    def monitor_tasks(self, tasks=None, status=None, age=None):
        '''Start monitoring specified or all tasks'''
        self.engine_ready.wait()

        if not tasks:
            tasks = self.task_status.keys()
        else:
            tasks = [x for x in tasks if x in self.task_status]

        # we only monitor running tasks
        with threading.Lock():
            for task in tasks:
                if self.task_status[task] in ('submitted', 'running') and not task in self.running_tasks:
                    # these tasks will be actively monitored
                    self.running_tasks.append(task)
        #
        if age is not None:
            age = expand_time(age, default_unit='d')
        return sorted([(x, self.task_status[x], self.task_date.get(x, time.time())) for x in tasks \
            if (status is None or self.task_status[x] in status) and (age is None or \
                ((age > 0 and time.time() - self.task_date.get(x, time.time()) > age)
                  or (age < 0 and time.time() - self.task_date.get(x, time.time()) < -age)))],
                key=lambda x: -x[2])

    def get_tasks(self):
        with threading.Lock():
            pending = copy.deepcopy(self.pending_tasks + list(self.submitting_tasks.keys()))
            running = copy.deepcopy(self.running_tasks)
        return pending, running

    def run(self):
        # get all system tasks that might have been running ...
        # this will be run only once when the task engine starts
        status_output = self.query_tasks([], verbosity=2, start_time=True)
        with threading.Lock():
            for line in status_output.split('\n'):
                if not line.strip():
                    continue
                try:
                    tid, _, ttm, tst = line.split('\t')
                    self.task_status[tid] = tst
                    self.task_date[tid] = float(ttm)
                except Exception as e:
                    env.logger.warning(f'Unrecognized response "{line}" ({e.__class__.__name__}): {e}')
        self._last_status_check = time.time()
        self.engine_ready.set()
        while True:
            # if no new task, does not do anything.
            if self.running_tasks and time.time() - self._last_status_check > self.status_check_interval:
                if self._status_checker is None:
                    self._status_checker = self._thread_workers.submit(self.query_tasks, self.running_tasks, 1)
                    continue
                elif self._status_checker.running():
                    time.sleep(0.01)
                    continue
                else:
                    status_output = self._status_checker.result()
                    self._status_checker = None
                #
                for line in status_output.splitlines():
                    if not line.strip():
                        continue
                    try:
                        tid, tst = line.split('\t')
                        if tid not in self.running_tasks:
                            env.logger.trace(f'Task {tid} removed since status check.')
                            continue
                        self.update_task_status(tid, tst)
                    except Exception as e:
                        env.logger.warning(f'Unrecognized response "{line}" ({e.__class__.__name__}): {e}')
                self.summarize_status()
                self._last_status_check = time.time()
            else:
                # if nothing to do, sleep to avoid empty loop. This will reduce CPU usage from 100% to 0.3%
                # sleeping longer will slightly reduce CPU usage but decrease responsiveness by a bit
                time.sleep(0.01)

            if self.submitting_tasks:
                with threading.Lock():
                    submitted = []
                    for k in self.submitting_tasks:
                        if not self.submitting_tasks[k].running():
                            submitted.append(k)
                            if self.submitting_tasks[k].result():
                                for tid in k:
                                    if tid in self.canceled_tasks:
                                        # task is canceled while being prepared
                                        self.notify(['change-status', self.agent.alias, tid, 'aborted'])
                                    else:
                                        self.running_tasks.append(tid)
                                        self.notify(['change-status', self.agent.alias, tid, 'submitted'])
                            else:
                                for tid in k:
                                    self.notify(['change-status', self.agent.alias, tid, 'failed'])
                                    self.task_status[tid] = 'failed'
                        #else:
                        #    env.logger.trace('{} is still being submitted.'.format(k))
                    for k in submitted:
                        self.submitting_tasks.pop(k)

            if self.pending_tasks:
                # check status
                num_active_tasks = len(self.submitting_tasks) + len(self.running_tasks)
                if num_active_tasks >= self.max_running_jobs:
                    continue

                # assign tasks to self.max_running_jobs workers
                slots = [[] for i in range(self.max_running_jobs)]
                sample_slots = list(range(self.max_running_jobs))
                random.shuffle(sample_slots)
                for i,tid in enumerate(self.pending_tasks[:self.batch_size * self.max_running_jobs]):
                    if self.task_status[tid] == 'running':
                        self.notify(f'{tid} ``runnng``')
                    elif tid in self.canceled_tasks:
                        # the job is canceled while being prepared to run
                        self.notify(f'{tid} ``canceled``')
                    else:
                        # randomly spread to tasks, but at most one.
                        slots[sample_slots[i % self.max_running_jobs]].append(tid)
                for slot in slots:
                    if not slot:
                        continue
                    for tid in slot:
                        env.logger.trace(f'Start submitting {tid} (status: {self.task_status.get(tid, "unknown")})')
                    self.submitting_tasks[tuple(slot)] = self._thread_workers.submit(self.execute_tasks, slot)
                #
                with threading.Lock():
                    for slot in slots:
                        for tid in slot:
                            self.pending_tasks.remove(tid)

    def submit_task(self, task_id):
        # we wait for the engine to start
        self.engine_ready.wait()

        # submit tasks simply add task_id to pending task list
        with threading.Lock():
            # if already in
            #if task_id in self.running_tasks or task_id in self.pending_tasks:
            #    self.notify('{} ``{}``'.format(task_id, self.task_status[task_id]))
            #    self.notify(['new-status', task_id, self.task_status[task_id]])
            #    return self.task_status[task_id]
            #
            if task_id in self.task_status and self.task_status[task_id]:
                if self.task_status[task_id] == 'running':
                    self.running_tasks.append(task_id)
                    self.notify(f'{task_id} ``already runnng``')
                    self.notify(['new-status', self.agent.alias, task_id, 'running',
                            self.task_date.get(task_id, time.time())])
                    return 'running'
                # there is a case when the job is already completed (complete-old), but
                # because we do not know if the user asks to rerun (-s force), we have to
                # resubmit the job. In the case of not-rerun, the task would be marked
                # completed very soon.
                elif self.task_status[task_id] == 'completed':
                    if env.config['sig_mode'] != 'force':
                        self.notify(f'{task_id} ``already completed``')
                        return 'completed'
                    elif task_id in env.config.get('resumed_tasks', []):
                        # force re-execution, but it is possible that this task has been
                        # executed but quit in no-wait mode (or canceled by user). More
                        # importantly, the Jupyter notebook would re-run complted workflow
                        # even if it has "-s force" signature.
                        #self.notify(['new-status', task_id, 'completed'])
                        self.notify(f'{task_id} ``resume with completed``')
                        return 'completed'
                    else:
                        self.notify(f'{task_id} ``re-execute completed``')
                else:
                    self.notify(f'{task_id} ``restart`` from status ``{self.task_status[task_id]}``')

            #self.notify('{} ``queued``'.format(task_id))
            self.pending_tasks.append(task_id)
            if task_id in self.canceled_tasks:
                self.canceled_tasks.remove(task_id)
            self.task_status[task_id] = 'pending'
            self.notify(['new-status', self.agent.alias, task_id, 'pending',
                    self.task_date.get(task_id, time.time())])
            return 'pending'

    def summarize_status(self):
        from collections import Counter
        statuses = Counter(self.task_status.values())
        env.logger.debug(
            ' '.join(f'{x}: {y}' for x, y in statuses.items()))

    def check_task_status(self, task_id, unknown='pending'):
        # we wait for the engine to start
        self.engine_ready.wait()
        try:
            with threading.Lock():
                return self.task_status[task_id]
        except Exception:
            # job not yet submitted
            return unknown

    def update_task_status(self, task_id, status):
        #
        env.logger.trace(f'STATUS {task_id}\t{status}\n')
        #
        with threading.Lock():
            if task_id in self.canceled_tasks and status != 'aborted':
                env.logger.debug(f'Task {task_id} is still not killed (status {status})')
                status = 'aborted'
            if status != 'missing':
                if task_id in self.task_status and self.task_status[task_id] == status:
                    self.notify(['pulse-status', self.agent.alias, task_id, status])
                else:
                    self.notify(['change-status', self.agent.alias, task_id, status])
            self.task_status[task_id] = status
            if status == 'pening' and task_id not in self.pending_tasks:
                self.pending_tasks.append(task_id)
            if status == 'running' and task_id not in self.running_tasks:
                self.running_tasks.append(task_id)
            # terminal states, remove tasks from task list
            if status in ('completed', 'failed', 'aborted', 'signature-mismatch') and task_id in self.running_tasks:
                self.running_tasks.remove(task_id)

    def remove_tasks(self, tasks):
        with threading.Lock():
            for task in tasks:
                self.notify(['remove-task', self.agent.alias, task])
                #if task in self.task_status:
                #    self.task_status.pop(task)
                #if task in self.running_tasks:
                #    self.running_tasks.remove(task)

    def query_tasks(self, tasks=None, verbosity=1, html=False, start_time=False, age=None, tags=None, status=None):
        try:
            return self.agent.check_output("sos status {} -v {} {} {} {} {} {}".format(
                '' if tasks is None else ' '.join(tasks), verbosity,
                '--html' if html else '',
                '--start-time' if start_time else '',
                f'--age {age}' if age else '',
                f'--tags {" ".join(tags)}' if tags else '',
                f'--status {" ".join(status)}' if status else '',
                ))
        except subprocess.CalledProcessError as e:
            env.logger.warning(f'Failed to query status of tasks on {self.alias}')
            return ''

    def kill_tasks(self, tasks, tags=None, all_tasks=False):
        # we wait for the engine to start
        self.engine_ready.wait()

        for task in tasks:
            with threading.Lock():
                self.task_status[task] = 'aborted'
        for task in tasks:
            with threading.Lock():
                if task in self.pending_tasks:
                    self.pending_tasks.remove(task)
                    env.logger.debug(f'Cancel pending task {task}')
                elif task in self.submitting_tasks:
                    # this is more troublesome because the task is being
                    # submitted at a new thread.
                    env.logger.debug(f'Cancel submission of task {task}')
                elif task in self.running_tasks:
                    env.logger.debug(f'Killing running task {task}')
                else:
                    # it is not in the system, so we need to know what the
                    # status of the task before we do anything...
                    pass

        self.canceled_tasks.extend(tasks)
        #
        cmd = "sos kill {} {} {}".format(' '.join(tasks),
                f'--tags {" ".join(tags)}' if tags else '',
                '-a' if all_tasks else '')

        try:
            ret = self.agent.check_output(cmd)
            env.logger.debug(f'"{cmd}" executed with response "{ret}"')
        except subprocess.CalledProcessError:
            env.logger.error(f'Failed to kill task {tasks}')
            return ''
        return ret

    def resume_task(self, task):
        # we wait for the engine to start
        self.engine_ready.wait()
        with threading.Lock():
            # it is possible that a task is aborted from an opened notebook with aborted status
            if task not in self.task_status or \
                    self.task_status[task] not in ('completed', 'failed', 'signature-mismatch', 'aborted'):
                env.logger.warning(f'Resume task called for non-canceled or non-completed/failed task {task}')
                return
            # the function might have been used multiple times (frontend multiple clicks)
            if task in self.canceled_tasks:
                self.canceled_tasks.remove(task)
            if task not in self.pending_tasks:
                self.pending_tasks.append(task)
            # tells the engine that preparation of task can fail
            self.resuming_tasks.add(task)
            self.task_status[task] = 'pending'

    def execute_tasks(self, task_ids):
        # we wait for the engine to start
        self.engine_ready.wait()
        # this is base class, the derived class will actually submit the tasks

        # if the task is being resumed, perhaps from another local host,
        # the preparation process can fail (e.g. no def file), but this
        # does not really matter. #587
        for task_id in task_ids:
            if task_id in self.resuming_tasks:
                self.resuming_tasks.remove(task_id)
                try:
                    self.agent.prepare_task(task_id)
                except Exception:
                    pass
            else:
                if not self.agent.prepare_task(task_id):
                    return False
        return True

    def purge_tasks(self, tasks, purge_all=False, age=None, status=None, tags=None, verbosity=2):
        try:
            return self.agent.check_output("sos purge {} {} {} {} {} -v {}".format(
                ' '.join(tasks), '--all' if purge_all else '',
                f'--age {age}' if age is not None else '',
                f'--status {" ".join(status)}' if status is not None else '',
                f'--tags {" ".join(tags)}' if tags is not None else '',
                verbosity))
        except subprocess.CalledProcessError:
            env.logger.error(f'Failed to purge tasks {tasks}')
            return ''


class BackgroundProcess_TaskEngine(TaskEngine):
    def __init__(self, agent):
        super(BackgroundProcess_TaskEngine, self).__init__(agent)
        if 'job_template' in self.config:
            self.job_template = self.config['job_template'].replace('\r\n', '\n')
        else:
            self.job_template = None
        #
        if 'batch_size' in self.config:
            self.batch_size = self.config['batch_size']
        else:
            # default allow stacking of up to 1000 jobs
            self.batch_size = 1000

    def execute_tasks(self, task_ids):
        if not super(BackgroundProcess_TaskEngine, self).execute_tasks(task_ids):
            env.logger.trace(f'Failed to prepare task {task_ids}')
            return False
        if self.job_template:
            if not self._submit_task_with_template(task_ids):
                return False
        else:
            if not self._submit_task(task_ids):
                return False
        return True

    def _submit_task(self, task_ids):
        # if no template, use a default command
        cmd = f"sos execute {' '.join(task_ids)} -v {env.verbosity} -s {env.config['sig_mode']} {'--dryrun' if env.config['run_mode'] == 'dryrun' else ''}"
        env.logger.trace(f'Execute "{cmd}" (waiting={self.wait_for_task})')
        self.agent.run_command(cmd, wait_for_task = self.wait_for_task)
        return True

    def _submit_task_with_template(self, task_ids):
        '''Submit tasks by interpolating a shell script defined in job_template'''
        runtime = self.config
        runtime.update({
            'cur_dir': os.getcwd(),
            'verbosity': env.verbosity,
            'sig_mode': env.config.get('sig_mode', 'default'),
            'run_mode': env.config.get('run_mode', 'run'),
            'home_dir': os.path.expanduser('~')})
        if '_runtime' in env.sos_dict:
            runtime.update({x:env.sos_dict['_runtime'][x] for x in ('nodes', 'cores', 'mem', 'walltime') if x in env.sos_dict['_runtime']})
        if 'nodes' not in runtime:
            runtime['nodes'] = 1
        if 'cores' not in runtime:
            runtime['cores'] = 1

        # let us first prepare a task file
        job_text = ''
        for task_id in task_ids:
            runtime['task'] = task_id
            try:
                job_text += cfg_interpolate(self.job_template, runtime)
                job_text += '\n'
            except Exception as e:
                raise ValueError(f'Failed to generate job file for task {task_id}: {e}')

        filename = task_ids[0] + ('.sh' if len(task_ids) == 1 else f'-{task_ids[-1]}.sh')
        # now we need to write a job file
        job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, filename)
        # do not translate newline under windows because the script will be executed
        # under linux/mac
        with open(job_file, 'w', newline='') as job:
            job.write(job_text)

        # then copy the job file to remote host if necessary
        self.agent.send_task_file(job_file)

        try:
            cmd = f'bash ~/.sos/tasks/{filename}'
            self.agent.run_command(cmd, wait_for_task = self.wait_for_task)
        except Exception as e:
            raise RuntimeError(f'Failed to submit task {task_ids}: {e}')
        return True
