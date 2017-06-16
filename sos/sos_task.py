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
from io import StringIO
from tokenize import generate_tokens
from collections.abc import Sequence, Mapping
import concurrent.futures

from sos.utils import env, short_repr, get_traceback, sample_of_file, tail_of_file, linecount_of_file, \
    format_HHMMSS, expand_time, expand_size
from sos.sos_eval import SoS_exec, SoS_eval

from .target import textMD5, RuntimeInfo, Undetermined, FileTarget, UnknownTarget, remote
from .sos_eval import interpolate
from .monitor import ProcessMonitor

from collections import OrderedDict
import subprocess


monitor_interval = 5
resource_monitor_interval = 60

class TaskParams(object):
    '''A parameter object that encaptulates parameters sending to
    task executors. This would makes the output of workers, especially
    in the web interface much cleaner (issue #259)'''
    def __init__(self, name, task, sos_dict, sigil):
        self.name = name
        self.task = task
        self.sos_dict = sos_dict
        self.sigil = sigil

    def save(self, job_file):
        with open(job_file, 'wb') as jf:
            try:
                pickle.dump(self, jf)
            except Exception as e:
                env.logger.warning(e)
                raise

    def __repr__(self):
        return self.name

class MasterTaskParams(TaskParams):
    def __init__(self, num_workers=0):
        self.ID = 'M_0'
        self.task = ''
        self.sos_dict = {'_runtime': {}, '_input': [], '_output': [], '_depends': []}
        self.sigil = None
        self.num_workers = num_workers
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
        # sigil
        if self.sigil is None:
            self.sigil = params.sigil
        elif self.sigil != params.sigil:
            raise ValueError('Cannot push a task with different sigil {}'.format(params.sigil))
        #
        # walltime
        if not self.task_stack:
            for key in ('walltime', 'max_walltime', 'cores', 'max_cores', 'mem', 'max_mem', 'preserved_vars', 'name'):
                if key in params.sos_dict['_runtime'] and params.sos_dict['_runtime'][key] is not None:
                    self.sos_dict['_runtime'][key] = params.sos_dict['_runtime'][key]
        else:
            for key in ('walltime', 'max_walltime', 'cores', 'max_cores', 'mem', 'max_mem', 'name'):
                val0 = self.task_stack[0][1].sos_dict['_runtime'].get(key, None)
                val = params.sos_dict['_runtime'].get(key, None)
                if val0 != val:
                    raise ValueError('All tasks should have the same resource {}'.format(key))
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
                    self.sos_dict['_runtime']['name'] = '{}_{}'.format(val0, len(self.task_stack) + 1)
        #
        # input, output, preserved vars etc
        for key in ['_input', '_output', '_depends']:
            if key in params.sos_dict and isinstance(params.sos_dict[key], list):
                self.sos_dict[key].extend(params.sos_dict[key])
        #
        self.task_stack.append((task_id, params))
        #
        self.ID = 'M{}_{}'.format(len(self.task_stack), self.task_stack[0][0])


def loadTask(filename):
    with open(filename, 'rb') as task:
        return pickle.load(task)

def collect_task_result(task_id, sigil, sos_dict):
    shared = {}
    if 'shared' in env.sos_dict['_runtime']:
        vars = env.sos_dict['_runtime']['shared']
        if isinstance(vars, str):
            if vars not in env.sos_dict:
                raise ValueError('Unavailable shared variable {} after the completion of task {}'.format(vars, task_id))
            shared[vars] = copy.deepcopy(env.sos_dict[vars])
        elif isinstance(vars, Mapping):
            for var, val in vars.items():
                if var != val:
                    env.sos_dict.set(var, SoS_eval(val, sigil))
                if var not in env.sos_dict:
                    raise ValueError('Unavailable shared variable {} after the completion of task {}'.format(var, task_id))
                shared[var] = copy.deepcopy(env.sos_dict[var])
        elif isinstance(vars, Sequence):
            # if there are dictionaries in the sequence, e.g.
            # shared=['A', 'B', {'C':'D"}]
            for item in vars:
                if isinstance(item, str):
                    if item not in env.sos_dict:
                        raise ValueError('Unavailable shared variable {} after the completion of task {}'.format(item, task_id))
                    shared[item] = copy.deepcopy(env.sos_dict[item])
                elif isinstance(item, Mapping):
                    for var, val in item.items():
                        if var != val:
                            env.sos_dict.set(var, SoS_eval(val, sigil))
                        if var not in env.sos_dict:
                            raise ValueError('Unavailable shared variable {} after the completion of task {}'.format(var, task_id))
                        shared[var] = copy.deepcopy(env.sos_dict[var])
                else:
                    raise ValueError('Option shared should be a string, a mapping of expression, or a list of string or mappings. {} provided'.format(vars))
        else:
            raise ValueError('Option shared should be a string, a mapping of expression, or a list of string or mappings. {} provided'.format(vars))
        env.logger.debug('task {} (index={}) return shared variable {}'.format(task_id, env.sos_dict['_index'], shared))
    # the difference between sos_dict and env.sos_dict is that sos_dict (the original version) can have remote() targets
    # which should not be reported.
    if env.sos_dict['_output'] is None:
        output = {}
    elif sos_dict['_output'] is None:
        output = {}
    else:
        output = {x:FileTarget(x).signature() for x in sos_dict['_output'] if isinstance(x, str)}
    return {'ret_code': 0, 'task': task_id,
            'output': output,
            'shared': {env.sos_dict['_index']: shared} }


def execute_task(task_id, verbosity=None, runmode='run', sigmode=None, monitor_interval=5,
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
        env.logger.trace('Executing subtask {}'.format(task_id))

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
        with open(master_out, 'ab') as out, open(master_err, 'ab') as err:
            def copy_out_and_err(result):
                tid = result['task']
                out.write('{}: {}\n'.format(tid, 'completed' if result['ret_code'] == 0 else 'failed').encode())
                out.write('output: {}\n'.format(result['output']).encode())
                sub_out = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', tid + '.out')
                if os.path.isfile(sub_out):
                    with open(sub_out, 'rb') as sout:
                        out.write(sout.read())

                sub_err = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', tid + '.err')
                err.write('{}: {}\n'.format(tid, 'completed' if result['ret_code'] == 0 else 'failed').encode())
                if os.path.isfile(sub_err):
                    with open(sub_err, 'rb') as serr:
                        err.write(serr.read())

            if params.num_workers > 1:
                from multiprocessing.pool import Pool
                p = Pool(params.num_workers)
                try:
                    results = []
                    for t in params.task_stack:
                        results.append(p.apply_async(execute_task, (t, verbosity, runmode,
                            sigmode, monitor_interval, resource_monitor_interval), callback=copy_out_and_err))
                    for idx,r in enumerate(results):
                        results[idx] = r.get()
                except Exception as e:
                    if env.verbosity > 2:
                        sys.stderr.write(get_traceback())
                    env.logger.error('{} ``failed`` with {} error {}'.format(task_id, e.__class__.__name__, e))
                    return {'ret_code': 1, 'exception': e}
            else:
                results = []
                for tid, tdef in params.task_stack:
                    try:
                        res = execute_task((tid, tdef), verbosity=verbosity, runmode=runmode,
                            sigmode=sigmode, monitor_interval=monitor_interval,
                            resource_monitor_interval=resource_monitor_interval)
                        copy_out_and_err(res)
                        results.append(res)
                    except Exception as e:
                        out.write('{}: failed\n'.format(tid).encode())
                        copy_out_and_err({'task': tid, 'ret_code': 1, 'output': []})
                        if env.verbosity > 2:
                            sys.stderr.write(get_traceback())
                        env.logger.error('{} ``failed`` with {} error {}'.format(task_id, e.__class__.__name__, e))
                        results.append({'ret_code': 1, 'exception': e})
        #
        # now we collect result
        all_res = {'ret_code': 0, 'output': {}, 'subtasks': {}, 'shared': {}}
        for tid, x in zip(params.task_stack, results):
            all_res['ret_code'] += x['ret_code']
            all_res['output'].update(x['output'])
            all_res['subtasks'][tid[0]] = x
            all_res['shared'].update(x['shared'])
        return all_res

    task, sos_dict, sigil = params.task, params.sos_dict, params.sigil

    if '_runtime' not in sos_dict:
        sos_dict['_runtime'] = {}

    # pulse thread
    m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
        resource_monitor_interval=resource_monitor_interval,
        max_walltime=sos_dict['_runtime'].get('max_walltime', None),
        max_mem=sos_dict['_runtime'].get('max_mem', None),
        max_procs=sos_dict['_runtime'].get('max_procs', None))

    m.start()
    if verbosity is not None:
        env.verbosity = verbosity
    if sigmode is not None:
        env.config['sig_mode'] = sigmode
    env.config['run_mode'] = runmode
    #
    # when the statements are executed in tasks, they are no longer remote_targets
    env.config['remote_targets'] = False

    if subtask:
        env.logger.debug('{} ``started``'.format(task_id))
    else:
        env.logger.info('{} ``started``'.format(task_id))

    env.sos_dict.quick_update(sos_dict)

    # if targets are defined as `remote`, they should be resolved during task execution
    def resolve_remote(x):
        if isinstance(x, remote):
            x = x.resolve()
            if isinstance(x, str):
                x = interpolate(x, sigil, env.sos_dict)
        return x

    for key in ['input', '_input',  'output', '_output', 'depends', '_depends']:
        if key in sos_dict and isinstance(sos_dict[key], list):
            # resolve remote() target
            env.sos_dict.set(key, [resolve_remote(x) for x in sos_dict[key]])

    skipped = False
    if env.config['sig_mode'] == 'ignore':
        sig = None
    else:
        tokens = [x[1] for x in generate_tokens(StringIO(task).readline)]
        # try to add #task so that the signature can be different from the step
        # if everything else is the same
        sig = RuntimeInfo(textMD5('#task\n' + ' '.join(tokens)), task,
            env.sos_dict['_input'], env.sos_dict['_output'], env.sos_dict['_depends'], env.sos_dict['__signature_vars__'])
        sig.lock()

        idx = env.sos_dict['_index']
        if env.config['sig_mode'] == 'default':
            matched = sig.validate()
            if isinstance(matched, dict):
                # in this case, an Undetermined output can get real output files
                # from a signature
                env.sos_dict.set('_input', matched['input'])
                env.sos_dict.set('_depends', matched['depends'])
                env.sos_dict.set('_output', matched['output'])
                env.sos_dict.set('_local_input', matched['local_output'])
                env.sos_dict.set('_local_output', matched['local_output'])
                env.sos_dict.set('local_input', env.sos_dict['_local_input'])
                env.sos_dict.set('local_output', env.sos_dict['_local_output'])
                env.sos_dict.update(matched['vars'])
                env.logger.info('Task ``{}`` (index={}) is ``ignored`` due to saved signature'.format(env.sos_dict['step_name'], idx))
                skipped = True
        elif env.config['sig_mode'] == 'assert':
            matched = sig.validate()
            if isinstance(matched, str):
                raise RuntimeError('Signature mismatch: {}'.format(matched))
            else:
                env.sos_dict.set('_input', matched['input'])
                env.sos_dict.set('_depends', matched['depends'])
                env.sos_dict.set('_output', matched['output'])
                env.sos_dict.set('_local_input', matched['local_output'])
                env.sos_dict.set('_local_output', matched['local_output'])
                env.sos_dict['local_input'].extend(env.sos_dict['_local_input'])
                env.sos_dict['local_output'].extend(env.sos_dict['_local_output'])
                env.sos_dict.update(matched['vars'])
                env.logger.info('Step ``{}`` (index={}) is ``ignored`` with matching signature'.format(env.sos_dict['step_name'], idx))
                skipped = True
        elif env.config['sig_mode'] == 'build':
            # build signature require existence of files
            if sig.write(
                env.sos_dict['_local_input_{}'.format(idx)],
                env.sos_dict['_local_output_{}'.format(idx)],
                rebuild=True):
                env.logger.info('Task ``{}`` (index={}) is ``ignored`` with signature constructed'.format(env.sos_dict['step_name'], idx))
                skipped = True
            else:
                env.logger.info('Task ``{}`` (index={}) is ``executed`` with failed signature constructed'.format(env.sos_dict['step_name'], idx))
        elif env.config['sig_mode'] == 'force':
            skipped = False
        else:
            raise RuntimeError('Unrecognized signature mode {}'.format(env.config['sig_mode']))

    if skipped:
        env.logger.info('{} ``skipped``'.format(task_id))
        return collect_task_result(task_id, sigil, sos_dict)

    try:
        # go to 'cur_dir'
        if '_runtime' in sos_dict and 'cur_dir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['cur_dir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
                except Exception as e:
                    raise RuntimeError('Failed to create cur_dir {}'.format(sos_dict['_runtime']['cur_dir']))
            os.chdir(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
        orig_dir = os.getcwd()

        if runmode != 'dryrun':
            # we will need to check existence of targets because the task might
            # be executed on a remote host where the targets are not available.
            for target in (sos_dict['_input'] if isinstance(sos_dict['_input'], list) else []) + \
                (sos_dict['_depends'] if isinstance(sos_dict['_depends'], list) else []):
                # if the file does not exist (although the signature exists)
                # request generation of files
                if isinstance(target, str):
                    if not FileTarget(target).exists('target'):
                        # remove the signature and regenerate the file
                        FileTarget(target).remove_sig()
                        raise UnknownTarget(target)
                elif not target.exists('target'):
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
                            env.logger.warning('Failed to create directory {}: {}'.format(parent_dir, e))


        # go to user specified workdir
        if '_runtime' in sos_dict and 'workdir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['workdir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['workdir']))
                except Exception as e:
                    raise RuntimeError('Failed to create workdir {}'.format(sos_dict['_runtime']['workdir']))
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
                    raise ValueError('Unacceptable input for option prepend_path: {}'.format(sos_dict['_runtime']['prepend_path']))

        # task output
        env.sos_dict.set('__std_out__', os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.out'))
        env.sos_dict.set('__std_err__', os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err'))
        env.logfile = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err')

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        # step process
        SoS_exec(task, sigil)

        if subtask:
            env.logger.debug('{} ``completed``'.format(task_id))
        else:
            env.logger.info('{} ``completed``'.format(task_id))

    except Exception as e:
        if env.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error('{} ``failed`` with {} error {}'.format(task_id, e.__class__.__name__, e))
        return {'ret_code': 1, 'exception': e, 'task': task_id, 'shared': {}}
    except KeyboardInterrupt:
        env.logger.error('{} ``interrupted``'.format(task_id))
        raise
    finally:
        env.sos_dict.set('__step_sig__', None)
        os.chdir(orig_dir)

    if sig:
        sig.write(env.sos_dict['_local_input_{}'.format(env.sos_dict['_index'])],
            env.sos_dict['_local_output_{}'.format(env.sos_dict['_index'])])
        sig.release()

    # the final result should be relative to cur_dir, not workdir
    # because output is defined outside of task
    return collect_task_result(task_id, sigil, sos_dict)

def check_task(task):
    #
    # status of the job, please refer to https://github.com/vatlab/SOS/issues/529
    # for details.
    #
    task_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.task')
    if not os.path.isfile(task_file):
        return 'non-exist'
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
            from .target import FileTarget
            with open(res_file, 'rb') as result:
                res = pickle.load(result)
            if ('ret_code' in res and res['ret_code'] == 0) or ('succ' in res and res['succ'] == 0):
                if isinstance(res['output'], dict):
                    for x,y in res['output'].items():
                        if not FileTarget(x).exists() or FileTarget(x).signature() != y:
                            env.logger.debug('{} not found or signature mismatch'.format(x))
                            return 'result-mismatch'
                            # otherwise, it can be submitted or pending...
                    # this is called "completed" remotely but will be
                    # translated to either completed or result-ready locally
                    return 'completed'
                else:
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
            env.logger.warning('{} is created in the future. Your system time might be problematic'.format(pulse_file))
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

def check_tasks(tasks, verbosity=1, html=False, start_time=False, age=None):
    # verbose is ignored for now
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    if not tasks:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getctime(x)) for x in tasks]
        if not all_tasks:
            return
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '{}*.task'.format(t)))
            matched = [(os.path.basename(x)[:-5], os.path.getctime(x)) for x in matched]
            if not matched:
                all_tasks.append((t, None))
            else:
                all_tasks.extend(matched)
    if age is not None:
        from sos.utils import expand_time
        age = expand_time(age, default_unit='d')
        if age > 0:
            all_tasks = [x for x in all_tasks if time.time() - x[1] >= age]
        else:
            all_tasks = [x for x in all_tasks if time.time() - x[1] <= -age]

    all_tasks = sorted(list(set(all_tasks)), key=lambda x: 0 if x[1] is None else x[1])
    if not all_tasks:
        env.logger.info('No matching tasks')
        return
    # at most 20 threads
    p = Pool(min(20, len(all_tasks)))
    status = p.map(check_task, [x[0] for x in all_tasks])
    #
    if html:
        verbosity = -1
    if verbosity == 0:
        print('\n'.join(status))
    elif verbosity == 1:
        for s, (t, d) in zip(status, all_tasks):
            print('{}\t{}'.format(t, s))
    elif verbosity == 2:
        from .utils import PrettyRelativeTime
        for s, (t, d) in zip(status, all_tasks):
            if start_time:
                if d is None:
                    print('{}\t{}\t{}'.format(t, time.time(), s))
                else:
                    print('{}\t{}\t{}'.format(t, d, s))
            else:
                if d is None:
                    print('{}\t{:>15}\t{}'.format(t, '', s))
                else:
                    print('{}\t{:>15} ago\t{}'.format(t, PrettyRelativeTime(time.time() - d), s))
    elif verbosity > 2:
        from .utils import PrettyRelativeTime
        import pprint
        import glob
        from .monitor import summarizeExecution

        for s, (t, d) in zip(status, all_tasks):
            print('{}\t{}\n'.format(t, s))
            if d is not None:
                print('Started {} ago'.format(PrettyRelativeTime(time.time() - d)))
            task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                continue
            params = loadTask(task_file)
            print('TASK:\n=====')
            print(params.task)
            print()
            print('ENVIRONMENT:\n============')
            job_vars = params.sos_dict
            for k in sorted(job_vars.keys()):
                v = job_vars[k]
                print('{:22}{}'.format(k, short_repr(v) if verbosity == 3 else pprint.pformat(v)))
            print()
            print('EXECUTION STATS:\n================')
            print(summarizeExecution(t, status=s))
            if verbosity == 4:
                # if there are other files such as job file, print them.
                files = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.*'))
                for f in sorted([x for x in files if os.path.splitext(x)[-1] not in ('.res',
                    '.task', '.pulse', '.status', '.def')]):
                    print('{}:\n{}'.format(os.path.basename(f), '='*(len(os.path.basename(f))+1)))
                    try:
                        with open(f) as fc:
                            print(fc.read())
                    except:
                        print('Binary file')
    else:
        # HTML output
        from .utils import PrettyRelativeTime, isPrimitive
        from .monitor import summarizeExecution
        import pprint
        import glob
        print('<table width="100%">')
        def row(th=None, td=None):
            if td is None:
                print('<tr><th align="right" width="30%"><font color="blue">{}</font></th><td></td></tr>'.format(th))
            elif th is None:
                print('<tr><td colspan="2" align="left"  width="30%">{}</td></tr>'.format(td))
            else:
                print('<tr><th align="right"  width="30%">{}</th><td align="left">{}</td></tr>'.format(th, td))
        for s, (t, d) in zip(status, all_tasks):
            row('ID', t)
            row('Status', s)
            if d is not None:
                row('Start', '{:>15} ago'.format(PrettyRelativeTime(time.time() - d)))
            task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                continue
            with open(task_file, 'rb') as task:
                params = pickle.load(task)
            row('Task')
            row(td='<pre style="text-align:left">{}</pre>'.format(params.task))
            row('Environment')
            job_vars = params.sos_dict
            for k in sorted(job_vars.keys()):
                v = job_vars[k]
                if not k.startswith('__') and not k == 'CONFIG':
                    if k == '_runtime':
                        for _k, _v in v.items():
                            if isPrimitive(_v) and _v not in (None, '', [], (), {}):
                                row(_k, _v)
                    elif isPrimitive(v) and v not in (None, '', [], (), {}):
                        row(k, '<pre style="text-align:left">{}</pre>'.format(pprint.pformat(v)))
            summary = summarizeExecution(t, status=s)
            if summary:
                row('Execution')
                for line in summary.split('\n'):
                    fields = line.split(None, 1)
                    if fields[0] == 'task':
                        continue
                    row(fields[0], '' if fields[1] is None else fields[1])
                # this is a placeholder for the frontend to draw figure
                row(td='<div id="res_{}"></div>'.format(t))
            #
            files = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', t + '.*'))
            for f in sorted([x for x in files if os.path.splitext(x)[-1] not in ('.def', '.res', '.task', '.pulse', '.status')]):
                numLines = linecount_of_file(f)
                row(os.path.basename(f), '(empty)' if numLines == 0 else '{} lines{}'.format(numLines, '' if numLines < 200 else ' (showing last 200)'))
                try:
                    row(td='<small><pre style="text-align:left">{}</pre></small>'.format(tail_of_file(f, 200, ansi2html=True)))
                except:
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
            except:
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

    var cpu = [''' + ','.join(['[{},{}]'.format(x*1000,y) for x,y in zip(etime, cpu)]) + '''];
    var mem = [''' + ','.join(['[{},{}]'.format(x*1000,y) for x,y in zip(etime, mem)]) + '''];

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


def kill_tasks(tasks):
    #
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    if not tasks:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [os.path.basename(x)[:-5] for x in tasks]
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '{}*.task'.format(t)))
            matched = [os.path.basename(x)[:-5] for x in matched]
            if not matched:
                env.logger.warning('{} does not match any existing task'.format(t))
            else:
                all_tasks.extend(matched)
    if not all_tasks:
        env.logger.warning('No task to kill')
        return
    all_tasks = sorted(list(set(all_tasks)))
    p = Pool(len(all_tasks))
    killed = p.map(kill_task, all_tasks)
    for s, t in zip(killed, all_tasks):
        print('{}\t{}'.format(t, s))

def kill_task(task):
    status = check_task(task)
    if status == 'pending':
        return 'cancelled'
    # remove job file as well
    job_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.sh')
    if os.path.isfile(job_file):
        try:
            os.remove(job_file)
        except:
            pass
    if status != 'running':
        return status
    # job is running
    pulse_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.pulse')
    from stat import S_IREAD, S_IRGRP, S_IROTH
    os.chmod(pulse_file, S_IREAD|S_IRGRP|S_IROTH)
    return 'killed'


def purge_tasks(tasks, purge_all=False, age=None, status=None, verbosity=2):
    # verbose is ignored for now
    import glob
    if tasks:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '{}*.task'.format(t)))
            matched = [(os.path.basename(x)[:-5], os.path.getctime(x)) for x in matched]
            all_tasks.extend(matched)
    elif purge_all or age or status:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getctime(x)) for x in tasks]
    else:
        sys.exit('Please specify a list of tasks and/or option -all, --age, or --status')
    if age is not None:
        from sos.utils import expand_time
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
    #
    # remoe all task files
    all_tasks = set([x[0] for x in all_tasks])
    if all_tasks:
        #
        # find all related files, including those in nested directories
        from collections import defaultdict
        to_be_removed = defaultdict(list)
        for dirname, dirlist, filelist in os.walk(os.path.join(os.path.expanduser('~'), '.sos', 'tasks')):
            for f in filelist:
                ID = os.path.basename(f).split('.', 1)[0]
                if ID in all_tasks:
                    to_be_removed[ID].append(os.path.join(dirname, f))
        #
        for task in all_tasks:
            for f in to_be_removed[task]:
                try:
                    env.logger.trace('Remove {}'.format(f))
                    os.remove(f)
                except Exception as e:
                    env.logger.warning('Failed to purge task {}'.format(task[0]))
            env.logger.info('Task ``{}`` removed.'.format(task))
    else:
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
                    env.logger.warning('Failed to remove {}'.format(f))
            else:
                try:
                    os.remove(f)
                    count += 1
                except Exception as e:
                    env.logger.warning('Failed to remove {}'.format(e))
        if count > 0:
            env.logger.info('{} other files and directories are removed.'.format(count))
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
            self.max_running_jobs = 10
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

    def monitor_tasks(self, tasks=[], status=None, age=None):
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
            from sos.utils import expand_time
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
                    tid, ttm, tst = line.split('\t')
                    self.task_status[tid] = tst
                    self.task_date[tid] = float(ttm)
                except Exception as e:
                    env.logger.warning('Unrecognized response "{}" ({}): {}'.format(line, e.__class__.__name__, e))
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
                            env.logger.trace('Task {} removed since status check.'.format(tid))
                            continue
                        self.update_task_status(tid, tst)
                    except Exception as e:
                        env.logger.warning('Unrecognized response "{}" ({}): {}'.format(line, e.__class__.__name__, e))
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
                                if k in self.canceled_tasks:
                                    # task is canceled while being prepared
                                    if hasattr(env, '__task_notifier__'):
                                        env.__task_notifier__(['change-status', self.agent.alias, k, 'aborted'])
                                else:
                                    self.running_tasks.append(k)
                                    if hasattr(env, '__task_notifier__'):
                                        env.__task_notifier__(['change-status', self.agent.alias, k, 'submitted'])
                            else:
                                if hasattr(env, '__task_notifier__'):
                                    env.__task_notifier__(['change-status', self.agent.alias, k, 'failed'])
                                self.task_status[k] = 'failed'
                        #else:
                        #    env.logger.trace('{} is still being submitted.'.format(k))
                    for k in submitted:
                        self.submitting_tasks.pop(k)

            if self.pending_tasks:
                to_run = []
                # check status
                num_active_tasks = len(self.submitting_tasks) + len(self.running_tasks)
                if num_active_tasks >= self.max_running_jobs:
                    continue

                to_run = self.pending_tasks[ : self.max_running_jobs - num_active_tasks]
                for tid in to_run:
                    if self.task_status[tid] == 'running':
                        env.logger.info('{} ``runnng``'.format(tid))
                    elif tid in self.canceled_tasks:
                        # the job is canceled while being prepared to run
                        env.logger.info('{} ``canceled``'.format(tid))
                    else:
                        env.logger.trace('Start submitting {} (status: {})'.format(tid, self.task_status.get(tid, 'unknown')))
                        self.submitting_tasks[tid] = self._thread_workers.submit(self.execute_task, tid)
                #
                with threading.Lock():
                    for tid in to_run:
                        self.pending_tasks.remove(tid)

    def submit_task(self, task_id):
        # we wait for the engine to start
        self.engine_ready.wait()

        # submit tasks simply add task_id to pending task list
        with threading.Lock():
            # if already in
            #if task_id in self.running_tasks or task_id in self.pending_tasks:
            #    env.logger.info('{} ``{}``'.format(task_id, self.task_status[task_id]))
            #    if hasattr(env, '__task_notifier__'):
            #        env.__task_notifier__(['new-status', task_id, self.task_status[task_id]])
            #    return self.task_status[task_id]
            #
            if task_id in self.task_status and self.task_status[task_id]:
                if self.task_status[task_id] == 'running':
                    self.running_tasks.append(task_id)
                    env.logger.info('{} ``already runnng``'.format(task_id))
                    if hasattr(env, '__task_notifier__'):
                        env.__task_notifier__(['new-status', self.agent.alias, task_id, 'running',
                            self.task_date.get(task_id, time.time())])
                    return 'running'
                # there is a case when the job is already completed (complete-old), but
                # because we do not know if the user asks to rerun (-s force), we have to
                # resubmit the job. In the case of not-rerun, the task would be marked
                # completed very soon.
                elif self.task_status[task_id] == 'completed':
                    if env.config['sig_mode'] != 'force':
                        env.logger.info('{} ``already completed``'.format(task_id))
                        return 'completed'
                    elif task_id in env.config.get('resumed_tasks', []):
                        # force re-execution, but it is possible that this task has been
                        # executed but quit in no-wait mode (or canceled by user). More
                        # importantly, the Jupyter notebook would re-run complted workflow
                        # even if it has "-s force" signature.
                        #if hasattr(env, '__task_notifier__'):
                        #    env.__task_notifier__(['new-status', task_id, 'completed'])
                        env.logger.info('{} ``resume with completed``'.format(task_id))
                        return 'completed'
                    else:
                        env.logger.info('{} ``re-execute completed``'.format(task_id))
                else:
                    env.logger.info('{} ``restart`` from status ``{}``'.format(task_id, self.task_status[task_id]))

            env.logger.info('{} ``queued``'.format(task_id))
            self.pending_tasks.append(task_id)
            if task_id in self.canceled_tasks:
                self.canceled_tasks.remove(task_id)
            self.task_status[task_id] = 'pending'
            if hasattr(env, '__task_notifier__'):
                env.__task_notifier__(['new-status', self.agent.alias, task_id, 'pending',
                    self.task_date.get(task_id, time.time())])
            return 'pending'

    def summarize_status(self):
        from collections import Counter
        statuses = Counter(self.task_status.values())
        env.logger.debug(
            ' '.join('{}: {}'.format(x, y) for x, y in statuses.items()))

    def check_task_status(self, task_id, unknown='pending'):
        # we wait for the engine to start
        self.engine_ready.wait()
        try:
            with threading.Lock():
                return self.task_status[task_id]
        except:
            # job not yet submitted
            return unknown

    def update_task_status(self, task_id, status):
        #
        env.logger.trace('STATUS {}\t{}\n'.format(task_id, status))
        #
        with threading.Lock():
            if task_id in self.canceled_tasks and status != 'aborted':
                env.logger.debug('Task {} is still not killed (status {})'.format(task_id, status))
                status = 'aborted'
            if hasattr(env, '__task_notifier__') and status != 'non-exist':
                if task_id in self.task_status and self.task_status[task_id] == status:
                    env.__task_notifier__(['pulse-status', self.agent.alias, task_id, status])
                else:
                    env.__task_notifier__(['change-status', self.agent.alias, task_id, status])
            self.task_status[task_id] = status
            if status == 'pening' and task_id not in self.pending_tasks:
                self.pending_tasks.append(task_id)
            if status == 'running' and task_id not in self.running_tasks:
                self.running_tasks.append(task_id)
            # terminal states, remove tasks from task list
            if status in ('completed', 'failed', 'aborted', 'result-mismatch') and task_id in self.running_tasks:
                self.running_tasks.remove(task_id)

    def remove_tasks(self, tasks):
        with threading.Lock():
            for task in tasks:
                if hasattr(env, '__task_notifier__'):
                    env.__task_notifier__(['remove-task', self.agent.alias, task])
                #if task in self.task_status:
                #    self.task_status.pop(task)
                #if task in self.running_tasks:
                #    self.running_tasks.remove(task)

    def pending_tasks(self):
        with threading.Lock():
            return self.pending_tasks()

    def query_tasks(self, tasks=None, verbosity=1, html=False, start_time=False, age=None):
        try:
            return self.agent.check_output("sos status {} -v {} {} {} {}".format(
                ' '.join(tasks), verbosity, '--html' if html else '',
                '--start-time' if start_time else '', '--age {}'.format(age) if age else ''))
        except subprocess.CalledProcessError as e:
            env.logger.error('Failed to query status of tasks {}: {}'.format(tasks, e))
            return ''

    def kill_tasks(self, tasks, all_tasks=False):
        # we wait for the engine to start
        self.engine_ready.wait()

        for task in tasks:
            with threading.Lock():
                self.task_status[task] = 'aborted'
        for task in tasks:
            with threading.Lock():
                if task in self.pending_tasks:
                    self.pending_tasks.remove(task)
                    env.logger.debug('Cancel pending task {}'.format(task))
                elif task in self.submitting_tasks:
                    # this is more troublesome because the task is being
                    # submitted at a new thread.
                    env.logger.debug('Cancel submission of task {}'.format(task))
                elif task in self.running_tasks:
                    env.logger.debug('Killing running task {}'.format(task))
                else:
                    # it is not in the system, so we need to know what the
                    # status of the task before we do anything...
                    pass

        self.canceled_tasks.extend(tasks)
        #
        cmd = "sos kill {} {}".format(' '.join(tasks), '-a' if all_tasks else '')
        try:
            ret = self.agent.check_output(cmd)
            env.logger.debug('"{}" executed with response "{}"'.format(cmd, ret))
        except subprocess.CalledProcessError as e:
            env.logger.error('Failed to kill task {}'.format(tasks))
            return ''
        return ret

    def resume_task(self, task):
        # we wait for the engine to start
        self.engine_ready.wait()
        with threading.Lock():
            # it is possible that a task is aborted from an opened notebook with aborted status
            if task not in self.task_status or \
                    self.task_status[task] not in ('completed', 'failed', 'result-mismatch', 'aborted'):
                env.logger.warning('Resume task called for non-canceled or non-completed/failed task {}'.format(task))
                return
            # the function might have been used multiple times (frontend multiple clicks)
            if task in self.canceled_tasks:
                self.canceled_tasks.remove(task)
            if task not in self.pending_tasks:
                self.pending_tasks.append(task)
            # tells the engine that preparation of task can fail
            self.resuming_tasks.add(task)
            self.task_status[task] = 'pending'

    def execute_task(self, task_id):
        # we wait for the engine to start
        self.engine_ready.wait()
        # this is base class

        # if the task is being resumed, perhaps from another local host,
        # the preparation process can fail (e.g. no def file), but this
        # does not really matter. #587
        if task_id in self.resuming_tasks:
            self.resuming_tasks.remove(task_id)
            try:
                self.agent.prepare_task(task_id)
            except:
                pass
            return True
        else:
            return self.agent.prepare_task(task_id)

    def purge_tasks(self, tasks, purge_all=False, age=None, status=None, verbosity=2):
        try:
            return self.agent.check_output("sos purge {} {} {} {} -v {}".format(
                ' '.join(tasks), '--all' if purge_all else '',
                '--age {}'.format(age) if age is not None else '',
                '--status {}'.format(' '.join(status)) if status is not None else '',
                verbosity))
        except subprocess.CalledProcessError as e:
            env.logger.error('Failed to purge tasks {}'.format(tasks))
            return ''


class BackgroundProcess_TaskEngine(TaskEngine):
    def __init__(self, agent):
        super(BackgroundProcess_TaskEngine, self).__init__(agent)

    def execute_task(self, task_id):
        if not super(BackgroundProcess_TaskEngine, self).execute_task(task_id):
            env.logger.trace('Failed to prepare task {}'.format(task_id))
            return False
        cmd = "sos execute {0} -v {1} -s {2} {3}".format(
            task_id, env.verbosity, env.config['sig_mode'], '--dryrun' if env.config['run_mode'] == 'dryrun' else '')
        env.logger.trace('Execute "{}" (waiting={})'.format(cmd, self.wait_for_task))
        self.agent.run_command(cmd, wait_for_task = self.wait_for_task)
        return True


