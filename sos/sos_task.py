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
from collections.abc import Sequence
import concurrent.futures

from sos.utils import env, short_repr, get_traceback
from sos.sos_eval import SoS_exec

from .target import textMD5, RuntimeInfo, Undetermined, FileTarget
from .monitor import ProcessMonitor

from collections import OrderedDict


monitor_interval = 5
resource_monitor_interval = 60

class TaskParams(object):
    '''A parameter object that encaptulates parameters sending to
    task executors. This would makes the output of workers, especially
    in the web interface much cleaner (issue #259)'''
    def __init__(self, name, data):
        self.name = name
        self.data = data

    def __repr__(self):
        return self.name

def execute_task(task_id, verbosity=None, runmode='run', sigmode=None, monitor_interval=5,
    resource_monitor_interval=60):
    '''A function that execute specified task within a local dictionary
    (from SoS env.sos_dict). This function should be self-contained in that
    it can be handled by a task manager, be executed locally in a separate
    process or remotely on a different machine.'''
    # start a monitoring file, which would be killed after the job
    # is done (killed etc)
    m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
        resource_monitor_interval=resource_monitor_interval)
    m.start()

    task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
    with open(task_file, 'rb') as task:
        params = pickle.load(task)

    task, sos_dict, sigil = params.data
    if verbosity is not None:
        env.verbosity = verbosity
    if sigmode is not None:
        env.config['sig_mode'] = sigmode
    env.config['run_mode'] = runmode
    env.register_process(os.getpid(), 'spawned_job with {} {}'
        .format(sos_dict['_input'], sos_dict['_output']))

    env.logger.info('{} ``started``'.format(task_id))
    env.sos_dict.quick_update(sos_dict)

    skipped = False
    if env.config['sig_mode'] == 'ignore' or env.sos_dict['_output'] is None:
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
        return {'succ': 0, 'output': env.sos_dict['_output'], 'path': os.environ['PATH']}

    try:
        # go to 'cur_dir'
        orig_dir = os.getcwd()
        if '_runtime' in sos_dict and 'cur_dir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['cur_dir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
                except Exception as e:
                    raise RuntimeError('Failed to create cur_dir {}'.format(sos_dict['_runtime']['cur_dir']))
            os.chdir(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
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

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        # step process
        SoS_exec(task, sigil)
        os.chdir(orig_dir)
    except Exception as e:
        if env.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error('{} ``failed`` with {} error {}'.format(task_id, e.__class__.__name__, e))
        return {'succ': 1, 'exception': e, 'path': os.environ['PATH']}
    except KeyboardInterrupt:
        env.logger.error('{} ``interrupted``'.format(task_id))
        raise
    finally:
        env.sos_dict.set('__step_sig__', None)

    if sig:
        sig.write(env.sos_dict['_local_input_{}'.format(env.sos_dict['_index'])],
            env.sos_dict['_local_output_{}'.format(env.sos_dict['_index'])])
        sig.release()
    env.deregister_process(os.getpid())
    env.logger.info('{} ``completed``'.format(task_id))
    return {'succ': 0, 'output': {} if env.sos_dict['_output'] is None else {x:FileTarget(x).signature() for x in env.sos_dict['_output'] if isinstance(x, str)}}


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
            if res['succ'] == 0:
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
            time.sleep(1)
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

def check_tasks(tasks, verbosity=1):
    # verbose is ignored for now
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    if not tasks:
        tasks = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getctime(x)) for x in tasks]
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', '{}*.task'.format(t)))
            matched = [(os.path.basename(x)[:-5], os.path.getctime(x)) for x in matched]
            if not matched:
                all_tasks.append((t, None))
            else:
                all_tasks.extend(matched)
    all_tasks = sorted(list(set(all_tasks)), key=lambda x: x[1])
    if not all_tasks:
        env.logger.warning('No matching tasks')
        return
    # at most 20 threads
    p = Pool(min(20, len(all_tasks)))
    status = p.map(check_task, [x[0] for x in all_tasks])
    if verbosity == 0:
        print('\n'.join(status))
    elif verbosity == 1:
        for s, (t, d) in zip(status, all_tasks):
            print('{}\t{}'.format(t, s))
    elif verbosity == 2:
        from .utils import PrettyRelativeTime
        for s, (t, d) in zip(status, all_tasks):
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
            with open(task_file, 'rb') as task:
                params = pickle.load(task)
            print('TASK:\n=====')
            print(params.data[0])
            print()
            print('ENVIRONMENT:\n============')
            job_vars = params.data[1]
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
                    '.task', '.pulse', '.status')]):
                    print('{}:\n{}'.format(os.path.basename(f), '='*(len(os.path.basename(f))+1)))
                    with open(f) as fc:
                        print(fc.read())

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
    elif status != 'running':
        return status
    # job is running
    pulse_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.pulse')
    from stat import S_IREAD, S_IRGRP, S_IROTH
    os.chmod(pulse_file, S_IREAD|S_IRGRP|S_IROTH)
    # remove job file as well
    job_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.sh')
    if os.path.isfile(job_file):
        try:
            os.remove(job_file)
        except:
            pass
    return 'killed'


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

        self.tasks = []
        self.pending_tasks = []
        self.submitting_tasks = {}

        self.task_status = {}
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
        if env.config['wait_for_task'] is not None:
            self.wait_for_task = env.config['wait_for_task']
        elif 'wait_for_task' in self.config:
            self.wait_for_task = self.config['wait_for_task']
        else:
            # default
            self.wait_for_task = True

    def reset(self):
        with threading.Lock():
            self.tasks = []
            self.pending_tasks = []
            self.submitting_tasks = {}
            self.task_status = {}

    def get_tasks(self):
        with threading.Lock():
            pending = copy.deepcopy(self.pending_tasks)
            running = copy.deepcopy(self.tasks)
        return pending, running

    def run(self):
        # get all system tasks that might have been running ...
        # this will be run only once when the task engine starts
        status_output = self.query_tasks([], verbosity=1)
        with threading.Lock():
            for line in status_output.split('\n'):
                if not line.strip():
                    continue
                try:
                    tid, tst = line.split('\t')
                    self.task_status[tid] = tst
                except Exception as e:
                    env.logger.warning('Unrecognized response "{}" ({}): {}'.format(line, e.__class__.__name__, e))
        self.engine_ready.set()
        while True:
            # if no new task, does not do anything.
            env.logger.trace('processing {}'.format(self.tasks))
            if self.tasks:
                status_output = self.query_tasks(self.tasks, verbosity=1)
                with threading.Lock():
                    for line in status_output.split('\n'):
                        if not line.strip():
                            continue
                        try:
                            tid, tst = line.split('\t')
                            env.logger.trace('STATUS {}\t{}\n'.format(tid, tst))
                            if hasattr(env, '__task_notifier__') and tst != 'non-exist':
                                if tid in self.task_status and self.task_status[tid] == tst:
                                    env.__task_notifier__(['pulse-status', tid, tst])
                                else:
                                    env.__task_notifier__(['change-status', tid, tst])
                            self.task_status[tid] = tst
                            # terminal states, remove tasks from task list
                            if tst in ('completed', 'failed', 'aborted', 'result-mismatch'):
                                self.tasks.remove(tid)
                        except Exception as e:
                            env.logger.warning('Unrecognized response "{}" ({}): {}'.format(line, e.__class__.__name__, e))
                    self.summarize_status()

            env.logger.trace('submitting {}'.format(self.submitting_tasks))
            if self.submitting_tasks:
                with threading.Lock():
                    submitted = []
                    for k in self.submitting_tasks:
                        if not self.submitting_tasks[k].running():
                            submitted.append(k)
                            if self.submitting_tasks[k].result():
                                self.tasks.append(k)
                            else:
                                if hasattr(env, '__task_notifier__'):
                                    env.__task_notifier__(['change-status', k, 'failed'])
                                self.task_status[k] = 'failed'
                    for k in submitted:
                        self.submitting_tasks.pop(k)

            env.logger.trace('pending {}'.format(self.pending_tasks))
            if self.pending_tasks:
                to_run = []
                # check status
                num_active_tasks = len([x for x in self.tasks if self.task_status[x] == 'running'])
                if num_active_tasks < self.max_running_jobs:
                    to_run = self.pending_tasks[ : self.max_running_jobs - num_active_tasks]
                env.logger.debug('{} pending {} running'.format(len(self.pending_tasks), num_active_tasks))
                for tid in to_run:
                    if self.task_status[tid] == 'running':
                        env.logger.info('{} ``runnng``'.format(tid))
                    else:
                        #env.logger.info('{} ``submitttt``'.format(tid))
                        #t = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                        #self.submitting_tasks[tid] = t.submit(self.execute_task, tid)
                        self.execute_task(tid)
                #
                with threading.Lock():
                    for tid in to_run:
                        self.pending_tasks.remove(tid)

            time.sleep(self.status_check_interval)

    def submit_task(self, task_id):
        # we wait for the engine to start
        self.engine_ready.wait()

        # submit tasks simply add task_id to pending task list
        with threading.Lock():
            # if already in
            #if task_id in self.tasks or task_id in self.pending_tasks:
            #    env.logger.info('{} ``{}``'.format(task_id, self.task_status[task_id]))
            #    if hasattr(env, '__task_notifier__'):
            #        env.__task_notifier__(['new-status', task_id, self.task_status[task_id]])
            #    return self.task_status[task_id]
            #
            if task_id in self.task_status and self.task_status[task_id]:
                if self.task_status[task_id] == 'running':
                    env.logger.info('{} ``already runnng``'.format(task_id))
                    if hasattr(env, '__task_notifier__'):
                        env.__task_notifier__(['new-status', task_id, 'running'])
                    return 'running'
                # there is a case when the job is already completed (complete-old), but
                # because we do not know if the user asks to rerun (-s force), we have to
                # resubmit the job. In the case of not-rerun, the task would be marked
                # completed very soon.
                elif self.task_status[task_id] == 'completed':
                    if env.config['sig_mode'] != 'force' or task_id in env.config['resumed_tasks']:
                        # force re-execution, but it is possible that this task has been
                        # executed but quit in no-wait mode (or canceled by user). More
                        # importantly, the Jupyter notebook would re-run complted workflow
                        # even if it has "-s force" signature.
                        env.logger.info('{} ``already completed``'.format(task_id))
                        #if hasattr(env, '__task_notifier__'):
                        #    env.__task_notifier__(['new-status', task_id, 'completed'])
                        return 'completed'
                    else:
                        env.logger.info('{} ``re-execute completed``'.format(task_id))
                else:
                    env.logger.info('{} ``restart`` from status ``{}``'.format(task_id, self.task_status[task_id]))

            env.logger.info('{} ``queued``'.format(task_id))
            self.pending_tasks.append(task_id)
            self.task_status[task_id] = 'pending'
            if hasattr(env, '__task_notifier__'):
                env.__task_notifier__(['new-status', task_id, 'pending'])
            return 'pending'

    def summarize_status(self):
        from collections import Counter
        statuses = Counter(self.task_status.values())
        env.logger.debug(
            ' '.join('{}: {}'.format(x, y) for x, y in statuses.items()))

    def check_task_status(self, task_id):
        try:
            with threading.Lock():
                return self.task_status[task_id]
        except:
            # job not yet submitted
            return 'pending'

    def remove_tasks(self, tasks):
        with threading.Lock():
            for task in tasks:
                if hasattr(env, '__task_notifier__'):
                    env.__task_notifier__(['remove-task', task])
                #if task in self.task_status:
                #    self.task_status.pop(task)
                #if task in self.tasks:
                #    self.tasks.remove(task)

    def pending_tasks(self):
        with threading.Lock():
            return self.pending_tasks()

    def query_tasks(self, tasks=None, verbosity=1):
        return self.agent.check_output("sos status {} -v {}".format(
                ' '.join(tasks), verbosity))

    def kill_tasks(self, tasks, all_tasks=False):
        return self.agent.check_output("sos kill {} {}".format(
            ' '.join(tasks), '-a' if all_tasks else ''))

    def execute_task(self, task_id):
        # this is base class
        return self.agent.prepare_task(task_id)


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


