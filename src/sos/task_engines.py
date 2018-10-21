#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import concurrent.futures
import copy
import os
import random
import subprocess
import threading
import time
from collections import OrderedDict

from .eval import cfg_interpolate
from .utils import env, expand_time


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
        self._thread_workers = concurrent.futures.ThreadPoolExecutor(
            max_workers=1)
        self._status_checker = None
        #
        if 'wait_for_task' in self.config:
            self.wait_for_task = self.config['wait_for_task']
        else:
            # default
            self.wait_for_task = True
        #
        # how to stack jobs ... currently backgroun process queue
        # allows stacking of up to 1000 tasks, but PBS queue does not
        # allow stacking.
        self.batch_size = 1

    def notify_controller(self, msg):
        if env.config['exec_mode']:
            # set cell_id to slave_id so that the frontend knows which
            # cell this task belong to
            msg['cell_id'] = env.config.get('slave_id', '')
            env.tapping_listener_socket.send_pyobj({
                'msg_type': 'task_status',
                'data': msg})

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
                if self.task_status[task] in ('submitted', 'running') and task not in self.running_tasks:
                    # these tasks will be actively monitored
                    self.running_tasks.append(task)
        #
        if age is not None:
            age = expand_time(age, default_unit='d')
        return sorted([(x, self.task_status[x], self.task_date.get(x, (time.time(), None, None))) for x in tasks
                       if (status is None or self.task_status[x] in status) and
                       (age is None or
                        ((age > 0 and time.time() - self.task_date.get(x, (time.time(), None, None))[0] > age)
                         or (age < 0 and time.time() - self.task_date.get(x, (time.time(), None, None))[0] < -age)))],
                      key=lambda x: -x[2][0])

    def get_tasks(self):
        with threading.Lock():
            pending = copy.deepcopy(
                self.pending_tasks + list(self.submitting_tasks.keys()))
            running = copy.deepcopy(self.running_tasks)
        return pending, running

    def run(self):
        # get all system tasks that might have been running ...
        # this will be run only once when the task engine starts
        status_output = self.query_tasks(
            check_all=True, verbosity=3, numeric_times=True)
        with threading.Lock():
            for line in status_output.split('\n'):
                if not line.strip():
                    continue
                try:
                    # return creation time, start time, and duration
                    tid, _, ct, st, dr, tst = line.split('\t')
                    # for some reason on windows there can be a \r at the end
                    self.task_status[tid] = tst.strip()
                    self.task_date[tid] = [float(ct), float(
                        st) if st.strip() else 0, float(dr) if dr.strip() else 0]
                except Exception as e:
                    env.logger.warning(
                        f'Unrecognized response "{line}" ({e.__class__.__name__}): {e}')
        self._last_status_check = time.time()
        self.engine_ready.set()
        while True:
            # if no new task, does not do anything.
            if self.running_tasks and time.time() - self._last_status_check > self.status_check_interval:
                if self._status_checker is None:
                    self._status_checker = self._thread_workers.submit(
                        self.query_tasks, self.running_tasks, check_all=False,
                        verbosity=3, numeric_times=True)
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
                        tid, _, ct, st, dr, tst = line.split('\t')
                        if tid not in self.running_tasks:
                            env.logger.trace(
                                f'Task {tid} removed since status check.')
                            continue
                        self.task_date[tid] = [float(ct), float(
                            st) if st.strip() else 0, float(dr) if dr.strip() else 0]
                        self.update_task_status(tid, tst)
                    except Exception as e:
                        env.logger.warning(
                            f'Unrecognized response "{line}" ({e.__class__.__name__}): {e}')
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
                                        self.notify_controller({
                                            'queue': self.agent.alias,
                                            'task_id': tid,
                                            'status': 'aborted'
                                        })
                                    else:
                                        self.running_tasks.append(tid)
                                        self.notify_controller(
                                            {
                                                'queue': self.agent.alias,
                                                'task_id': tid,
                                                'status': 'submitted'
                                            })
                            else:
                                for tid in k:
                                    self.notify_controller(
                                        {
                                            'queue': self.agent.alias,
                                            'task_id': tid,
                                            'status': 'failed'
                                        })
                                    self.task_status[tid] = 'failed'
                        # else:
                        #    env.logger.trace('{} is still being submitted.'.format(k))
                    for k in submitted:
                        self.submitting_tasks.pop(k)

            if self.pending_tasks:
                # check status
                num_active_tasks = len(
                    self.submitting_tasks) + len(self.running_tasks)
                if num_active_tasks >= self.max_running_jobs:
                    continue

                # assign tasks to self.max_running_jobs workers
                slots = [[] for i in range(self.max_running_jobs)]
                sample_slots = list(range(self.max_running_jobs))
                random.shuffle(sample_slots)
                for i, tid in enumerate(self.pending_tasks[:self.batch_size * self.max_running_jobs]):
                    if self.task_status[tid] == 'running':
                        env.logger.info(f'{tid} ``runnng``')
                    elif tid in self.canceled_tasks:
                        # the job is canceled while being prepared to run
                        env.logger.info(f'{tid} ``canceled``')
                    else:
                        # randomly spread to tasks, but at most one.
                        slots[sample_slots[i %
                                           self.max_running_jobs]].append(tid)
                for slot in slots:
                    if not slot:
                        continue
                    for tid in slot:
                        env.logger.trace(
                            f'Start submitting {tid} (status: {self.task_status.get(tid, "unknown")})')
                    self.submitting_tasks[tuple(slot)] = self._thread_workers.submit(
                        self.execute_tasks, slot)
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
            # if task_id in self.running_tasks or task_id in self.pending_tasks:
            #    self.notify_controller('{} ``{}``'.format(task_id, self.task_status[task_id]))
            #    self.notify_controller(['new-status', task_id, self.task_status[task_id]])
            #    return self.task_status[task_id]
            #
            if task_id in self.task_status and self.task_status[task_id]:
                if self.task_status[task_id] == 'running':
                    self.running_tasks.append(task_id)
                    env.logger.info(f'{task_id} ``already runnng``')
                    self.notify_controller({
                        'queue': self.agent.alias,
                        'task_id': task_id,
                        'status': 'running'
                    })
                    return 'running'
                # there is a case when the job is already completed (complete-old), but
                # because we do not know if the user asks to rerun (-s force), we have to
                # resubmit the job. In the case of not-rerun, the task would be marked
                # completed very soon.
                elif self.task_status[task_id] == 'completed':
                    if task_id in env.config.get('resumed_tasks', []):
                        # force re-execution, but it is possible that this task has been
                        # executed but quit in no-wait mode (or canceled by user). More
                        # importantly, the Jupyter notebook would re-run complted workflow
                        # even if it has "-s force" signature.
                        # self.notify_controller(['new-status', task_id, 'completed'])
                        env.logger.info(f'{task_id} ``resume with completed``')
                        return 'completed'
                    else:
                        env.logger.info(f'{task_id} ``re-execute completed``')
                elif self.task_status[task_id] != 'new':
                    env.logger.info(
                        f'{task_id} ``restart`` from status ``{self.task_status[task_id]}``')

            # self.notify_controller('{} ``queued``'.format(task_id))
            self.pending_tasks.append(task_id)
            if task_id in self.canceled_tasks:
                self.canceled_tasks.remove(task_id)
            self.task_status[task_id] = 'pending'
            self.notify_controller({
                'queue': self.agent.alias,
                'task_id': task_id,
                'status': 'pending'
            })
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
                env.logger.debug(
                    f'Task {task_id} is still not killed (status {status})')
                status = 'aborted'
            if status != 'missing':
                if task_id not in self.task_date:
                    self.task_date[task_id] = [None, None, None]
                if task_id in self.task_status and self.task_status[task_id] == status:
                    self.notify_controller(
                        {
                            'queue': self.agent.alias,
                            'task_id': task_id,
                            'status': status,
                            'start_time': self.task_date[task_id][1]
                        })
                else:
                    if status == 'running':
                        if task_id not in self.task_date:
                            self.task_date[task_id] = [
                                time.time(), time.time(), 0]
                        elif not self.task_date[task_id][1]:
                            self.task_date[task_id][1] = time.time()
                    self.notify_controller(
                        {
                            'queue': self.agent.alias,
                            'task_id': task_id,
                            'status': status,
                            'start_time': self.task_date[task_id][1]
                        })
            self.task_status[task_id] = status
            if status == 'pening' and task_id not in self.pending_tasks:
                self.pending_tasks.append(task_id)
            if status == 'running' and task_id not in self.running_tasks:
                self.running_tasks.append(task_id)
            # terminal states, remove tasks from task list
            if status in ('completed', 'failed', 'aborted') and task_id in self.running_tasks:
                self.running_tasks.remove(task_id)

    def query_tasks(self, tasks=None, check_all=False, verbosity=1, html=False, numeric_times=False, age=None, tags=None, status=None):
        try:
            return self.agent.check_output("sos status {} -v {} {} {} {} {} {} {}".format(
                '' if tasks is None else ' '.join(tasks), verbosity,
                '--all' if check_all else '',
                '--html' if html else '',
                '--numeric-times' if numeric_times else '',
                f'--age {age}' if age else '',
                f'--tags {" ".join(tags)}' if tags else '',
                f'--status {" ".join(status)}' if status else '',
            ))
        except subprocess.CalledProcessError as e:
            if verbosity >= 3:
                env.logger.warning(
                    f'Failed to query status of tasks on {self.alias}: {e.output}')
            return ''

    def kill_tasks(self, tasks, tags=None, all_tasks=False):
        # we wait for the engine to start
        self.engine_ready.wait()
        if all_tasks:
            tasks = self.pending_tasks + \
                list(self.submitting_tasks.keys()) + self.running_tasks

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
                    # the job will be killed from command line with
                    # status updated here
                    env.logger.debug(f'Killing running task {task}')
                else:
                    # it is not in the system, so we need to know what the
                    # status of the task before we do anything...
                    pass

        self.canceled_tasks.extend(tasks)
        #
        # verbosity cannot be send to underlying command because task engines
        # rely on the output of certain verbosity (-v1) to post kill the jobs
        cmd = "sos kill {} {} {}".format('' if all_tasks else ' '.join(tasks),
                                         f'--tags {" ".join(tags)}' if tags else '',
                                         '-a' if all_tasks else '')

        try:
            ret = self.agent.check_output(cmd)
            env.logger.debug(f'"{cmd}" executed with response "{ret}"')
        except subprocess.CalledProcessError:
            env.logger.error(
                f'Failed to kill all tasks' if all_tasks else f'Failed to kill tasks {" ".join(tasks)}')
            return ''
        return ret

    def resume_task(self, task):
        # we wait for the engine to start
        self.engine_ready.wait()
        with threading.Lock():
            # it is possible that a task is aborted from an opened notebook with aborted status
            if task not in self.task_status or \
                    self.task_status[task] not in ('completed', 'failed', 'aborted'):
                env.logger.warning(
                    f'Resume task called for non-canceled or non-completed/failed task {task}')
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
            self.job_template = self.config['job_template'].replace(
                '\r\n', '\n')
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
        self.agent.run_command(cmd, wait_for_task=self.wait_for_task)
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
            runtime.update({x: env.sos_dict['_runtime'][x] for x in (
                'nodes', 'cores', 'mem', 'walltime') if x in env.sos_dict['_runtime']})
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
                raise ValueError(
                    f'Failed to generate job file for task {task_id}: {e}')

        filename = task_ids[0] + ('.sh' if len(task_ids)
                                  == 1 else f'-{task_ids[-1]}.sh')
        # now we need to write a job file
        job_file = os.path.join(os.path.expanduser(
            '~'), '.sos', 'tasks', filename)
        # do not translate newline under windows because the script will be executed
        # under linux/mac
        with open(job_file, 'w', newline='') as job:
            job.write(job_text)

        # then copy the job file to remote host if necessary
        self.agent.send_task_file(job_file)

        try:
            cmd = f'bash ~/.sos/tasks/{filename}'
            self.agent.run_command(cmd, wait_for_task=self.wait_for_task)
        except Exception as e:
            raise RuntimeError(f'Failed to submit task {task_ids}: {e}')
        return True
