#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
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
import psutil
import threading
import time
from datetime import datetime
import stat
from .utils import env, expand_time, format_HHMMSS

class ProcessMonitor(threading.Thread):
    def __init__(self, task_id, monitor_interval, resource_monitor_interval, max_walltime=None, max_mem=None, max_procs=None):
        threading.Thread.__init__(self)
        self.task_id = task_id
        self.pid = os.getpid()
        self.monitor_interval = monitor_interval
        self.resource_monitor_interval = max(resource_monitor_interval // monitor_interval, 1)
        self.daemon = True
        self.max_walltime = max_walltime
        if self.max_walltime is not None:
            self.max_walltime = expand_time(self.max_walltime)
        self.max_mem = max_mem
        self.max_procs = max_procs
        self.pulse_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.pulse')
        # remove previous status file, which could be readonly if the job is killed
        if os.path.isfile(self.pulse_file):
            if not os.access(self.pulse_file, os.W_OK):
                os.chmod(self.pulse_file, stat.S_IREAD | stat.S_IWRITE)
            os.remove(self.pulse_file)
        with open(self.pulse_file, 'w') as pd:
            pd.write('#task: {}\n'.format(task_id))
            pd.write('#started at {}\n#\n'.format(datetime.now().strftime("%A, %d. %B %Y %I:%M%p")))
            pd.write('#time\tproc_cpu\tproc_mem\tchildren\tchildren_cpu\tchildren_mem\n')

    def _check(self):
        current_process = psutil.Process(self.pid)
        par_cpu = current_process.cpu_percent()
        par_mem = current_process.memory_info()[0]
        ch_cpu = 0
        ch_mem = 0
        children = current_process.children(recursive=True)
        n_children = len(children)
        for child in children:
            ch_cpu += child.cpu_percent()
            ch_mem += child.memory_info()[0]
        return par_cpu, par_mem, n_children, ch_cpu, ch_mem

    def _exceed_resource(self, msg):
        err_file =  os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.task_id + '.err')
        with open(err_file, 'a') as err:
            err.write(msg + '\n')
        # kill the task
        from stat import S_IREAD, S_IRGRP, S_IROTH
        os.chmod(self.pulse_file, S_IREAD|S_IRGRP|S_IROTH)
        p = psutil.Process(self.pid)
        p.kill()

    def run(self):
        counter = 0
        start_time = time.time()
        while True:
            try:
                if not os.access(self.pulse_file, os.W_OK):
                    # the job should be killed
                    p = psutil.Process(self.pid)
                    p.kill()
                # most of the time we only update 
                if counter % self.resource_monitor_interval:
                    os.utime(self.pulse_file, None)
                else:
                    cpu, mem, nch, ch_cpu, ch_mem = self._check()
                    with open(self.pulse_file, 'a') as pd:
                        pd.write('{}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(time.time(), cpu, mem, nch, ch_cpu, ch_mem))
                    if self.max_procs is not None and cpu > self.max_procs:
                        self._exceed_resource('Task {} exits because of excessive use of procs (used {}, limit {})'.format(
                            self.task_id, cpu, self.max_procs))
                    if self.max_mem is not None and mem > self.max_mem:
                        self._exceed_resource('Task {} exits because of excessive use of max_mem (used {}, limit {})'.format(
                            self.task_id, mem, self.max_mem))
                # walltime can be checked more frequently and does not have to wait for resource option
                elapsed = time.time() - start_time
                if self.max_walltime is not None and elapsed > self.max_walltime:
                    self._exceed_resource('Task {} exits because of excessive run time (used {}, limit {})'.format(
                        self.task_id, format_HHMMSS(int(elapsed)), format_HHMMSS(self.max_walltime)))
                time.sleep(self.monitor_interval)
                counter += 1
            except Exception as e:
                # if the process died, exit the thread
                # the warning message is usually:
                # WARNING: psutil.NoSuchProcess no process found with pid XXXXX
                #env.logger.warning(e)
                env.logger.debug('Monitor of {} failed with message {}'.format(self.task_id, e))
                break

def summarizeExecution(task_id, status='Unknown'):
    pulse_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.pulse')
    if not os.path.isfile(pulse_file):
        pulse_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.status')
        if not os.path.isfile(pulse_file):
            return
    peak_cpu = 0
    accu_cpu = 0
    peak_mem = 0
    accu_mem = 0
    peak_nch = 0
    start_time = None
    end_time = None
    count = 0
    with open(pulse_file) as proc:
        for line in proc:
            if line.startswith('#'):
                continue
            try:
                t, c, m, nch, cc, cm = line.split()
            except Exception as e:
                env.logger.warning('Unrecognized resource line "{}": {}'.format(line.strip(), e))
            if start_time is None:
                start_time = float(t)
                end_time = float(t)
            else:
                end_time = float(t)
            accu_cpu += float(c) + float(cc)
            accu_mem += float(m) + float(cm)
            count += 1
            if float(c) + float(cc) > peak_cpu:
                peak_cpu = float(c) + float(cc)
            if float(m) + float(cm) > peak_mem:
                peak_mem = float(m) + float(cm)
            if int(nch) > peak_nch:
                peak_nch = int(nch)
    try:
        second_elapsed = end_time - start_time
    except:
        second_elapsed = 0
    result = [
        ('status', status),
        ('task', task_id),
        ('nproc', str(peak_nch)),
        ('start', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))),
        ('end', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))),
        ('duration', ('' if second_elapsed < 86400 else '{} day{} '.format(int(second_elapsed/86400), 's' if second_elapsed > 172800 else '')) + \
                time.strftime('%H:%M:%S', time.gmtime(second_elapsed))),
        ('cpu_peak', '{:.1f}'.format(peak_cpu)),
        ('cpu_avg', '{:.1f}'.format(0 if count == 0 else accu_cpu/count)),
        ('mem_peak', '{:.1f}Mb'.format(peak_mem/1024/1024)),
        ('mem_avg', '{:.1f}Mb'.format(0 if count == 0 else accu_mem/1024/1024/count))
        ]
    return '\n'.join('{:20s} {}'.format(x,y) for x,y in result)

        

