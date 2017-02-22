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
from .utils import env

class ProcessMonitor(threading.Thread):
    def __init__(self, task_id, interval=10):
        threading.Thread.__init__(self)
        self.pid = os.getpid()
        self.interval = interval
        self.status_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.status')
        with open(self.status_file, 'w') as pd:
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

    def run(self):
        while True:
            try:
                cpu, mem, nch, ch_cpu, ch_mem = self._check()
                time.sleep(self.interval)
                with open(self.status_file, 'a') as pd:
                    pd.write('{}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(time.time(), cpu, mem, nch, ch_cpu, ch_mem))
            except Exception:
                # if the process died, exit the thread
                # the warning message is usually:
                # WARNING: psutil.NoSuchProcess no process found with pid XXXXX
                #env.logger.warning(e)
                break

def summarizeExecution(task_id, status='Unknown'):
    status_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.status')
    if not os.path.isfile(status_file):
        return
    peak_cpu = 0
    accu_cpu = 0
    peak_mem = 0
    accu_mem = 0
    peak_nch = 0
    start_time = None
    end_time = None
    count = 0
    with open(status_file) as proc:
        for line in proc:
            if line.startswith('#'):
                continue
            try:
                t, c, m, nch, cc, cm = line.split()
            except Exception as e:
                env.logger.warning('Unrecognized resource line "{}": {}'.format(line.strip(), e))
            if start_time is None:
                start_time = float(t)
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
    second_elapsed = end_time - start_time
    result = [
        ('status', status),
        ('task', task_id),
        ('nproc', str(peak_nch)),
        ('start', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))),
        ('end', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))),
        ('duration', ('' if second_elapsed < 86400 else '{} day{} '.format(int(second_elapsed/86400), 's' if second_elapsed > 172800 else '')) + \
                time.strftime('%H:%M:%S', time.gmtime(second_elapsed))),
        ('cpu_peak', '{:.1f}'.format(peak_cpu)),
        ('cpu_avg', '{:.1f}'.format(accu_cpu/count)),
        ('mem_peak', '{:.1f}Mb'.format(peak_mem/1024/1024)),
        ('mem_avg', '{:.1f}Mb'.format(accu_mem/1024/1024/count))
        ]
    max_width = [max(len(x) for x in col) for col in result]
    return ' '.join(s.ljust(l) for s,l in zip([x[0] for x in result], max_width)) + '\n' + \
        ' '.join(s.ljust(l) for s,l in zip([x[1] for x in result], max_width))

        

