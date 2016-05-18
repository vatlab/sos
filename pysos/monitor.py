#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
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
from .utils import env

class ProcessMonitor(threading.Thread):
    def __init__(self, pid, interval=1):
        threading.Thread.__init__(self)
        self.pid = pid
        self.interval = interval
        self.start_time = time.time()
        self.proc_file = os.path.join(env.exec_dir, '.sos/{}.proc'.format(self.pid))
                 
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
                time.sleep(self.interval)
                cpu, mem, nch, ch_cpu, ch_mem = self._check()
                with open(self.proc_file, 'a') as pd:
                    pd.write('{:.4f}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(time.time() - self.start_time, cpu, mem, nch, ch_cpu, ch_mem))
            except Exception as e:
                # if the process died, exit the thread
                #with open(self.proc_file, 'a') as pd:
                #    pd.write('Proc exited after {} seconds.'.format(time.time() - self.start_time))
                env.logger.warning(e)
                break

        
def summarizeExecution(pid):
    proc_file = os.path.join(env.exec_dir, '.sos/{}.proc'.format(pid))
    if not os.path.isfile(proc_file):
        return ''
    peak_cpu = 0
    peak_mem = 0
    peak_nch = 0
    run_time = 0
    with open(proc_file) as proc:
        for line in proc:
            t, c, m, nch, cc, cm = line.split()
            run_time = float(t)
            if float(c) > peak_cpu:
                peak_cpu = float(c) + float(cc)
            if float(m) > peak_mem:
                peak_mem = float(m) + float(cm)
            if int(nch) > peak_nch:
                peak_nch = int(nch)
    return ('Completed in {:.1f} seconds with {} children and {:.2f} % peak CPU and {:.1f} Mb peak memory usage'
        .format(run_time, peak_nch, peak_cpu * 100, peak_mem/1024/1024))
        

