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
import time

from sos.utils import env
from sos.sos_step import SP_Step_Executor, TaskParams, execute_task

class RQ_Step_Executor(SP_Step_Executor):
    #
    # This is not working yet
    def __init__(self, step, queue, redis_queue):
        SP_Step_Executor.__init__(self, step, queue)
        self.redis_queue = redis_queue

    def submit_task(self, signature):
        if 'walltime' in env.sos_dict['_runtime']:
            walltime = env.sos_dict['_runtime']
            if isinstance(walltime, str):
                if walltime.count(':') > 2:
                    raise ValueError('Incorrect format.')
                try:
                    walltime = sum([int(val)*60**idx  for idx, val in enumerate(walltime.split(':')[-1::-1])])
                except Exception:
                    raise ValueError('Unacceptable walltime {} (can be "HH:MM:SS" or a number (seconds))'.format(walltime))
        else:
            # bioinformatics can be running for long time...
            # let me assume a longest running time of 1 month
            walltime = 60*60*24*30 
        #
        # tell subprocess where pysos.runtime is
        param = TaskParams(
            name = '{} (index={})'.format(self.step.step_name(), env.sos_dict['_index']),
            data = (
                self.step.task,          # task
                self.step.global_def,    # global process
                self.step.global_sigil,
                # if pool, it must not be in prepare mode and have
                # __signature_vars__
                env.sos_dict.clone_selected_vars(env.sos_dict['__signature_vars__'] \
                    | {'_input', '_output', '_depends', 'input', 'output',
                        'depends', '_index', '__args__', 'step_name', '_runtime',
                        '__workflow_sig__', '__report_output__',
                        '_local_input_{}'.format(env.sos_dict['_index']),
                        '_local_output_{}'.format(env.sos_dict['_index'])
                        }),
                signature,
                self.step.sigil
            ))

        self.proc_results.append(
            self.redis_queue.enqueue(
            execute_task,            # function
            args=(param,),
            timeout=walltime))

    def wait_for_results(self):
        while True:
            # wait for results
            try:
                if any(x.result is None for x in self.proc_results):
                    time.sleep(1)
                else:
                    # the executor interface expect results ...
                    self.proc_results = [x.result for x in self.proc_results]
                    return
            except KeyboardInterrupt:
                # if keyboard interrupt
                raise RuntimeError('KeyboardInterrupt fro m {} (master)'.format(os.getpid()))
            except Exception as e:
                # if keyboard interrupt etc
                env.logger.error('Caught {}'.format(e))
                raise

