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
import pickle
from sos.utils import env
from rq import Queue as rqQueue
from redis import Redis
from sos.sos_task import TaskEngine, execute_task

class RQ_TaskEngine(TaskEngine):
    def __init__(self, agent):
        super(RQ_TaskEngine, self).__init__(agent)
        # we have self.config for configurations
        #
        # redis_host
        # redis_port
        #
        self.redis_host = self.config['redis_host'] if 'redis_host' in self.config else 'localhost'
        self.redis_port = self.config['redis_port'] if 'redis_port' in self.config else 6379

        try:
            redis_conn = Redis(host=self.redis_host, port=self.redis_port)
        except Exception as e:
            env.logger.error('Failed to connect to redis server with host {} and port {}: {}'.format(
                self.redis_server, self.redis.port, e))

        self.redis_queue = rqQueue(connection=redis_conn)

    def execute_task(self, task_id):
        # read the task file and look for runtime info
        # 
        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.task')
        with open(task_file, 'rb') as task:
            params = pickle.load(task)
            task, sos_dict, sigil = params.data
        # bioinformatics can be running for long time...
        # let me assume a longest running time of 1 month
        walltime = sos_dict['_runtime']['walltime'] if 'walltime' in sos_dict['_runtime'] else 60*60*24*30 

        if isinstance(walltime, str):
            if walltime.count(':') > 2:
                raise ValueError('Incorrect format.')
            try:
                walltime = sum([int(val)*60**idx  for idx, val in enumerate(walltime.split(':')[-1::-1])])
            except Exception:
                raise ValueError('Unacceptable walltime {} (can be "HH:MM:SS" or a number (seconds))'.format(walltime))

        # tell subprocess where pysos.runtime is
        self.proc_results.append(
            self.redis_queue.enqueue(
            execute_task,
            args=(task_id, env.verbosity, env.sig_mode),
            timeout=walltime))


