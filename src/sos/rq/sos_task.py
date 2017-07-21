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
from sos.utils import env
from rq import Queue as rqQueue
from redis import Redis
from sos.sos_task import TaskEngine, execute_task, loadTask
from sos.sos_eval import cfg_interpolate

class RQ_TaskEngine(TaskEngine):

    def __init__(self, agent):
        super(RQ_TaskEngine, self).__init__(agent)
        #
        self.redis_host = cfg_interpolate(self.config.get('redis_host', self.config.get('address', 'localhost')))
        self.redis_port = self.config.get('redis_port', 6379)
        self.redis_queue = cfg_interpolate(self.config.get('queue', 'default'))

        try:
            env.logger.debug('Connecting to redis server {} at port {}'.format(self.redis_host, self.redis_port))
            redis_conn = Redis(host=self.redis_host, port=self.redis_port)
        except Exception as e:
            raise RuntimeError('Failed to connect to redis server with host {} and port {}: {}'.format(
                self.redis_server, self.redis.port, e))

        self.redis_queue = rqQueue(self.redis_queue, connection=redis_conn)
        self.jobs = {}

    def execute_task(self, task_id):
        #
        if not super(RQ_TaskEngine, self).execute_task(task_id):
            return False

        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.task')
        sos_dict = loadTask(task_file).sos_dict

        # however, these could be fixed in the job template and we do not need to have them all in the runtime
        runtime = self.config
        runtime.update({x:sos_dict['_runtime'][x] for x in ('walltime', 'cur_dir', 'home_dir', 'name') if x in sos_dict['_runtime']})
        runtime['task'] = task_id
        runtime['verbosity'] = env.verbosity
        runtime['sig_mode'] = env.config['sig_mode']
        runtime['run_mode'] = env.config['run_mode']
        if 'name' in runtime:
            runtime['job_name'] = cfg_interpolate(runtime['name'], sos_dict)
        else:
            runtime['job_name'] = cfg_interpolate('${step_name}_${_index}', sos_dict)
        if 'nodes' not in runtime:
            runtime['nodes'] = 1
        if 'cores' not in runtime:
            runtime['cores'] = 1
    
        # tell subprocess where pysos.runtime is
        self.jobs[task_id] = self.redis_queue.enqueue(
            execute_task, args=(task_id, env.verbosity, runtime['run_mode'],
                runtime['sig_mode'], 5, 60),
            job_id = runtime['job_name'],
            # result expire after one day
            result_ttl=86400,
            timeout=runtime.get('walltime', 86400*30))
        return True

