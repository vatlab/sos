#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
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
from sos.utils import env
from sos.sos_executor import MP_Executor

from .sos_step import RQ_Step_Executor

class RQ_Executor(MP_Executor):
    def __init__(self, workflow, args=[], nested=False, config={}):
        MP_Executor.__init__(self, workflow, args, nested=nested, config=config)
        env.__task_engine__ = 'rq'

        from rq import Queue as rqQueue
        from redis import Redis
        import yaml

        connection_file = os.path.join(env.exec_dir, '.sos', 'redis_connection.yaml')
        if os.path.isfile(connection_file):
            try:
                with open(connection_file) as conn:
                    cfg = yaml.safe_load(conn)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to parse redis connection file {}. ({}}'.format(connection_file, e))
                sys.exit(1)
            try:
                redis_conn = Redis(host=cfg['host'], port=cfg.get('port', 6379))
            except Exception as e:
                env.logger.error('Failed to connect to redis server: {}'.format(e))
                sys.exit(1)
        else:
            redis_conn = Redis()
        self.redis_queue = rqQueue(connection=redis_conn)

    def step_executor(self, section, queue):
        return RQ_Step_Executor(section, queue, self.redis_queue)

