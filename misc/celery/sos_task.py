#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import pickle
from sos.utils import env
from sos.tasks import TaskEngine, execute_task
from celery import Celery
from .app import app

class Celery_TaskEngine(TaskEngine):
    def __init__(self, agent):
        super(Celery_TaskEngine, self).__init__(agent)
        # we have self.config for configurations
        #
        # broker
        # backend
        #
        self.broker = self.config['broker'] if 'broker' in self.config else 'redis://localhost'
        self.backend = self.config['backend'] if 'backend' in self.config else 'redis://localhost'

        try:
            self.app = Celery('sos.celery.celery',
                broker = self.broker,
                backend  = self.backend)
        except Exception as e:
            env.logger.error('Failed to start to Celery app broker {} and backend {}: {}'.format(
                self.broker, self.backend, e))

        @self.app.task
        def celery_execute_task(task, global_def, sos_dict):
            return execute_task(task, global_def, sos_dict)

        # Optional configuration, see the application user guide.
        self.app.conf.update(
            CELERY_TASK_RESULT_EXPIRES=3600,
        )


    def execute_task(self, task_id):
        # read the task file and look for runtime info
        #
        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.task')
        with open(task_file, 'rb') as task:
            params = pickle.load(task)
            task, sos_dict  = params.data
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
            self.celery_execute_task.apply_async(
                (task, env.verbosity, env.sig_mode)
            ))


