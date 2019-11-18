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
from collections import OrderedDict, defaultdict

from .eval import cfg_interpolate
from .utils import env, expand_time
from .workflows import TaskFile
from .messages import encode_msg


class WorkflowEngine:
    def __init__(self):
        pass

    def submit_workflow(self, workflow_id):
        pass

class BackgroundProcess_WorkflowEngine(WorkflowEngine):
    def __init__(self, agent):
        super(BackgroundProcess_WorkflowEngine, self).__init__(agent)
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

    def execute_workflow(self, workflow_ids):
        if not super(BackgroundProcess_WorkflowEngine,
                     self).execute_workflows(workflow_ids):
            env.log_to_file('TASK', f'Failed to prepare workflow {workflow_ids}')
            return False
        if self.job_template:
            if not self._submit_workflow_with_template(workflow_ids):
                return False
        else:
            if not self._submit_workflow(workflow_ids):
                return False
        return True

    def _submit_workflow(self, workflow_ids):
        # if no template, use a default command
        cmd = f"sos execute {' '.join(workflow_ids)} -v {env.verbosity} -s {env.config['sig_mode']} -m {env.config['run_mode']}"
        env.log_to_file('TASK',
                        f'Execute "{cmd}" (waiting={self.wait_for_workflow})')
        self.agent.run_command(cmd, wait_for_workflow=self.wait_for_workflow)
        return True

    def _submit_workflow_with_template(self, workflow_ids):
        '''Submit workflows by interpolating a shell script defined in job_template'''
        runtime = self.config
        runtime.update({
            'workdir': os.getcwd(),
            'cur_dir': os.getcwd(),  # for backward compatibility
            'verbosity': env.verbosity,
            'sig_mode': env.config.get('sig_mode', 'default'),
            'run_mode': env.config.get('run_mode', 'run'),
            'home_dir': os.path.expanduser('~')
        })
        if '_runtime' in env.sos_dict:
            runtime.update({
                x: env.sos_dict['_runtime'][x]
                for x in ('nodes', 'cores', 'workdir', 'mem', 'walltime')
                if x in env.sos_dict['_runtime']
            })
        if 'nodes' not in runtime:
            runtime['nodes'] = 1
        if 'cores' not in runtime:
            runtime['cores'] = 1

        # let us first prepare a workflow file
        job_text = ''
        for workflow_id in workflow_ids:
            runtime['workflow'] = workflow_id
            try:
                job_text += cfg_interpolate(self.job_template, runtime)
                job_text += '\n'
            except Exception as e:
                raise ValueError(
                    f'Failed to generate job file for workflow {workflow_id}: {e}')

        filename = workflow_ids[0] + ('.sh' if len(workflow_ids) == 1 else
                                  f'-{workflow_ids[-1]}.sh')
        # now we need to write a job file
        job_file = os.path.join(
            os.path.expanduser('~'), '.sos', 'workflows', filename)
        # do not translate newline under windows because the script will be executed
        # under linux/mac
        with open(job_file, 'w', newline='') as job:
            job.write(job_text)

        # then copy the job file to remote host if necessary
        self.agent.send_workflow_file(job_file)

        try:
            cmd = f'bash ~/.sos/workflows/{filename}'
            env.log_to_file('TASK', f'Execute "{cmd}" with script {job_text}')
            self.agent.run_command(cmd, wait_for_workflow=self.wait_for_workflow)
        except Exception as e:
            raise RuntimeError(f'Failed to submit workflow {workflow_ids}: {e}')
        return True
