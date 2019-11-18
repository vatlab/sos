#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import os

from .eval import cfg_interpolate
from .utils import env
from .messages import encode_msg


class WorkflowEngine:

    def __init__(self, agent):
        self.agent = agent
        self.config = agent.config
        self.alias = self.config['alias']

    def submit_workflow(self, cmd, **kwargs):
        pass


class BackgroundProcess_WorkflowEngine(WorkflowEngine):

    def __init__(self, agent):
        super(BackgroundProcess_WorkflowEngine, self).__init__(agent)
        if 'job_template' in self.config:
            self.job_template = self.config['job_template'].replace(
                '\r\n', '\n')
        else:
            self.job_template = None

    def submit_workflow(self, cmd, **kwargs):
        if not super(BackgroundProcess_WorkflowEngine,
                     self).submit_workflow(cmd, **kwargs):
            env.log_to_file('WORKFLOW',
                            f'Failed to prepare workflow with command "{cmd}"')
            return False
        if self.job_template:
            if not self._submit_workflow_with_template(cmd, **kwargs):
                return False
        else:
            if not self._submit_workflow(cmd, **kwargs):
                return False
        return True

    def _submit_workflow(self, cmd, **kwargs):
        # if no template, use a default command
        env.log_to_file('WORKDLOW', f'Execute "{cmd}"')
        self.agent.run_command(cmd)
        return True

    def _submit_workflow_with_template(self, cmd, **kwargs):
        '''Submit workflows by interpolating a shell script defined in job_template'''
        runtime = self.config
        runtime.update({'workdir': os.getcwd(), 'cmd': cmd})
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
        for workflow_id in cmd:
            runtime['workflow'] = workflow_id
            try:
                job_text += cfg_interpolate(self.job_template, runtime)
                job_text += '\n'
            except Exception as e:
                raise ValueError(
                    f'Failed to generate job file for workflow {workflow_id}: {e}'
                )

        filename = cmd[0] + ('.sh' if len(cmd) == 1 else f'-{cmd[-1]}.sh')
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
            self.agent.run_command(cmd)
        except Exception as e:
            raise RuntimeError(f'Failed to submit workflow {cmd}: {e}')
        return True
