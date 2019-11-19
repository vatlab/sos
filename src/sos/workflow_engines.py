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

    def remove_arg(self, argv, arg):
        r_idx = [idx for idx, x in enumerate(argv) if x.startswith(arg)]
        if not r_idx:
            return argv
        else:
            r_idx = r_idx[0]
        # find next option
        r_next = [
            idx for idx, x in enumerate(argv[r_idx + 1:]) if x.startswith('-')
        ]
        if r_next:
            argv = argv[:r_idx] + argv[r_idx + 1 + r_next[0]:]
        else:
            argv = argv[:r_idx]
        return argv


    def execute_workflow(self, script, cmd, **kwargs):
        # there is no need to prepare workflow (e.g. copy over)
        if os.path.isfile(script):
            dest = self.agent.send_to_host(script)
        elif os.path.isfile(script + '.sos'):
            dest = self.agent.send_to_host(script + '.sos')
        elif os.path.isfile(script + '.ipynb'):
            dest = self.agent.send_to_host(script + '.ipynb')
        else:
            raise RuntimeError(f'Failed to locate script {script}')

        # perhaps a different filename is presented
        self.script = dest.values()[0]
                
        self.cmd = self.remove_arg(cmd, '-r')
        # -c only point to local config file.
        self.cmd = self.remove_arg(self.cmd, '-c')
        # remove --slave mode because the master cannot reach remote slave
        self.cmd = self.remove_arg(self.cmd, '-m')
        # replace absolute path with relative one because remote sos might have
        # a different path.
        if os.path.basename(argv[0]) == 'sos':
            self.cmd[0] = 'sos'
            self.cmd[2] = self.script
        elif os.path.basename(argv[0]) == 'sos-runner':
            self.cmd[0] = 'sos-runner'
            self.cmd[1] = self.script

        return True


class BackgroundProcess_WorkflowEngine(WorkflowEngine):

    def __init__(self, agent):
        super(BackgroundProcess_WorkflowEngine, self).__init__(agent)
        if 'workflow_template' in self.config:
            self.workflow_template = self.config['workflow_template'].replace(
                '\r\n', '\n')
        else:
            self.workflow_template = None

    def execute_workflow(self, script, cmd, **kwargs):
        if not super(BackgroundProcess_WorkflowEngine,
                     self).execute_workflow(script, cmd, **kwargs):
            env.log_to_file('WORKFLOW',
                            f'Failed to prepare workflow with command "{cmd}"')
            return False
            
        if self.workflow_template:
            if not self._execute_workflow_with_template(self.script, self.cmd, **kwargs):
                return False
        else:
            if not self._execute_workflow(self.script, self.cmd, **kwargs):
                return False
        return True

    def _execute_workflow(self, script, cmd, **kwargs):
        # if no template, use a default command
        env.log_to_file('WORKDLOW', f'Execute "{cmd}"')
        self.agent.run_command(cmd)
        return True

    def _execute_workflow_with_template(self, script, cmd, **kwargs):
        '''Submit workflows by interpolating a shell script defined in workflow_template'''
        template_args = kwargs
        template_args['script'] = script
        template_args['cmd'] = cmd

        try:
            job_text += cfg_interpolate(self.workflow_template, template_args)
            job_text += '\n'
        except Exception as e:
            raise ValueError(
                f'Failed to generate job file for the execution of workflow {script}: {e}'
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
