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

    def execute_workflow(self, filename, command, **template_args):
        # there is no need to prepare workflow (e.g. copy over)
        if os.path.isfile(filename):
            dest = self.agent.send_to_host(filename)
            self.filename = filename
        elif os.path.isfile(filename + '.sos'):
            dest = self.agent.send_to_host(filename + '.sos')
            self.filename = filename + '.sos'
        elif os.path.isfile(filename + '.ipynb'):
            dest = self.agent.send_to_host(filename + '.ipynb')
            self.filename = filename + '.ipynb'
        else:
            raise RuntimeError(f'Failed to locate script {script}')

        self.command = self.remove_arg(command, '-r')
        # -c only point to local config file.
        self.command = self.remove_arg(self.command, '-c')
        # remove --slave mode because the master cannot reach remote slave
        self.command = self.remove_arg(self.command, '-m')
        # replace absolute path with relative one because remote sos might have
        # a different path.
        if os.path.basename(command[0]) == 'sos':
            self.command[0] = 'sos'
            self.command[2] = dest.values()[0]
        elif os.path.basename(command[0]) == 'sos-runner':
            self.command[0] = 'sos-runner'
            self.command[1] = dest.values()[0]

        self.template_args = template_args
        return True


class BackgroundProcess_WorkflowEngine(WorkflowEngine):

    def __init__(self, agent):
        super(BackgroundProcess_WorkflowEngine, self).__init__(agent)
        if 'workflow_template' in self.config:
            self.workflow_template = self.config['workflow_template'].replace(
                '\r\n', '\n')
        else:
            self.workflow_template = None

    def execute_workflow(self, filename, command, **template_args):
        #
        # calling super execute_workflow would set cleaned versions
        # of self.filename, self.command, and self.template_args
        if not super(BackgroundProcess_WorkflowEngine, self).execute_workflow(
                filename, command, **template_args):
            env.log_to_file(
                'WORKFLOW',
                f'Failed to prepare workflow with command "{command}"')
            return False

        if self.workflow_template:
            if not self._execute_workflow_with_template():
                return False
        else:
            if not self._execute_workflow():
                return False
        return True

    def _execute_workflow(self):
        # if no template, use a default command
        env.log_to_file('WORKDLOW', f'Execute "{self.command}"')
        self.agent.run_command(self.command)
        return True

    def _execute_workflow_with_template(self):
        '''Submit workflows by interpolating a shell script defined in workflow_template'''

        try:
            if self.filename.lower().endswith('.ipynb'):
                from .converter import extract_workflow
                script = extract_workflow(self.filename)
            else:
                with open(self.filename) as script_file:
                    script = script_file.read()

            self.template_args['filename'] = self.filename
            self.template_args['script'] = script
            self.template_args['cmd'] = self.command
            job_text = cfg_interpolate(self.workflow_template,
                                       self.template_args) + '\n'

            filename = cmd[0] + ('.sh' if len(cmd) == 1 else f'-{cmd[-1]}.sh')

        except Exception as e:
            raise ValueError(
                f'Failed to generate job file for the execution of workflow {script}: {e}'
            )

        try:
            wf_dir = os.path.join(os.path.expanduser('~'), '.sos', 'workflows')
            if not os.path.isdir(wf_dir):
                os.makedirs(wf_dir)

            job_file = tempfile.TemporaryDirectory(
                dir=wf_dir, prefix='tmp_wf', suffix='.sh', delete=False).name
            # do not translate newline under windows because the script will be executed
            # under linux/mac
            with open(job_file, 'w', newline='') as job:
                job.write(job_text)

            # then copy the job file to remote host if necessary
            self.agent.send_job_file(job_file, dir='workflows')

            cmd = f'bash ~/.sos/workflows/{filename}'
            env.log_to_file('TASK', f'Execute "{cmd}" with script {job_text}')
            self.agent.run_command(cmd)
        except Exception as e:
            raise RuntimeError(f'Failed to submit workflow {cmd}: {e}')
        finally:
            try:
                os.remove(job_file)
            except Exception as e:
                env.logger.debug(
                    f'Failed to remove temporary workflow {job_file}')
        return True
