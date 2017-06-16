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
from sos.sos_eval import interpolate
from sos.sos_task import TaskEngine, loadTask
from sos.pattern import extract_pattern

class PBS_TaskEngine(TaskEngine):
    def __init__(self, agent):
        super(PBS_TaskEngine, self).__init__(agent)
        # we have self.config for configurations
        #
        # job_template
        # submit_cmd
        #
        # => status_cmd (perhaps not needed)
        # => kill_cmd (perhaps not needed)
        if 'job_template' in self.config:
            self.job_template = self.config['job_template'].replace('\r\n', '\n')
        elif 'template_file' in self.config:
            if not os.path.isfile(os.path.expanduser(self.config['template_file'])):
                raise ValueError('Missing job_template file {} for queue {}'.format(self.config['job_template'], self.alias))
            else:
                with open(os.path.expanduser(self.config['template_file'])) as tmpl:
                    self.job_template = tmpl.read()
                env.logger.warning('Option template_file is deprecated and will be removed from next formal release of SoS.')
        else:
            raise ValueError('A job_template is required for queue {}'.format(self.alias))

        if 'submit_cmd' not in self.config:
            raise ValueError('Missing configuration submit_cmd for queue {}'.format(self.alias))
        else:
            self.submit_cmd = self.config['submit_cmd']

        if 'status_cmd' not in self.config:
            raise ValueError('Missing configuration status_cmd for queue {}'.format(self.alias))
        else:
            self.status_cmd = self.config['status_cmd']

        if 'kill_cmd' not in self.config:
            raise ValueError('Missing configuration kill_cmd for queue {}'.format(self.alias))
        else:
            self.kill_cmd = self.config['kill_cmd']

    def execute_task(self, task_id):
        #
        if not super(PBS_TaskEngine, self).execute_task(task_id):
            return False

        try:
            return self._prepare_script(task_id)
        except Exception as e:
            env.logger.error(e)
            return False

    def _prepare_script(self, task_id):
        # read the task file and look for runtime info
        # 
        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.task')
        params = loadTask(task_file)
        sos_dict = params.sos_dict

        # for this task, we will need walltime, nodes, cores, mem
        # however, these could be fixed in the job template and we do not need to have them all in the runtime
        runtime = self.config
        runtime.update({x:sos_dict['_runtime'][x] for x in ('nodes', 'cores', 'ppn', 'mem', 'walltime', 'cur_dir', 'home_dir', 'name') if x in sos_dict['_runtime']})
        runtime['task'] = task_id
        runtime['verbosity'] = env.verbosity
        runtime['sig_mode'] = env.config['sig_mode']
        runtime['run_mode'] = env.config['run_mode']
        if 'name' in runtime:
            runtime['job_name'] = interpolate(runtime['name'], '${ }', sos_dict)
        else:
            runtime['job_name'] = interpolate('${step_name}_${_index}', '${ }', sos_dict)
        if 'nodes' not in runtime:
            runtime['nodes'] = 1
        if 'cores' not in runtime:
            if 'ppn' in runtime:
                env.logger.warning('Option ppn is deprecated and will be removed from a formal release of SoS')
                runtime['cores'] = runtime['ppn']
            else:
                runtime['cores'] = 1
        # for backward compatibility
        runtime['ppn'] = runtime['cores']
        runtime['job_file'] = '~/.sos/tasks/{}.sh'.format(task_id)

        # let us first prepare a task file
        try:
            job_text = interpolate(self.job_template, '${ }', runtime)
        except Exception as e:
            raise ValueError('Failed to generate job file for task {}: {}'.format(task_id, e))

        # now we need to write a job file
        job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.sh')
        # do not translate newline under windows because the script will be executed
        # under linux/mac
        with open(job_file, 'w', newline='') as job:
            job.write(job_text)

        # then copy the job file to remote host if necessary
        self.agent.send_task_file(task_id + '.sh')

        if env.config['run_mode'] == 'dryrun':
            try:
                cmd = 'bash ~/.sos/tasks/{}.sh'.format(task_id)
                print(self.agent.check_output(cmd))
            except Exception as e:
                raise RuntimeError('Failed to submit task {}: {}'.format(task_id, e))
        else:
            #
            # now we need to figure out a command to submit the task
            try:
                cmd = interpolate(self.submit_cmd, '${ }', runtime)
            except Exception as e:
                raise ValueError('Failed to generate job submission command from template "{}": {}'.format(
                    self.submit_cmd, e))
            env.logger.debug('submit {}: {}'.format(task_id, cmd))
            try:
                try:
                    cmd_output = self.agent.check_output(cmd).strip()
                except Exception as e:
                    raise RuntimeError('Failed to submit task {}: {}'.format(task_id, e))

                if 'submit_cmd_output' not in self.config:
                    submit_cmd_output = '{job_id}'
                else:
                    submit_cmd_output = self.config['submit_cmd_output']
                #
                if not '{job_id}' in submit_cmd_output:
                    raise ValueError('Option submit_cmd_output should have at least a pattern for job_id, "{}" specified.'.format(submit_cmd_output))
                #
                # try to extract job_id from command output
                # let us write an job_id file so that we can check status of tasks more easily
                job_id_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.job_id')
                with open(job_id_file, 'w') as job:
                    try:
                        res = extract_pattern(submit_cmd_output, [cmd_output.strip()])
                        if 'job_id' not in res or len(res['job_id']) != 1:
                            env.logger.warning('Failed to extract job_id from "{}" using pattern "{}"'.format(
                                cmd_output.strip(), submit_cmd_output))
                            job_id = '000000'
                            job.write('job_id: {}\n'.format(job_id))
                        else:
                            job_id = res['job_id'][0]
                            # other variables
                            for k,v in res.items():
                                job.write('{}: {}\n'.format(k, v[0]))
                    except Exception as e:
                        env.logger.warning('Failed to extract job_id from "{}" using pattern "{}"'.format(
                            cmd_output.strip(), submit_cmd_output))
                        job_id = '000000'
                        job.write('job_id: {}\n'.format(job_id))
                # output job id to stdout
                env.logger.info('{} ``submitted`` to {} with job id {}'.format(task_id, self.alias, job_id))
                return True
            except Exception as e:
                raise RuntimeError('Failed to submit task {}: {}'.format(task_id, e))

    def _get_job_id(self, task_id):
        job_id_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias, task_id + '.job_id')
        if not os.path.isfile(job_id_file):
            return {}
        with open(job_id_file) as job:
            result = {}
            for line in job:
                k, v = line.split(':', 1)
                result[k.strip()] = v.strip()
            return result

    def _query_job_status(self, job_id, task_id):
        job_id.update({'task': task_id, 'verbosity': 1})
        cmd = interpolate(self.status_cmd, '${ }', job_id)
        return self.agent.check_output(cmd)

    def query_tasks(self, tasks, verbosity=1, html=False, start_time=True, age=None):
        if verbosity <= 2:
            status_lines = super(PBS_TaskEngine, self).query_tasks(tasks, verbosity, html, start_time, age=age)
            # there is a change that a job is submitted, but failed before the sos command is executed
            # so we will have to ask the task engine about the submitted jobs #608
            if not html:
                res = ''
                for line in status_lines.split('\n'):
                    if not line.strip():
                        continue
                    fields = line.split('\t')
                    if len(fields) < verbosity:
                        env.logger.error(fields)
                        env.logger.warning('Suspicious status line {}'.format(line))
                        continue
                    task_id = fields[0]
                    if fields[verbosity] == 'submitted':
                        try:
                            job_id = self._get_job_id(task_id)
                            if not job_id:
                                raise RuntimeError('failed to obtain job id for task {}'.format(task_id))
                            self._query_job_status(job_id, task_id)
                        except Exception as e:
                            env.logger.trace('Failed to query status for task {}: {}'.format(task_id, e))
                            fields[verbosity] = 'failed'
                    res += '\t'.join(fields) + '\n'
                return res
            else:
                # ID line: <tr><th align="right"  width="30%">ID</th><td align="left">5173b80bf85d3d03153b96f9a5b4d6cc</td></tr>
                task_id = status_lines.split('>ID<', 1)[-1].split('</td',1)[0].split('>')[-1]
                status = status_lines.split('>Status<', 1)[-1].split('</td',1)[0].split('>')[-1]
                if status == 'submitted':
                    try:
                        job_id = self._get_job_id(task_id)
                        if not job_id:
                            raise RuntimeError('failed to obtain job id for task {}'.format(task_id))
                        self._query_job_status(job_id, task_id)
                    except Exception as e:
                        env.logger.trace('Failed to query status for task {}: {}'.format(task_id, e))
                        status_lines = status_lines.replace('submitted', 'failed', 1)
                return status_lines

        # for more verbose case, we will call pbs's status_cmd to get more accurate information
        status_lines = super(PBS_TaskEngine, self).query_tasks(tasks, 1)
        res = ''
        for line in status_lines.split('\n'):
            if not line.strip():
                continue
            task_id, status = line.split('\t')
            # call query_tasks again for more verbose output
            res += super(PBS_TaskEngine, self).query_tasks([task_id], verbosity, html) + '\n'
            #
            try:
                job_id = self._get_job_id(task_id)
                if not job_id:
                    # no job id file
                    raise RuntimeError('failed to obtain job id for task {}'.format(task_id))
                res += self._query_job_status(job_id, task_id)
            except Exception as e:
                env.logger.debug('Failed to get status of task {} (job_id: {}) from template "{}": {}'.format(
                    task_id, job_id, self.status_cmd, e))
        return res

    def kill_tasks(self, tasks, all_tasks=False):
        # remove the task from SoS task queue, this would also give us a list of
        # tasks on the remote server
        output = super(PBS_TaskEngine, self).kill_tasks(tasks, all_tasks)
        env.logger.trace('Output of local kill: {}'.format(output))
        # then we call the real PBS commands to kill tasks
        res = ''
        for line in output.split('\n'):
            if not line.strip():
                continue
            task_id, status = line.split('\t')
            res += '{}\t{} (old status)\t'.format(task_id, status)

            job_id = self._get_job_id(task_id)
            if not job_id:
                env.logger.debug('No job_id for task {}'.format(task_id))
                continue
            try:
                job_id.update({'task': task_id})
                cmd = interpolate(self.kill_cmd, '${ }', job_id)
                env.logger.debug('Running {}'.format(cmd))
                res += self.agent.check_output(cmd) + '\n'
            except Exception as e:
                env.logger.debug('Failed to kill job {} (job_id: {}) from template "{}": {}'.format(
                    task_id, job_id, self.kill_cmd, e))
        return res
