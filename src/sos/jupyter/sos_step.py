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


from sos.sos_step import Step_Executor, Base_Step_Executor, PendingTasks
from sos.hosts import Host
from sos.utils import env, short_repr
import time

class Interactive_Step_Executor(Step_Executor):
    def __init__(self, step):
        # This is the only interesting part of this executor. Basically
        # it derives everything from SP_Step_Executor but does not
        # use the Queue mechanism, so the __init__ and the run
        # functions are copied from Base_Step_Executor
        Base_Step_Executor.__init__(self, step)
        self.run_mode='interactive'

    def pending_tasks(self, tasks):
        if not tasks:
            return
        if 'queue' in env.sos_dict['_runtime']:
            queue = env.sos_dict['_runtime']['queue']
        elif env.config['default_queue']:
            queue = env.config['default_queue']
        else:
            queue = 'localhost'

        host = Host(queue)
        res = [host.submit_task(task) for task in tasks]
        if all(x == 'completed' for x in host.check_status(tasks)):
            if len(tasks) > 4:
                print('!sos_hint: {} task{} completed: {}, {}, ..., {}'.format(len(tasks), 's' if len(tasks) > 1 else '',
                    """<a onclick="task_info('{}', '{}')">{}</a>""".format(tasks[0], queue, tasks[0][:4]),
                    """<a onclick="task_info('{}', '{}')">{}</a>""".format(tasks[1], queue, tasks[1][:4]),
                    """<a onclick="task_info('{}', '{}')">{}</a>""".format(tasks[-1], queue, tasks[-1][:4])))
            else:
                print('!sos_hint: {} task{} completed: {}'.format(len(tasks), 's' if len(tasks) > 1 else '',
                    ','.join(["""<a onclick="task_info('{}', '{}')">{}</a>""".format(x, queue, x[:4]) for x in tasks])))
            host._task_engine.remove_tasks(tasks)
            return host.retrieve_results(tasks)
        while True:
            res = host.check_status(tasks)
            if all(x not in ('submitted', 'pending', 'running') for x in res):
                #completed = [task for task, status in zip(tasks, res) if status == 'completed']
                host._task_engine.remove_tasks(tasks)
                return host.retrieve_results(tasks)
            # no pending
            elif not env.config['wait_for_task']:
                raise PendingTasks([(queue, x) for x,y in zip(tasks, res) if y in ('pending', 'submitted', 'running')])
            time.sleep(1)


    def run(self):
        return Base_Step_Executor.run(self)

    def log(self, stage=None, msg=None):
        if stage == 'start':
            env.logger.debug('{} ``{}``: {}'.format('Checking' if self.run_mode == 'dryrun' else 'Executing',
                self.step.step_name(), self.step.comment.strip()))
        elif stage == 'input':
            if env.sos_dict['input'] is not None:
                env.logger.debug('input:    ``{}``'.format(short_repr(env.sos_dict['input'])))
        elif stage == 'output':
            if env.sos_dict['output'] is not None:
                env.logger.debug('output:   ``{}``'.format(short_repr(env.sos_dict['output'])))

