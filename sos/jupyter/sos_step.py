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

from sos.utils import env

from sos.sos_step import Step_Executor, Base_Step_Executor
from sos.hosts import Host
import time

class Interactive_Step_Executor(Step_Executor):
    def __init__(self, step):
        # This is the only interesting part of this executor. Basically
        # it derives everything from SP_Step_Executor but does not
        # use the Queue mechanism, so the __init__ and the run
        # functions are copied from Base_Step_Executor
        Base_Step_Executor.__init__(self, step)

    def pending_tasks(self, tasks):
        host = Host()
        for task in tasks:
            host.submit_task(task)
        while True:
             res = host.check_status(tasks)
             if any(x in  ('pending', 'running', 'completed-old', 'failed-old', 'failed-missing-output', 'failed-old-missing-output') for x in res):
                continue
             elif all(x == 'completed' for x in res):
                return host.retrieve_results(tasks)
             time.sleep(1)

    def run(self):
        return Base_Step_Executor.run(self)

    def log(self, stage=None, msg=None):
        if stage == 'start':
            env.logger.info('Running ``{}``: {}'.format(self.step.step_name(), self.step.comment.strip()))

