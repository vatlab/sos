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

from sos.utils import env
from sos.sos_executor import MP_Executor

from .sos_step import Celery_Step_Executor

class Celery_Executor(MP_Executor):
    def __init__(self, workflow, args=[], nested=False, config={}):
        MP_Executor.__init__(self, workflow, args, nested=nested, config=config)
        env.__task_engine__ = 'Celery'

    def step_executor(self, section, queue):
        # pass celery_app if needed
        return Celery_Step_Executor(section, queue)

