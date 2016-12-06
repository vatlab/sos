#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
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

from sos.sos_step import Base_Step_Executor, _expand_file_list

class Interactive_Step_Executor(Base_Step_Executor):
    def __init__(self, step):
        Base_Step_Executor.__init__(self, step)

    def expand_input_files(self, value, *args):
        # We ignore 'dynamic' option in run mode
        # if unspecified, use __step_output__ as input (default)
        if not args:
            return env.sos_dict['input']
        else:
            return _expand_file_list(False, *args)

    def expand_depends_files(self, *args):
        '''handle directive depends'''
        return _expand_file_list(False, *args)

    def expand_output_files(self, value, *args):
        return _expand_file_list(True, *args)

    def log(self, stage=None, msg=None):
        if stage == 'start':
            env.logger.info('Running ``{}``: {}'.format(self.step.step_name(), self.step.comment.strip()))

    def collect_result(self):
        return self.last_res
