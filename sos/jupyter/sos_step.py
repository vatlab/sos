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
from sos.sos_eval import Undetermined, SoS_eval
from sos.target import FileTarget, dynamic, RemovedTarget

class Interactive_Step_Executor(Base_Step_Executor):
    def __init__(self, step):
        Base_Step_Executor.__init__(self, step)

    def log(self, stage=None, msg=None):
        if stage == 'start':
            env.logger.info('Running ``{}``: {}'.format(self.step.step_name(), self.step.comment.strip()))

    def collect_result(self):
        return self.last_res

    def verify_input(self):
        # now, if we are actually going to run the script, we
        # need to check the input files actually exists, not just the signatures
        if isinstance(env.sos_dict['_input'], list):
            for target in env.sos_dict['_input']:
                # if the file does not exist (although the signature exists)
                # request generation of files
                if isinstance(target, str) and not FileTarget(target).exists('target'):
                    # remove the signature and regenerate the file
                    FileTarget(target).remove('signature')
                    raise RemovedTarget(target)

    def expand_input_files(self, value, *args):
        # if unspecified, use __step_output__ as input (default)
        # resolve dynamic input.
        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        if not args:
            return env.sos_dict['input']
        else:
            return _expand_file_list(False, *args)

    def expand_depends_files(self, *args, **kwargs):
        '''handle directive depends'''
        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        return _expand_file_list(False, *args)

    def reevaluate_output(self):
        # re-process the output statement to determine output files
        args, kwargs = SoS_eval('__null_func__({})'.format(env.sos_dict['output'].expr), self.step.sigil)
        # handle dynamic args
        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        env.sos_dict.set('output', self.expand_output_files('', *args))

    def verify_output(self):
        if env.sos_dict['output'] is None:
            return
        if isinstance(env.sos_dict['output'], Undetermined):
            raise RuntimeError('Output of a completed step cannot be undetermined.')
        for target in env.sos_dict['output']:
            if isinstance(target, str):
                if not FileTarget(target).exists('target' if '__hard_target__' in env.sos_dict else 'any'):
                    raise RuntimeError('Output target {} does not exist after the completion of step {}'
                            .format(target, env.sos_dict['step_name']))
            elif not target.exists('any'):
                raise RuntimeError('Output target {} does not exist after the completion of step {}'
                            .format(target, env.sos_dict['step_name']))


