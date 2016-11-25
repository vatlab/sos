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
import os
import yaml
import shlex
import argparse
from sos.utils import env, frozendict, dict_merge
from sos.sos_eval import get_default_global_sigil
from sos._version import __version__
from sos.sos_script import SoS_Script
from sos.sos_executor import Base_Executor, __null_func__
from sos.sos_syntax import SOS_SECTION_HEADER

from .sos_step import Interactive_Step_Executor

class Interactive_Executor(Base_Executor):
    '''Interactive executor called from by iPython Jupyter or Spyder'''
    def __init__(self):
        # we actually do not have our own workflow, everything is passed from ipython
        # by nested = True we actually mean no new dictionary
        Base_Executor.__init__(self, nested=True)

    def parse_command_line(self, command_line):
        parser = argparse.ArgumentParser()
        # no default workflow so it will execute any workflow if the code piece
        # defines only one workflow
        # 
        # parser.add_argument('-j', type=int, metavar='JOBS', default=1, dest='__max_jobs__')
        parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE')
        #parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE',
        #    default=os.path.join('.sos', '__step_report.md'))
        runmode = parser.add_argument_group(title='Run mode options')
        runmode.add_argument('-f', action='store_true', dest='__rerun__')
        runmode.add_argument('-F', action='store_true', dest='__construct__')
        # default to 1 to avoid output env.logger.info to notebook
        parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=1)
        #
        args, workflow_args = parser.parse_known_args(shlex.split(command_line))
        return args, workflow_args

    def parse_script(self, code):
        '''Used by the kernel to judge if the code is complete'''
        return SoS_Script(content=code, global_sigl=get_default_global_sigil())

    def set_dict(self, args):
        env.sos_dict.set('__null_func__', __null_func__)
        env.sos_dict.set('SOS_VERSION', __version__)

        # load configuration files
        cfg = {}
        sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yaml')
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        # local config file
        sos_config_file = 'config.yaml'
        if os.path.isfile(sos_config_file):
            try:
                with open(sos_config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse local sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
        if args.__config__ is not None:
            # user-specified configuration file.
            if not os.path.isfile(args.__config__):
                raise RuntimeError('Config file {} not found'.format(args.__config__))
            try:
                with open(args.__config__) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format? ({})'.format(self.config_file, e))
        # set config to CONFIG
        env.sos_dict.set('CONFIG', frozendict(cfg))

    def run(self, block, command_line=''):
        '''Execute a block of SoS script that is sent by iPython/Jupyer/Spyer
        The code can be simple SoS/Python statements, one SoS step, or more
        or more SoS workflows with multiple steps. This executor,
        1. adds a section header to the script if there is no section head
        2. execute the workflow in interactive mode, which is different from
           batch mode in a number of ways, which most notably without support
           for nested workflow.
        3. Optionally execute the workflow in preparation mode for debugging purposes.
        '''
        # if there is no valid code do nothing
        if not block.strip():
            return
        # if there is no section header, add a header so that the block
        # appears to be a SoS script with one section
        if not any([SOS_SECTION_HEADER.match(line) for line in block.split()]):
            block = '[interactive_0]\n' + block

        script = SoS_Script(content=block, global_sigil=get_default_global_sigil())
        env.run_mode = 'interactive'
        try:
            args, workflow_args = self.parse_command_line(command_line)
            env.sos_dict.set('__args__', workflow_args)
            env.sos_dict.set('__unknown_args__', workflow_args)
            self.set_dict(args)
            self.workflow = script.workflow()

            if args.__rerun__:
                env.sig_mode = 'ignore'
            elif args.__construct__:
                env.sig_mode = 'construct'
            else:
                env.sig_mode = 'default'

            #if os.path.isfile(args.__report__):
            #    os.remove(args.__report__)

            # this is the result returned by the workflow, if the
            # last stement is an expression.
            last_res = None
            #
            # clear __step_input__, __step_output__ etc because there is
            # no concept of passing input/outputs across cells.
            env.sos_dict.set('__step_output__', [])
            for k in ['__step_input__', '__default_output__', 'input', 'output', \
                'depends', '_input', '_output', '_depends']:
                env.sos_dict.pop(k, None)

            for idx, section in enumerate(self.workflow.sections):
                if 'skip' in section.options:
                    val_skip = section.options['skip']
                    if val_skip is None or val_skip is True:
                        continue
                    elif val_skip is not False:
                        raise RuntimeError('The value of section option skip can only be None, True or False, {} provided'.format(val_skip))
                #
                last_res = Interactive_Step_Executor(section).run()
                # if the step is failed
                if isinstance(last_res, Exception):
                    raise RuntimeError(last_res)
            return last_res
        finally:
            env.verbosity = 1
            env.sig_mode = 'default'



