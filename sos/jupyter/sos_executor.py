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
import sys
import os
import yaml
import shlex
import argparse
from sos.utils import env, frozendict, dict_merge, _parse_error, get_traceback
from sos.sos_eval import SoS_exec
from sos._version import __version__
from sos.__main__ import add_run_arguments
from sos.sos_script import SoS_Script
from sos.sos_executor import Base_Executor, __null_func__
from sos.sos_syntax import SOS_SECTION_HEADER
from sos.target import FileTarget, UnknownTarget, RemovedTarget, UnavailableLock
from .sos_step import Interactive_Step_Executor

class Interactive_Executor(Base_Executor):
    '''Interactive executor called from by iPython Jupyter or Spyder'''
    def __init__(self, workflow=None, args=[], config_file=None, **kwargs):
        # we actually do not have our own workflow, everything is passed from ipython
        # by nested = True we actually mean no new dictionary
        Base_Executor.__init__(self, workflow=workflow, args=args, nested=True)
        self.__config__ = config_file
        env.__task_engine__ = 'interactive'

    def set_dict(self):
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
        if self.__config__ is not None:
            # user-specified configuration file.
            if not os.path.isfile(self.__config__):
                raise RuntimeError('Config file {} not found'.format(self.__config__))
            try:
                with open(self.__config__) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format? ({})'.format(self.config_file, e))
        # set config to CONFIG
        env.sos_dict.set('CONFIG', frozendict(cfg))

    def run(self, targets=None):
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
        self.set_dict()

        # this is the result returned by the workflow, if the
        # last stement is an expression.
        last_res = None

        # process step of the pipelinp
        dag = self.initialize_dag(targets=targets)
        #
        # if targets are specified and there are only signatures for them, we need
        # to remove the signature and really generate them
        if targets:
            for t in targets:
                if not FileTarget(t).exists('target'):
                    FileTarget(t).remove('signature')
        #
        self.set_dict()
        while True:
            # find any step that can be executed and run it, and update the DAT
            # with status.
            runnable = dag.find_executable()
            if runnable is None:
                # no runnable
                #dag.show_nodes()
                break
            # find the section from runnable
            section = self.workflow.section_by_id(runnable._step_uuid)
            #
            # this is to keep compatibility of dag run with sequential run because
            # in sequential run, we evaluate global section of each step in
            # order to determine values of options such as skip.
            # The consequence is that global definitions are available in
            # SoS namespace.
            try:
                SoS_exec(section.global_def, section.global_sigil)
            except Exception as e:
                if env.verbosity > 2:
                    sys.stderr.write(get_traceback())
                raise RuntimeError('Failed to execute statements\n"{}"\n{}'.format(
                    section.global_def, e))

            # clear existing keys, otherwise the results from some random result
            # might mess with the execution of another step that does not define input
            for k in ['__step_input__', '__default_output__', '__step_output__']:
                if k in env.sos_dict:
                    env.sos_dict.pop(k)
            # if the step has its own context
            env.sos_dict.quick_update(runnable._context)
            # execute section with specified input
            runnable._status = 'running'
            try:
                executor = Interactive_Step_Executor(section)
                last_res = executor.run()
                # set context to the next logic step.
                for edge in dag.out_edges(runnable):
                    node = edge[1]
                    # if node is the logical next step...
                    if node._node_index is not None and runnable._node_index is not None:
                        #and node._node_index == runnable._node_index + 1:
                        node._context.update(env.sos_dict.clone_selected_vars(
                            node._context['__signature_vars__'] | node._context['__environ_vars__'] \
                            | {'_input', '__step_output__', '__default_output__', '__args__'}))
                runnable._status = 'completed'
            except UnknownTarget as e:
                runnable._status = None
                target = e.target
                if self.resolve_dangling_targets(dag, [target]) == 0:
                    raise RuntimeError('Failed to resolve {}{}.'
                        .format(target, dag.steps_depending_on(target, self.workflow)))
                # now, there should be no dangling targets, let us connect nodes
                # this can be done more efficiently
                runnable._depends_targets.append(target)
                dag._all_dependent_files[target].append(runnable)
                #
                dag.build(self.workflow.auxiliary_sections)
                #dag.show_nodes()
                cycle = dag.circular_dependencies()
                if cycle:
                    raise RuntimeError('Circular dependency detected {}. It is likely a later step produces input of a previous step.'.format(cycle))
            except RemovedTarget as e:
                runnable._status = None
                dag.regenerate_target(e.target)
            except UnavailableLock:
                runnable._status = 'pending'
                runnable._signature = (e.output, e.sig_file)
                env.logger.info('Waiting on another process for step {}'.format(section.step_name()))
            # if the job is failed
            except Exception as e:
                runnable._status = 'failed'
                raise

        return last_res

#
# function runfile that is used by spyder to execute complete script
#
def runfile(script=None, args='', wdir='.', code=None, **kwargs):
    parser = argparse.ArgumentParser(description='''Execute a sos script''')
    add_run_arguments(parser, interactive=True)
    parser.error = _parse_error
    env.sos_dict.set('__args__', args)
    if isinstance(args, str):
        args = shlex.split(args)
    if (script is None and code is None) or '-h' in args:
        parser.print_help()
        return
    args, workflow_args = parser.parse_known_args(args)
    # calling the associated functions
    
    env.max_jobs = args.__max_jobs__
    env.verbosity = args.verbosity
    env.__task_engine__ = 'interactive'

    #
    if args.__rerun__:
        env.sig_mode = 'ignore'
    elif args.__construct__:
        env.sig_mode = 'construct'
    else:
        env.sig_mode = 'default'

    if args.__bin_dirs__:
        import fasteners
        for d in args.__bin_dirs__:
            if d == '~/.sos/bin' and not os.path.isdir(os.path.expanduser(d)):
                with fasteners.InterProcessLock('/tmp/sos_lock_bin'):
                    os.makedirs(os.path.expanduser(d))
            elif not os.path.isdir(os.path.expanduser(d)):
                raise ValueError('directory does not exist: {}'.format(d))
        os.environ['PATH'] = os.pathsep.join([os.path.expanduser(x) for x in args.__bin_dirs__]) + os.pathsep + os.environ['PATH']

    # clear __step_input__, __step_output__ etc because there is
    # no concept of passing input/outputs across cells.
    env.sos_dict.set('__step_output__', [])
    for k in ['__step_input__', '__default_output__', 'input', 'output', \
        'depends', '_input', '_output', '_depends']:
        env.sos_dict.pop(k, None)

    try:
        if script is None:
            if not code.strip():
                return
            # if there is no section header, add a header so that the block
            # appears to be a SoS script with one section
            if not any([SOS_SECTION_HEADER.match(line) for line in code.splitlines()]):
                code = '[interactive_0]\n' + code
            script = SoS_Script(content=code)
        else:
            script = SoS_Script(filename=script)
        workflow = script.workflow(args.workflow)
        executor = Interactive_Executor(workflow, args=workflow_args, config_file=args.__config__)
        #
        # if dag is None, the script will be run sequentially and cannot handle
        # make-style steps.
        return executor.run(args.__targets__)
    except Exception:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        raise
    finally:
        env.sig_mode = 'default'
        env.verbosity = 1

