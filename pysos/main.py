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

import os
import sys
import yaml
import atexit
import fnmatch
import tempfile

from .utils import env, get_traceback, dict_merge
from .sos_script import SoS_Script
from .sos_executor import Base_Executor, Sequential_Executor
from .converter import script_to_html, workflow_to_html, script_to_markdown, \
    workflow_to_markdown, script_to_notebook, workflow_to_notebook, \
    script_to_term, workflow_to_term, notebook_to_script
from .sos_syntax import CONFIG_NAME
from .target import FileTarget
#
# subcommand convert
#
def cmd_convert(args, style_args):
    env.verbosity = args.verbosity
    # convert from ...
    try:
        if args.sos:
            if not args.from_file.endswith('ipynb'):
                raise RuntimeError('Can only convert from a .ipynb file to SoS script: {} specified.'.format(args.from_file))
            notebook_to_script(args.from_file, args.sos, style_args)
        elif args.notebook and args.from_file.lower().endswith('.ipynb'):
            try:
                sos_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.sos', delete=False).name
                notebook_to_script(args.from_file, sos_file, style_args)
                transcript_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.transcript', delete=False).name
                script_to_notebook(sos_file, args.notebook)
                sys.exit(0)
            finally:
                os.remove(sos_file)
                os.remove(transcript_file)
        else:
            transcript_file = os.path.join('.sos', '{}.transcript'.format(os.path.basename(args.from_file)))
            with open(transcript_file, 'w') as transcript:
                try:
                    script = SoS_Script(filename=args.from_file, transcript=transcript)
                except Exception as e:
                    script = None
                    env.logger.warning(e)
            if args.workflow:
                if not script:
                    raise RuntimeError('workflow {} is not available due to syntax error in script {}'.format(args.workflow, args.from_file))
                workflow = script.workflow(args.workflow)
                if args.html is not None:
                    workflow_to_html(workflow, args.from_file, args.html, style_args)
                elif args.markdown is not None:
                    workflow_to_markdown(workflow, args.from_file, args.markdown, style_args)
                elif args.notebook is not None:
                    workflow_to_notebook(workflow, args.from_file, args.notebook)
                elif args.term:
                    workflow_to_term(workflow, args.from_file, style_args)
                else:
                    workflow.show()
            else:
                if args.html is not None:
                    script_to_html(transcript_file, args.from_file, args.html, style_args)
                elif args.markdown is not None:
                    script_to_markdown(transcript_file, args.from_file, args.markdown)
                elif args.notebook is not None:
                    script_to_notebook(args.from_file, args.notebook)
                elif args.term:
                    script_to_term(transcript_file, args.from_file, style_args)
                else:
                    script.show()
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand inspect
#
def cmd_inspect(args, workflow_args):
    env.verbosity = args.verbosity
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        executor = Base_Executor(workflow, args=workflow_args, config_file=args.__config__)
        executor.inspect()
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand prepare
#
def cmd_prepare(args, workflow_args):
    env.verbosity = args.verbosity

    if args.__rerun__:
        env.sig_mode = 'ignore'
    elif args.__construct__:
        env.sig_mode = 'construct'
    else:
        env.sig_mode = 'default'

    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        executor = Base_Executor(workflow, args=workflow_args, config_file=args.__config__)
        executor.inspect()
        executor.prepare()
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand run
#
def cmd_run(args, workflow_args):
    env.max_jobs = args.__max_jobs__
    env.verbosity = args.verbosity
    # kill all remainging processes when the master process is killed.
    atexit.register(env.cleanup)
    #
    if args.__rerun__:
        env.sig_mode = 'ignore'
    elif args.__construct__:
        env.sig_mode = 'construct'
    else:
        env.sig_mode = 'default'

    for filename in args.__discard__:
        # remove also signature of file if it exists
        FileTarget(filename).remove('both')
            
    run_all = not (args.__inspect__ or args.__prepare__ or args.__run__)
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        executor = Sequential_Executor(workflow, args=workflow_args, config_file=args.__config__)
        if run_all or args.__inspect__:
            executor.inspect()
        if run_all or args.__prepare__:
            dag = executor.prepare(args.__targets__)
        else:
            dag = None
        if run_all or args.__run__:
            # if dag is None, the script will be run sequentially and cannot handle
            # make-style steps.
            executor.run(dag)
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand config
#
def cmd_config(args, workflow_args):
    if workflow_args:
        raise RuntimeError('Unrecognized arguments {}'.format(' '.join(workflow_args)))
    #
    if args.__global_config__:
        config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yaml')
    elif args.__config_file__:
        config_file = os.path.expanduser(args.__config_file__)
    else:
        config_file = 'config.yaml'
    if args.__get_config__ is not None:
        if os.path.isfile(config_file):
            try:
                with open(config_file) as config:
                    cfg = yaml.safe_load(config)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to parse sos config file {}, is it in YAML/JSON format? ({}}'.format(config_file, e))
                sys.exit(1)
            for option in (args.__get_config__ if args.__get_config__ else ['*']):
                for k, v in cfg.items():
                    if fnmatch.fnmatch(k, option):
                        print('{}\t{!r}'.format(k, v))
    elif args.__unset_config__:
        if os.path.isfile(config_file):
            try:
                with open(config_file) as config:
                    cfg = yaml.safe_load(config)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to parse sos config file {}, is it in YAML/JSON format? ({})'.format(config_file, e))
                sys.exit(1)
        else:
            env.logger.error('Config file {} does not exist'.format(config_file))
        #
        unset = []
        for option in args.__unset_config__:
            for k in cfg.keys():
                if fnmatch.fnmatch(k, option):
                    unset.append(k)
                    print('Unset {}'.format(k))
        #
        if unset:
            for k in set(unset):
                cfg.pop(k)
            #
            if unset:
                with open(config_file, 'w') as config:
                    config.write(yaml.safe_dump(cfg, default_flow_style=False))
    elif args.__set_config__:
        if os.path.isfile(config_file):
            try:
                with open(config_file) as config:
                    cfg = yaml.safe_load(config)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to sos config file {}, is it in YAML/JSON format? ({})'.format(config_file, e))
                sys.exit(1)
        else:
            cfg = {}
        #
        if len(args.__set_config__) == 1:
            env.logger.error('Please specify a value for key {}'.format(args.__set_config__[0]))
            sys.exit(1)
        #
        k = args.__set_config__[0]
        if not CONFIG_NAME.match(k):
            env.logger.error('Unacceptable variable name for config file: {}'.format(k))
            sys.exit(1)
        #
        values = []
        for v in args.__set_config__[1:]:
            try:
                v_val = eval(v)
                # test if the value can be saved by yaml
                yaml.safe_dump(v_val)
                v = v_val
            except Exception:
                env.logger.warning('Value "{}" is an invalid expression and is treated as a string.'.format(v))
            values.append(v)
        #
        if len(values) == 1:
            values = values[0]
        #
        # say v   = 1
        #     key = a.b.c
        #
        # this gives:
        #     {'a': {'b': {'c': 1}}}
        dv = values
        for key in reversed(k.split('.')):
            new_dv = {}
            new_dv[key] = dv
            dv = new_dv
        # however we can not update directly and has to merge two
        # dictionaries. For example, if
        #
        # cfg = {'a': {'b': {'d': 2}}, {'c': 1}}
        #
        # we need to get
        #
        # cfg = {'a': {'b': {'d': 2, 'c': 1}}, {'c': 1}}
        #
        dict_merge(cfg, dv)
        # reporting assignment with existing values
        print('Set {} to {!r}'.format(k.split('.')[0], cfg[k.split('.')[0]]))
        #
        with open(config_file, 'w') as config:
            config.write(yaml.safe_dump(cfg, default_flow_style=False))



