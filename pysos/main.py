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

from .utils import env, get_traceback
from .sos_script import SoS_Script
#
# subcommmand show
#
def sos_show(args, workflow_args):
    try:
        script = SoS_Script(filename=args.script)
        if args.workflow:
            workflow = script.workflow(args.workflow)
            workflow.show()
        else:
            script.show()
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand dryrun
#
def sos_dryrun(args, workflow_args):
    args.__max_jobs__ = 1
    args.__dryrun__ = True
    args.__prepare__ = False
    args.__run__ = False
    args.__rerun__ = False
    args.__config__ = None
    sos_run(args, workflow_args)

#
# subcommand prepare
#
def sos_prepare(args, workflow_args):
    args.__max_jobs__ = 1
    args.__dryrun__ = True
    args.__prepare__ = True
    args.__run__ = False
    args.__rerun__ = False
    args.__config__ = None
    sos_run(args, workflow_args)

#
# subcommand run
#
def sos_run(args, workflow_args):
    env.max_jobs = args.__max_jobs__
    env.verbosity = args.verbosity
    # kill all remainging processes when the master process is killed.
    atexit.register(env.cleanup)
    # default mode: run in dryrun mode
    args.__run__ = not (args.__rerun__ or args.__prepare__ or args.__dryrun__)
    #
    if args.__run__ or args.__rerun__:
        args.__prepare__ = True
    #
    # always run in dryrun mode
    env.run_mode = 'dryrun'
    # if this is not the last step, use verbosity 1 (warning)
    #if args.__prepare__:
    #    env.verbosity = min(args.verbosity, 1)
    #else:
    #
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        workflow.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)
    # then prepare mode
    if args.__prepare__:
        # if this is not the last step, use verbosity 1 (warning)
        #if args.__run__ or args.__rerun__:
        #    env.verbosity = min(args.verbosity, 1)
        #else:
        #    env.verbosity = args.verbosity
        #
        env.run_mode = 'prepare'
        try:
            script = SoS_Script(filename=args.script)
            workflow = script.workflow(args.workflow)
            workflow.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
        except Exception as e:
            if args.verbosity and args.verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error(e)
            sys.exit(1)
    # then run mode
    if args.__run__ or args.__rerun__:
        env.run_mode = 'run'
        # env.verbosity = args.verbosity
        if args.__rerun__:
            env.sig_mode = 'ignore'
        try:
            script = SoS_Script(filename=args.script)
            workflow = script.workflow(args.workflow)
            workflow.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
        except Exception as e:
            if args.verbosity and args.verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error(e)
            sys.exit(1)

#
# subcommand config
#
def sos_config(args, workflow_args):
    if workflow_args:
        raise RuntimeError('Unrecognized arguments {}'.format(' '.join(workflow_args)))
    #
    if args.__global_config__:
        config_file = os.path.expanduser('~/.sos/config.json')
    elif args.__config_file__:
        config_file = os.path.expanduser(args.__config_file__)
    else:
        config_file = os.path.expanduser('.sos/config.json')
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
        for option in args.__set_config__:
            k, v = option.split('=', 1)
            try:
                v = eval(v)
            except Exception as e:
                env.logger.error('Cannot interpret option {}. Please quote the string if it is a string option. ({})'.format(option, e))
                sys.exit(1)
            cfg[k] = v
            print('Set {} to {!r}'.format(k, v))
        #
        with open(config_file, 'w') as config:
            config.write(yaml.safe_dump(cfg, default_flow_style=False))
        

        
