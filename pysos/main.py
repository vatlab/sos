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
#
import os
import sys
import fasteners

#
# subcommand convert
#
def cmd_convert(args, style_args):
    import tempfile
    from .utils import env, get_traceback
    from .sos_script import SoS_Script
    from .converter import script_to_html, workflow_to_html, script_to_markdown, \
        workflow_to_markdown, script_to_notebook, workflow_to_notebook, \
        script_to_term, workflow_to_term, notebook_to_script
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
# subcommand run
#
def cmd_run(args, workflow_args):
    import atexit
    from .utils import env, get_traceback
    from .sos_script import SoS_Script
    from .target import FileTarget
    from .sos_executor import Base_Executor, MP_Executor, RQ_Executor, Celery_Executor
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

    if args.__bin_dirs__:
        for d in args.__bin_dirs__:
            with fasteners.InterProcessLock('/tmp/sos_lock_bin'):
                if d == '~/.sos/bin' and not os.path.isdir(os.path.expanduser(d)):
                    os.makedirs(os.path.expanduser(d))
                elif not os.path.isdir(os.path.expanduser(d)):
                    raise ValueError('directory does not exist: {}'.format(d))
        os.environ['PATH'] = os.pathsep.join([os.path.expanduser(x) for x in args.__bin_dirs__]) + os.pathsep + os.environ['PATH']
            
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        if args.__queue__ is None:
            if args.__max_jobs__ == 1:
                # single process executor
                executor = Base_Executor(workflow, args=workflow_args, config_file=args.__config__)
            else:
                executor = MP_Executor(workflow, args=workflow_args, config_file=args.__config__)
        elif args.__queue__ == 'rq':
            executor = RQ_Executor(workflow, args=workflow_args, config_file=args.__config__)
        elif args.__queue__ == 'celery':
            executor = Celery_Executor(workflow, args=workflow_args, config_file=args.__config__)
        else:
            raise ValueError('Only the default multiprocessing and a rq engine is allowed')
        #
        if args.__dryrun__:
            executor.dryrun(args.__targets__)
        elif args.__prepare__:
            executor.prepare(args.__targets__)
        else:
            # if dag is None, the script will be run sequentially and cannot handle
            # make-style steps.
            executor.run(args.__targets__)
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand dryrun
#
def cmd_dryrun(args, workflow_args):
    args.__rerun__ = False
    args.__construct__ = False
    args.__queue__ = None
    args.__max_jobs__ = 1
    args.__dryrun__ = True
    args.__prepare__ = True
    cmd_run(args, workflow_args)

#
# subcommand prepare
#
def cmd_prepare(args, workflow_args):
    args.__rerun__ = False
    args.__construct__ = False
    args.__queue__ = None
    args.__max_jobs__ = 1
    args.__dryrun__ = False
    args.__prepare__ = True
    cmd_run(args, workflow_args)

#
# command clean
#
def get_tracked_files(sig_file):
    with open(sig_file) as sig:
        start = False
        files = []
        for line in sig:
            if line.startswith('# input and dependent files'):
                start = True
                continue
            if line.startswith('#'):
                continue
            if start:
                files.append(line.strip())
        return set(files)

def cmd_clean(args, unknown_args):
    import glob
    from .utils import env
    import shutil
    from collections import OrderedDict
    sig_files = glob.glob('.sos/*.sig')
    if not sig_files:
        sys.exit('No executed workflows.')
    tracked_files = set()
    for sig_file in sig_files:
        tracked_files |= get_tracked_files(sig_file)
    #
    tracked_files = {os.path.abspath(os.path.expanduser(x)) for x in tracked_files}
    tracked_dirs = {os.path.dirname(x) for x in tracked_files}
    removed_files = []
    removed_dirs = []
    if tracked_files:
        env.logger.info('{} tracked files from {} run{} are identified.'
                .format(len(tracked_files), len(sig_files), 's' if len(sig_files) > 1 else ''))
    else:
        env.logger.info('No tracked file from {} run{} are identified.'
                .format(len(sig_files), 's' if len(sig_files) > 1 else ''))
    for dir in args.dirs:
        if not os.path.isdir(dir):
            sys.exit('Invalid directory to be cleaned: {}'.format(dir))
        relpath = os.path.relpath(dir, '.')
        if relpath.startswith('..'):
            sys.exit('Only subdirectories of the current directory can be cleaned. {} specified.'.format(dir))
        for dirname, subdir, filelist in os.walk(dir):
            # we do not remove files under the current directory
            if dirname != '.':
                removed_files.extend([os.path.join(dirname, x) for x in filelist if os.path.abspath(os.path.join(dirname, x)) \
                    not in tracked_files and not x.startswith('.')])
            # we do not track dot directories 
            untracked_dirs = [x for x in subdir if os.path.abspath(os.path.join(dirname, x)) not in tracked_dirs and not x.startswith('.')]
            removed_dirs.extend([os.path.join(dirname, x) for x in untracked_dirs])
            # do not scan the directory if it does not contain any tracked files
            subdir[:] = [d for d in subdir if d not in untracked_dirs and not d.startswith('.')]
    #
    def dedup(_list):
        return OrderedDict((item, None) for item in _list).keys()

    if args.__dryrun__:
        if not args.__files__ and not args.__dirs__:
            print('\n'.join(tracked_files))
            print('{} files are tracked.'.format(len(tracked_files)))
        if args.__files__:
            rf = dedup(removed_files)
            for r in rf:
                print('Would remove {}'.format(r))
            print('{} files to be removed.'.format(len(rf)))
        if args.__dirs__:
            rd = dedup(removed_dirs)
            for d in rd:
                print('Would remove {}'.format(d))
            print('{} directories to be removed.'.format(len(rd)))
    else:
        if not args.__files__ and not args.__dirs__:
            sys.exit('One of options -n -f -d needs to be specified.')
        if args.__files__:
            rf = dedup(removed_files)
            for f in rf:
                print('Removing {}'.format(f))
                os.remove(f)
            print('{} files are removed'.format(len(rf)))
        if args.__dirs__:
            rd = dedup(removed_dirs)
            for d in rd:
                print('Removing {}'.format(rd))
                if os.path.isdir(d):
                    shutil.rmtree(d)
                else:
                    os.unlink(d)
            print('{} directories are removed'.format(len(rd)))

#
# command start
#
def cmd_start(args, unknown_args):
    import subprocess
    if args.server_type == 'server':
        # TODO: run it in background so that sos would quit
        # TODO: write .sos/redis_connection.yaml
        subprocess.call('redis-server')
    elif args.server_type == 'worker':
        # read .sos/redis_connection.yaml
        # test redis connection???
        # write .sos/rq_worker_settings.py
        subprocess.call('rq worker -c .sos/rq_worker_settings')

#
# subcommand config
#
def cmd_config(args, workflow_args):
    import fnmatch
    import yaml
    from .utils import env, dict_merge
    from .sos_syntax import CONFIG_NAME
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



