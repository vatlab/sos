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
# command remove
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

def cmd_remove(args, unknown_args):
    import glob
    from .utils import env
    import shutil
    from collections import OrderedDict
    from .target import FileTarget

    sig_files = glob.glob('.sos/*.sig')
    if not sig_files:
        sys.exit('No executed workflows.')
    tracked_files = set()
    for sig_file in sig_files:
        tracked_files |= get_tracked_files(sig_file)
    #
    tracked_files = {os.path.abspath(os.path.expanduser(x)) for x in tracked_files}
    tracked_dirs = set()
    # need to get all directories along the path
    for x in tracked_files:
        relpath = os.path.relpath(x, '.')
        if relpath.startswith('..'):
            # we do not care about tracked files outside of current directory
            continue
        # add all path to tracked file as tracked directories
        tmp = relpath
        while os.sep in tmp:
            tmp = os.path.dirname(tmp)
            tracked_dirs.add(os.path.abspath(tmp))

    if tracked_files:
        env.logger.info('{} tracked files from {} run{} are identified.'
                .format(len(tracked_files), len(sig_files), 's' if len(sig_files) > 1 else ''))
    else:
        env.logger.info('No tracked file from {} run{} are identified.'
                .format(len(sig_files), 's' if len(sig_files) > 1 else ''))
    #
    specified_tracked_files = []
    specified_tracked_dirs = []
    specified_untracked_files = []
    specified_untracked_dirs = []
    #
    if not args.targets:
        sys.exit('No files or directories to be removed.')
    # 
    for target in args.targets:
        target = os.path.expanduser(target)
        if not os.path.exists(target):
            sys.exit('Invalid file or directory to be removed: {}'.format(target))
        relpath = os.path.relpath(target, '.')
        if relpath.startswith('..'):
            # we do not care about tracked files outside of current directory
            sys.exit('Only subdirectories of the current directory can be removed. {} specified.'.format(target))
        # file
        if os.path.isfile(target):
            if os.path.abspath(target) in tracked_files:
                specified_tracked_files.append(target)
            else:
                specified_untracked_files.append(target)
            continue
        # we will not remove . itself
        elif os.path.isdir(target) and os.path.abspath(target) != os.path.abspath('.'):
            if os.path.abspath(target) in tracked_dirs:
                specified_tracked_dirs.append(target)
            else:
                specified_untracked_dirs.append(target)
        # directory
        for dirname, dirlist, filelist in os.walk(target):
            # we do not remove files under the current directory
            if dirname != '.':
                for x in filelist:
                    # ignore hidden file
                    if x.startswith('.'):
                        continue
                    if os.path.abspath(os.path.join(dirname, x)) in tracked_files:
                        specified_tracked_files.append(os.path.join(dirname, x))
                    else:
                        specified_untracked_files.append(os.path.join(dirname, x))
            # we do not track dot directories 
            dir_with_tracked_files = []
            for x in dirlist:
                # ignore hidden directories such as .git
                if x.startswith('.'):
                    continue
                if any(y.startswith(os.path.abspath(os.path.join(dirname, x))) for y in tracked_dirs):
                    dir_with_tracked_files.append(x)
                    specified_tracked_dirs.append(os.path.join(dirname, x))
                else:
                    specified_untracked_dirs.append(os.path.join(dirname, x))
            # do not scan the directory if it does not contain any tracked files because
            # they will be handled as a total directory
            dirlist[:] = dir_with_tracked_files
    #
    def dedup(_list):
        return OrderedDict((item, None) for item in _list).keys()

    def get_response(msg):
        if args.__confirm__:
            print(msg)
            return True
        while True:
            res = input('{} (y/n/a)? '.format(msg))
            if res == 'a':
                args.__confirm__ = True
                return True
            elif res == 'y':
                return True
            elif res == 'n':
                return False

    # in case of tracked or all, we need to remove signature
    if args.__tracked__ or not args.__untracked__:
        for f in specified_tracked_files:
            if args.__dryrun__:
                if args.__tracked__:
                    print('Would remove tracked file {} and its signature'.format(f))
            else:
                if get_response('Remove tracked file {} and its signature'.format(f)):
                    FileTarget(f).remove('both')
        # note: signatures of tracked files under
        # these directories should have been removed.
        for d in sorted(specified_tracked_dirs, key=len, reverse=True):
            if args.__dryrun__:
                if args.__tracked__:
                    print('Would remove {} with tracked files if empty'.format(d))
            else:
                #if os.listdir(d):
                if not os.listdir(d):
                    if get_response('Remove {} with tracked files'.format(d)):
                        if os.path.isdir(d):
                            shutil.rmtree(d)
                        else:
                            os.unlink(d)
                else:
                    print('Do not remove {} with tracked file because it is not empty'.format(d))
    elif args.__untracked__:
        for f in specified_untracked_files:
            if args.__dryrun__:
                print('Would remove untracked file {}'.format(f))
            else:
                if get_response('Remove untracked file {}'.format(f)):
                    os.remove(f)
        # note: signatures of tracked files under
        # these directories should have been removed.
        for d in specified_untracked_dirs:
            if args.__dryrun__:
                print('Would remove untracked directory {}'.format(d))
            else:
                if get_response('Remove untracked directory {}'.format(d)):
                    if os.path.isdir(d):
                        shutil.rmtree(d)
                    else:
                        os.unlink(d)
    # in case of all, we need to remove everything
    if not args.__tracked__ and not args.__untracked__:
        for target in args.targets:
            target = os.path.expanduser(target)
            if args.__dryrun__:
                print('Would remove {}'.format(target))
            elif os.path.exists(target):
                if get_response('Remove {}'.format(target)):
                    if os.path.isfile(target):
                        os.remove(target)
                    elif os.path.isdir(target):
                        if os.path.abspath(target) != os.path.abspath('.'):
                            shutil.rmtree(target)
                    else:
                        os.unlink(target)
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



