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
#
import os
import sys
import argparse
import pkg_resources

script_help = '''A SoS script that defines one or more workflows, in format
    .sos or .ipynb. The script can be a filename or a URL from which the
    content of a SoS will be read. If a valid file cannot be located or
    downloaded, SoS will search for the script, with specified name, and
    with added extension .sos and .ipynb, in a search path specified by
    variable `sos_path` defined in the global SoS configuration file
    (~/.sos/config.yml).'''
workflow_spec =  '''Name of the workflow to execute. This option can be
    ignored if the script defines a default workflow (with no name or with
    name `default`) or defines only a single workflow. A subworkflow or a
    combined workflow can also be specified, where a subworkflow executes a
    subset of workflow (`name_steps` where `steps` can be `n` (a step `n`),
    `:n` (up to step `n`), `n:m` (from step `n` to `m`), and `n:` (from step
    `n`)), and a combined workflow executes to multiple (sub)workflows
    combined by `+` (e.g. `A_0+B+C`).'''
workflow_options = '''Arbitrary parameters defined by the [parameters] step
    of the script, and [parameters] steps of other scripts if nested workflows
    are defined in other SoS files (option `source`). The name, default and
    type of the parameters are specified in the script. Single value parameters
    should be passed using option `--name value` and multi-value parameters
    should be passed using option `--name value1 value2`. '''
transcript_help = '''Name of a file that records the execution transcript of
    the script. The transcript will be recorded in dryrun and run mode but
    might differ in content because of dynamic input and output of scripts.
    If the option is specified wiht no value, the transcript will be written
    to standard error output.'''
bindir_help = '''Extra directories in which SoS will look for executables before
    standard $PATH. This option essentially prefix $PATH with these directories.
    Note that the default value '~/.sos/bin' is by convention a default
    directory for commands that are installed by SoS. You can use option '-b'
    without value to disallow commands under ~/.sos/bin.'''

#
# subcommand install
#
# def get_install_parser(desc_only=False):
#     parser = argparse.ArgumentParser('install',
#         description='''Install various component to the system''')
#     parser.short_description = '''Install various components to the system'''
#     if desc_only:
#         return parser
#     parser.set_defaults(func=cmd_install)
#     subparsers = parser.add_subparsers(title='installers', dest='installer_name')
#     for entrypoint in pkg_resources.iter_entry_points(group='sos_installers'):
#         try:
#             name = entrypoint.name
#             if not name.endswith('.parser'):
#                 continue
#             item = name.rsplit('.',1)[0]
#             subparser = add_sub_parser(subparsers, entrypoint.load()(), name=item)
#         except Exception as e:
#             print('Failed to load installer {}: {}'.format(entrypoint.name, e))
#     return parser
# 
# 
# def cmd_install(args, unknown_args):
#     from .utils import env, get_traceback
#     for entrypoint in pkg_resources.iter_entry_points(group='sos_installers'):
#         try:
#             if entrypoint.name == args.installer_name + '.func':
#                 func = entrypoint.load()
#                 func(args)
#         except Exception as e:
#             # if no other parameter, with option list all
#             if args.verbosity and args.verbosity > 2:
#                 sys.stderr.write(get_traceback())
#             env.logger.error('Failed to execute installer {}: {}'.format(entrypoint.name.rsplit('.', 1)[0], e))
#             sys.exit(1)
# 
#
# subcommand convert
#
def get_convert_parser(desc_only=False):
    parser = argparse.ArgumentParser('convert',
        description='''Converts .sos to various formats including
            .html for web display, to jupyter notebook (.ipynb), and to terminal
            for syntax highlighted viewing on terminal. It also allows converting
            from jupyter notebook (.ipynb) to sos script (.sos).''',
        epilog='''Extra command line argument could be specified to customize
            the style of html, markdown, and terminal output. ''',
        )
    parser.short_description = '''Convert between .sos, .ipynb and other formats'''
    if desc_only:
        return parser
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_convert)
    subparsers = parser.add_subparsers(title='converters (name of converter is not needed from command line)',
        dest='converter_name')
    for entrypoint in pkg_resources.iter_entry_points(group='sos_converters'):
        try:
            name = entrypoint.name
            if not name.endswith('.parser'):
                continue
            f_format, t_format = name.rsplit('.',1)[0].split('-')
            subparser = add_sub_parser(subparsers, entrypoint.load()(), name='{}-{}'.format(f_format, t_format))
            subparser.add_argument('from_file', metavar='FROM', nargs='?',
                help='''File to be converted.''')
            subparser.add_argument('to_file', metavar='TO', nargs='?',
                help='''File to convert to, default to standard output if no file
                    name is given.''')
            subparser.add_argument('-t', '--to', dest='__to_format__', metavar='TO_FORMAT',
                help='''Destination file format, which is usually determined from
                    extension of `to_file` filename, but is needed if `to_file` is
                    unspecified.''')
        except Exception as e:
            print('Failed to load converter {}: {}'.format(entrypoint.name, e))
    return parser


def get_converter_formats(argv):
    if len(argv) == 1 and '-' in argv[0] and '.' not in argv[0]:
        return argv[0].split('-', 1)
    parser = argparse.ArgumentParser('convert')
    parser.add_argument('from_file', nargs='?')
    parser.add_argument('to_file', nargs='?')
    parser.add_argument('-t', '--to', dest='__to_format__')
    args, _ = parser.parse_known_args(argv)
    from_format = args.from_file[1:] if args.from_file.startswith('.') and args.from_file.count('.') == 1 \
        else os.path.splitext(args.from_file)[-1][1:]
    if not from_format:
        from_format = 'sos'
    if args.__to_format__:
        to_format = args.__to_format__
    elif args.to_file:
        to_format = args.to_file[1:] if args.to_file.startswith('.') and args.to_file.count('.') == 1 \
            else os.path.splitext(args.to_file)[-1][1:]
    else:
        return None, None
    return from_format, to_format


def print_converter_help():
    from_format, to_format = get_converter_formats([x for x in sys.argv[2:] if x != '-h'])
    if from_format is None or to_format is None:
        return
    for entrypoint in pkg_resources.iter_entry_points(group='sos_converters'):
        try:
            name = entrypoint.name
            if not name.endswith('.parser'):
                continue
            f_format, t_format = name.rsplit('.',1)[0].split('-')
            if from_format != f_format or to_format != t_format:
                continue
            parser = entrypoint.load()()
            sys.exit(parser.print_help())
        except Exception as e:
            sys.exit('Failed to load converter {}: {}'.format(entrypoint.name, e))


def cmd_convert(args, unknown_args):
    from .utils import env, get_traceback
    for entrypoint in pkg_resources.iter_entry_points(group='sos_converters'):
        try:
            if entrypoint.name == args.converter_name + '.func':
                func = entrypoint.load()
                func(args.from_file, args.to_file, args, unknown_args)
        except Exception as e:
            # if no other parameter, with option list all
            if args.verbosity and args.verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error('Failed to execute converter {}: {}'.format(entrypoint.name.rsplit('.', 1)[0], e))
            sys.exit(1)

#
# subcommand run
#
def get_run_parser(interactive=False, with_workflow=True, desc_only=False):
    parser = argparse.ArgumentParser(prog='run',
        description='Execute default or specified workflow defined in script',
        epilog=workflow_options)
    parser.short_description = 'Execute default or specified workflow in script'
    if desc_only:
        return parser
    if not interactive:
        parser.add_argument('script', metavar='SCRIPT', help=script_help)
    if with_workflow:
        parser.add_argument('workflow', metavar='WORKFLOW', nargs='?',
            help=workflow_spec)
    parser.add_argument('-j', type=int, metavar='JOBS',
        dest='__max_procs__', default=max(os.cpu_count() // 2, 1),
        help='''Maximum number of worker processes for the execution of steps in
            a workflow and input groups in a step (with input option concurrent),
            default to half of number of CPUs''')
    parser.add_argument('-J', type=int, metavar='EXTERNAL_JOBS',
        dest='__max_running_jobs__',
        help='''Maximum number of externally running tasks. This option
            overrides option "max_running_jobs" of a task queue (option -q)
            so that you can, for example, submit one job at a time (with
            -J 1) to test the task queue.''')
    parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE',
        help='''A configuration file in the format of YAML/JSON. The content
            of the configuration file will be available as a dictionary
            CONF in the SoS script being executed.''')
    parser.add_argument('-t', dest='__targets__', metavar='FILE', default=[],
        nargs='+', help='''One of more files or alias of other targets that
            will be the target of execution. If specified, SoS will execute
            only part of a workflow or multiple workflows or auxiliary steps
            to generate specified targets.''')
    parser.add_argument('-b', dest='__bin_dirs__', nargs='*', metavar='BIN_DIR',
        default=['~/.sos/bin'], help=bindir_help)
    parser.add_argument('-q', dest='__queue__', nargs='?', const='', metavar='QUEUE',
        help='''host (server) or job queues to execute all tasks in the
            workflow. The queue can be defined in global or local sos
            configuration file, or a file specified by option  --config. A host is
            assumed to be a remote machine with process type if no configuration
            is found. If this option is specified without value, SoS will list
            all configured queues and exit.''')
    parser.add_argument('-w', dest='__wait__', action='store_true',
        help='''Wait for the completion of external tasks regardless of the
            setting of individual task queue.''')
    parser.add_argument('-W', dest='__no_wait__', action='store_true',
        help='''Do not wait for the completion of external tasks and quit SoS
            if all tasks are being executed by external task queues. This option
            overrides the default wait setting of task queues.''')
    parser.add_argument('-r', dest='__remote__', metavar='HOST', nargs='?', const='',
        help='''Execute the workflow in specified remote host, which should
            be defined under key host of sos configuration files (preferrably
            in ~/.sos/hosts.yml). This option basically copy the workflow
            to remote host and invoke sos command there. No path translation
            and input/output file synchronization will be performed before or
            after the execution of the workflow. If this option is specified without
            value, SoS will list all configured queues and exit.''')
    #parser.add_argument('-r', dest='__remote__', action='store_true',
    #    help='''Forcing all targets specified in input, output, and
    #        depends are remote targets so that they are not synchronized
    #        between local and remote hosts.''')
    #parser.add_argument('-t', dest='__transcript__', nargs='?',
    #    metavar='TRANSCRIPT', const='__STDERR__', help=transcript_help)
    runmode = parser.add_argument_group(title='Run mode options',
        description='''Control how sos scirpt is executed.''')
    runmode.add_argument('-n', action='store_true', dest='dryrun',
        help='''Execute a workflow in dryrun mode. Please check command
        sos dryrun for details of the dryrun mode.''')
    runmode.add_argument('-s', choices=['default', 'ignore', 'force', 'build', 'assert'],
        default='default', metavar='SIGMODE',
        dest='__sig_mode__',
        help='''How runtime signature would be handled, which can be "default"
            (save and use signature, default mode in batch mode), "ignore"
            (ignore runtime signature, default mode in interactive mode),
            "force" (ignore existing signature and overwrite them while
            executing the workflow), "build" (build new or overwrite
            existing signature from existing environment and output files), and
            "assert" for validating existing files against their signatures.
            Please refer to online documentation for details about the
            use of runtime signatures.''')
    output = parser.add_argument_group(title='Output options',
        description='''Output of workflow''')
    output.add_argument('-d', nargs='?', default='', metavar='DAG', dest='__dag__',
        help='''Output Direct Acyclic Graph (DAGs) in graphiviz .dot format. An
            exntesion of ".dot" would be added automatically. Because DAG could
            change during the execution of workflow, multiple DAGs could be
            outputed with names $FILE_1.dot, $FILE_2.dot. If this option is
            specified without a name, the DAG would be wrritten to the standard
            output.''')
    output.add_argument('-v', dest='verbosity', type=int, choices=range(5),
        default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_run)
    return parser


def cmd_run(args, workflow_args):
    #import multiprocessing as mp
    # #562, #558, #493
    #
    #if sys.platform != 'win32':
    #    mp.set_start_method('forkserver')

    from .utils import env, get_traceback, load_config_files
    from .parser import SoS_Script

    if args.__queue__ == '':
        cfg = load_config_files(args.__config__)
        from .hosts import list_queues
        list_queues(cfg, args.verbosity)
        return

    if args.__wait__ and args.__no_wait__:
        sys.exit('Please specify only one of -w (wait) or -W (no-wait)')

    if args.__remote__ is not None:
        if args.__remote__ == '':
            cfg = load_config_files(args.config)
            from .hosts import list_queues
            list_queues(cfg, args.verbosity)
            return

        # if executing on a remote host...
        from .hosts import Host
        cfg = load_config_files(args.__config__)
        host = Host(args.__remote__)
        from .utils import remove_arg
        #
        # copy script to remote host...
        host.send_to_host(args.script)

        argv = remove_arg(sys.argv, '-r')
        # -c only point to local config file.
        argv = remove_arg(argv, '-c')
        # replace absolute path with relative one because remote sos might have
        # a different path.
        if os.path.basename(argv[0]) == 'sos':
            argv[0] = 'sos'
        elif os.path.basename(argv[0]) == 'sos-runner':
            argv[0] = 'sos-runner'
        # execute the command on remote host
        sys.exit(host._host_agent.check_call(argv))


    # '' means no -d
    if args.__dag__ is None:
        args.__dag__ = '-'
    elif args.__dag__ == '':
        args.__dag__ = None
    env.verbosity = args.verbosity

    from .workflow_executor import Base_Executor

    if args.__bin_dirs__:
        for d in args.__bin_dirs__:
            if d == '~/.sos/bin' and not os.path.isdir(os.path.expanduser(d)):
                os.makedirs(os.path.expanduser(d), exist_ok=True)
            elif not os.path.isdir(os.path.expanduser(d)):
                raise ValueError('directory does not exist: {}'.format(d))
        os.environ['PATH'] = os.pathsep.join([os.path.expanduser(x) for x in args.__bin_dirs__]) + os.pathsep + os.environ['PATH']

    try:
        # workflow args has to be in the format of --arg, not positional, not -a
        if workflow_args and not workflow_args[0].startswith('--'):
            raise ValueError("Unrecognized command line option {}".format(' '.join(workflow_args)))
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow, use_default=not args.__targets__)
        executor = Base_Executor(workflow, args=workflow_args, config={
                'config_file': args.__config__,
                'output_dag': args.__dag__,
                # wait if -w or in dryrun mode, not wait if -W, otherwise use queue default
                'wait_for_task': True if args.__wait__ is True or args.dryrun else (False if args.__no_wait__ else None),
                'default_queue': '' if args.__queue__ is None else args.__queue__,
                'max_procs': args.__max_procs__,
                'max_running_jobs': args.__max_running_jobs__,
                'sig_mode': args.__sig_mode__,
                'run_mode': 'dryrun' if args.dryrun else 'run',
                'resume_mode': getattr(args, '__resume__', False),
                'verbosity': args.verbosity,
                # for infomration and resume only
                'workdir': os.getcwd(),
                'script': args.script,
                'workflow': args.workflow,
                'targets': args.__targets__,
                'bin_dirs': args.__bin_dirs__,
                'workflow_args': workflow_args
                })
        executor.run(args.__targets__, mode='dryrun' if args.dryrun else 'run')
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)


def get_resume_parser(interactive=False, with_workflow=True, desc_only=False):
    parser = argparse.ArgumentParser(prog='resume',
        description='''Resume the execution or list status of suspended workflows''')
    parser.short_description = 'Resume the execution of a suspended workflow'
    if desc_only:
        return parser
    parser.add_argument('workflow_id', nargs='?',
        help='''First few characters of the ID of a workflow as long as it
            uniquely identifies the workflow. The last executed workflow will
            be resumed if no workflow is specified.''')
    parser.add_argument('-s', '--status', action='store_true',
        help='''Return the status of all or specified workflows without
            resuming them. This option will list status of all pending
            tasks.''')
    parser.add_argument('-w', dest='__wait__', action='store_true',
        help='''Wait for the completion of external tasks regardless of the
            setting of individual task queue.''')
    parser.add_argument('-W', dest='__no_wait__', action='store_true',
        help='''Do not wait for the completion of external tasks and quit SoS
            if all tasks are being executed by external task queues. This option
            overrides the default wait setting of task queues.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.add_argument('-j', type=int, metavar='JOBS',
        default=4, dest='__max_procs__',
        help='''Maximum number of worker processes for the execution of the
            workflow if the workflow can be executed in parallel (namely
            having multiple starting points or execution branches).''')
    parser.add_argument('-J', type=int, metavar='EXTERNAL_JOBS',
        dest='__max_running_jobs__',
        help='''Maximum number of externally running tasks. This option
            overrides option "max_running_jobs" of a task queue (option -q)
            so that you can, for example, submit one job at a time (with
            -J 1) to test the task queue.''')
    parser.add_argument('-r', dest='__remote__', nargs='?', const='',
        help='''Resume workflow that was executed on remote host (sos run with
            option -r)''')
    parser.set_defaults(func=cmd_resume)
    return parser

def workflow_status(workflow):
    import time
    from .utils import env, load_var, load_config_files, PrettyRelativeTime
    from .hosts import Host
    from .tasks import check_tasks
    from .eval import interpolate
    from io import StringIO
    from contextlib import redirect_stdout
    from collections import defaultdict
    import re
    pending_tasks = defaultdict(list)
    res = {'task_status': []}
    with open(workflow) as wf:
        config_file = None
        for line in wf:
            if line.startswith('#') or not line.strip():
                continue
            try:
                k, v = load_var(line)
                res[k] = v
                if k == 'config_file':
                    config_file = v
                elif k == 'pending_task':
                    pending_tasks[v[0]].append(v[1])
            except Exception as e:
                raise ValueError('Unrecognizable status line {}: {}'.format(line, e))
    if 'script' not in res:
        env.logger.error('Cannot resume a workflow with script file (it must have been started programmatically with content of a script).')
        sys.exit(1)
    #
    env.logger.info('{:15s} \t{}'.format('Workflow ID:', os.path.basename(workflow)[:-7]))
    env.logger.info('{:15s} \t{}'.format('Command:', re.sub(r'\s+', ' ', interpolate(
        'sos run {script} {workflow if workflow else ""} '
        '{("-c " + config_file) if config_file else ""} '
        '{("-s " + sig_mode) if sig_mode not in ("", "default") else ""} '
        '{("-q " + default_queue) if default_queue else ""} '
        '{("-d " + output_dag) if output_dag not in ("", None) else ""} '
        '{("-b " + " ".join(bin_dirs)) if bin_dirs and bin_dirs != ["~/.sos/bin"] else ""} '
        '{("-j " + str(max_procs)) if max_procs != 4 else ""} '
        '{("-J " + str(max_running_jobs)) if max_running_jobs else ""} '
        '{("-t " + " ".join(targets)) if targets else ""} '
        '{" ".join(workflow_args)} '
        , res))))
    env.logger.info('{:15s} \t{} ago'.format('Started:',
                PrettyRelativeTime(time.time() - os.path.getmtime(workflow))))
    env.logger.info('{:15s} \t{}'.format('Working dir:', res['workdir']))
    #
    for k,v in pending_tasks.items():
        if k in ('', 'localhost'):
            with StringIO() as buf, redirect_stdout(buf):
                check_tasks(v, 0, False, False, None)
                status = buf.getvalue().strip().split('\n')
        else:
            # remote host?
            load_config_files(config_file)
            try:
                host = Host(k)
                status = host._task_engine.query_tasks(v, 0, False, False, None).strip().split('\n')
            except Exception as e:
                env.logger.warning('Failed to check status of task {} at host {}'.format(v, k))
                status = ['unknown'] * len(v)
        for v,s in zip(v, status):
            env.logger.info('{:15s} \t{} at {}, currently ``{}``'.format('Pending task:', v, k, s))
        res['task_status'].extend(status)
    return res

def cmd_resume(args, workflow_args):
    if workflow_args:
        sys.exit('No additional parameter is allowed for command resume: {} provided'.format(workflow_args))

    if args.__remote__ is not None:
        if args.__remote__ == '':
            from .utils import load_config_files
            from .hosts import list_queues
            cfg = load_config_files()
            list_queues(cfg, args.verbosity)
            return

        # resume executing on a remote host...
        from .hosts import Host
        host = Host(args.__remote__)
        #
        r_idx = [idx for idx, x in enumerate(sys.argv) if x.startswith('-r')][0]
        if sys.argv[r_idx] == '-r':
            # in case of -r host
            argv = sys.argv[:r_idx] + sys.argv[r_idx+2:]
        else:
            # in case of -r=host...
            argv = sys.argv[:r_idx] + sys.argv[r_idx+1:]
        # replace absolute path with relative one because remote sos might have
        # a different path.
        if os.path.basename(argv[0]) == 'sos':
            argv[0] = 'sos'
        # execute the command on remote host
        sys.exit(host._host_agent.check_call(argv))

    import glob
    import time
    from .utils import env, PrettyRelativeTime
    env.verbosity = args.verbosity

    workflows = glob.glob(os.path.join(os.path.expanduser('~'), '.sos', '{}*.status').format(args.workflow_id if args.workflow_id else ''))
    if not workflows:
        env.logger.info('No resumable workflow')
        sys.exit(0)
    #
    if args.status:
        for wf in workflows:
            workflow_status(wf)
        sys.exit(0)
    elif len(workflows) > 1:
        workflows = sorted(workflows, key=os.path.getmtime)
        for wf in workflows:
            env.logger.info('{}\tstarted {} ago'.format(os.path.basename(wf)[:-7],
                PrettyRelativeTime(time.time() - os.path.getmtime(wf))))
    #
    # resume execution...
    if args.workflow_id and len(workflows) > 1:
        env.logger.error('{} matches more than one resumable workflows {}'.format(args.workflow_id,
            ', '.join([os.path.basename(x)[:-4] for x in workflows])))
        sys.exit(1)
    else:
        workflow = workflows[-1]
    #
    res = workflow_status(workflow)
    if not res['task_status']:
        env.logger.warn('Removing workflow {} because it does not have any pending task. The workflow might have been interrupted.'.format(os.path.basename(workflow)[:-4]))
        os.remove(workflow)
        sys.exit(1)
    if all(x == 'running' for x in res['task_status']) and args.__wait__ is not True:
        env.logger.info('Cannot resume workflow {} because all tasks are still running'.format(
            os.path.basename(workflow)[:-7]))
        sys.exit(0)
    #
    args.__config__ = res['config_file']
    args.__sig_mode__ = res['sig_mode']
    args.__max_procs__ = args.__max_procs__ if args.__max_procs__ != 4 else res['max_procs']
    args.__resume__ = True
    args.__max_running_jobs__ = args.__max_running_jobs__ if args.__max_running_jobs__ is not None else res['max_running_jobs']
    args.dryrun = False
    args.__wait__ = args.__wait__ if args.__wait__ is True else None
    args.__no_wait__ = args.__no_wait__ if args.__no_wait__ is True else None
    args.__bin_dirs__ = res['bin_dirs']
    args.__queue__ = None if res['default_queue'] == '' else res['default_queue']
    args.__dag__ = None if res['output_dag'] == '-' else ('' if res['output_dag'] is None else res['output_dag'])
    args.__targets__ = res['targets']
    args.script = res['script']
    args.workflow = res['workflow']
    if 'workdir' in res:
        os.chdir(res['workdir'])
    cmd_run(args, res['workflow_args'])

#
# subcommand dryrun
#
def get_dryrun_parser(desc_only=False):
    parser = argparse.ArgumentParser('dryrun',
        description='''Execute workflow in dryrun mode. This mode is identical
        to run mode except that 1). Actions might behavior differently. In
        particular, script-running steps would print instead of execute script.
        2). Steps will generate empty output files if specified output do not
        exist after execution. 3). Signature mode is set to ignore. 4). Option
        -q is ignored so all tasks are executed locally. 5). Tasks are generated
        but not executed.''',
        epilog=workflow_options)
    parser.short_description = '''Execute workflow in dryrun mode'''
    if desc_only:
        return parser
    parser.add_argument('script', metavar='SCRIPT', help=script_help)
    parser.add_argument('workflow', metavar='WORKFLOW', nargs='?',
        help=workflow_spec)
    parser.add_argument('-c', dest='__config__', metavar='CONFIG_FILE',
        help='''A configuration file in the format of YAML/JSON. The content
            of the configuration file will be available as a dictionary
            CONF in the SoS script being executed.''')
    parser.add_argument('-t', dest='__targets__', metavar='FILES', default=[],
        nargs='+', help='''One of more files or alias of other targets that
            will be the target of execution. If specified, SoS will execute
            only part of a workflow or multiple workflows or auxiliary steps
            to generate specified targets. ''')
    parser.add_argument('-q', dest='__queue__', nargs='?', const='', metavar='QUEUE',
        help='''host (server) or job queues to execute all tasks in the
            workflow. The queue can be defined in global or local sos
            configuration file, or a file specified by option  --config. A host is
            assumed to be a remote machine with process type if no configuration
            is found. If this option is specified without value, SoS will list all
            configured queues and exit.''')
    output = parser.add_argument_group(title='Output options',
        description='''Output of workflow''')
    output.add_argument('-d', nargs='?', default='', metavar='DAG', dest='__dag__',
        help='''Output Direct Acyclic Graph (DAGs) in graphiviz .dot format. An
            exntesion of ".dot" would be added automatically. Because DAG could
            change during the execution of workflow, multiple DAGs could be
            outputed with names $FILE_1.dot, $FILE_2.dot. If this option is
            specified without a name, the DAG would be wrritten to the standard
            output.''')
    output.add_argument('-v', dest='verbosity', type=int, choices=range(5),
        default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_dryrun)
    return parser

def cmd_dryrun(args, workflow_args):
    args.__sig_mode__ = 'ignore'
    args.__max_procs__ = 1
    args.__max_running_jobs__ = 1
    args.dryrun = True
    args.__wait__ = True
    args.__no_wait__ = False
    args.__bin_dirs__ = []
    args.__remote__ = None
    cmd_run(args, workflow_args)

#
# subcommand push
#
def get_push_parser(desc_only=False):
    parser = argparse.ArgumentParser('push',
        description='''Push local files or directory to a remote host''')
    if desc_only:
        return parser
    parser.add_argument('items', nargs='+', help='''Files or directories to be sent
        to remote host. The location of remote files are determined by "path_map"
        determined by "paths" definitions of local and remote hosts.''')
    parser.add_argument('-t', '--to', dest='host', nargs='?', const='',
        help='''Remote host to which the files will be sent. SoS will list all configured
        queues and exit''')
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_push)
    return parser


def cmd_push(args, workflow_args):
    from .utils import env, load_config_files
    from .hosts import Host
    env.verbosity = args.verbosity
    cfg = load_config_files(args.config)
    if args.host == '':
        from .hosts import list_queues
        list_queues(cfg, args.verbosity)
        return
    try:
        host = Host(args.host)
        #
        sent = host.send_to_host(args.items)
        #
        print('{} item{} sent:\n{}'.format(len(sent),
            ' is' if len(sent) <= 1 else 's are',
            '\n'.join(['{} => {}'.format(x, sent[x]) for x in sorted(sent.keys())])))
    except Exception as e:
        from .utils import get_traceback
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)
#
# subcommand pull
#
def get_pull_parser(desc_only=False):
    parser = argparse.ArgumentParser('pull',
        description='''Pull files or directories from remote host to local host''')
    if desc_only:
        return parser
    parser.add_argument('items', nargs='+', help='''Files or directories to be
        retrieved from remote host. The files should be relative to local file
        system. The files to retrieve are determined by "path_map"
        determined by "paths" definitions of local and remote hosts.''')
    parser.add_argument('-f', '--from', dest='host', nargs='?', const='',
        help='''Remote host to which the files will be sent. SoS will list all configured
        queues and exit''')
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_pull)
    return parser


def cmd_pull(args, workflow_args):
    from .utils import env, load_config_files
    from .hosts import Host
    env.verbosity = args.verbosity
    cfg = load_config_files(args.config)
    if args.host == '':
        from .hosts import list_queues
        list_queues(cfg, args.verbosity)
        return
    try:
        host = Host(args.host)
        #
        received = host.receive_from_host(args.items)
        #
        print('{} item{} received:\n{}'.format(len(received),
            ' is' if len(received) <= 1 else 's are',
            '\n'.join(['{} <= {}'.format(x, received[x]) for x in sorted(received.keys())])))
    except Exception as e:
        from .utils import get_traceback
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# command preview
#
def get_preview_parser(desc_only=False):
    parser = argparse.ArgumentParser(prog='preview',
        description='''Preview files, sos variables, or expressions in the
            side panel, or notebook if side panel is not opened, unless
            options --panel or --notebook is specified.''')
    parser.short_description = '''Preview files on local or remote host'''
    if desc_only:
        return parser
    parser.add_argument('items', nargs='*',
        help='''Filename, variable name, or expression. Wildcard characters
            such as '*' and '?' are allowed for filenames.''')
    # this option is currently hidden
    parser.add_argument('-s', '--style', choices=['table', 'scatterplot', 'png'],
        help='''Option for preview file or variable, which by default is "table"
        for Pandas DataFrame. The %%preview magic also accepts arbitrary additional
        keyword arguments, which would be interpreted by individual style. Passing
        '-h' with '--style' would display the usage information of particular
        style.''')
    parser.add_argument('-r', '--host', dest='host', metavar='HOST', nargs='?', const='',
        help='''Preview files on specified remote host, which should
        be defined under key host of sos configuration files (preferrably
        in ~/.sos/hosts.yml). If this option is specified without
        value, SoS will list all configured queues and exit.''')
    parser.add_argument('--html', action='store_true',
        help=argparse.SUPPRESS)
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_preview)
    return parser


def preview_file(previewers, filename, style=None):
    from .utils import pretty_size
    from IPython.core.display import HTML
    msg = []
    if not os.path.isfile(filename):
        msg.append(['stream', {
                        'name': 'stderr',
                        'text': '\n> ' + filename + ' does not exist'}])
        return msg
    msg.append(['display_data',
         {'metadata': {},
         'data': {
             'text/plain': '\n> {} ({}):'.format(filename, pretty_size(os.path.getsize(filename))),
             'text/html': HTML('<div class="sos_hint">> {} ({}):</div>'.format(filename, pretty_size(os.path.getsize(filename)))).data,
            }
         }])
    previewer_func = None
    # lazy import of previewers
    import fnmatch
    for x, y, _ in previewers:
        if isinstance(x, str):
            if fnmatch.fnmatch(os.path.basename(filename), x):
                # we load entrypoint only before it is used. This is to avoid
                # loading previewers that require additional external modules
                # we can cache the loaded function but there does not seem to be
                # a strong performance need for this.
                previewer_func = y.load()
                break
        else:
            # it should be a function
            try:
                if x(filename):
                    try:
                        previewer_func = y.load()
                    except Exception as e:
                        msg.append(['stream', {
                            'name': 'stderr',
                            'text': 'Failed to load previewer {}: {}'.format(y, e) }])
                        continue
                    break
            except Exception as e:
                msg.append(['stream', {
                            'name': 'stderr',
                            'text': str(e)}])
                continue
    #
    # if no previewer can be found
    if previewer_func is None:
        return msg
    try:
        result = previewer_func(filename, None, style)
        if not result:
            return msg
        if isinstance(result, str):
            msg.append(['stream',
                {'name': 'stdout', 'text': result}])
        elif isinstance(result, dict):
            msg.append(['display_data',
                {'source': filename, 'data': result, 'metadata': {}}])
        elif isinstance(result, [list, tuple]) and len(result) == 2:
            msg.append(['display_data',
                {'source': filename, 'data': result[0], 'metadata': result[1]}])
        else:
            msg.append(['stream', {
                'name': 'stderr',
                'text': 'Unrecognized preview content: {}'.format(result)}])
    except Exception as e:
        msg.append(['stream', {
            'name': 'stderr',
            'text': 'Failed to preview {}: {}'.format(filename, e)}])
    return msg


def cmd_preview(args, unknown_args):
    from .utils import env, load_config_files
    from .hosts import Host
    cfg = load_config_files(args.config)
    env.verbosity = args.verbosity
    if args.host == '':
        from .hosts import list_queues
        list_queues(cfg, args.verbosity)
        return
    if args.host:
        # remote host?
        host = Host(args.host)
        rargs = ['sos', 'preview'] + args.items + ['--html']
        if args.style:
            rargs += ['-s', args.style] + unknown_args
        env.logger.debug('Running "{}"'.format(' '.join(rargs)))
        msgs = eval(host._host_agent.check_output(rargs))
    else:
        from .preview import get_previewers
        previewers = get_previewers()
        msgs = []
        style = {'style': args.style, 'options': unknown_args } if args.style else None
        for filename in args.items:
            msgs.extend(preview_file(previewers, filename, style))
    if args.html:
        print(msgs)
    else:
        from .utils import colorstr, dehtml
        for msg in msgs:
            if msg[0] == 'stream':
                if msg[1]['name'] == 'stdout':
                    print(msg[1]['text'])
                else:
                    print(colorstr(msg[1]['text'], 'PURPLE'))
            elif msg[0] == 'display_data':
                if 'text/plain' in msg[1]['data']:
                    print(msg[1]['data']['text/plain'])
                elif 'text/HTML' in msg[1]['data']:
                    print(dehtml(msg[1]['data']['text/html']))
                else:
                    print('BINARY DATA of type {}'.format(', '.join(msg[1]['data'].keys())))
            else:
                raise RuntimeError('Unrecognized preview output: {}'.format(msg))
    # exit with code 1 if error happens
    sys.exit(1 if any(msg[1]['name'] == 'stderr' for msg in msgs if msg[0] == 'stream') else 0)

#
# subcommand execute
#
def get_execute_parser(desc_only=False):
    parser = argparse.ArgumentParser('execute',
        description='''Execute a packages task''')
    if desc_only:
        return parser
    parser.add_argument('tasks', nargs='+', help='''IDs of the task.''')
    parser.add_argument('-s', choices=['default', 'ignore', 'force', 'build', 'assert'],
        default='default', metavar='SIGMODE',
        dest='__sig_mode__',
        help='''How runtime signature would be handled, which can be "default"
            (save and use signature, default mode in batch mode), "ignore"
            (ignore runtime signature, default mode in interactive mode),
            "force" (ignore existing signature and overwrite them while
            executing the workflow), "build" (build new or overwrite
            existing signature from existing environment and output files), and
            "assert" for validating existing files against their signatures.
            Please refer to online documentation for details about the
            use of runtime signatures.''')
    parser.add_argument('-v', dest='verbosity', type=int, choices=range(5),
        default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.add_argument('-n', '--dryrun', action='store_true', dest='dryrun',
        help='''Dryrun mode, which will cause actions to print scripts instead
            of executing them.''')
    parser.add_argument('-q', '--queue', nargs='?', const='',
        help='''Check the status of job on specified tasks queue or remote host
        if the tasks . The queue can be defined in global or local sos
        configuration file, or a file specified by option  --config. A host is
        assumed to be a remote machine with process type if no configuration
        is found. If this option is specified without value, SoS will list all
        configured queues and exit.''')
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-w', '--wait', action='store_true', help='''Wait for the
        completion of the task, and retrieve job results if needed after the
        completion of the task. This option is only valid with the specification
        of the -q option.''')
    parser.set_defaults(func=cmd_execute)
    return parser


def cmd_execute(args, workflow_args):
    from .tasks import execute_task, check_task, monitor_interval, resource_monitor_interval
    from .monitor import summarizeExecution
    from .utils import env, load_config_files
    import glob
    if args.queue is None:
        # local machine ...
        exit_code = []
        for task in args.tasks:
            #
            matched = [os.path.basename(x)[:-5] for x in glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '*.task'))]
            if not matched:
                env.logger.error('{} does not match any existing task'.format(task))
                exit_code.append(1)
                continue
            elif len(matched) > 1:
                env.logger.error('"{}" matches more than one task ID {}'.format(task, ', '.join(matched)))
                exit_code.append(1)
                continue
            else:
                task = matched[0]
            # this is for local execution, perhaps on a remote host, and
            # there is no daemon process etc. It also does not handle job
            # preparation.
            status = check_task(task)
            res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.res')
            if status == 'running':
                if args.verbosity <= 1:
                    print(status)
                else:
                    print(summarizeExecution(task, status=status))
                exit_code.append(1)
                continue
            if status == 'completed' and args.__sig_mode__ != 'force':
                # touch the result file, this will effective change task
                # status from completed-old to completed
                os.utime(res_file, None)
                #if args.verbosity <= 1:
                env.logger.info('{} ``already completed``'.format(task))
                with open(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task + '.err'), 'a') as err:
                    err.write('{} already completed'.format(task))
                #else:
                #    print(summarizeExecution(task, status=status))
                exit_code.append(0)
                continue
            #
            if os.path.isfile(res_file):
                os.remove(res_file)
            exit_code.append(execute_task(task, verbosity=args.verbosity, runmode='dryrun' if args.dryrun else 'run',
                sigmode=args.__sig_mode__,
                monitor_interval=monitor_interval, resource_monitor_interval=resource_monitor_interval))
        sys.exit(sum(exit_code))
    elif args.queue == '':
        cfg = load_config_files(args.config)
        from .hosts import list_queues
        list_queues(cfg, args.verbosity)
        return
    # with queue definition
    from .hosts import Host
    import time
    # this is for local execution using a task queue. The task queue
    # will prepare the task, sync files, and execute this command remotely
    # if needed.
    cfg = load_config_files(args.config)
    env.verbosity = args.verbosity
    env.config['sig_mode'] = args.__sig_mode__
    env.config['run_mode'] = 'dryrun' if args.dryrun else 'run'
    host = Host(args.queue)
    for task in args.tasks:
        host.submit_task(task)
    failed_tasks = set()
    while True:
        res = host.check_status(args.tasks)
        if any(x in ('failed', 'aborted', 'signature-mismatch') for x in res):
            for t, s in zip(args.tasks, res):
                if s in ('failed', 'aborted', 'signature-mismatch') and t not in failed_tasks:
                    env.logger.warning('{} ``{}``'.format(t, s))
                    failed_tasks.add(t)
            if all(x in ('completed', 'failed', 'aborted', 'signature-mismatch') for x in res):
                raise RuntimeError('{} completed, {} failed, {} aborted, {} signature-mismatch)'.format(
                    len([x for x in res if x == 'completed']), len([x for x in res if x=='failed']),
                    len([x for x in res if x.startswith('aborted')]), len([x for x in res if x=='signature-mismatch'])))
        if all(x == 'completed' for x in res):
            env.logger.debug('Put results for {}'.format(args.tasks))
            res = host.retrieve_results(args.tasks)
            return
        elif all(x != 'pending' for x in res) and not args.wait:
            return
        elif any(x in ('pending', 'running', 'submitted') for x in res):
            continue
        else:
            raise RuntimeError('Job returned with status {}'.format(res))
        time.sleep(0.01)

#
# command status
#
def get_status_parser(desc_only=False):
    parser = argparse.ArgumentParser('status',
        description='''Check the status of specified tasks''')
    if desc_only:
        return parser
    parser.add_argument('tasks', nargs='*', help='''ID of the task. All tasks
        will be checked if unspecified. There is no need to specify compelete
        task IDs because SoS will match specified name with tasks starting with
        these names.''')
    parser.add_argument('-q', '--queue', nargs='?', const='',
        help='''Check the status of job on specified tasks queue or remote host
        if the tasks . The queue can be defined in global or local sos
        configuration file, or a file specified by option  --config. A host is
        assumed to be a remote machine with process type if no configuration
        is found. If this option is specified without value, SoS will list all
        configured queues and exit.''')
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-v', dest='verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.add_argument('-t', '--tags', nargs='*', help='''Only list tasks with
        one of the specified tags.''')
    parser.add_argument('-s', '--status', nargs='*', help='''Display tasks with
        one of the specified status.''')
    parser.add_argument('--age', help='''Limit to tasks that are created more than
        (default) or within specified age. Value of this parameter can be in units
        s (second), m (minute), h (hour), or d (day, default), or in the foramt of
        HH:MM:SS, with optional prefix + for older (default) and - for newer than
        specified age.''')
    parser.add_argument('--html', action='store_true',
        help='''Output results in HTML format. This option will override option
            verbosity and output detailed status information in HTML tables and
            figures.''')
    parser.add_argument('--start-time', action='store_true',
        help=argparse.SUPPRESS)
    parser.set_defaults(func=cmd_status)
    return parser


def cmd_status(args, workflow_args):
    from .tasks import check_tasks
    from .utils import env, load_config_files, get_traceback
    from .hosts import Host
    #from .monitor import summarizeExecution
    env.verbosity = args.verbosity
    try:
        cfg = load_config_files(args.config)
        if args.queue == '':
            from .hosts import list_queues
            list_queues(cfg, args.verbosity, check_status=True)
            return
        if not args.queue:
            check_tasks(tasks=args.tasks, verbosity=args.verbosity, html=args.html, start_time=args.start_time,
                    age=args.age, tags=args.tags, status=args.status)
        else:
            # remote host?
            host = Host(args.queue)
            print(host._task_engine.query_tasks(tasks=args.tasks, verbosity=args.verbosity, html=args.html,
                start_time=args.start_time, age=args.age, tags=args.tags, status=args.status))
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# command purge
#
def get_purge_parser(desc_only=False):
    parser = argparse.ArgumentParser('purge',
        description='''Remove local or remote tasks''')
    if desc_only:
        return parser
    parser.add_argument('tasks', nargs='*', help='''ID of the tasks to be removed.
        There is no need to specify compelete task IDs because SoS will match specified
        name with tasks starting with these names. If no task ID is specified,
        all tasks related to specified workflows (option -w) will be removed.''')
    parser.add_argument('-a', '--all', action='store_true',
        help='''Clear all task information on local or specified remote task queue,
        including tasks created by other workflows.''')
    parser.add_argument('--age', help='''Limit to tasks that are created more than
        (default) or within specified age. Value of this parameter can be in units
        s (second), m (minute), h (hour), or d (day, default), or in the foramt of
        HH:MM:SS, with optional prefix + for older (default) and - for newer than
        specified age.''')
    parser.add_argument('-s', '--status', nargs='+', help='''Only remove tasks with
        specified status, which can be pending, submitted, running, completed, failed,
        aborted, and signature-mismatch. One of more status can be specified.''')
    parser.add_argument('-t', '--tags', nargs='*', help='''Only remove tasks with
        one of the specified tags.''')
    parser.add_argument('-q', '--queue', nargs='?', const='',
        help='''Remove tasks on specified tasks queue or remote host
        if the tasks . The queue can be defined in global or local sos
        configuration file, or a file specified by option  --config. A host is
        assumed to be a remote machine with process type if no configuration
        is found. If this option is specified without value, SoS will list all
        configured queues and exit.''')
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-v', dest='verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_purge)
    return parser


def cmd_purge(args, workflow_args):
    from .tasks import purge_tasks
    from .utils import env, load_config_files, get_traceback
    from .hosts import Host
    #from .monitor import summarizeExecution
    env.verbosity = args.verbosity
    try:
        if args.queue == '':
            cfg = load_config_files(args.config)
            from .hosts import list_queues
            list_queues(cfg, args.verbosity)
            return
        if not args.queue:
            purge_tasks(args.tasks, args.all, args.age, args.status, args.tags, args.verbosity)
        else:
            # remote host?
            cfg = load_config_files(args.config)
            host = Host(args.queue)
            print(host._task_engine.purge_tasks(args.tasks, args.all, args.age, args.status, args.tags, args.verbosity))
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# command kill
#
#
def get_kill_parser(desc_only=False):
    parser = argparse.ArgumentParser('kill',
        description='''Stop the execution of running task''')
    if desc_only:
        return parser
    parser.add_argument('tasks', nargs='*', help='''IDs of the tasks
        that will be killed. There is no need to specify compelete task IDs because
        SoS will match specified name with tasks starting with these names.''')
    parser.add_argument('-a', '--all', action='store_true',
        help='''Kill all tasks in local or specified remote task queue''')
    parser.add_argument('-q', '--queue', nargs='?', const='',
        help='''Kill jobs on specified tasks queue or remote host
        if the tasks . The queue can be defined in global or local sos
        configuration file, or a file specified by option  --config. A host is
        assumed to be a remote machine with process type if no configuration
        is found. If this option is specified without value, SoS will list all
        configured queues and exit.''')
    parser.add_argument('-t', '--tags', nargs='*', help='''Only kill tasks with
        one of the specified tags.''')
    parser.add_argument('-c', '--config', help='''A configuration file with host
        definitions, in case the definitions are not defined in global sos config.yml files.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_kill)
    return parser


def cmd_kill(args, workflow_args):
    from .tasks import kill_tasks
    from .utils import env, load_config_files
    from .hosts import Host
    env.verbosity = args.verbosity
    if args.queue == '':
        cfg = load_config_files(args.config)
        from .hosts import list_queues
        list_queues(cfg, args.verbosity)
        return
    if not args.queue:
        if args.all:
            if args.tasks:
                env.logger.warning('Task ids "{}" are ignored with option --all'.format(' '.join(args.tasks)))
            if args.tags:
                env.logger.warning('Option tags is ignored with option --all')
            kill_tasks([])
        else:
            if not args.tasks and not args.tags:
                env.logger.warning('Please specify task id, or one of options --all and --tags')
            else:
                kill_tasks(tasks=args.tasks, tags=args.tags)
    else:
        # remote host?
        cfg = load_config_files(args.config)
        host = Host(args.queue)
        print(host._task_engine.kill_tasks(tasks=args.tasks, tags=args.tags))

#
# command remove
#
def get_remove_parser(desc_only=False):
    parser = argparse.ArgumentParser('remove',
        description='''Remove specified files and/or their signatures''')
    if desc_only:
        return parser
    parser.add_argument('targets', nargs='*', metavar='FILE_OR_DIR',
        help='''Files and directories to be removed. Directories will be
        scanned for files to removed but no directory will be removed.''')
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-t', '--tracked', action='store_true', default=False,
        help='''Limit files to only files tracked by SoS, namely files that are
            input, output, or dependent files of steps.''')
    group.add_argument('-u', '--untracked', action='store_true', default=False,
        help='''Limit files to untracked files, namely files that are not
            tracked by SoS steps.''')
    group.add_argument('-s', '--signature', action='store_true', default=False,
        help='''Remove signatures of specified files (not files themselves).
        As a special case, all local signatures will be removed if this option
        is specified without target.''')
    group.add_argument('-z', '--zap', action='store_true', default=False,
        help='''Replace files with their signatures. The file will not be
        regenerated by SoS unless is it actually needed by other steps. This
        option is usually used to remove large intermediate files from
        completed workflows while allowing relevant steps to be skipped
        during re-execution of the workflow.''')
    parser.add_argument('-e', '--external', action='store_true', default=False,
        help='''By default the remove command will only remove files and
        signatures under the current project directory. This option allows
        sos to remove files and/or signature of external files.''')
    parser.add_argument('--size',
        help='''Limit to files that exceed or smaller than specified size.
        Value of option should be in unit K, M, KB, MB, MiB, GB, etc, with
        optional prefix + for larger than (default), or - for smaller than
        specified size.''')
    parser.add_argument('--age', help='''Limit to files that are modified more than
        (default) or within specified age. Value of this parameter can be in units
        s (second), m (minute), h (hour), or d (day, default), or in the foramt of
        HH:MM:SS, with optional prefix + for older (default) and - for newer than
        specified age.''')
    parser.add_argument('-n', '--dryrun', action='store_true',
        help='''List files or directories to be removed, without actually
            removing them.''')
    parser.add_argument('-y', '--yes', action='store_true', dest='__confirm__',
        help='''Remove files without confirmation, suitable for batch removal
            of files.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_remove)
    return parser

class AnswerMachine:
    def __init__(self, always_yes=False, confirmed=False):
        self._always_yes = always_yes
        self._confirmed = confirmed

    def get(self, msg):
        if self._always_yes:
            return True
        if self._confirmed:
            print(msg)
            return True
        while True:
            res = input('{} (y/n/a)? '.format(msg))
            if res == 'a':
                self._confirmed = True
                return True
            elif res == 'y':
                return True
            elif res == 'n':
                return False

def get_tracked_files(sig_file):
    from .targets import file_target
    with open(sig_file) as sig:
        tracked_files = []
        script_files = []
        runtime_files = [sig_file]
        for line in sig:
            if line.startswith('IN_FILE') or line.startswith('OUT_FILE'):
                # format is something like IN_FILE\tfilename=xxxx\tsession=...
                tracked_files.append(line.rsplit('\t', 4)[1][9:])
                t = file_target(tracked_files[-1])
                if t.target_exists('signature'):
                    runtime_files.append(t.sig_file())
            elif line.startswith('EXE_SIG'):
                runtime_files.append('.sos/.runtime/{}.exe_info'.format(line.split('session=', 1)[1].strip()))
            elif line.startswith('# script:'):
                script_files.append(line.split(':', 1)[1].strip())
            elif line.startswith('# included:'):
                script_files.extend(line.split(':', 1)[-1].strip().split(','))
    return set([x for x in script_files if x.strip()]), set(tracked_files), set(runtime_files)


def cmd_remove(args, unknown_args):
    import glob
    from .utils import env
    import shutil
    from .targets import file_target
    env.verbosity = args.verbosity

    # what about global signature?
    sig_files = glob.glob('.sos/*.sig')
    tracked_files = set()
    runtime_files = set()
    for sig_file in sig_files:
        s, t, r = get_tracked_files(sig_file)
        tracked_files |= t
        runtime_files |= r
    #
    if args.signature and not args.targets:
        # a special case where all file and runtime signatures are removed.
        # no other options are allowed.
        removed_cnt = 0
        for s in glob.glob(os.path.join('.sos', '.runtime', '*.file_info')):
            try:
                if args.dryrun:
                    print('Would remove {}'.format(os.path.basename(s)))
                else:
                    env.logger.debug('Remove {}'.format(s))
                    os.remove(s)
                removed_cnt += 1
            except Exception as e:
                env.logger.warning('Failed to remove signature {}: {}'.format(s, e))
        if args.dryrun:
            env.logger.info('Would remove {} runtime signatures'.format(removed_cnt))
        elif removed_cnt:
            env.logger.info('{} runtime signatures removed'.format(removed_cnt))
        else:
            env.logger.info('No runtime signatures removed')
        return
    #
    tracked_files = {os.path.abspath(os.path.expanduser(x)) for x in tracked_files}

    if tracked_files:
        env.logger.info('{} tracked files from {} run{} are identified.'
                .format(len(tracked_files), len(sig_files), 's' if len(sig_files) > 1 else ''))
    else:
        env.logger.info('No tracked file from {} run{} are identified.'
                .format(len(sig_files), 's' if len(sig_files) > 1 else ''))
    #
    if not args.targets:
        args.targets = ['.']
    #
    if args.size:
        from .utils import expand_size
        args.size = expand_size(args.size)
    if args.age:
        import time
        from .utils import expand_time
        args.age = expand_time(args.age, default_unit='d')
    if args.signature:
        def func(filename, resp):
            if os.path.abspath(filename) not in tracked_files:
                return False
            target = file_target(filename)
            if target.is_external() and not args.external():
                env.logger.debug('Ignore external file {}'.format(filename))
                return False
            if not target.target_exists('signature'):
                return False
            if args.size:
                if (args.size > 0 and os.path.getsize(filename) < args.size) or \
                    (args.size < 0 and os.path.getsize(filename) > -args.size):
                    env.logger.debug('{} ignored due to size limit {}'.format(filename, args.size))
                    return False
            if args.age:
                if (args.age > 0 and time.time() - os.path.getmtime(filename) < args.age) or \
                    (args.age < 0 and time.time() - os.path.getmtime(filename) > -args.age):
                    env.logger.debug('{} ignored due to age limit {}'.format(filename, args.age))
                    return False
            if not args.dryrun:
                env.logger.debug('Remove {}'.format(s))
                try:
                    os.remove(target.sig_file())
                except Exception as e:
                    env.logger.warning('Failed to remove signature of {}: {}'.format(filename, e))
                return True
            return False
    elif args.tracked:
        def func(filename, resp):
            if os.path.abspath(filename) not in tracked_files:
                return False
            target = file_target(filename)
            if target.is_external() and not args.external():
                env.logger.debug('Ignore external file {}'.format(filename))
                return False
            if args.size:
                if (args.size > 0 and os.path.getsize(filename) < args.size) or \
                    (args.size < 0 and os.path.getsize(filename) > -args.size):
                    env.logger.debug('{} ignored due to size limit {}'.format(filename, args.size))
                    return False
            if args.age:
                if (args.age > 0 and time.time() - os.path.getmtime(filename) < args.age) or \
                    (args.age < 0 and time.time() - os.path.getmtime(filename) > -args.age):
                    env.logger.debug('{} ignored due to age limit {}'.format(filename, args.age))
                    return False
            if resp.get('{} tracked file {}'.format('Would remove' if args.dryrun else 'Remove', filename)):
                if not args.dryrun:
                    env.logger.debug('Remove {}'.format(s))
                    try:
                        target.remove('both')
                    except Exception as e:
                        env.logger.warning('Failed to remove {}: {}'.format(filename, e))
                    return True
            else:
                env.logger.debug('No signature exists for tracked file {}'.format(filename))
            return False
    elif args.untracked:
        def func(filename, resp):
            if os.path.abspath(filename) in tracked_files:
                return False
            target = file_target(filename)
            if target.is_external() and not args.external():
                env.logger.debug('Ignore external file {}'.format(filename))
                return False
            if args.size:
                if (args.size > 0 and os.path.getsize(filename) < args.size) or \
                    (args.size < 0 and os.path.getsize(filename) > -args.size):
                    env.logger.debug('{} ignored due to size limit {}'.format(filename, args.size))
                    return False
            if args.age:
                if (args.age > 0 and time.time() - os.path.getmtime(filename) < args.age) or \
                    (args.age < 0 and time.time() - os.path.getmtime(filename) > -args.age):
                    env.logger.debug('{} ignored due to age limit {}'.format(filename, args.age))
                    return False
            if resp.get('{} untracked file {}'.format('Would remove' if args.dryrun else 'Remove', filename)):
                if not args.dryrun:
                    env.logger.debug('Remove {}'.format(s))
                    try:
                        target.remove('both')
                    except Exception as e:
                        env.logger.warning('Failed to remove {}: {}'.format(filename, e))
                    return True
            else:
                env.logger.debug('No signature exists for tracked file {}'.format(filename))
            return False
    elif args.zap:
        def func(filename, resp):
            if os.path.abspath(filename) not in tracked_files:
                return False
            target = file_target(filename)
            if target.is_external() and not args.external():
                env.logger.debug('Ignore external file {}'.format(filename))
                return False
            if args.size:
                if (args.size > 0 and os.path.getsize(filename) < args.size) or \
                    (args.size < 0 and os.path.getsize(filename) > -args.size):
                    env.logger.debug('{} ignored due to size limit {}'.format(filename, args.size))
                    return False
            if args.age:
                if (args.age > 0 and time.time() - os.path.getmtime(filename) < args.age) or \
                    (args.age < 0 and time.time() - os.path.getmtime(filename) > -args.age):
                    env.logger.debug('{} ignored due to age limit {}'.format(filename, args.age))
                    return False
            if resp.get('{} tracked file {}'.format('Would zap' if args.dryrun else 'Zap', filename)):
                if not args.dryrun:
                    env.logger.debug('Zap {}'.format(s))
                    try:
                        target.write_sig()
                        shutil.copy(target.sig_file(), filename + '.zapped')
                        os.remove(filename)
                    except Exception as e:
                        env.logger.warning('Failed to zap {}: {}'.format(filename, e))
                    return True
            else:
                env.logger.debug('No signature exists for tracked file {}'.format(filename))
            return False
    else:
        # default behavior
        def func(filename, resp):
            target = file_target(filename)
            if target.is_external() and not args.external():
                env.logger.debug('Ignore external file {}'.format(filename))
                return False
            if args.size:
                if (args.size > 0 and os.path.getsize(filename) < args.size) or \
                    (args.size < 0 and os.path.getsize(filename) > -args.size):
                    env.logger.debug('{} ignored due to size limit {}'.format(filename, args.size))
                    return False
            if args.age:
                if (args.age > 0 and time.time() - os.path.getmtime(filename) < args.age) or \
                    (args.age < 0 and time.time() - os.path.getmtime(filename) > -args.age):
                    env.logger.debug('{} ignored due to age limit {}'.format(filename, args.age))
                    return False
            if resp.get('{} file {}'.format('Would remove' if args.dryrun else 'Remove', filename)):
                if not args.dryrun:
                    env.logger.debug('Remove {}'.format(s))
                    try:
                        target.remove('both')
                    except Exception as e:
                        env.logger.warning('Failed to remove {}: {}'.format(filename, e))
                    return True
            else:
                env.logger.debug('No signature exists for tracked file {}'.format(filename))
            return False

    removed = 0
    resp = AnswerMachine(always_yes = args.dryrun, confirmed = args.__confirm__)
    for target in args.targets:
        target = os.path.expanduser(target)
        if file_target(target).is_external():
            if os.path.isdir(target):
                sys.exit('Canot remove external directory {}'.format(target))
            elif not args.external:
                sys.exit('Only subdirectories of the current directory can be removed unless option --external is specified. {} specified.'.format(target))
        if os.path.isfile(target):
            removed += func(target, resp)
            continue
        # directory
        for dirname, dirlist, filelist in os.walk(target):
            # we do not remove files under the current directory
            if dirname != '.':
                for x in filelist:
                    # ignore hidden file
                    if x.startswith('.'):
                        continue
                    removed += func(os.path.join(dirname, x), resp)
            dirlist[:] = [x for x in dirlist if not x.startswith('.')]
    env.logger.info('{}{} file{} {}'.format('Signagure of ' if args.signature else '', removed,
        's' if removed > 1 else '', 'zapped' if args.zap else 'removed'))

#
# subcommand config
#
def get_config_parser(desc_only=False):
    parser = argparse.ArgumentParser('config',
        description='''Displays configurations in host, global, local, and user specified
            configuration files. ''')
    parser.short_description = '''Read and write sos configuration files'''
    if desc_only:
        return parser
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--site', action='store_true', dest='__site_config__',
        help='''Set (--set) or unset (--unset) options in system site configuration file
            (${SOS}/site_config.yml).''')
    group.add_argument('-g', '--global', action='store_true', dest='__global_config__',
        help=argparse.SUPPRESS)
    group.add_argument('--hosts', action='store_true', dest='__hosts_config__',
        help='''Set (--set) or unset (--unset) options in hosts (~/.sos/hosts.yml)''')
    group.add_argument('-c', '--config', dest='__config_file__', metavar='CONFIG_FILE',
        help='''Set (--set) or unset (--unset) options in user specified configuration file,
            or display options (--get) also in this file.''')
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--get', nargs='*', metavar='OPTION', dest='__get_config__',
        help='''Display values of options that contain one of the specified words
            from all configuration files.''')
    group.add_argument('--unset', nargs='+', metavar='OPTION',  dest='__unset_config__',
        help='''Unset (remove) settings for specified options. The arguments of this
        option can be a single configuration option or a list of option. Wildcard
        characters are allowed to match more options (e.g. '*timeout', or '*' for
        all options, quotation is needed to avoid shell expansion).''')
    group.add_argument('--set', nargs='+', metavar='KEY VALUE', dest='__set_config__',
        help='''--set KEY VALUE sets VALUE to variable KEY. The value can be any valid
        python expression (e.g. 5 for integer 5 and '{"c": 2, "d": 1}' for a dictionary)
        with invalid expression (e.g. val without quote) considered as string. Syntax
        'A.B=v' can be used to add {'B': v} to dictionary 'A', and --set KEY VALUE1 VALUE2 ...
        will create a list with multiple values.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_config)
    return parser

def cmd_config(args, workflow_args):
    import fnmatch
    import yaml
    from .utils import env, dict_merge, load_config_files
    from .syntax import CONFIG_NAME
    if workflow_args:
        raise RuntimeError('Unrecognized arguments {}'.format(' '.join(workflow_args)))
    #
    if args.__unset_config__:
        if args.__site_config__:
            config_file = os.path.join(os.path.split(__file__)[0], 'site_config.yml')
        elif args.__hosts_config__:
            config_file = os.path.join(os.path.expanduser('~'), '.sos', 'hosts.yml')
        elif args.__config_file__:
            config_file = os.path.expanduser(args.__config_file__)
        else:
            config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yml')

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
            sys.exit(1)
        #
        def unset_key(prefix, cfg):
            changed = 0
            unset = []
            for option in args.__unset_config__:
                for k in cfg.keys():
                    if fnmatch.fnmatch('.'.join(prefix + [k]), option):
                        unset.append(k)
                        print('Unset {}'.format('.'.join(prefix + [k])))
            #
            if unset:
                changed += len(unset)
                for k in set(unset):
                    cfg.pop(k)

            for k,v in cfg.items():
                if isinstance(v, dict):
                    changed += unset_key(prefix + [k], v)
            return changed
        #
        if unset_key([], cfg):
            with open(config_file, 'w') as config:
                config.write(yaml.safe_dump(cfg, default_flow_style=False))
        else:
            env.logger.warning('{} does not match any configuration key'.format(', '.join(args.__unset_config__)))
    elif args.__set_config__:
        if args.__site_config__:
            config_file = os.path.join(os.path.split(__file__)[0], 'site_config.yml')
        elif args.__hosts_config__:
            config_file = os.path.join(os.path.expanduser('~'), '.sos', 'hosts.yml')
        elif args.__config_file__:
            config_file = os.path.expanduser(args.__config_file__)
        else:
            config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yml')

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
            env.logger.error('Unacceptable variable name {}.'.format(k))
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
                env.logger.debug('Value "{}" is an invalid expression and is treated as a string.'.format(v))
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
    elif args.__get_config__ is not None:
        cfg = load_config_files(args.__config_file__)
        def disp_matched(obj, options, prefix=[]):
            for k, v in obj.items():
                key = '.'.join(prefix + [k])
                if isinstance(v, dict):
                    disp_matched(v, options, prefix + [k])
                elif not options or any(option in key for option in options):
                    print('{}\t{!r}'.format(key, v))

        disp_matched(cfg, args.__get_config__)
    else:
        from pprint import PrettyPrinter
        pp = PrettyPrinter(indent=2)
        cfg = load_config_files(args.__config_file__)
        pp.pprint(cfg)


#
# command pack
#
def get_pack_parser(desc_only=False):
    parser = argparse.ArgumentParser('pack',
        description='''Collect sos scripts, all input, output, and tracked intermediate
        files related to a workflow run and bundle them into a single archive.
        The archive can be examined (without unpacking) with command "sos
        show" and be unpacked with command "sos unpack". This command does not
        include files outside of the current working directory unless they
        are specified by option --include, or --all.''')
    parser.short_description = '''Pack workflow related files into an archive'''
    if desc_only:
        return parser
    parser.add_argument('session', nargs='?',
        help='''ID of the session to be saved, which can be any number of
        digits as long as it can uniquely determine a workflow session. This
        parameter can be ignored if only one session is available.''')
    parser.add_argument('-o', '--output', default="-",
        help='''Output file, which can be a file with extension ".sar" (the
        extension will be be automatically appended if needed), or "-" for
        standard output (default).''')
    parser.add_argument('-i', '--include', nargs='*', default=[],
        help='''Additional files or directories to be incldued in the archive.
        SoS will archive all files under specified directories, including
        hidden directories such as ".git". Option --exclude could be used
        to exclude these files.''')
    parser.add_argument('-e', '--exclude', nargs='*', default=[],
        help='''Files that should be excluded from archive. The parameter
        should be one or more patterns that match the whole path (e.g.
        "output/*.log" or file or directory names such as "tmp" or "*.bam".
        ''')
    parser.add_argument('-a', '--all', action='store_true', dest='__all__',
        help='''Include all tracked files even if they reside outside of the
        current working directory.''')
    parser.add_argument('-m', '--message',
        help='''A short message to be included into the archive. Because the
        message would be lost during unpacking, it is highly recommended that
        you create a README file and include it with option --include.''')
    parser.add_argument('-n', '--dryrun', action='store_true',
        help='''List files to be included and total file size without actually
        archiving them''')
    parser.add_argument('-y', '--yes', action='store_true', dest='__confirm__',
        help='''Overwrite output file if it already exists''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_pack)
    return parser

def locate_files(session, include, exclude, all_files):
    import fnmatch
    import glob
    from .utils import env
    from .targets import file_target
    sig_files = glob.glob('.sos/*.sig')
    if not sig_files:
        raise ValueError('No executed workflow is identified.')
    if not session:
        if len(sig_files) == 1:
            sig_file = sig_files[0]
        else:
            raise ValueError('More than one sessions have been executed. '
                'Please specify one of the sessions to save.\n'
                'Available sessions are:\n' +
                '\n'.join(os.path.basename(x)[:-4] for x in sig_files))
    else:
        matched = [x for x in sig_files if os.path.basename(x).startswith(session)]
        if len(matched) == 1:
            sig_file = matched[0]
        elif not matched:
            raise ValueError('No session matches specified session ID ({}). '.format(session) +
                'Available sessions are:\n' +
                '\n'.join(os.path.basename(x)[:-4] for x in sig_files))
        else:
            raise ValueError('More than one matching sessions have been located. '
                'Please specify one of the sessions to save.\n '
                'Available sessions are:\n' +
                '\n'.join(os.path.basename(x)[:-4] for x in sig_files))
    #
    script_files, tracked_files, runtime_files = get_tracked_files(sig_file)
    # all
    if not all_files:
        external_files = []
        for x in tracked_files:
            if file_target(x).is_external():
                env.logger.info('{} is excluded. Use option --all to include tracked files outside of current directory.'.format(x))
                external_files.append(x)
        tracked_files -= set(external_files)
    # include
    for inc in include:
        if os.path.isfile(os.path.expanduser(inc)):
            tracked_files.add(os.path.expanduser(inc))
        elif os.path.isdir(os.path.expanduser(inc)):
            for dirname, dirlist, filelist in os.walk(os.path.expanduser(inc)):
                tracked_files.update(os.path.join(dirname, x) for x in filelist if not x.startswith('.'))
                dirlist[:] = [x for x in dirlist if not x.startswith('.')]
        else:
            raise ValueError('Extra include file {} does not exist'.format(inc))
    # excludle
    for ex in exclude:
        tracked_files = [x for x in tracked_files if not fnmatch.fnmatch(x, ex)]
    #
    return script_files, tracked_files, runtime_files

def cmd_pack(args, unknown_args):
    import tarfile
    import tempfile
    from tqdm import tqdm as ProgressBar
    from .utils import pretty_size, env, ProgressFileObj
    from .targets import file_target
    #
    env.verbosity = args.verbosity
    try:
        script_files, tracked_files, runtime_files = locate_files(args.session, args.include, args.exclude, args.__all__)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
    #
    # get information about files
    file_sizes = {x: file_target(x).size() for x in tracked_files}
    # getting file size to create progress bar
    total_size = sum(file_sizes.values())

    if args.output == '-':
        tar_args = {'fileobj': sys.stdout.buffer, 'mode': 'w:gz'}
    elif not args.output.endswith('.sar'):
        tar_args = {'name': args.output + '.sar', 'mode': 'w:gz'}
    else:
        tar_args = {'name': args.output, 'mode': 'w:gz'}

    resp = AnswerMachine(always_yes=False, confirmed=args.__confirm__)
    if os.path.isfile(args.output) and not args.__confirm__ and not resp.get('Overwrite {}'.format(args.output)):
        env.logger.info('Operation aborted due to existing output file')
        sys.exit(0)

    prog = ProgressBar(desc='Checking', total=total_size, disable=args.verbosity != 1)
    manifest_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(manifest_file, 'w') as manifest:
        # write message in repr format (with "\n") to keep it in the same line
        manifest.write('# {!r}\n'.format(args.message if args.message else ''))
        # add .archive.info file
        for f in script_files:
            if f == 'None':
                continue
            ft = file_target(f)
            if not ft.target_exists():
                env.logger.warning('Missing script file {}'.format(ft.target_name()))
            else:
                manifest.write('SCRIPTS\t{}\t{}\t{}\t{}\n'.format(os.path.basename(f), ft.mtime(), ft.size(), ft.target_signature()))
        for f in tracked_files:
            env.logger.info('Checking {}'.format(f))
            ft = file_target(f)
            if not ft.target_exists():
                env.logger.warning('Missing tracked file {}'.format(ft.target_name()))
            elif ft.is_external():
                manifest.write('EXTERNAL\t{}\t{}\t{}\t{}\n'.format(f.replace('\\', '/'), ft.mtime(), ft.size(), ft.target_signature()))
            else:
                manifest.write('TRACKED\t{}\t{}\t{}\t{}\n'.format(f.replace('\\', '/'), ft.mtime(), ft.size(), ft.target_signature()))
        for f in runtime_files:
            ft = file_target(f)
            if not ft.target_exists():
                env.logger.warning('Missing runtime file {}'.format(ft.target_name()))
            else:
                manifest.write('RUNTIME\t{}\t{}\t{}\t{}\n'.format(os.path.basename(f), ft.mtime(), ft.size(), ft.target_signature()))
    prog.close()
    #
    if args.dryrun:
        print('A total of {} files ({}) with additional scripts and runtime files would be archived.'.
            format(len(tracked_files), pretty_size(total_size)))
        sys.exit(0)
    else:
        env.logger.info('Archiving {} files ({})...'.format(len(tracked_files), pretty_size(total_size)))
    #
    prog = ProgressBar(desc=args.output, total=total_size, disable=args.verbosity != 1)
    with tarfile.open(**tar_args) as archive:
        # add manifest
        archive.add(manifest_file, arcname='MANIFEST.txt')
        # add .archive.info file
        for f in script_files:
            if not os.path.isfile(f):
                continue
            env.logger.info('Adding {}'.format(os.path.basename(f)))
            archive.add(f, arcname='scripts/' + os.path.basename(f))
        for f in tracked_files:
            if not os.path.isfile(f):
                if os.path.isfile(f + '.zapped'):
                    f = f + '.zapped'
                else:
                    continue
            env.logger.info('Adding {}'.format(f))
            if file_target(f).is_external():
                # external files
                if args.verbosity == 1:
                    tarinfo = archive.gettarinfo(f, arcname='external/' + f)
                    archive.addfile(tarinfo, fileobj=ProgressFileObj(prog, f, 'rb'))
                else:
                    archive.add(f, arcname='external/' + f)
            else:
                if args.verbosity == 1:
                    tarinfo = archive.gettarinfo(f, arcname='tracked/' + f)
                    archive.addfile(tarinfo, fileobj=ProgressFileObj(prog, f, 'rb'))
                else:
                    archive.add(f, arcname='tracked/' + f)
        env.logger.info('Adding runtime files')
        for f in runtime_files:
            if not os.path.isfile(f):
                continue
            env.logger.trace('Adding {}'.format(os.path.basename(f)))
            archive.add(f, arcname='runtime/' + os.path.basename(f))
    prog.close()

#
# command unpack
#
def get_unpack_parser(desc_only=False):
    parser = argparse.ArgumentParser('unpack',
        description='''Unpack a sos archive to a specified directory. For security
        reasons, files that were outside of the project directory would be
        extracted in this directory unless option -e is specified.''')
    parser.short_description = '''Unpack workflow related files from an SoS archive'''
    if desc_only:
        return parser
    parser.add_argument('archive',
        help='''SoS archive saved by command sos pack''')
    parser.add_argument('files', nargs='*',
        help='''An optional list of files to be processed, which can be exact
        filenames or patterns (e.g. "*.bam"). No runtime information will be
        extracted if this option is specified.''')
    parser.add_argument('-d', '--dest', default='.',
        help='''Directory where a sos archive would be unpacked. Default to
        current directory.''')
    parser.add_argument('-s', '--script', action='store_true',
        help='''If specified, extract sos script(s) related to the workflow
        to the current or specified directory (option --dest), regardless of
        their original locations in the filesystem. Note that sos scripts are
        not extracted by default because they are usually external and it is
        dangerous to overwrite existing scripts with archived ones.''')
    parser.add_argument('-l', '--list', action='store_true', dest='__list__',
        help='''List content of the archive instead of extracting it. The names,
        uncompressed file sizes  and  modification  dates and times of the
        specified files are printed, along with totals for all files specified.''')
    parser.add_argument('-e', '--external', action='store_true',
        help='''Extract files outside of the project to their external destinations.
        This option can be dangerous because it can overwrite system files silently
        if accompanied with option -y.''')
    parser.add_argument('-n', '--no', action='store_true', dest='__no_overwrite__',
        help='''Do not overwrite existing files without promoting users.''')
    parser.add_argument('-y', '--yes', action='store_true', dest='__confirm__',
        help='''Overwrite existing files without promoting users. This option
        can be dangerous to use. Note that SoS checks file signature and
        ignores existing files that are identical to those in the archive.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_unpack)
    return parser

def cmd_unpack(args, unknown_args):
    import tarfile
    from tqdm import tqdm as ProgressBar
    from .utils import env, pretty_size, ProgressFileObj
    from .targets import fileMD5
    import fnmatch
    import time

    resp = AnswerMachine(always_yes=False, confirmed=False)

    env.verbosity = args.verbosity
    prog = ProgressBar(desc='Extracting {}'.format(args.archive), total=os.path.getsize(args.archive),
        disable=args.verbosity!=1 or args.__list__)
    try:
        with tarfile.open(fileobj=ProgressFileObj(prog, args.archive, 'rb')) as archive:
            manifest_file = archive.next()
            if manifest_file.name != 'MANIFEST.txt':
                raise ValueError('Manifest not found in SoS archive {}.'.format(args.archive))
            archive.extract(manifest_file)
            md5 = {}
            if args.__list__:
                print('    Length  Date     Time  Type    Name')
                print('    ------  ----     ----  ----    ----')
            with open(manifest_file.name) as manifest:
                line = manifest.readline()
                total_size = 0
                total_files = 0
                for line in manifest:
                    fields = line.split()
                    # type, name, mtime, size, md5
                    if args.__list__:
                        if fields[0] == 'RUNTIME':
                            continue
                        print(' {:>9s}  {:>12s} {:>6s} {}'.format(pretty_size(int(fields[3])),
                            time.strftime('%m-%d-%y %H:%M', time.gmtime(float(fields[2]))),
                            fields[0], fields[1]))
                        total_size += int(fields[3])
                        total_files += 1
                    md5[fields[1]] = fields[4].strip()
            os.remove(manifest_file.name)
            if args.__list__:
                print('   ------                  ----')
                print('{:>9s}                  {} files'.format(pretty_size(total_size), total_files) )
                return
            while True:
                f = archive.next()
                if f is None:
                    break
                if args.files:
                    # see if filename matches specified name
                    selected = False
                    for pattern in args.files:
                        # runtime information is not extracted with the specification of any file
                        if f.name.startswith('runtime/'):
                            continue
                        m_name = f.name.split('/', 1)[-1]
                        if fnmatch.fnmatch(m_name, pattern) or fnmatch.fnmatch(os.path.basename(m_name), pattern):
                            selected=True
                            break
                    if not selected:
                        env.logger.debug('Ignore {}'.format(m_name))
                        continue
                dest = args.dest
                is_runtime = False
                # hacking f.name to correct destination
                if f.name.startswith('external/'):
                    if not resp.get('Extract {} to outside of current directory'.format(f.name[9:])):
                        continue
                    f.name = f.name[9:]
                elif f.name.startswith('tracked/'):
                    f.name = f.name[8:]
                elif f.name.startswith('runtime/'):
                    is_runtime = True
                    if f.name.endswith('.sig'):
                        # this goes to local directory
                        dest = os.path.join(args.dest, '.sos')
                    else:
                        # this goes to global signature directory
                        dest = '.sos/.runtime'
                    f.name = f.name[8:]
                elif f.name.startswith('scripts/'):
                    if not args.script:
                        env.logger.debug('Ignore {}'.format(f.name[8:]))
                        continue
                    f.name = f.name[8:]
                else:
                    raise RuntimeError('Unexpected file {} from SoS archive'.format(f.name))

                dest_file = os.path.join(dest, f.name)
                # runtime file?
                if os.path.isfile(dest_file):
                    # signature files should not have md5
                    if dest_file.endswith('.zapped'):
                        continue
                    if fileMD5(dest_file) == md5[f.name]:
                        if not is_runtime:
                            env.logger.info('Ignore identical {}'.format(f.name))
                        continue
                    if not resp.get('Overwrite existing file {}'.format(f.name)):
                        continue
                if not f.name.startswith('.sos'):
                    env.logger.info('Extracting {}'.format(f.name))
                else:
                    env.logger.debug('Extracting {}'.format(f.name))
                archive.extract(f, path=dest)
    except Exception as e:
        raise ValueError('Failed to unpack SoS archive: {}'.format(e))
    prog.close()


#
# Handling addon commands
#
def handle_addon(args, unknown_args):
    for entrypoint in pkg_resources.iter_entry_points(group='sos_addons'):
        name = entrypoint.name.strip()
        if name.endswith('.func') and name.rsplit('.', 1)[0] == args.addon_name:
            func = entrypoint.load()
            func(args, unknown_args)

#
# this is the sos-runner command
#
def sosrunner():
    parser = get_run_parser()
    parser.prog = 'sos-runner'
    if len(sys.argv) > 2 and '-h' in sys.argv:
        try:
            from .parser import SoS_Script
            script = SoS_Script(filename=sys.argv[1])
            script.print_help()
            sys.exit(0)
        except Exception as e:
            sys.exit('No help information is available for script {}: {}'.format(sys.argv[1], e))
    args, workflow_args = parser.parse_known_args()
    cmd_run(args, workflow_args)


# add another ArgumentParser to an existing ArgumentParser as
# a subparser
#
def add_sub_parser(subparsers, parser, name=None, hidden=False):
    if hidden:
        return subparsers.add_parser(parser.prog if name is None else name,
            description=parser.description,
            epilog=parser.epilog,
            parents=[parser],
            add_help=False)
    else:
        return subparsers.add_parser(parser.prog if name is None else name,
            description=parser.description,
            epilog=parser.epilog,
        help=parser.short_description if hasattr(parser, 'short_description') else parser.description,
            parents=[parser],
            add_help=False)


def main():
    from ._version import SOS_FULL_VERSION

    # only load specific subparser to save on start-up time
    if len(sys.argv) == 1 or sys.argv[1] == '-h':
        subcommand = None
    else:
        subcommand = sys.argv[1]

    master_parser = argparse.ArgumentParser(description='''A workflow system
            for the execution of commands and scripts in different languages.''',
        prog='sos',
        fromfile_prefix_chars='@',
        epilog='''Use 'sos cmd -h' for details about each subcommand. Please
            contact Bo Peng (bpeng at mdanderson.org) if you have any question.''')

    try:
        master_parser.add_argument('--version', action='version',
            version='%(prog)s {}'.format(SOS_FULL_VERSION))
        subparsers = master_parser.add_subparsers(title='subcommands',
            # hide pack and unpack
            metavar = '{install,run,resume,dryrun,status,push,pull,execute,kill,purge,config,convert,remove}')

        # command install
        # add_sub_parser(subparsers, get_install_parser(desc_only='install'!=subcommand))
        #
        # command run
        add_sub_parser(subparsers, get_run_parser(desc_only='run'!=subcommand))
        #
        # command resume
        add_sub_parser(subparsers, get_resume_parser(desc_only='resume'!=subcommand))
        #
        # command dryrun
        add_sub_parser(subparsers, get_dryrun_parser(desc_only='dryrun'!=subcommand))

        #
        # command status
        add_sub_parser(subparsers, get_status_parser(desc_only='status'!=subcommand))
        #
        # command push
        add_sub_parser(subparsers, get_push_parser(desc_only='push'!=subcommand))
        #
        # command pull
        add_sub_parser(subparsers, get_pull_parser(desc_only='pull'!=subcommand))
        #
        # command preview
        add_sub_parser(subparsers, get_preview_parser(desc_only='preview'!=subcommand), hidden=True)
        #
        # command execute
        add_sub_parser(subparsers, get_execute_parser(desc_only='execute'!=subcommand))
        #
        # command kill
        add_sub_parser(subparsers, get_kill_parser(desc_only='kill'!=subcommand))
        #
        # command purge
        add_sub_parser(subparsers, get_purge_parser(desc_only='purge'!=subcommand))

        #
        # command config
        add_sub_parser(subparsers, get_config_parser(desc_only='config'!=subcommand))
        #
        # command convert
        add_sub_parser(subparsers, get_convert_parser(desc_only='convert'!=subcommand))
        #
        # command remove
        add_sub_parser(subparsers, get_remove_parser(desc_only='remove'!=subcommand))
        #
        # command pack
        add_sub_parser(subparsers, get_pack_parser(desc_only='pack'!=subcommand), hidden=True)
        #
        # command unpack
        add_sub_parser(subparsers, get_unpack_parser(desc_only='unpack'!=subcommand), hidden=True)
        #
        # addon packages
        if subcommand is None or subcommand not in ['install', 'run', 'dryrun', 'convert', 'push', 'pull',
                'remove', 'config', 'pack', 'unpack']:
            for entrypoint in pkg_resources.iter_entry_points(group='sos_addons'):
                if entrypoint.name.strip().endswith('.parser'):
                    name = entrypoint.name.rsplit('.', 1)[0]
                    func = entrypoint.load()
                    parser = add_sub_parser(subparsers, func(), name=name)
                    parser.add_argument('--addon-name', help=argparse.SUPPRESS,
                            default=name)
                    parser.set_defaults(func=handle_addon)
        #
        if len(sys.argv) == 1 or sys.argv[1] == '-h':
            master_parser.print_help()
            sys.exit(0)
        if '-h' in sys.argv:
            if len(sys.argv) > 3 and sys.argv[1] == 'run' and not sys.argv[2].startswith('-'):
                try:
                    from .parser import SoS_Script
                    script = SoS_Script(filename=sys.argv[2])
                    script.print_help()
                    sys.exit(0)
                except Exception as e:
                    sys.exit('No help information is available for script {}: {}'.format(sys.argv[1], e))
            if len(sys.argv) > 3 and sys.argv[1] == 'convert':
                print_converter_help()
        elif sys.argv[1] == 'convert':
            # this command has to be processed separately because I hat to use
            # sos convert sos-html FROM TO etc
            from_format, to_format = get_converter_formats(sys.argv[2:])
            if from_format is None or to_format is None:
                sys.exit('Cannot determine from or to format')
            sys.argv.insert(2, '{}-{}'.format(from_format, to_format))
        args, workflow_args = master_parser.parse_known_args()
        if not hasattr(args, 'func'):
            # in case of sos -v etc that no subcommand is specified
            master_parser.print_help()
            sys.exit(0)
        # calling the associated functions
        args.func(args, workflow_args)
    except KeyboardInterrupt:
        sys.exit('KeyboardInterrupt')
    except Exception as e:
        from .utils import env, get_traceback
        if env.verbosity and env.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

if __name__ == '__main__':
    main()
