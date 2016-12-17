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
import argparse
import pkg_resources

script_help = '''A SoS script that defines one or more workflows. The
    script can be a filename or a URL from which the content of a SoS will
    be read. If a valid file cannot be located or downloaded, SoS will
    search for the script in a search path specified by variable `sos_path`
    defined in the global SoS configuration file (~/.sos/config.yml).'''
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
# subcommand convert
#
def get_convert_parser():
    parser = argparse.ArgumentParser('convert',
        description='''Converts .sos to various formats including
            .html for web display, to jupyter notebook (.ipynb), and to terminal
            for syntax highlighted viewing on terminal. It also allows converting
            from jupyter notebook (.ipynb) to sos script (.sos).''',
        epilog='''Extra command line argument could be specified to customize
            the style of html, markdown, and terminal output. ''',
        )
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).'''),
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
    parser = argparse.ArgumentParser('convert')
    parser.add_argument('from_file', nargs='?')
    parser.add_argument('to_file', nargs='?')
    parser.add_argument('-t', '--to', dest='__to_format__')
    args, unknown_args = parser.parse_known_args(argv)
    from_format = args.from_file[1:] if args.from_file.startswith('.') and args.from_file.count('.') == 1 \
        else os.path.splitext(args.from_file)[-1][1:]
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
    try:
        for entrypoint in pkg_resources.iter_entry_points(group='sos_converters'):
            try:
                if entrypoint.name == args.converter_name + '.func':
                    func = entrypoint.load()
                    func(args.from_file, args.to_file, args, unknown_args)
            except Exception as e:
                raise RuntimeError('Failed to load converter {}: {}'.format(entrypoint.name, e))
    except Exception as e:
        # if no other parameter, with option list all
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand run
#
def get_run_parser(interactive=False, with_workflow=True):
    parser = argparse.ArgumentParser(prog='run',
        description='Execute default or specified workflow defined in script',
        epilog=workflow_options)
    if not interactive:
        parser.add_argument('script', metavar='SCRIPT', help=script_help)
    if with_workflow:
        parser.add_argument('workflow', metavar='WORKFLOW', nargs='?',
            help=workflow_spec)
    if not interactive:
        parser.add_argument('-j', type=int, metavar='JOBS', default=4, dest='__max_jobs__',
        help='''Number of concurrent process allowed. A workflow is by default
            executed sequentially (-j 1). If a greater than 1 number is specified
            SoS will execute the workflow in parallel mode and execute up to
            specified processes concurrently. These include looped processes
            within a step (with runtime option `concurrent=True`) and steps with
            non-missing required files.''')
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
    if not interactive:
        parser.add_argument('-q', dest='__queue__', metavar='QUEUE',
        help='''Task-processing queue. SoS by default uses a local multiprocessing
            queue where tasks are executed by different processes. Supported task
            queues include a 'rq' engine where tasks will be distributed to one or
            more rq-workers with assistance from a redis server, and a 'celery'
            quque where tasks will be distributed to celery workers.''')
    parser.add_argument('-r', dest='__report__', metavar='REPORT_FILE', nargs='?',
         help='''Default output of action report, which is by default the
            standard output but you can redirect it to another file with this
            option. Files specified by this option will be interpolated so you
            can use names such as 'out_${step_name}.md' to out to step-depedent
            files.''')
    #parser.add_argument('-t', dest='__transcript__', nargs='?',
    #    metavar='TRANSCRIPT', const='__STDERR__', help=transcript_help)
    runmode = parser.add_argument_group(title='Run mode options',
        description='''Control how sos scirpt is executed.''')
    runmode.add_argument('-n', action='store_true', dest='__dryrun__',
        help='''Execute a workflow without executing any actions. This can be
            used to check the syntax of a SoS file.''')
    runmode.add_argument('-s', choices=['default', 'ignore', 'force', 'build', 'assert'],
        default='ignore' if interactive else 'default', metavar='SIGMODE',
        dest='__sigmode__',
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
        default=1 if interactive else 2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).''')
    parser.set_defaults(func=cmd_run)
    return parser

def cmd_run(args, workflow_args):
    import atexit
    from .utils import env, get_traceback
    from .sos_script import SoS_Script
    
    # '' means no -d
    if args.__dag__ is None:
        args.__dag__ = '-'
    elif args.__dag__ == '':
        args.__dag__ = None
    env.max_jobs = args.__max_jobs__
    env.verbosity = args.verbosity

    if args.__queue__:
        # import all executors
        executor_class = None
        for entrypoint in pkg_resources.iter_entry_points(group='sos_executors'):
            # Grab the function that is the actual plugin.
            name = entrypoint.name
            if name == args.__queue:
                try:
                    executor_class = entrypoint.load()
                except Exception as e:
                    print('Failed to load queue executor {}: {}'.format(entrypoint.name, e))

        if not executor_class:
            sys.exit('Could not locate specified queue executor {}'.format(args.__queue__))
    else:
        from .sos_executor import Base_Executor, MP_Executor
        if args.__max_jobs__ == 1:
            executor_class = Base_Executor
        else:
            executor_class = MP_Executor

    # kill all remainging processes when the master process is killed.
    atexit.register(env.cleanup)
    #
    env.sig_mode = args.__sigmode__

    if args.__bin_dirs__:
        import fasteners
        for d in args.__bin_dirs__:
            if d == '~/.sos/bin' and not os.path.isdir(os.path.expanduser(d)):
                with fasteners.InterProcessLock('/tmp/sos_lock_bin'):
                    os.makedirs(os.path.expanduser(d))
            elif not os.path.isdir(os.path.expanduser(d)):
                raise ValueError('directory does not exist: {}'.format(d))
        os.environ['PATH'] = os.pathsep.join([os.path.expanduser(x) for x in args.__bin_dirs__]) + os.pathsep + os.environ['PATH']

    try:
        # workflow args has to be in the format of --arg, not positional, not -a
        if workflow_args and not workflow_args[0].startswith('--'):
            raise ValueError("Unrecognized command line option {}".format(' '.join(workflow_args)))
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow, use_default=not args.__targets__)
        executor = executor_class(workflow, args=workflow_args, config={
                'config_file': args.__config__,
                'output_dag': args.__dag__,
                'report_output': args.__report__})
        #
        if args.__dryrun__:
            executor.dryrun(args.__targets__)
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
def get_dryrun_parser():
    parser = argparse.ArgumentParser('dryrun',
        description='''Inspect specified script for syntax errors''',
        epilog=workflow_options)
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
    args.__sigmode__ = 'ignore'
    args.__queue__ = None
    args.__max_jobs__ = 1
    args.__report__ = None
    args.__dryrun__ = True
    args.__bin_dirs__ = []
    cmd_run(args, workflow_args)

#
# command remove
#
def get_remove_parser():
    parser = argparse.ArgumentParser('remove',
        description='''Remove specified files and directories and their
            signatures (if available). Optionally, you can remove only
            tracked files (input, output and intermediate files of executed
            workflows) or untracked file from specified files and/or
            directories.''')
    parser.add_argument('targets', nargs='*', metavar='FILE_OR_DIR',
        help='''Files and directories to be removed, which should be under the
            current directory (default). All, tracked, or untracked files
            will be removed depending on other options ('-t' or '-u').
            For safety reasons, files under the current directory have to be
            listed (not as files under .) to be removed.''')
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('-t', action='store_true', dest='__tracked__', default=False,
        help='''Remove tracked files and their signatures from specified files
            and directories.''')
    group.add_argument('-u', action='store_true', dest='__untracked__', default=False,
        help='''Remove untracked files from specified files and directories.''')
    parser.add_argument('-n', action='store_true', dest='__dryrun__',
        help='''List files or directories to be removed, without actually
            removing them.''')
    parser.add_argument('-y', '--yes', action='store_true', dest='__confirm__',
        help='''Remove files without confirmation, suitable for batch removal
            of files.''')
    parser.add_argument('-v', '--verbosity', type=int, choices=range(5), default=2,
        help='''Output error (0), warning (1), info (2), debug (3) and trace (4)
            information to standard output (default to 2).'''),
    parser.set_defaults(func=cmd_remove)
    return parser

def get_tracked_files(sig_file):
    from .target import FileTarget
    with open(sig_file) as sig:
        tracked_files = []
        script_files = []
        runtime_files = [sig_file]
        for line in sig:
            if line.startswith('IN_FILE') or line.startswith('OUT_FILE'):
                # format is something like IN_FILE\tfilename=xxxx\tsession=...
                tracked_files.append(line.rsplit('\t', 4)[1][9:])
                t = FileTarget(tracked_files[-1])
                if t.exists('signature'):
                    runtime_files.append(t.sig_file())
            elif line.startswith('EXE_SIG'):
                runtime_files.append('.sos/.runtime/{}.exe_info'.format(line.split('session=', 1)[1].strip()))
            elif line.startswith('# script:'):
                script_files.append(line.split(':', 1)[1].strip())
            elif line.startswith('# included:'):
                script_files.extend(line.split(':', 1)[-1].strip().split(','))
    return set(script_files), set(tracked_files), set(runtime_files)

def cmd_remove(args, unknown_args):
    import glob
    from .utils import env
    import shutil
    from collections import OrderedDict
    from .target import FileTarget

    sig_files = glob.glob('.sos/*.sig')
    if not sig_files:
        raise sys.exit('No executed workflow is identified.')
    tracked_files = set()
    for sig_file in sig_files:
        tracked_files |= get_tracked_files(sig_file)[1]
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
        relpath = os.path.relpath(target, '.')
        if relpath.startswith('..'):
            # we do not care about tracked files outside of current directory
            sys.exit('Only subdirectories of the current directory can be removed. {} specified.'.format(target))
        # file
        if os.path.abspath(target) in tracked_files:
            specified_tracked_files.append(target)
        elif os.path.isfile(target):
            specified_untracked_files.append(target)
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
            elif os.path.exists(target) and target != '.':
                if get_response('Remove {}'.format(target)):
                    if os.path.isfile(target):
                        os.remove(target)
                    elif os.path.isdir(target):
                        if os.path.abspath(target) != os.path.abspath('.'):
                            shutil.rmtree(target)
                    else:
                        os.unlink(target)

#
# subcommand config
#
def get_config_parser():
    parser = argparse.ArgumentParser('config',
        description='''Displays, set, and unset configuration
            variables defined in global or local configuration files.''')        
    parser.add_argument('-g', '--global', action='store_true', dest='__global_config__',
        help='''If set, change global (~/.sos/config.yml) instead of local
        (.sos/config.yml) configuration''')
    parser.add_argument('-c', '--config', dest='__config_file__', metavar='CONFIG_FILE',
        help='''User specified configuration file in YAML format. This file will not be
        automatically loaded by SoS but can be specified using option `-c`''')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--get', nargs='*', metavar='OPTION', dest='__get_config__',
        help='''Display values of specified configuration. The arguments of this
        option can be a single configuration option or a list of option. Wildcard
        characters are allowed to match more options (e.g. '*timeout', quotation
        is needed to avoid shell expansion). If no option is given, all options
        will be outputted.''')
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
            information to standard output (default to 2).'''),
    parser.set_defaults(func=cmd_config)
    return parser

def cmd_config(args, workflow_args):
    import fnmatch
    import yaml
    from .utils import env, dict_merge
    from .sos_syntax import CONFIG_NAME
    if workflow_args:
        raise RuntimeError('Unrecognized arguments {}'.format(' '.join(workflow_args)))
    #
    if args.__global_config__:
        config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yml')
    elif args.__config_file__:
        config_file = os.path.expanduser(args.__config_file__)
    else:
        config_file = 'config.yml'
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


#
# command pack
#
def get_pack_parser():
    parser = argparse.ArgumentParser('pack',
        description='''Collect sos scripts, all input, output, and tracked intermediate
        files related to a workflow run and bundle them into a single archive.
        The archive can be examined (without unpacking) with command "sos
        show" and be unpacked with command "sos unpack". This command does not
        include files outside of the current working directory unless they
        are specified by option --include, or --all.''')
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
            information to standard output (default to 2).'''),
    parser.set_defaults(func=cmd_pack)
    return parser

def locate_files(session, include, exclude, all_files):
    import fnmatch
    import glob
    from .utils import env
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
            relpath = os.path.relpath(x, '.')
            if relpath.startswith('..'):
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
    from .utils import pretty_size, env, ProgressBar, ProgressFileObj
    from .target import FileTarget
    #
    env.verbosity = args.verbosity
    try:
        script_files, tracked_files, runtime_files = locate_files(args.session, args.include, args.exclude, args.__all__)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
    #
    # get information about files
    file_sizes = {x: os.path.getsize(x) for x in tracked_files}
    # getting file size to create progress bar
    total_size = sum(file_sizes.values())

    if args.output == '-':
        tar_args = {'fileobj': sys.stdout.buffer, 'mode': 'w:gz'}
    elif not args.output.endswith('.sar'):
        tar_args = {'name': args.output + '.sar', 'mode': 'w:gz'}
    else:
        tar_args = {'name': args.output, 'mode': 'w:gz'}

    def get_response(msg):
        if args.__confirm__:
            print(msg)
            return True
        while True:
            res = input('{} (y/n)? '.format(msg))
            if res == 'y':
                return True
            elif res == 'n':
                return False

    if os.path.isfile(args.output) and not args.__confirm__ and not get_response('Overwrite {}'.format(args.output)):
        env.logger.info('Operation aborted due to existing output file')
        sys.exit(0)

    prog = ProgressBar('Checking', total_size, disp=args.verbosity == 1)
    manifest_file = tempfile.NamedTemporaryFile(delete=False).name
    with open(manifest_file, 'w') as manifest:
        # write message in repr format (with "\n") to keep it in the same line
        manifest.write('# {!r}\n'.format(args.message if args.message else ''))
        # add .archive.info file
        for f in script_files:
            if f == 'None':
                continue
            ft = FileTarget(f)
            if not ft.exists():
                env.logger.warning('Missing script file {}'.format(ft.name()))
            else:
                manifest.write('{}\t{}\t{}\t{}\n'.format(os.path.basename(f), ft.mtime(), ft.size(), ft.signature()))
        for f in tracked_files:
            env.logger.info('Checking {}'.format(f))
            ft = FileTarget(f)
            if not ft.exists():
                env.logger.warning('Missing tracked file {}'.format(ft.name()))
            else:
                manifest.write('{}\t{}\t{}\t{}\n'.format(f, ft.mtime(), ft.size(), ft.signature()))
        for f in runtime_files:
            ft = FileTarget(f)
            if not ft.exists():
                env.logger.warning('Missing runtime file {}'.format(ft.name()))
            else:
                manifest.write('{}\t{}\t{}\t{}\n'.format(f, ft.mtime(), ft.size(), ft.signature()))
    prog.done()
    #
    if args.dryrun:
        print('A total of {} files ({}) with additional scripts and runtime files would be archived.'.
            format(len(tracked_files), pretty_size(total_size)))
        sys.exit(0)
    else:
        env.logger.info('Archiving {} files ({})...'.format(len(tracked_files), pretty_size(total_size)))
    #
    prog = ProgressBar(args.output, total_size, disp=args.verbosity == 1)
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
                continue
            env.logger.info('Adding {}'.format(f))
            relpath = os.path.relpath(f, '.')
            if relpath.startswith('..'):
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
            env.logger.trace('Adding {}'.format(f))
            archive.add(f, arcname='runtime/' + f[5:])
    prog.done()

#
# command unpack
#
def get_unpack_parser():
    parser = argparse.ArgumentParser('unpack',
        description='''Unpack a sos archive to a specified directory. For security
        reasons, files that were outside of the project directory would be
        extracted in this directory unless option -e is specified.''')
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
            information to standard output (default to 2).'''),
    parser.set_defaults(func=cmd_unpack)
    return parser

def cmd_unpack(args, unknown_args):
    import tarfile
    from .utils import env, ProgressBar, pretty_size, ProgressFileObj
    from .target import fileMD5
    import fnmatch
    import time

    env.verbosity = args.verbosity
    def get_response(msg):
        if args.__confirm__:
            print(msg + '? y')
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

    prog = ProgressBar('Extracting {}'.format(args.archive), os.path.getsize(args.archive),
        disp=args.verbosity==1 and not args.__list__)
    try:
        with tarfile.open(fileobj=ProgressFileObj(prog, args.archive, 'rb')) as archive:
            manifest_file = archive.next()
            if manifest_file.name != 'MANIFEST.txt':
                raise ValueError('Manifest not found in SoS archive {}.'.format(args.archive))
            archive.extract(manifest_file)
            md5 = {}
            if args.__list__:
                print('   Length    Date   Time   Name')
                print('   ------    ----   ----   ----')
            with open(manifest_file.name) as manifest:
                line = manifest.readline()
                total_size = 0
                total_files = 0
                for line in manifest:
                    fields = line.split()
                    # name, mtime, size, md5
                    if args.__list__:
                        if fields[0].startswith('.sos'):
                            continue
                        print('{:>9s}  {:>12s}  {}'.format(pretty_size(int(fields[2])),
                            time.strftime('%m-%d-%y %H:%M', time.gmtime(float(fields[1]))), fields[0]))
                        total_size += int(fields[2])
                        total_files += 1
                    md5[fields[0]] = fields[3].strip()
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
                # hacking f.name to correct destination
                if f.name.startswith('external/'):
                    if not get_response('Extract {} to outside of current directory'.format(f.name[9:])):
                        continue
                    f.name = f.name[9:]
                elif f.name.startswith('tracked/'):
                    f.name = f.name[8:]
                elif f.name.startswith('runtime/'):
                    f.name = os.path.join('.sos', f.name[8:])
                elif f.name.startswith('scripts/'):
                    if not args.script:
                        env.logger.debug('Ignore {}'.format(f.name[8:]))
                        continue
                    f.name = f.name[8:]
                else:
                    raise RuntimeError('Unexpected file {} from SoS archive'.format(f.name))

                dest_file = os.path.join(args.dest, f.name)
                # runtime file?
                if os.path.isfile(dest_file):
                    # signature files should not have md5
                    if fileMD5(dest_file) == md5[f.name]:
                        if not f.name.startswith('.sos'):
                            env.logger.info('Ignore identical {}'.format(f.name))
                        continue
                    if not get_response('Overwrite existing file {}'.format(f.name)):
                        continue
                if not f.name.startswith('.sos'):
                    env.logger.info('Extracting {}'.format(f.name))
                else:
                    env.logger.debug('Extracting {}'.format(f.name))
                archive.extract(f, path=args.dest)
    except Exception as e:
        raise ValueError('Failed to unpack SoS archive: {}'.format(e))
    prog.done()


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
            from .sos_script import SoS_Script
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
def add_sub_parser(subparsers, parser, name=None):
    return subparsers.add_parser(parser.prog if name is None else name,
        description=parser.description,
        epilog=parser.epilog,
        help=parser.description, parents=[parser],
        add_help=False)


def main():
    from ._version import SOS_FULL_VERSION
    master_parser = argparse.ArgumentParser(description='''A workflow system
            for the execution of commands and scripts in different languages.''',
        prog='sos',
        fromfile_prefix_chars='@',
        epilog='''Use 'sos cmd -h' for details about each subcommand. Please
            contact Bo Peng (bpeng at mdanderson.org) if you have any question.''')

    master_parser.add_argument('--version', action='version',
        version='%(prog)s {}'.format(SOS_FULL_VERSION))
    subparsers = master_parser.add_subparsers(title='subcommands')
    #
    # command run
    add_sub_parser(subparsers, get_run_parser())
    #
    # command dryrun
    add_sub_parser(subparsers, get_dryrun_parser())
    #
    # command convert
    add_sub_parser(subparsers, get_convert_parser())
    #
    # command remove
    add_sub_parser(subparsers, get_remove_parser())
    #
    # command config
    add_sub_parser(subparsers, get_config_parser())
    #
    # command pack
    add_sub_parser(subparsers, get_pack_parser())
    #
    # command unpack
    add_sub_parser(subparsers, get_unpack_parser())
    #
    # addon packages
    for entrypoint in pkg_resources.iter_entry_points(group='sos_addons'):
        if entrypoint.name.strip().endswith('.parser'):
            name = entrypoint.name.rsplit('.', 1)[0]
            func = entrypoint.load()
            parser = add_sub_parser(subparsers, func(), name=name)
            parser.add_argument('--addon-name', help=argparse.SUPPRESS,
                    default=name)
            parser.set_defaults(func=handle_addon)
    #
    if len(sys.argv) == 1:
        master_parser.print_help()
        sys.exit(0)
    if '-h' in sys.argv:
        if len(sys.argv) > 3 and sys.argv[1] == 'run' and not sys.argv[2].startswith('-'):
            try:
                from .sos_script import SoS_Script
                script = SoS_Script(filename=sys.argv[2])
                script.print_help()
                sys.exit(0)
            except Exception as e:
                sys.exit('No help information is available for script {}: {}'.format(sys.argv[1], e))
        if len(sys.argv) > 4 and sys.argv[1] == 'convert':
            print_converter_help()
    elif sys.argv[1] == 'convert':
        # this command has to be processed separately because I hat to use
        # sos convert sos-html FROM TO etc
        from_format, to_format = get_converter_formats(sys.argv[2:])
        if from_format is None or to_format is None:
            sys.exit('Cannot determine from or to format')
        sys.argv.insert(2, '{}-{}'.format(from_format, to_format))
    args, workflow_args = master_parser.parse_known_args()
    # calling the associated functions
    args.func(args, workflow_args)

if __name__ == '__main__':
    main()
