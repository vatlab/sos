#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
from ._version import __version__

from .workflow_executor import Base_Executor
from sos.parser import SoS_Script

assert __version__


def execute_workflow(script: str,
                     workflow=None,
                     targets=None,
                     args=[],
                     options={},
                     config={}):
    '''
    Execute a SoS workflow with the following parameters:

    script:
        A (multi-line) string that defines one or more workflows.

    workflow (string, optional):
        Name of workflow to execute. This option can be ignored if
        ignored if the script defines a default workflow (with no name
        or with name `default`) or defines only a single workflow, or
        if the workflow is triggered by option "targets".

    targets (string or list, or sos_targets, optional):
        One (string) or more (list) of filenames as targets to be generated
        by the workflow, equivalent to option "-t" from command line.

    args (list or dict, optional):
        Command line arguments as a list in the format of "[`--cutoff', '0.5']"
        or a dictionary in the format of {"cutoff": 0.5} for workflow parameters
        defined in "parameter" statements.

    options (dict, optional):
        Dictionary with the following configuration options. Please
        refer to output of "sos run -h" for details about each option.
            config_file: configuration file ("-c"). The content of the config
                file can also be specified with option config.
            output_dag: option "-d"
            output_report: option "-p"
            default_queue: option "-q"
            worker_procs: option "-j"
            max_running_jobs: option "-J"
            sig_mode: option "-s"
            run_mode: "run" or "dryrun"
            verbosity: option "-v"
            trace_existing: option "-T"

    config (dict, optional):
        config as if loaded from "options['config_file']"

    Note: executing on specified host (option "-r") is not supported by this function.
    '''
    try:
        script = SoS_Script(script)
    except Exception as e:
        # show script with error
        raise ValueError(f'Failed to parse script {script}: {e}')
    #
    if workflow and targets:
        raise ValueError(
            "Only one of parameters workflow and targets should be specified.")
    wf = script.workflow(workflow, use_default=not targets)

    if not isinstance(config, dict):
        raise ValueError(f'Option config should be a dictionary.')
    run_options = {
        'config_file': None,
        'extra_config': config,
        'output_dag': None,
        'output_report': None,
        'default_queue': None,
        'worker_procs': None,
        'max_running_jobs': None,
        'sig_mode': 'default',
        'run_mode': 'run',
        'verbosity': 2,
        'workdir': os.getcwd(),
        'script': script,
        'workflow': workflow,
        'targets': targets,
        'workflow_args': args,
        'trace_existing': False,
        'error_mode': 'default',
        'exec_mode': None
    }
    # a convenience feature
    if isinstance(args, dict):
        run_options['workflow_vars'] = args
        args = []
    run_options.update(options)

    from .utils import env
    env.verbosity = run_options['verbosity']

    executor = Base_Executor(wf, args=args, config=run_options)
    if isinstance(targets, str):
        targets = [targets]
    return executor.run(targets=targets)
