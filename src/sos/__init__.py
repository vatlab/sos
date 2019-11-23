#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
from ._version import __version__

assert __version__

from .workflow_executor import Base_Executor
from sos.parser import SoS_Script


def execute_workflow(script: str,
                     workflow=None,
                     targets=None,
                     args=[],
                     config={}):
    '''
    Execute a SoS workflow with the following parameters:

    script:
        A (multi-line) string that defines one or more workflows.

    workflow (optional):
        Name of workflow to execute. This option can be ignored if
        ignored if the script defines a default workflow (with no name
        or with name `default`) or defines only a single workflow, or
        if the workflow is triggered by option "targets".

    targets (string or list, or sos_targets, optional):

    args (list or dict, optional):
        Command line arguments (a list in the format of "[`--name', 'value']")
        for workflow parameters defined in "parameter" statements. A dictionary
        (e.g. {"name", "value"}) is also acceptable and is equivalent to setting
        `config["workflow_vars"]'.

    config (dict, optional):
        Dictionary with the following configuration options. Please
        refer to output of "sos run -h" for details about each option.
            config_file: configuration file ("-c")
            output_dag: option "-d"
            output_report: option "-p"
            default_queue: option "-q"
            worker_procs: option "-j"
            max_running_jobs: option "-J"
            sig_mode: option "-s"
            run_mode: "run" or "dryrun"
            verbosity: option "-v"
            trace_existing: option "-T"
    '''
    try:
        script = SoS_Script(script)
    except Exception as e:
        # show script with error
        raise ValueError(f'Failed to parse script {script}: {e}')
    #
    wf = script.workflow(workflow, use_default=not targets)
    run_config = {
        'config_file': None,
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
        'keep_going': False,
        'exec_mode': None
    }
    # a convenience feature
    if isinstance(args, dict):
        run_config['workflow_vars'] = args
        args = []
    run_config.update(config)
    executor = Base_Executor(wf, args=args, config=run_config)
    if isinstance(targets, str):
        targets = [targets]
    return executor.run(targets=targets)