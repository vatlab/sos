#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import subprocess
import sys
import contextlib
import zmq

from io import StringIO

from .eval import SoS_exec
from .targets import (RemovedTarget, RuntimeInfo, UnavailableLock,
                      UnknownTarget)
from .executor_utils import (prepare_env, clear_output, verify_input, kill_all_subprocesses,
        reevaluate_output, validate_step_sig, create_task, get_traceback_msg, statementMD5)

from .utils import (StopInputGroup, TerminateExecution, ArgumentError, env)


@contextlib.contextmanager
def stdoutIO():
    oldout = sys.stdout
    olderr = sys.stderr
    stdout = StringIO()
    stderr = StringIO()
    sys.stdout = stdout
    sys.stderr = stderr
    yield stdout, stderr
    sys.stdout = oldout
    sys.stderr = olderr


def execute_substep(stmt, global_def='', task='', proc_vars={}, shared_vars=[], config={}):
    '''Execute a substep with specific input etc

    Substep executed by this function should be self-contained. It can contain
    tasks (which will be sent to the master process) but not nested workflows.

    The executor checks step signatures and might skip the substep if it has
    been executed and the signature matches.

    The executor accepts connections to the controller, and a socket using
    which the results will be returned. However, the calling process should
    take care of the connection and disconnection of controller sockets and
    this function only takes care of the connection and disconnection of
    result socket.

    stmt:
        Main statement of the substep

    global_def:
        Global definitions, might define functions useful to the substep

    task:
        External task

    proc_vars:
        Environmental variables, signature variables etc

    shared_vars:
        Variables that should be returned after the execution

    config:
        Runmode, signature mode, verbosity, etc.

    The return value should be a dictionary with the following keys:

    index: index of the substep within the step
    ret_code: (all) return code, 0 for successful
    sig_skipped: (optional) return if the step is skipped due to signature
    shared: (optional) shared variable as specified by 'shared_vars'
    stdout: (optional) if in interactive mode
    stderr: (optional) if in interactive mode
    exception: (optional) if an exception occures
    '''
    assert not env.zmq_context.closed
    assert not env.controller_push_socket.closed
    assert not env.controller_req_socket.closed
    assert not env.signature_push_socket.closed
    assert not env.signature_req_socket.closed
    assert 'workflow_id' in proc_vars
    assert 'step_id' in proc_vars
    assert '_input' in proc_vars
    assert '_output' in proc_vars
    assert '_depends' in proc_vars
    assert 'step_output' in proc_vars
    assert '_index' in proc_vars
    assert 'result_push_socket' in config["sockets"]

    try:
        res_socket = env.zmq_context.socket(zmq.PUSH)
        res_socket.connect(f'tcp://127.0.0.1:{config["sockets"]["result_push_socket"]}')
        res = _execute_substep(stmt=stmt, global_def=global_def, task=task, proc_vars=proc_vars,
            shared_vars=shared_vars, config=config)
        res_socket.send_pyobj(res)
    finally:
        res_socket.close()

def _execute_substep(stmt, global_def, task, proc_vars, shared_vars, config):
    # passing configuration and port numbers to the subprocess
    env.config.update(config)
    # prepare a working environment with sos symbols and functions
    prepare_env(global_def)
    # update it with variables passed from master process
    env.sos_dict.quick_update(proc_vars)
    if env.config['sig_mode'] == 'ignore' or env.sos_dict['_output'].unspecified():
        sig = None
    else:
        sig = RuntimeInfo(
            statementMD5([stmt, task]),
            env.sos_dict['_input'],
            env.sos_dict['_output'],
            env.sos_dict['_depends'],
            env.sos_dict['__signature_vars__'],
            shared_vars=shared_vars)
    outmsg = ''
    errmsg = ''
    capture_output = env.config['run_mode'] == 'interactive'
    try:
        if sig:
            matched = validate_step_sig(sig)
            if matched:
                # avoid sig being released in the final statement
                sig = None
                # complete case: concurrent ignore without task
                env.controller_push_socket.send_pyobj(['progress', 'substep_ignored', env.sos_dict['step_id']])
                return {'index': env.sos_dict['_index'], 'ret_code': 0, 'sig_skipped': 1, 'output': matched['output'],
                    'shared': matched['vars']}
            sig.lock()

        # check if input and depends targets actually exist
        #
        # if depends on a sos_variable but the variable is not actually used in
        # the substep, it is ok to ignore it. If the variable is used in the substep
        # it should have been included as part of the signature variables.
        verify_input(ignore_internal_targets=True)

        if stmt:
            # statement can be empty for task only substep
            if capture_output:
                with stdoutIO() as (out, err):
                    SoS_exec(stmt, return_result=False)
                    outmsg = out.getvalue()
                    errmsg = err.getvalue()
            else:
                SoS_exec(stmt, return_result=False)
        if task:
            task_id, taskdef, task_vars = create_task(global_def, task)
            res = {'index': env.sos_dict['_index'], 'task_id': task_id, 'task_def': taskdef, 'task_vars': task_vars}
        else:
            if env.sos_dict['step_output'].undetermined():
                # the pool worker does not have __null_func__ defined
                env.sos_dict.set('_output', reevaluate_output())
            res = {'index': env.sos_dict['_index'], 'ret_code': 0}
            if sig:
                sig.set_output(env.sos_dict['_output'])
                # sig.write will use env.signature_push_socket
                if sig.write():
                    res.update({'output': list(sig.content['output'].keys()), 'shared': sig.content['end_context']})
            if capture_output:
                res.update({'stdout': outmsg, 'stderr': errmsg})
            # complete case: concurrent execution without task
            env.controller_push_socket.send_pyobj(['progress', 'substep_completed', env.sos_dict['step_id']])
        return res
    except (StopInputGroup, TerminateExecution, UnknownTarget, RemovedTarget, UnavailableLock) as e:
        clear_output()
        res = {'index': env.sos_dict['_index'], 'ret_code': 1, 'exception': e}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    except (KeyboardInterrupt, SystemExit) as e:
        clear_output()
        kill_all_subprocesses()
        raise e
    except subprocess.CalledProcessError as e:
        clear_output()
        # cannot pass CalledProcessError back because it is not pickleable
        res = {'index': env.sos_dict['_index'], 'ret_code': e.returncode, 'exception': RuntimeError(e.stderr)}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    except ArgumentError as e:
        clear_output()
        return {'index': env.sos_dict['_index'], 'ret_code': 1, 'exception': e}
    except Exception as e:
        clear_output()
        res = {'index': env.sos_dict['_index'], 'ret_code': 1,
            'exception': RuntimeError(get_traceback_msg(e))}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    finally:
        # release the lock even if the process becomes zombie? #871
        if sig:
            sig.release(quiet=True)
