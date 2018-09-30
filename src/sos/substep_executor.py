#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.


import os
import subprocess
import sys
import traceback
import contextlib
import psutil

from io import StringIO

from .eval import SoS_exec, stmtHash
from .controller import connect_controllers
from .targets import (RemovedTarget, RuntimeInfo, UnavailableLock,
                      UnknownTarget)
from .executor_utils import (__null_func__, clear_output,
                      verify_input, reevaluate_output,
                     validate_step_sig)

from .utils import (StopInputGroup, TerminateExecution, ArgumentError, env)


# overwrite concurrent_execute defined in Base_Step_Executor because sos notebook
# can only handle stdout/stderr from the master process
#
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

def execute_substep(stmt, proc_vars={}, step_md5=None, step_tokens=[],
    shared_vars=[], config={}, capture_output=False):
    '''Execute statements in the passed dictionary'''
    # passing configuration and port numbers to the subprocess
    env.config.update(config)
    connect_controllers()
    # prepare a working environment with sos symbols and functions
    from ._version import __version__
    env.sos_dict.set('__null_func__', __null_func__)
    # initial values
    env.sos_dict.set('SOS_VERSION', __version__)
    SoS_exec('import os, sys, glob', None)
    SoS_exec('from sos.runtime import *', None)
    # update it with variables passed from master process
    env.sos_dict.quick_update(proc_vars)
    sig = None if env.config['sig_mode'] == 'ignore' or env.sos_dict['_output'].unspecified() else RuntimeInfo(
        step_md5, step_tokens,
        env.sos_dict['_input'],
        env.sos_dict['_output'],
        env.sos_dict['_depends'],
        env.sos_dict['__signature_vars__'],
        shared_vars=shared_vars)
    outmsg = ''
    errmsg = ''
    try:
        if sig:
            matched = validate_step_sig(sig)
            if matched:
                # avoid sig being released in the final statement
                sig = None
                # complete case: concurrent ignore without task
                env.controller_push_socket.send_pyobj(['progress', 'substep_ignored', env.sos_dict['step_id']])
                return {'ret_code': 0, 'sig_skipped': 1, 'output': matched['output'], 'shared': matched['vars']}
            sig.lock()
        verify_input()

        if capture_output:
            with stdoutIO() as (out, err):
                SoS_exec(stmt, return_result=False)
                outmsg = out.getvalue()
                errmsg = err.getvalue()
        else:
            SoS_exec(stmt, return_result=False)
        if env.sos_dict['step_output'].undetermined():
            # the pool worker does not have __null_func__ defined
            env.sos_dict.set('_output', reevaluate_output())
        res = {'ret_code': 0}
        if sig:
            sig.set_output(env.sos_dict['_output'])
            if sig.write():
                res.update({'output': sig.content['output'], 'shared': sig.content['end_context']})
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        # complete case: concurrent execution without task
        env.controller_push_socket.send_pyobj(['progress', 'substep_completed', env.sos_dict['step_id']])
        return res
    except (StopInputGroup, TerminateExecution, UnknownTarget, RemovedTarget, UnavailableLock) as e:
        clear_output()
        res = {'ret_code': 1, 'exception': e}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    except (KeyboardInterrupt, SystemExit) as e:
        clear_output()
        # Note that KeyboardInterrupt is not an instance of Exception so this piece is needed for
        # the subprocesses to handle keyboard interrupt. We do not pass the exception
        # back to the master process because the master process would handle KeyboardInterrupt
        # as well and has no chance to handle the returned code.
        procs = psutil.Process().children(recursive=True)
        if procs:
            if env.verbosity > 2:
                env.logger.info(
                    f'{os.getpid()} interrupted. Killing subprocesses {" ".join(str(x.pid) for x in procs)}')
            for p in procs:
                p.terminate()
            gone, alive = psutil.wait_procs(procs, timeout=3)
            if alive:
                for p in alive:
                    p.kill()
            gone, alive = psutil.wait_procs(procs, timeout=3)
            if alive:
                for p in alive:
                    env.logger.warning(f'Failed to kill subprocess {p.pid}')
        elif env.verbosity > 2:
            env.logger.info(f'{os.getpid()} interrupted. No subprocess.')
        raise e
    except subprocess.CalledProcessError as e:
        clear_output()
        # cannot pass CalledProcessError back because it is not pickleable
        res = {'ret_code': e.returncode, 'exception': RuntimeError(e.stderr)}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    except ArgumentError as e:
        clear_output()
        return {'ret_code': 1, 'exception': e}
    except Exception as e:
        clear_output()
        error_class = e.__class__.__name__
        cl, exc, tb = sys.exc_info()
        msg = ''
        for st in reversed(traceback.extract_tb(tb)):
            if st.filename.startswith('script_'):
                code = stmtHash.script(st.filename)
                line_number = st.lineno
                code = '\n'.join([f'{"---->" if i+1 == line_number else "     "} {x.rstrip()}' for i,
                                  x in enumerate(code.splitlines())][max(line_number - 3, 0):line_number + 3])
                msg += f'''\
{st.filename} in {st.name}
{code}
'''
        detail = e.args[0] if e.args else ''
        res = {'ret_code': 1, 'exception': RuntimeError(f'''
---------------------------------------------------------------------------
{error_class:42}Traceback (most recent call last)
{msg}
{error_class}: {detail}''') if msg else RuntimeError(f'{error_class}: {detail}')}
        if capture_output:
            res.update({'stdout': outmsg, 'stderr': errmsg})
        return res
    finally:
        # release the lock even if the process becomes zombie? #871
        if sig:
            sig.release(quiet=True)
