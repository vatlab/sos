#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import copy
import os
import subprocess
import time
from collections import OrderedDict, Mapping, Sequence

from .eval import SoS_eval, SoS_exec, interpolate
from .monitor import ProcessMonitor
from .targets import (InMemorySignature, UnknownTarget, file_target,
                      remote, sos_step, dynamic, sos_targets)
from .utils import StopInputGroup, env, short_repr, pickleable
from .tasks import TaskFile, remove_task_files
from .step_executor import parse_shared_vars
from .executor_utils import __null_func__, get_traceback_msg

def collect_task_result(task_id, sos_dict, skipped=False, signature=None):
    shared = {}
    if 'shared' in env.sos_dict['_runtime']:
        svars = env.sos_dict['_runtime']['shared']
        if isinstance(svars, str):
            if svars not in env.sos_dict:
                raise ValueError(
                    f'Unavailable shared variable {svars} after the completion of task {task_id}')
            if not pickleable(env.sos_dict[svars], svars):
                env.logger.warning(
                    f'{svars} of type {type(env.sos_dict[svars])} is not sharable')
            else:
                shared[svars] = copy.deepcopy(env.sos_dict[svars])
        elif isinstance(svars, Mapping):
            for var, val in svars.items():
                if var != val:
                    env.sos_dict.set(var, SoS_eval(val))
                if var not in env.sos_dict:
                    raise ValueError(
                        f'Unavailable shared variable {var} after the completion of task {task_id}')
                if not pickleable(env.sos_dict[var], var):
                    env.logger.warning(
                        f'{var} of type {type(env.sos_dict[var])} is not sharable')
                else:
                    shared[var] = copy.deepcopy(env.sos_dict[var])
        elif isinstance(svars, Sequence):
            # if there are dictionaries in the sequence, e.g.
            # shared=['A', 'B', {'C':'D"}]
            for item in svars:
                if isinstance(item, str):
                    if item not in env.sos_dict:
                        raise ValueError(
                            f'Unavailable shared variable {item} after the completion of task {task_id}')
                    if not pickleable(env.sos_dict[item], item):
                        env.logger.warning(
                            f'{item} of type {type(env.sos_dict[item])} is not sharable')
                    else:
                        shared[item] = copy.deepcopy(env.sos_dict[item])
                elif isinstance(item, Mapping):
                    for var, val in item.items():
                        if var != val:
                            env.sos_dict.set(var, SoS_eval(val))
                        if var not in env.sos_dict:
                            raise ValueError(
                                f'Unavailable shared variable {var} after the completion of task {task_id}')
                        if not pickleable(env.sos_dict[var], var):
                            env.logger.warning(
                                f'{var} of type {type(env.sos_dict[var])} is not sharable')
                        else:
                            shared[var] = copy.deepcopy(env.sos_dict[var])
                else:
                    raise ValueError(
                        f'Option shared should be a string, a mapping of expression, or a list of string or mappings. {svars} provided')
        else:
            raise ValueError(
                f'Option shared should be a string, a mapping of expression, or a list of string or mappings. {svars} provided')
        env.logger.debug(
            f'task {task_id} (index={env.sos_dict["_index"]}) return shared variable {shared}')
    # the difference between sos_dict and env.sos_dict is that sos_dict (the original version) can have remote() targets
    # which should not be reported.
    if env.sos_dict['_output'] is None:
        output = {}
    elif env.sos_dict['_output'].undetermined():
        env.sos_dict.set('__null_func__', __null_func__)
        # re-process the output statement to determine output files
        args, _ = SoS_eval(
            f'__null_func__({env.sos_dict["_output"]._undetermined})')
        # handle dynamic args
        args = [x.resolve() if isinstance(x, dynamic) else x for x in args]
        output = {x: file_target(x).target_signature()
                  for x in sos_targets(*args)}
    elif sos_dict['_output'] is None:
        output = {}
    else:
        output = {x: file_target(x).target_signature()
                  for x in sos_dict['_output'] if isinstance(x, (str, file_target))}

    input = {} if env.sos_dict['_input'] is None or sos_dict['_input'] is None else {x: file_target(
        x).target_signature() for x in sos_dict['_input'] if isinstance(x, (str, file_target))}
    depends = {} if env.sos_dict['_depends'] is None or sos_dict['_depends'] is None else {
        x: file_target(x).target_signature() for x in sos_dict['_depends'] if isinstance(x, (str, file_target))}
    return {'ret_code': 0, 'task': task_id, 'input': input, 'output': output, 'depends': depends,
            'shared': shared, 'skipped': skipped,
            'start_time': sos_dict.get('start_time', ''),
            'peak_cpu': sos_dict.get('peak_cpu', 0),
            'peak_mem': sos_dict.get('peak_mem', 0),
            'end_time': time.time(),
            'signature': {task_id: signature.write()} if signature else {}
            }


def execute_task(task_id, verbosity=None, runmode='run', sigmode=None, monitor_interval=5,
                 resource_monitor_interval=60):
    tf = TaskFile(task_id)

    tf.status = 'running'
    # write result file
    res = _execute_task(task_id, verbosity, runmode, sigmode,
                        monitor_interval, resource_monitor_interval)

    if res['ret_code'] != 0 and 'exception' in res:
        with open(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err'), 'a') as err:
            err.write(f'Task {task_id} exits with code {res["ret_code"]}')

    if res.get('skipped', False):
        tf.status = 'skipped'
    else:
        tf.status = 'completed' if res['ret_code'] == 0 else 'failed'
        tf.add_outputs()
        sig = res.get('signature', {})
        res.pop('signature', None)
        tf.add_result(res)
        if sig:
            tf.add_signature(sig)

        # **after** result file is created, remove other files
        remove_task_files(task_id, ['.pulse', '.out', '.err', '.job_id', '.sh'])

    return res['ret_code']


def _validate_task_signature(sig, saved_sig, task_id):
    idx = env.sos_dict['_index']
    if env.config['sig_mode'] == 'default':
        matched = sig.validate(saved_sig)
        if isinstance(matched, dict):
            # in this case, an Undetermined output can get real output files
            # from a signature
            env.sos_dict.set('_input', sos_targets(matched['input']))
            env.sos_dict.set('_depends', sos_targets(matched['depends']))
            env.sos_dict.set('_output', sos_targets(matched['output']))
            env.sos_dict.update(matched['vars'])
            env.logger.info(
                f'Task ``{task_id}`` for substep ``{env.sos_dict["step_name"]}`` (index={idx}) is ``ignored`` due to saved signature')
            return True
    elif env.config['sig_mode'] == 'assert':
        matched = sig.validate(saved_sig)
        if isinstance(matched, str):
            raise RuntimeError(f'Signature mismatch: {matched}')
        else:
            env.sos_dict.set('_input', sos_targets(matched['input']))
            env.sos_dict.set('_depends', sos_targets(matched['depends']))
            env.sos_dict.set('_output', sos_targets(matched['output']))
            env.sos_dict.update(matched['vars'])
            env.logger.info(
                f'Task ``{task_id}`` for substep ``{env.sos_dict["step_name"]}`` (index={idx}) is ``ignored`` with matching signature')
            sig.content = saved_sig
            return True
    elif env.config['sig_mode'] == 'build':
        # The signature will be write twice
        if sig.write(rebuild=True):
            env.logger.info(
                f'Task ``{task_id}`` for substep ``{env.sos_dict["step_name"]}`` (index={idx}) is ``ignored`` with signature constructed')
            return True
        else:
            env.logger.info(
                f'Task ``{task_id}`` for substep ``{env.sos_dict["step_name"]}`` (index={idx}) is ``executed`` with failed signature constructed')
            return False
    elif env.config['sig_mode'] == 'force':
        return False
    else:
        raise RuntimeError(
            f'Unrecognized signature mode {env.config["sig_mode"]}')


def _execute_sub_tasks(task_id, params, sig_content, verbosity, runmode, sigmode,
                       monitor_interval, resource_monitor_interval):
    '''If this is a master task, execute as individual tasks'''
    m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
                       resource_monitor_interval=resource_monitor_interval,
                       max_walltime=params.sos_dict['_runtime'].get(
                           'max_walltime', None),
                       max_mem=params.sos_dict['_runtime'].get(
                           'max_mem', None),
                       max_procs=params.sos_dict['_runtime'].get(
                           'max_procs', None),
                       sos_dict=params.sos_dict)
    m.start()

    master_out = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task_id + '.out')
    master_err = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task_id + '.err')
    # if this is a master task, calling each sub task
    with open(master_out, 'wb') as out, open(master_err, 'wb') as err:
        def copy_out_and_err(result):
            tid = result['task']
            out.write(
                f'{tid}: {"completed" if result["ret_code"] == 0 else "failed"}\n'.encode())
            if 'output' in result:
                out.write(f'output: {result["output"]}\n'.encode())
            sub_out = os.path.join(os.path.expanduser(
                '~'), '.sos', 'tasks', tid + '.out')
            if os.path.isfile(sub_out):
                with open(sub_out, 'rb') as sout:
                    out.write(sout.read())
                try:
                    os.remove(sub_out)
                except Exception as e:
                    env.logger.warning(f'Failed to remove {sub_out}: {e}')

            sub_err = os.path.join(os.path.expanduser(
                '~'), '.sos', 'tasks', tid + '.err')
            err.write(
                f'{tid}: {"completed" if result["ret_code"] == 0 else "failed"}\n'.encode())
            if os.path.isfile(sub_err):
                with open(sub_err, 'rb') as serr:
                    err.write(serr.read())
                try:
                    os.remove(sub_err)
                except Exception as e:
                    env.logger.warning(f'Failed to remove {sub_err}: {e}')

            # remove other files as well
            try:
                remove_task_files(tid, ['.out', '.err'])
            except Exception as e:
                env.logger.debug(f'Failed to remove files {tid}: {e}')

        if params.num_workers > 1:
            from multiprocessing.pool import Pool
            p = Pool(params.num_workers)
            results = []
            for t in params.task_stack:
                results.append(p.apply_async(_execute_task,
                                             ((*t, {t[0]: sig_content.get(t[0], {})}), verbosity, runmode,
                                              sigmode, None, None), callback=copy_out_and_err))
            for idx, r in enumerate(results):
                results[idx] = r.get()
            p.close()
            p.join()
            # we wait for all results to be ready to return or raise
            # but we only raise exception for one of the subtasks
            for res in results:
                if 'exception' in res:
                    failed = [x.get("task", "")
                              for x in results if "exception" in x]
                    env.logger.error(
                        f'{task_id} ``failed`` due to failure of subtask{"s" if len(failed) > 1 else ""} {", ".join(failed)}')
                    return {'ret_code': 1, 'exception': res['exception'], 'task': task_id}
        else:
            results = []
            for tid, tdef in params.task_stack:
                # no monitor process for subtasks
                res = _execute_task((tid, tdef, {tid: sig_content.get(tid, {})}), verbosity=verbosity, runmode=runmode,
                                    sigmode=sigmode, monitor_interval=None, resource_monitor_interval=None)
                copy_out_and_err(res)
                results.append(res)
            for res in results:
                if 'exception' in res:
                    failed = [x.get("task", "")
                              for x in results if "exception" in x]
                    env.logger.error(
                        f'{task_id} ``failed`` due to failure of subtask{"s" if len(failed) > 1 else ""} {", ".join(failed)}')
                    return {'ret_code': 1, 'exception': res['exception'], 'task': task_id}
    #
    # now we collect result
    all_res = {'ret_code': 0, 'output': {},
               'subtasks': {}, 'shared': {}, 'skipped': False, 'signature': {}}
    for tid, x in zip(params.task_stack, results):
        all_res['ret_code'] += x['ret_code']
        all_res['output'].update(x['output'])
        all_res['subtasks'][tid[0]] = x
        all_res['shared'].update(x['shared'])
        # does not care if one or all subtasks are executed or skipped.
        all_res['skipped'] = x['skipped']
        all_res['signature'].update(x['signature'])

    return all_res


def _execute_task(task_id, verbosity=None, runmode='run', sigmode=None, monitor_interval=5,
                  resource_monitor_interval=60):
    '''A function that execute specified task within a local dictionary
    (from SoS env.sos_dict). This function should be self-contained in that
    it can be handled by a task manager, be executed locally in a separate
    process or remotely on a different machine.'''
    # start a monitoring file, which would be killed after the job
    # is done (killed etc)
    if isinstance(task_id, str):
        params = TaskFile(task_id).params
        sig_content = TaskFile(task_id).signature
        subtask = False
    else:
        # subtask
        subtask = True
        task_id, params, sig_content = task_id
        env.logger.trace(f'Executing subtask {task_id}')

    if hasattr(params, 'task_stack'):
        return _execute_sub_tasks(task_id, params, sig_content, verbosity, runmode, sigmode,
                                  monitor_interval, resource_monitor_interval)

    global_def, task, sos_dict = params.global_def, params.task, params.sos_dict

    # task output
    env.sos_dict.set('__std_out__', os.path.join(
        os.path.expanduser('~'), '.sos', 'tasks', task_id + '.out'))
    env.sos_dict.set('__std_err__', os.path.join(
        os.path.expanduser('~'), '.sos', 'tasks', task_id + '.err'))
    env.logfile = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task_id + '.err')
    # clear the content of existing .out and .err file if exists, but do not create one if it does not exist
    if os.path.exists(env.sos_dict['__std_out__']):
        open(env.sos_dict['__std_out__'], 'w').close()
    if os.path.exists(env.sos_dict['__std_err__']):
        open(env.sos_dict['__std_err__'], 'w').close()

    if verbosity is not None:
        env.verbosity = verbosity
    try:
        # global def could fail due to execution on remote host...
        # we also execute global_def way before others and allows variables set by
        # global_def be overwritten by other passed variables
        #
        # note that we do not handle parameter in tasks because values should already be
        # in sos_task dictionary
        SoS_exec('''\
import os, sys, glob
from sos.runtime import *
CONFIG = {}
del sos_handle_parameter_
''' + global_def)
    except Exception as e:
        env.logger.trace(
            f'Failed to execute global definition {short_repr(global_def)}: {e}')

    if '_runtime' not in sos_dict:
        sos_dict['_runtime'] = {}

    # pulse thread
    if monitor_interval is not None:
        m = ProcessMonitor(task_id, monitor_interval=monitor_interval,
                           resource_monitor_interval=resource_monitor_interval,
                           max_walltime=sos_dict['_runtime'].get(
                               'max_walltime', None),
                           max_mem=sos_dict['_runtime'].get('max_mem', None),
                           max_procs=sos_dict['_runtime'].get(
                               'max_procs', None),
                           sos_dict=sos_dict)

        m.start()
    env.config['run_mode'] = runmode
    if runmode == 'dryrun':
        env.config['sig_mode'] = 'ignore'
    elif sigmode is not None:
        env.config['sig_mode'] = sigmode
    #
    if subtask:
        env.logger.debug(f'{task_id} ``started``')
    else:
        env.logger.info(f'{task_id} ``started``')

    env.sos_dict.quick_update(sos_dict)

    # if targets are defined as `remote`, they should be resolved during task execution
    def resolve_remote(x):
        if isinstance(x, remote):
            x = x.resolve()
            if isinstance(x, str):
                x = interpolate(x, env.sos_dict._dict)
        return x

    for key in ['step_input', '_input',  'step_output', '_output', 'step_depends', '_depends']:
        if key in sos_dict and isinstance(sos_dict[key], (list, sos_targets)):
            # resolve remote() target
            env.sos_dict.set(key, sos_targets(resolve_remote(x)
                                              for x in sos_dict[key] if not isinstance(x, sos_step)))

    sig = None if env.config['sig_mode'] == 'ignore' or env.sos_dict['_output'].unspecified() else InMemorySignature(
        env.sos_dict['_input'], env.sos_dict['_output'],
        env.sos_dict['_depends'], env.sos_dict['__signature_vars__'],
        shared_vars=parse_shared_vars(env.sos_dict['_runtime'].get('shared', None)))

    if sig and _validate_task_signature(sig, sig_content.get(task_id, {}), task_id):
        #env.logger.info(f'{task_id} ``skipped``')
        return collect_task_result(task_id, sos_dict, skipped=True, signature=sig)

    # if we are to really execute the task, touch the task file so that sos status shows correct
    # execution duration.
    if not subtask:
        sos_dict['start_time'] = time.time()

    try:
        # go to 'cur_dir'
        if '_runtime' in sos_dict and 'cur_dir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['cur_dir'])):
                try:
                    os.makedirs(os.path.expanduser(
                        sos_dict['_runtime']['cur_dir']))
                    os.chdir(os.path.expanduser(
                        sos_dict['_runtime']['cur_dir']))
                except Exception as e:
                    # sometimes it is not possible to go to a "cur_dir" because of
                    # file system differences, but this should be ok if a work_dir
                    # has been specified.
                    env.logger.debug(
                        f'Failed to create cur_dir {sos_dict["_runtime"]["cur_dir"]}: {e}')
            else:
                os.chdir(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
        #
        orig_dir = os.getcwd()

        # we will need to check existence of targets because the task might
        # be executed on a remote host where the targets are not available.
        for target in (sos_dict['_input'] if isinstance(sos_dict['_input'], list) else []) + \
                (sos_dict['_depends'] if isinstance(sos_dict['_depends'], list) else []):
            # if the file does not exist (although the signature exists)
            # request generation of files
            if isinstance(target, str):
                if not file_target(target).target_exists('target'):
                    # remove the signature and regenerate the file
                    raise UnknownTarget(target)
            # the sos_step target should not be checked in tasks because tasks are
            # independently executable units.
            elif not isinstance(target, sos_step) and not target.target_exists('target'):
                raise UnknownTarget(target)

        # create directory. This usually has been done at the step level but the task can be executed
        # on a remote host where the directory does not yet exist.
        ofiles = env.sos_dict['_output']
        if ofiles.valid():
            for ofile in ofiles:
                parent_dir = ofile.parent
                if not parent_dir.is_dir():
                    parent_dir.mkdir(parents=True, exist_ok=True)

                    # go to user specified workdir
        if '_runtime' in sos_dict and 'workdir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['workdir'])):
                try:
                    os.makedirs(os.path.expanduser(
                        sos_dict['_runtime']['workdir']))
                except Exception as e:
                    raise RuntimeError(
                        f'Failed to create workdir {sos_dict["_runtime"]["workdir"]}: {e}')
            os.chdir(os.path.expanduser(sos_dict['_runtime']['workdir']))
        # set environ ...
        # we join PATH because the task might be executed on a different machine
        if '_runtime' in sos_dict:
            if 'env' in sos_dict['_runtime']:
                for key, value in sos_dict['_runtime']['env'].items():
                    if 'PATH' in key and key in os.environ:
                        new_path = OrderedDict()
                        for p in value.split(os.pathsep):
                            new_path[p] = 1
                        for p in value.split(os.environ[key]):
                            new_path[p] = 1
                        os.environ[key] = os.pathsep.join(new_path.keys())
                    else:
                        os.environ[key] = value
            if 'prepend_path' in sos_dict['_runtime']:
                if isinstance(sos_dict['_runtime']['prepend_path'], str):
                    os.environ['PATH'] = sos_dict['_runtime']['prepend_path'] + \
                        os.pathsep + os.environ['PATH']
                elif isinstance(env.sos_dict['_runtime']['prepend_path'], Sequence):
                    os.environ['PATH'] = os.pathsep.join(
                        sos_dict['_runtime']['prepend_path']) + os.pathsep + os.environ['PATH']
                else:
                    raise ValueError(
                        f'Unacceptable input for option prepend_path: {sos_dict["_runtime"]["prepend_path"]}')

        # step process
        SoS_exec(task)

        if subtask:
            env.logger.debug(f'{task_id} ``completed``')
        else:
            env.logger.info(f'{task_id} ``completed``')

    except StopInputGroup as e:
        # task ignored with stop_if exception
        if e.message:
            env.logger.warning(f'{task_id} ``stopped``: {e.message}')
        return {'ret_code': 0, 'task': task_id, 'input': [],
                'output': [], 'depends': [], 'shared': {}}
    except KeyboardInterrupt:
        env.logger.error(f'{task_id} ``interrupted``')
        raise
    except subprocess.CalledProcessError as e:
        return {'ret_code': e.returncode, 'task': task_id, 'shared': {},
                'exception': RuntimeError(e.stderr)}
    except Exception as e:
        msg = get_traceback_msg(e)
        env.logger.error(f'{task_id} ``failed``: {msg}')
        return {'ret_code': 1, 'exception': RuntimeError(msg), 'task': task_id, 'shared': {}}
    finally:
        os.chdir(orig_dir)

    # the final result should be relative to cur_dir, not workdir
    # because output is defined outside of task. This is the only place
    # when signature matters (worth saving)
    return collect_task_result(task_id, sos_dict, signature=sig)
