#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import concurrent.futures
import copy
import os
import fasteners
import pickle
import random
import subprocess
import threading
import time
from collections import OrderedDict
from stat import S_IREAD, S_IRGRP, S_IROTH
from typing import Union, Dict

from .eval import cfg_interpolate
from .utils import (env, expand_time, linecount_of_file, sample_of_file,
                    short_repr, tail_of_file)
from .task_executor import monitor_interval, taskTags, taskDuration, loadTask

# return {} if result is unchanged
# otherwise return a dictionary of with keys 'status' and 'files'


def check_task(task, hint={}) -> Dict[str, Union[str, Dict[str, float]]]:
    #
    # for pending status, a new file might have been created to change its status
    # so we cannot use cached signature.
    #
    # for running status, static pulse file actually means something wrong.
    #
    # when testing. if the timestamp is 0, the file does not exist originally, it should
    # still does not exist. Otherwise the file should exist and has the same timestamp
    if hint and hint['status'] not in ('pending', 'running') and \
            all((os.path.isfile(f) and os.stat(f).st_mtime == v) if v else (not os.path.isfile(f)) for f, v in hint['files'].items()):
        return {}

    # status of the job, please refer to https://github.com/vatlab/SOS/issues/529
    # for details.
    #
    task_file = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task + '.task')
    if not os.path.isfile(task_file):
        return dict(status='missing', files={task_file: 0})
    pulse_file = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task + '.pulse')
    if not os.path.isfile(pulse_file):
        pulse_file = os.path.join(os.path.expanduser(
            '~'), '.sos', 'tasks', task + '.status')

    def has_pulse():
        # for whatever reason, sometimes the pulse file might appear to be slightly
        # before the task file, and if the task is very short so the pulse file is
        # not updated, the task will appear to be in pending mode forever
        return os.path.isfile(pulse_file) and os.stat(pulse_file).st_mtime >= os.stat(task_file).st_mtime - 1

    res_file = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task + '.res')

    def has_res():
        return os.path.isfile(res_file) and os.stat(res_file).st_mtime >= os.stat(task_file).st_mtime

    job_file = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task + '.sh')

    def has_job():
        job_id_file = os.path.join(os.path.expanduser(
            '~'), '.sos', 'tasks', task + '.job_id')
        return os.path.isfile(job_file) and os.stat(job_file).st_mtime >= os.stat(task_file).st_mtime \
            and os.path.isfile(job_id_file) and os.stat(job_id_file).st_mtime >= os.stat(job_file).st_mtime

    if has_res():
        try:
            from .targets import file_target
            with open(res_file, 'rb') as result:
                res = pickle.load(result)
            status_files = {task_file: os.stat(task_file).st_mtime,
                            res_file: os.stat(res_file).st_mtime,
                            pulse_file: os.stat(pulse_file).st_mtime if os.path.isfile(
                                pulse_file) else 0
                            }

            if ('ret_code' in res and res['ret_code'] == 0) or ('succ' in res and res['succ'] == 0):
                for var in ('input', 'output', 'depends'):
                    if var not in res or not isinstance(res[var], dict):
                        continue
                    for x, y in res[var].items():
                        if not file_target(x).target_exists() or file_target(x).target_signature() != y:
                            env.logger.debug(
                                f'{x} not found or signature mismatch')
                            return dict(status='signature-mismatch',
                                        files=status_files)
                return dict(status='completed', files=status_files)
            else:
                return dict(status='failed', files=status_files)
        except Exception as e:
            # sometimes the resfile is changed while we are reading it
            # so we wait a bit and try again.
            env.logger.warning(e)
            time.sleep(.5)
            return check_task(task)
    #
    if has_pulse():
        status_files = {task_file: os.stat(task_file).st_mtime,
                        pulse_file: os.stat(pulse_file).st_mtime}
        # dead?
        # if the status file is readonly
        if not os.access(pulse_file, os.W_OK):
            return dict(status='aborted', files={task_file: os.stat(task_file).st_mtime,
                                                 pulse_file: os.stat(pulse_file).st_mtime})
        start_stamp = os.stat(pulse_file).st_mtime
        elapsed = time.time() - start_stamp
        if elapsed < 0:
            env.logger.debug(
                f'{pulse_file} is created in the future. Your system time might be problematic')
        # if the file is within 5 seconds
        if elapsed < monitor_interval:
            # if running, we return old hint files even if the timestamp has been changed
            # because we will check the status of running jobs anyway.
            if hint and hint['status'] == 'running':
                return {}
            else:
                return dict(status='running', files=status_files)
        elif elapsed > 2 * monitor_interval:
            if has_res():
                # result file appears during sos tatus run
                return check_task(task)
            else:
                return dict(status='aborted', files=status_files)
        # otherwise, let us be patient ... perhaps there is some problem with the filesystem etc
        time.sleep(2 * monitor_interval)
        end_stamp = os.stat(pulse_file).st_mtime
        # the process is still alive
        if has_res():
            return check_task(task)
        elif start_stamp != end_stamp:
            if hint and hint['status'] == 'running':
                return {}
            else:
                return dict(status='running', files=status_files)
        else:
            return dict(status='aborted', files=status_files)
    # if there is no status file
    if has_job():
        return dict(status='submitted', files={task_file: os.stat(task_file).st_mtime,
                                               job_file: os.stat(job_file).st_mtime,
                                               pulse_file: 0})
    else:
        # status not changed
        if hint and hint['status'] == 'pending' and hint['files'][task_file] == os.stat(task_file).st_mtime:
            return {}
        else:
            return dict(status='pending', files={task_file: os.stat(task_file).st_mtime,
                                                 job_file: 0})


def print_task_status(tasks, verbosity: int=1, html: bool=False, start_time=False, age=None, tags=None, status=None):
    # verbose is ignored for now
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    check_all: bool = not tasks
    if check_all:
        tasks = glob.glob(os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getmtime(x))
                     for x in tasks]
        if not all_tasks:
            return
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'),
                                             '.sos', 'tasks', f'{t}*.task'))
            matched = [(os.path.basename(x)[:-5], os.path.getmtime(x))
                       for x in matched]
            if not matched:
                all_tasks.append((t, None))
            else:
                all_tasks.extend(matched)
    if age is not None:
        age = expand_time(age, default_unit='d')
        if age > 0:
            all_tasks = [x for x in all_tasks if time.time() - x[1] >= age]
        else:
            all_tasks = [x for x in all_tasks if time.time() - x[1] <= -age]

    all_tasks = sorted(list(set(all_tasks)),
                       key=lambda x: 0 if x[1] is None else x[1])

    if tags:
        all_tasks = [x for x in all_tasks if any(
            x in tags for x in taskTags(x[0]).split(' '))]

    if not all_tasks:
        env.logger.info('No matching tasks')
        return
    cache_file: str = os.path.join(
        os.path.expanduser('~'), '.sos', 'tasks', 'status_cache.pickle')
    #
    status_cache = {}
    if os.path.isfile(cache_file):
        with fasteners.InterProcessLock(cache_file + '_'):
            with open(cache_file, 'rb') as cache:
                status_cache = pickle.load(cache)
    #
    # at most 20 threads
    p = Pool(min(20, len(all_tasks)))
    # the result can be {} for unchanged, or real results
    raw_status = p.starmap(
        check_task, [(x[0], status_cache.get(x[0], {})) for x in all_tasks])
    # if check all, we clear the cache and record all existing tasks
    has_changes: bool = any(x for x in raw_status)
    if has_changes:
        if check_all:
            status_cache = {k[0]: v if v else status_cache[k[0]]
                            for k, v in zip(all_tasks, raw_status)}
        else:
            status_cache.update(
                {k[0]: v for k, v in zip(all_tasks, raw_status) if v})
    obtained_status = [status_cache[x[0]]['status'] for x in all_tasks]
    if status:
        all_tasks = [x for x, s in zip(
            all_tasks, obtained_status) if s in status]
        obtained_status = [x for x in obtained_status if x in status]
    #
    # automatically remove non-running tasks that are more than 30 days old
    to_be_removed = []
    #
    if html:
        verbosity = -1
    if verbosity == 0:
        print('\n'.join(obtained_status))
    elif verbosity == 1:
        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30 * 24 * 60 * 60 and s != 'running':
                to_be_removed.append(t)
                continue
            print(f'{t}\t{s}')
    elif verbosity == 2:
        from .utils import PrettyRelativeTime
        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30 * 24 * 60 * 60 and s != 'running':
                to_be_removed.append(t)
                continue
            if start_time:
                if d is None:
                    print(f'{t}\t{taskTags(t)}\t{time.time()}\t{s}')
                else:
                    print(f'{t}\t{taskTags(t)}\t{d}\t{s}')
            else:
                if d is None:
                    print(f'{t}\t{taskTags(t)}\t{"":>15}\t{s}')
                elif s in ('pending', 'submitted', 'running'):
                    print(
                        f'{t}\t{taskTags(t)}\t{PrettyRelativeTime(time.time() - d):>15}\t{s}')
                else:
                    # completed or failed
                    print(
                        f'{t}\t{taskTags(t)}\t{PrettyRelativeTime(taskDuration(t)):>15}\t{s}')
    elif verbosity > 2:
        from .utils import PrettyRelativeTime
        import pprint
        from .monitor import summarizeExecution

        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30 * 24 * 60 * 60 and s != 'running':
                to_be_removed.append(t)
                continue
            print(f'{t}\t{s}\n')
            task_file = os.path.join(os.path.expanduser(
                '~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                continue
            if d is not None:
                print(f'Started {PrettyRelativeTime(time.time() - d)} ago')
            if s not in ('pending', 'submitted', 'running'):
                print(f'Duration {PrettyRelativeTime(taskDuration(t))}')
            params = loadTask(task_file)
            print('TASK:\n=====')
            print(params.task)
            print('TAGS:\n=====')
            print(' '.join(params.tags))
            print()
            if params.global_def:
                print('GLOBAL:\n=======')
                print(params.global_def)
                print()
            print('ENVIRONMENT:\n============')
            job_vars = params.sos_dict
            for k in sorted(job_vars.keys()):
                v = job_vars[k]
                print(
                    f'{k:22}{short_repr(v) if verbosity == 3 else pprint.pformat(v)}')
            print()
            print('EXECUTION STATS:\n================')
            print(summarizeExecution(t, status=s))
            if verbosity == 4:
                # if there are other files such as job file, print them.
                files = glob.glob(os.path.join(
                    os.path.expanduser('~'), '.sos', 'tasks', t + '.*'))
                for f in sorted([x for x in files if os.path.splitext(x)[-1] not in ('.res',
                                                                                     '.task', '.pulse', '.status', '.def')]):
                    print(
                        f'{os.path.basename(f)}:\n{"="*(len(os.path.basename(f))+1)}')
                    try:
                        with open(f) as fc:
                            print(fc.read())
                    except Exception:
                        print('Binary file')
    else:
        # HTML output
        from .utils import PrettyRelativeTime, isPrimitive
        from .monitor import summarizeExecution
        import pprint
        print('<table width="100%" class="resource_table">')

        def row(th=None, td=None):
            if td is None:
                print(
                    f'<tr><th align="right" width="30%">{th}</th><td></td></tr>')
            elif th is None:
                print(
                    f'<tr><td colspan="2" align="left"  width="30%">{td}</td></tr>')
            else:
                print(
                    f'<tr><th align="right"  width="30%">{th}</th><td align="left"><div class="one_liner">{td}</div></td></tr>')
        for s, (t, d) in zip(obtained_status, all_tasks):
            if d is not None and time.time() - d > 30 * 24 * 60 * 60:
                to_be_removed.append(t)
                continue
            row('ID', t)
            row('Status', s)
            task_file = os.path.join(os.path.expanduser(
                '~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                print('</table>')
                continue
            if d is not None:
                row('Start', f'{PrettyRelativeTime(time.time() - d):>15} ago')
            if s not in ('pending', 'submitted', 'running'):
                row('Duration', f'{PrettyRelativeTime(taskDuration(t)):>15}')
            task_file = os.path.join(os.path.expanduser(
                '~'), '.sos', 'tasks', t + '.task')
            if not os.path.isfile(task_file):
                continue
            params = loadTask(task_file)
            row('Task')
            row(td=f'<pre style="text-align:left">{params.task}</pre>')
            row('Tags')
            row(td=f'<pre style="text-align:left">{" ".join(params.tags)}</pre>')
            if params.global_def:
                row('Global')
                row(
                    td=f'<pre style="text-align:left">{params.global_def}</pre>')
            # row('Environment')
            job_vars = params.sos_dict
            for k in sorted(job_vars.keys()):
                v = job_vars[k]
                if not k.startswith('__') and not k == 'CONFIG':
                    if k == '_runtime':
                        for _k, _v in v.items():
                            if isPrimitive(_v) and _v not in (None, '', [], (), {}):
                                row(_k, _v)
                    elif isPrimitive(v) and v not in (None, '', [], (), {}):
                        row(k,
                            f'<pre style="text-align:left">{pprint.pformat(v)}</pre>')
            summary = summarizeExecution(t, status=s)
            if summary:
                # row('Execution')
                for line in summary.split('\n'):
                    fields = line.split(None, 1)
                    if fields[0] == 'task':
                        continue
                    row(fields[0], '' if fields[1] is None else fields[1])
                # this is a placeholder for the frontend to draw figure
                row(td=f'<div id="res_{t}"></div>')
            #
            files = glob.glob(os.path.join(
                os.path.expanduser('~'), '.sos', 'tasks', t + '.*'))
            for f in sorted([x for x in files if os.path.splitext(x)[-1] not in ('.def', '.res', '.task', '.pulse', '.status')]):
                numLines = linecount_of_file(f)
                row(os.path.splitext(f)[-1], '(empty)' if numLines ==
                    0 else f'{numLines} lines{"" if numLines < 200 else " (showing last 200)"}')
                try:
                    row(
                        td=f'<small><pre style="text-align:left">{tail_of_file(f, 200, ansi2html=True)}</pre></small>')
                except Exception:
                    row(td='<small><pre style="text-align:left">ignored.</pre><small>')
            print('</table>')
            #
            # supplement run time information
            pulse_file = os.path.join(os.path.expanduser(
                '~'), '.sos', 'tasks', t + '.pulse')
            if not os.path.isfile(pulse_file):
                return
            else:
                # A sample of 400 point should be enough to show the change of resources
                lines = sample_of_file(pulse_file, 400).splitlines()
                if len(lines) <= 2:
                    return
            # read the pulse file and plot it
            # time   proc_cpu        proc_mem        children        children_cpu    children_mem
            try:
                etime = []
                cpu = []
                mem = []
                for line in lines:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.split()
                    etime.append(float(fields[0]))
                    cpu.append(float(fields[1]) + float(fields[4]))
                    mem.append(float(fields[2]) / 1e6 + float(fields[5]) / 1e6)
                if not etime:
                    return
            except Exception:
                return
            #
            print('''
<script>
    function loadFiles(files, fn) {
        if (!files.length) {
            files = [];
        }
        var head = document.head || document.getElementsByTagName('head')[0];

        function loadFile(index) {
            if (files.length > index) {
                if (files[index].endsWith('.css')) {
                    var fileref = document.createElement('link');
                    fileref.setAttribute("rel", "stylesheet");
                    fileref.setAttribute("type", "text/css");
                    fileref.setAttribute("href", files[index]);
                } else {
                    var fileref = document.createElement('script');
                    fileref.setAttribute("type", "text/javascript");
                    fileref.setAttribute("src", files[index]);
                }
                console.log('Load ' + files[index]);
                head.appendChild(fileref);
                index = index + 1;
                // Used to call a callback function
                fileref.onload = function() {
                    loadFile(index);
                }
            } else if (fn) {
                fn();
            }
        }
        loadFile(0);
    }

function plotResourcePlot_''' + t + '''() {
    // get the item
    // parent element is a table cell, needs enlarge
    document.getElementById(
        "res_''' + t + '''").parentElement.setAttribute("height", "300px;");
    $("#res_''' + t + '''").css("height", "300px");
    $("#res_''' + t + '''").css("width", "100%");
    $("#res_''' + t + '''").css("min-height", "300px");

    var cpu = [''' + ','.join([f'[{x*1000},{y}]' for x, y in zip(etime, cpu)]) + '''];
    var mem = [''' + ','.join([f'[{x*1000},{y}]' for x, y in zip(etime, mem)]) + '''];

    $.plot('#res_''' + t + '''', [{
                data: cpu,
                label: "CPU (%)"
            },
            {
                data: mem,
                label: "mem (M)",
                yaxis: 2
            }
        ], {
            xaxes: [{
                mode: "time"
            }],
            yaxes: [{
                min: 0
            }, {
                position: "right",
                tickFormatter: function(v, axis) {
                    return v.toFixed(1) + 'M';
                }
            }],
            legend: {
                position: "nw"
            }
        });
}

var dt = 100;
// the frontend might be notified before the table is inserted as results.
function showResourceFigure_''' + t + '''() {
    if ( $("#res_''' + t + '''").length === 0) {
          dt = dt * 1.5; // slow-down checks for datatable as time goes on;
          setTimeout(showResourceFigure_''' + t + ''', dt);
          return;
    } else {
        $("#res_''' + t + '''").css('width', "100%").css('height', "300px");
        loadFiles(["http://www.flotcharts.org/flot/jquery.flot.js",
             "http://www.flotcharts.org/flot/jquery.flot.time.js"
            ], plotResourcePlot_''' + t + ''');
    }
}
showResourceFigure_''' + t + '''()
</script>
''')
    if has_changes:
        with fasteners.InterProcessLock(cache_file + '_'):
            with open(cache_file, 'wb') as cache:
                pickle.dump(status_cache, cache)
    else:
        env.logger.debug('No new status detected')
    # remove jobs that are older than 1 month
    if to_be_removed:
        purge_tasks(to_be_removed, verbosity=0)


def kill_tasks(tasks, tags=None):
    #
    import glob
    from multiprocessing.pool import ThreadPool as Pool
    if not tasks:
        tasks = glob.glob(os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [os.path.basename(x)[:-5] for x in tasks]
    else:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'),
                                             '.sos', 'tasks', f'{t}*.task'))
            matched = [os.path.basename(x)[:-5] for x in matched]
            if not matched:
                env.logger.warning(f'{t} does not match any existing task')
            else:
                all_tasks.extend(matched)
    if tags:
        all_tasks = [x for x in all_tasks if any(
            x in tags for x in taskTags(x).split(' '))]

    if not all_tasks:
        env.logger.warning('No task to kill')
        return
    all_tasks = sorted(list(set(all_tasks)))
    # at most 20 threads
    p = Pool(min(20, len(all_tasks)))
    killed = p.map(kill_task, all_tasks)
    for s, t in zip(killed, all_tasks):
        print(f'{t}\t{s}')


def kill_task(task):
    status = check_task(task)['status']
    if status == 'pending':
        # the task engine will remove it from lists if killed through a task engine
        # but it can also be killed through command line externally, so we will need
        # to mark the task as killed
        pulse_file = os.path.join(os.path.expanduser(
            '~'), '.sos', 'tasks', task + '.pulse')
        os.open(pulse_file, flags=os.O_CREAT | os.O_RDONLY)
        return 'aborted'
    # remove job file as well
    job_file = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task + '.sh')
    if os.path.isfile(job_file):
        try:
            os.remove(job_file)
        except Exception:
            pass
    if status != 'running':
        return status
    # job is running
    pulse_file = os.path.join(os.path.expanduser(
        '~'), '.sos', 'tasks', task + '.pulse')
    os.chmod(pulse_file, S_IREAD | S_IRGRP | S_IROTH)
    return 'killed'


def purge_tasks(tasks, purge_all=False, age=None, status=None, tags=None, verbosity=2):
    # verbose is ignored for now
    import glob
    if tasks:
        all_tasks = []
        for t in tasks:
            matched = glob.glob(os.path.join(os.path.expanduser('~'),
                                             '.sos', 'tasks', f'{t}*.task'))
            matched = [(os.path.basename(x)[:-5], os.path.getmtime(x))
                       for x in matched]
            all_tasks.extend(matched)
    else:
        tasks = glob.glob(os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', '*.task'))
        all_tasks = [(os.path.basename(x)[:-5], os.path.getmtime(x))
                     for x in tasks]
    #
    if age is not None:
        age = expand_time(age, default_unit='d')
        if age > 0:
            all_tasks = [x for x in all_tasks if time.time() - x[1] >= age]
        else:
            all_tasks = [x for x in all_tasks if time.time() - x[1] <= -age]

    if status:
        # at most 20 threads
        from multiprocessing.pool import ThreadPool as Pool
        p = Pool(min(20, len(all_tasks)))
        task_status = p.map(check_task, [x[0] for x in all_tasks])
        all_tasks = [x for x, s in zip(
            all_tasks, task_status) if s['status'] in status]

    if tags:
        all_tasks = [x for x in all_tasks if any(
            x in tags for x in taskTags(x[0]).split(' '))]
    #
    # remoe all task files
    all_tasks = set([x[0] for x in all_tasks])
    if all_tasks:
        #
        # find all related files, including those in nested directories
        from collections import defaultdict
        to_be_removed = defaultdict(list)
        for dirname, _, filelist in os.walk(os.path.join(os.path.expanduser('~'), '.sos', 'tasks')):
            for f in filelist:
                ID = os.path.basename(f).split('.', 1)[0]
                if ID in all_tasks:
                    to_be_removed[ID].append(os.path.join(dirname, f))
        #
        cache_file: str = os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', 'status_cache.pickle')

        if os.path.isfile(cache_file):
            with fasteners.InterProcessLock(cache_file + '_'):
                with open(cache_file, 'rb') as cache:
                    status_cache = pickle.load(cache)
        else:
            status_cache = {}
        for task in all_tasks:
            removed = True
            for f in to_be_removed[task]:
                try:
                    if verbosity > 3:
                        env.logger.trace(f'Remove {f}')
                    os.remove(f)
                except Exception as e:
                    removed = False
                    if verbosity > 0:
                        env.logger.warning(f'Failed to purge task {task[0]}')
            status_cache.pop(task, None)
            if removed and verbosity > 1:
                env.logger.info(f'Task ``{task}`` removed.')
        with fasteners.InterProcessLock(cache_file + '_'):
            with open(cache_file, 'wb') as cache:
                pickle.dump(status_cache, cache)
    elif verbosity > 1:
        env.logger.info('No matching tasks')
    if purge_all:
        matched = glob.glob(os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', '*'))
        count = 0
        for f in matched:
            if os.path.isdir(f):
                import shutil
                try:
                    shutil.rmtree(f)
                    count += 1
                except Exception as e:
                    if verbosity > 0:
                        env.logger.warning(f'Failed to remove {f}')
            else:
                try:
                    os.remove(f)
                    count += 1
                except Exception as e:
                    if verbosity > 0:
                        env.logger.warning(f'Failed to remove {e}')
        if count > 0 and verbosity > 1:
            env.logger.info(
                f'{count} other files and directories are removed.')
    return ''


class TaskEngine(threading.Thread):
    def __init__(self, agent):
        threading.Thread.__init__(self)
        self.daemon = True
        #
        # agent is the agent that provides function
        #
        #    run_command
        #
        # to submit command, which can be a direct process call, or a call
        # on the remote server.
        #
        self.agent = agent
        self.config = agent.config
        self.alias = self.config['alias']
        self.engine_ready = threading.Event()

        self.running_tasks = []
        self.pending_tasks = []
        self.submitting_tasks = {}
        self.canceled_tasks = []
        self.resuming_tasks = set()

        self.task_status = OrderedDict()
        self.task_date = {}
        self.last_checked = None
        if 'status_check_interval' not in self.config:
            self.status_check_interval = 10
        else:
            self.status_check_interval = self.config['status_check_interval']
        #
        if env.config['max_running_jobs'] is not None:
            # override from command line
            self.max_running_jobs = env.config['max_running_jobs']
        elif 'max_running_jobs' in self.config:
            # queue setting
            self.max_running_jobs = self.config['max_running_jobs']
        else:
            # default
            self.max_running_jobs = max(os.cpu_count() // 2, 1)
        #
        # multiple thread job submission does not work because the threads share the
        # same namespace with variables such as sos_dict. Variables changed by
        # one thread can be changed again by another thread, which makes it unsafe
        # to submit jobs this way. Multi-processing is possible but that can be done
        # later.
        #
        self._thread_workers = concurrent.futures.ThreadPoolExecutor(
            max_workers=1)
        self._status_checker = None
        #
        if env.config['wait_for_task'] is not None:
            self.wait_for_task = env.config['wait_for_task']
        elif 'wait_for_task' in self.config:
            self.wait_for_task = self.config['wait_for_task']
        else:
            # default
            self.wait_for_task = True
        #
        # how to stack jobs ... currently backgroun process queue
        # allows stacking of up to 1000 tasks, but PBS queue does not
        # allow stacking.
        self.batch_size = 1

    def notify(self, msg):
        # GUI ...
        if hasattr(env, '__task_notifier__'):
            if not isinstance(msg, str):
                env.__task_notifier__(msg)
        elif isinstance(msg, str):
            env.logger.info(msg)
        # text mode does not provide detailed message change information

    def monitor_tasks(self, tasks=None, status=None, age=None):
        '''Start monitoring specified or all tasks'''
        self.engine_ready.wait()

        if not tasks:
            tasks = self.task_status.keys()
        else:
            tasks = [x for x in tasks if x in self.task_status]

        # we only monitor running tasks
        with threading.Lock():
            for task in tasks:
                if self.task_status[task] in ('submitted', 'running') and task not in self.running_tasks:
                    # these tasks will be actively monitored
                    self.running_tasks.append(task)
        #
        if age is not None:
            age = expand_time(age, default_unit='d')
        return sorted([(x, self.task_status[x], self.task_date.get(x, time.time())) for x in tasks
                       if (status is None or self.task_status[x] in status) and (age is None or
                                                                                 ((age > 0 and time.time() - self.task_date.get(x, time.time()) > age)
                                                                                  or (age < 0 and time.time() - self.task_date.get(x, time.time()) < -age)))],
                      key=lambda x: -x[2])

    def get_tasks(self):
        with threading.Lock():
            pending = copy.deepcopy(
                self.pending_tasks + list(self.submitting_tasks.keys()))
            running = copy.deepcopy(self.running_tasks)
        return pending, running

    def run(self):
        # get all system tasks that might have been running ...
        # this will be run only once when the task engine starts
        status_output = self.query_tasks([], verbosity=2, start_time=True)
        with threading.Lock():
            for line in status_output.split('\n'):
                if not line.strip():
                    continue
                try:
                    tid, _, ttm, tst = line.split('\t')
                    self.task_status[tid] = tst
                    self.task_date[tid] = float(ttm)
                except Exception as e:
                    env.logger.warning(
                        f'Unrecognized response "{line}" ({e.__class__.__name__}): {e}')
        self._last_status_check = time.time()
        self.engine_ready.set()
        while True:
            # if no new task, does not do anything.
            if self.running_tasks and time.time() - self._last_status_check > self.status_check_interval:
                if self._status_checker is None:
                    self._status_checker = self._thread_workers.submit(
                        self.query_tasks, self.running_tasks, 1)
                    continue
                elif self._status_checker.running():
                    time.sleep(0.01)
                    continue
                else:
                    status_output = self._status_checker.result()
                    self._status_checker = None
                #
                for line in status_output.splitlines():
                    if not line.strip():
                        continue
                    try:
                        tid, tst = line.split('\t')
                        if tid not in self.running_tasks:
                            env.logger.trace(
                                f'Task {tid} removed since status check.')
                            continue
                        self.update_task_status(tid, tst)
                    except Exception as e:
                        env.logger.warning(
                            f'Unrecognized response "{line}" ({e.__class__.__name__}): {e}')
                self.summarize_status()
                self._last_status_check = time.time()
            else:
                # if nothing to do, sleep to avoid empty loop. This will reduce CPU usage from 100% to 0.3%
                # sleeping longer will slightly reduce CPU usage but decrease responsiveness by a bit
                time.sleep(0.01)

            if self.submitting_tasks:
                with threading.Lock():
                    submitted = []
                    for k in self.submitting_tasks:
                        if not self.submitting_tasks[k].running():
                            submitted.append(k)
                            if self.submitting_tasks[k].result():
                                for tid in k:
                                    if tid in self.canceled_tasks:
                                        # task is canceled while being prepared
                                        self.notify(
                                            ['change-status', self.agent.alias, tid, 'aborted'])
                                    else:
                                        self.running_tasks.append(tid)
                                        self.notify(
                                            ['change-status', self.agent.alias, tid, 'submitted'])
                            else:
                                for tid in k:
                                    self.notify(
                                        ['change-status', self.agent.alias, tid, 'failed'])
                                    self.task_status[tid] = 'failed'
                        # else:
                        #    env.logger.trace('{} is still being submitted.'.format(k))
                    for k in submitted:
                        self.submitting_tasks.pop(k)

            if self.pending_tasks:
                # check status
                num_active_tasks = len(
                    self.submitting_tasks) + len(self.running_tasks)
                if num_active_tasks >= self.max_running_jobs:
                    continue

                # assign tasks to self.max_running_jobs workers
                slots = [[] for i in range(self.max_running_jobs)]
                sample_slots = list(range(self.max_running_jobs))
                random.shuffle(sample_slots)
                for i, tid in enumerate(self.pending_tasks[:self.batch_size * self.max_running_jobs]):
                    if self.task_status[tid] == 'running':
                        self.notify(f'{tid} ``runnng``')
                    elif tid in self.canceled_tasks:
                        # the job is canceled while being prepared to run
                        self.notify(f'{tid} ``canceled``')
                    else:
                        # randomly spread to tasks, but at most one.
                        slots[sample_slots[i %
                                           self.max_running_jobs]].append(tid)
                for slot in slots:
                    if not slot:
                        continue
                    for tid in slot:
                        env.logger.trace(
                            f'Start submitting {tid} (status: {self.task_status.get(tid, "unknown")})')
                    self.submitting_tasks[tuple(slot)] = self._thread_workers.submit(
                        self.execute_tasks, slot)
                #
                with threading.Lock():
                    for slot in slots:
                        for tid in slot:
                            self.pending_tasks.remove(tid)

    def submit_task(self, task_id):
        # we wait for the engine to start
        self.engine_ready.wait()

        # submit tasks simply add task_id to pending task list
        with threading.Lock():
            # if already in
            # if task_id in self.running_tasks or task_id in self.pending_tasks:
            #    self.notify('{} ``{}``'.format(task_id, self.task_status[task_id]))
            #    self.notify(['new-status', task_id, self.task_status[task_id]])
            #    return self.task_status[task_id]
            #
            if task_id in self.task_status and self.task_status[task_id]:
                if self.task_status[task_id] == 'running':
                    self.running_tasks.append(task_id)
                    self.notify(f'{task_id} ``already runnng``')
                    self.notify(['new-status', self.agent.alias, task_id, 'running',
                                 self.task_date.get(task_id, time.time())])
                    return 'running'
                # there is a case when the job is already completed (complete-old), but
                # because we do not know if the user asks to rerun (-s force), we have to
                # resubmit the job. In the case of not-rerun, the task would be marked
                # completed very soon.
                elif self.task_status[task_id] == 'completed':
                    if env.config['sig_mode'] != 'force':
                        self.notify(f'{task_id} ``already completed``')
                        return 'completed'
                    elif task_id in env.config.get('resumed_tasks', []):
                        # force re-execution, but it is possible that this task has been
                        # executed but quit in no-wait mode (or canceled by user). More
                        # importantly, the Jupyter notebook would re-run complted workflow
                        # even if it has "-s force" signature.
                        # self.notify(['new-status', task_id, 'completed'])
                        self.notify(f'{task_id} ``resume with completed``')
                        return 'completed'
                    else:
                        self.notify(f'{task_id} ``re-execute completed``')
                else:
                    self.notify(
                        f'{task_id} ``restart`` from status ``{self.task_status[task_id]}``')

            # self.notify('{} ``queued``'.format(task_id))
            self.pending_tasks.append(task_id)
            if task_id in self.canceled_tasks:
                self.canceled_tasks.remove(task_id)
            self.task_status[task_id] = 'pending'
            self.notify(['new-status', self.agent.alias, task_id, 'pending',
                         self.task_date.get(task_id, time.time())])
            return 'pending'

    def summarize_status(self):
        from collections import Counter
        statuses = Counter(self.task_status.values())
        env.logger.debug(
            ' '.join(f'{x}: {y}' for x, y in statuses.items()))

    def check_task_status(self, task_id, unknown='pending'):
        # we wait for the engine to start
        self.engine_ready.wait()
        try:
            with threading.Lock():
                return self.task_status[task_id]
        except Exception:
            # job not yet submitted
            return unknown

    def update_task_status(self, task_id, status):
        #
        env.logger.trace(f'STATUS {task_id}\t{status}\n')
        #
        with threading.Lock():
            if task_id in self.canceled_tasks and status != 'aborted':
                env.logger.debug(
                    f'Task {task_id} is still not killed (status {status})')
                status = 'aborted'
            if status != 'missing':
                if task_id in self.task_status and self.task_status[task_id] == status:
                    self.notify(
                        ['pulse-status', self.agent.alias, task_id, status])
                else:
                    self.notify(
                        ['change-status', self.agent.alias, task_id, status])
            self.task_status[task_id] = status
            if status == 'pening' and task_id not in self.pending_tasks:
                self.pending_tasks.append(task_id)
            if status == 'running' and task_id not in self.running_tasks:
                self.running_tasks.append(task_id)
            # terminal states, remove tasks from task list
            if status in ('completed', 'failed', 'aborted', 'signature-mismatch') and task_id in self.running_tasks:
                self.running_tasks.remove(task_id)

    def remove_tasks(self, tasks):
        with threading.Lock():
            for task in tasks:
                self.notify(['remove-task', self.agent.alias, task])
                # if task in self.task_status:
                #    self.task_status.pop(task)
                # if task in self.running_tasks:
                #    self.running_tasks.remove(task)

    def query_tasks(self, tasks=None, verbosity=1, html=False, start_time=False, age=None, tags=None, status=None):
        try:
            return self.agent.check_output("sos status {} -v {} {} {} {} {} {}".format(
                '' if tasks is None else ' '.join(tasks), verbosity,
                '--html' if html else '',
                '--start-time' if start_time else '',
                f'--age {age}' if age else '',
                f'--tags {" ".join(tags)}' if tags else '',
                f'--status {" ".join(status)}' if status else '',
            ))
        except subprocess.CalledProcessError as e:
            if verbosity >= 3:
                env.logger.warning(
                    f'Failed to query status of tasks on {self.alias}')
            return ''

    def kill_tasks(self, tasks, tags=None, all_tasks=False):
        # we wait for the engine to start
        self.engine_ready.wait()
        if all_tasks:
            tasks = self.pending_tasks + \
                list(self.submitting_tasks.keys()) + self.running_tasks

        for task in tasks:
            with threading.Lock():
                self.task_status[task] = 'aborted'
        for task in tasks:
            with threading.Lock():
                if task in self.pending_tasks:
                    self.pending_tasks.remove(task)
                    env.logger.debug(f'Cancel pending task {task}')
                elif task in self.submitting_tasks:
                    # this is more troublesome because the task is being
                    # submitted at a new thread.
                    env.logger.debug(f'Cancel submission of task {task}')
                elif task in self.running_tasks:
                    # the job will be killed from command line with
                    # status updated here
                    env.logger.debug(f'Killing running task {task}')
                else:
                    # it is not in the system, so we need to know what the
                    # status of the task before we do anything...
                    pass

        self.canceled_tasks.extend(tasks)
        #
        # verbosity cannot be send to underlying command because task engines
        # rely on the output of certain verbosity (-v1) to post kill the jobs
        cmd = "sos kill {} {} {}".format('' if all_tasks else ' '.join(tasks),
                                         f'--tags {" ".join(tags)}' if tags else '',
                                         '-a' if all_tasks else '')

        try:
            ret = self.agent.check_output(cmd)
            env.logger.debug(f'"{cmd}" executed with response "{ret}"')
        except subprocess.CalledProcessError:
            env.logger.error(
                f'Failed to kill all tasks' if all_tasks else f'Failed to kill tasks {" ".join(tasks)}')
            return ''
        return ret

    def resume_task(self, task):
        # we wait for the engine to start
        self.engine_ready.wait()
        with threading.Lock():
            # it is possible that a task is aborted from an opened notebook with aborted status
            if task not in self.task_status or \
                    self.task_status[task] not in ('completed', 'failed', 'signature-mismatch', 'aborted'):
                env.logger.warning(
                    f'Resume task called for non-canceled or non-completed/failed task {task}')
                return
            # the function might have been used multiple times (frontend multiple clicks)
            if task in self.canceled_tasks:
                self.canceled_tasks.remove(task)
            if task not in self.pending_tasks:
                self.pending_tasks.append(task)
            # tells the engine that preparation of task can fail
            self.resuming_tasks.add(task)
            self.task_status[task] = 'pending'

    def execute_tasks(self, task_ids):
        # we wait for the engine to start
        self.engine_ready.wait()
        # this is base class, the derived class will actually submit the tasks

        # if the task is being resumed, perhaps from another local host,
        # the preparation process can fail (e.g. no def file), but this
        # does not really matter. #587
        for task_id in task_ids:
            if task_id in self.resuming_tasks:
                self.resuming_tasks.remove(task_id)
                try:
                    self.agent.prepare_task(task_id)
                except Exception:
                    pass
            else:
                if not self.agent.prepare_task(task_id):
                    return False
        return True

    def purge_tasks(self, tasks, purge_all=False, age=None, status=None, tags=None, verbosity=2):
        try:
            return self.agent.check_output("sos purge {} {} {} {} {} -v {}".format(
                ' '.join(tasks), '--all' if purge_all else '',
                f'--age {age}' if age is not None else '',
                f'--status {" ".join(status)}' if status is not None else '',
                f'--tags {" ".join(tags)}' if tags is not None else '',
                verbosity))
        except subprocess.CalledProcessError:
            env.logger.error(f'Failed to purge tasks {tasks}')
            return ''


class BackgroundProcess_TaskEngine(TaskEngine):
    def __init__(self, agent):
        super(BackgroundProcess_TaskEngine, self).__init__(agent)
        if 'job_template' in self.config:
            self.job_template = self.config['job_template'].replace(
                '\r\n', '\n')
        else:
            self.job_template = None
        #
        if 'batch_size' in self.config:
            self.batch_size = self.config['batch_size']
        else:
            # default allow stacking of up to 1000 jobs
            self.batch_size = 1000

    def execute_tasks(self, task_ids):
        if not super(BackgroundProcess_TaskEngine, self).execute_tasks(task_ids):
            env.logger.trace(f'Failed to prepare task {task_ids}')
            return False
        if self.job_template:
            if not self._submit_task_with_template(task_ids):
                return False
        else:
            if not self._submit_task(task_ids):
                return False
        return True

    def _submit_task(self, task_ids):
        # if no template, use a default command
        cmd = f"sos execute {' '.join(task_ids)} -v {env.verbosity} -s {env.config['sig_mode']} {'--dryrun' if env.config['run_mode'] == 'dryrun' else ''}"
        env.logger.trace(f'Execute "{cmd}" (waiting={self.wait_for_task})')
        self.agent.run_command(cmd, wait_for_task=self.wait_for_task)
        return True

    def _submit_task_with_template(self, task_ids):
        '''Submit tasks by interpolating a shell script defined in job_template'''
        runtime = self.config
        runtime.update({
            'cur_dir': os.getcwd(),
            'verbosity': env.verbosity,
            'sig_mode': env.config.get('sig_mode', 'default'),
            'run_mode': env.config.get('run_mode', 'run'),
            'home_dir': os.path.expanduser('~')})
        if '_runtime' in env.sos_dict:
            runtime.update({x: env.sos_dict['_runtime'][x] for x in (
                'nodes', 'cores', 'mem', 'walltime') if x in env.sos_dict['_runtime']})
        if 'nodes' not in runtime:
            runtime['nodes'] = 1
        if 'cores' not in runtime:
            runtime['cores'] = 1

        # let us first prepare a task file
        job_text = ''
        for task_id in task_ids:
            runtime['task'] = task_id
            try:
                job_text += cfg_interpolate(self.job_template, runtime)
                job_text += '\n'
            except Exception as e:
                raise ValueError(
                    f'Failed to generate job file for task {task_id}: {e}')

        filename = task_ids[0] + ('.sh' if len(task_ids)
                                  == 1 else f'-{task_ids[-1]}.sh')
        # now we need to write a job file
        job_file = os.path.join(os.path.expanduser(
            '~'), '.sos', 'tasks', self.alias, filename)
        # do not translate newline under windows because the script will be executed
        # under linux/mac
        with open(job_file, 'w', newline='') as job:
            job.write(job_text)

        # then copy the job file to remote host if necessary
        self.agent.send_task_file(job_file)

        try:
            cmd = f'bash ~/.sos/tasks/{filename}'
            self.agent.run_command(cmd, wait_for_task=self.wait_for_task)
        except Exception as e:
            raise RuntimeError(f'Failed to submit task {task_ids}: {e}')
        return True
