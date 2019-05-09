#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import glob
import multiprocessing as mp
import os
import shutil
import socket
import stat
import subprocess
import sys
from collections import Sequence

import pexpect
import pkg_resources

from .eval import Undetermined, cfg_interpolate
from .syntax import SOS_LOGLINE
from .targets import path, sos_targets
from .task_engines import BackgroundProcess_TaskEngine
from .tasks import TaskFile
from .utils import (env, expand_size, expand_time, format_HHMMSS, short_repr)

#
# A 'queue' is defined by queue configurations in SoS configuration files.
# It encapsulate properties of a queue and tells sos how to interact with
# the queue. A queue can be a local host, a remote host without queue, or
# a remote host with a task queue, or a RQ or Celery server. There are
# two main categories of properties of a host.
#
# 1. host properties, namely how to copy files and how to execute comamnds
#   on the host. Note that even for queues that communicate with sos
#   through network, the workers might be on a different host with different
#   file systems.
#
#   Keys for host configuration include:
#   * path_map: path map between local and remote hosts
#   * shared: paths that are shared between local and remote hosts
#   * send_cmd (optional): alternative command to send files
#   * receive_cmd (optional): alternative command to receive files
#   * execute_cmd (optional): alternative command to execute commands
#
# 2. task properties, namely how to manage running jobs. These include
#   direct execution, PBS and various cluster systems, and various task
#   queues.
#
#   Keys for task configuration depend largely one queue type.
#
#   * task_engine: type of task engine
#   * max_jobs: maximum number of concurrent jobs on the host.
#
#
# Implementation wise, a queue instance is created for each queue.
#

from sos.targets import file_target
from typing import Any, Dict, List, Optional, Union


class DaemonizedProcess(mp.Process):

    def __init__(self, cmd, *args, **kwargs):
        super(DaemonizedProcess, self).__init__(*args, **kwargs)
        self.cmd = cmd

    def run(self):
        try:
            pid = os.fork()
            if pid > 0:
                # exit from second parent
                sys.exit(0)
        except OSError as err:
            env.logger.error(f'_Fork #1 failed: {err}\n')
            sys.exit(1)

        os.setsid()
        os.umask(0)
        # do second fork
        try:
            pid = os.fork()
            if pid > 0:
                # exit from second parent
                sys.exit(0)
        except OSError as err:
            env.logger.error(f'_Fork #2 failed: {err}\n')
            sys.exit(1)
        # the following is also need to properly daemonize the process
        # redirect standard file descriptors
        sys.stdout.flush()
        sys.stderr.flush()
        try:
            si = open(os.devnull, 'r')
            so = open(os.devnull, 'w')
            se = open(os.devnull, 'w')
            os.dup2(si.fileno(), sys.stdin.fileno())
            os.dup2(so.fileno(), sys.stdout.fileno())
            os.dup2(se.fileno(), sys.stderr.fileno())
        except Exception:
            # #493
            pass

        # fork a new process
        subprocess.Popen(self.cmd, shell=True, close_fds=True)
        return


def _show_err_and_out(task_id, res) -> None:
    if 'stdout' in res:
        sys.stderr.write(f'\n{task_id}.out:\n')
        ends_with_newline = False
        for line in res['stdout'].splitlines():
            if not SOS_LOGLINE.match(line):
                sys.stderr.write(line)
                ends_with_newline = line.endswith('\n')
        if not ends_with_newline:
            sys.stderr.write('\n')
    if 'stderr' in res:
        sys.stderr.write(f'\n{task_id}.err:\n')
        ends_with_newline = False
        for line in res['stderr'].splitlines():
            if not SOS_LOGLINE.match(line):
                sys.stderr.write(line)
                ends_with_newline = line.endswith('\n')
        if not ends_with_newline:
            sys.stderr.write('\n')


class LocalHost:
    '''For local host, no path map, send and receive ...'''

    def __init__(self, config: Dict[str, Union[str, int, List[str]]]) -> None:
        # even if the config has an alias, we use localhost to make it clear that the host is localhost
        self.alias = config.get('alias', 'localhost')
        self.address = 'localhost'
        # we checkk local jobs more aggressively
        self.config = {'alias': self.alias, 'status_check_interval': 2}
        self.config.update(config)

    def send_to_host(self, items):
        return {x: x for x in items}

    def receive_from_host(self, items):
        return {x: x for x in items}

    def prepare_task(self, task_id):
        task_file = os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
        # add server restriction on task file
        if not os.path.isfile(task_file):
            raise ValueError(f'Missing task definition {task_file}')

        tf = TaskFile(task_id)
        params, old_runtime = tf.get_params_and_runtime()
        # clear possible previous result
        task_vars = params.sos_dict
        runtime = {
            '_runtime': {
                x: task_vars['_runtime'][x]
                for x in ('verbosity', 'sig_mode', 'run_mode', 'walltime',
                          'cores', 'mem')
                if x in task_vars['_runtime']
            }
        }
        runtime['_runtime']['workdir'] = task_vars['_runtime'][
            'workdir'] if 'workdir' in task_vars['_runtime'] else os.getcwd()

        if 'max_mem' in self.config or 'max_cores' in self.config or 'max_walltime' in self.config:
            for key in ('max_mem', 'max_cores', 'max_walltime'):
                if key in self.config:
                    runtime['_runtime'][key] = format_HHMMSS(
                        self.config[key]
                    ) if key == 'max_walltime' else self.config[key]

            if self.config.get('max_mem', None) is not None and task_vars['_runtime'].get('mem', None) is not None \
                    and self.config['max_mem'] < task_vars['_runtime']['mem']:
                env.logger.error(
                    f'Task {task_id} requested more mem ({task_vars["_runtime"]["mem"]}) than allowed max_mem ({self.config["max_mem"]})'
                )
                return False
            if self.config.get('max_cores', None) is not None and task_vars['_runtime'].get('cores', None) is not None \
                    and self.config['max_cores'] < task_vars['_runtime']['cores']:
                env.logger.error(
                    f'Task {task_id} requested more cores ({task_vars["_runtime"]["cores"]}) than allowed max_cores ({self.config["max_cores"]})'
                )
                return False
            if self.config.get('max_walltime', None) is not None and task_vars['_runtime'].get('walltime', None) is not None \
                    and expand_time(self.config['max_walltime']) < expand_time(task_vars['_runtime']['walltime']):
                env.logger.error(
                    f'Task {task_id} requested more walltime ({task_vars["_runtime"]["walltime"]}) than allowed max_walltime ({self.config["max_walltime"]})'
                )
                return False

        # if the task has been running remotely, we need to reset runtime for local execution
        if len(runtime) > 1 or runtime['_runtime'] or runtime != old_runtime:
            tf.runtime = runtime
        tf.status = 'pending'
        #
        if 'to_host' in task_vars['_runtime'] and isinstance(
                task_vars['_runtime']['to_host'], dict):
            for l, r in task_vars['_runtime']['to_host'].items():
                if l != r:
                    shutil.copy(l, r)
        self.send_task_file(task_file)
        return True

    def send_task_file(self, task_file):
        # on the same file system, no action is needed.
        dest_task_file = os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks',
            os.path.basename(task_file))
        if task_file != dest_task_file:
            shutil.copyfile(task_file, dest_task_file)

    def check_output(self, cmd, under_workdir=False):
        # get the output of command
        if isinstance(cmd, list):
            cmd = subprocess.list2cmdline(cmd)
        try:
            cmd = cfg_interpolate(cmd)
            return subprocess.check_output(
                cmd, shell=isinstance(cmd, str)).decode()
        except Exception as e:
            env.logger.warning(f'Check output of {cmd} failed: {e}')
            raise

    def check_call(self, cmd, under_workdir=False, **kwargs):
        # get the output of command
        try:
            return subprocess.check_call(
                cmd, shell=isinstance(cmd, str), **kwargs)
        except Exception as e:
            env.logger.warning(f'Check output of {cmd} failed: {e}')
            raise

    def run_command(self, cmd, wait_for_task, realtime=False, **kwargs):
        # run command but does not wait for result.
        if realtime:
            from .utils import pexpect_run
            return pexpect_run(cmd)
        elif wait_for_task or sys.platform == 'win32':
            return subprocess.Popen(cmd, shell=True, **kwargs)
        else:
            p = DaemonizedProcess(cmd, **kwargs)
            p.start()
            p.join()

    def receive_result(self, task_id: str) -> Dict[str, Any]:
        tf = TaskFile(task_id)
        params = tf.params
        job_dict = params.sos_dict

        if 'from_host' in job_dict['_runtime'] and isinstance(
                job_dict['_runtime']['from_host'], dict):
            for l, r in job_dict['_runtime']['from_host'].items():
                if l != r:
                    shutil.copy(r, l)

        res = tf.result
        if not res or 'ret_code' not in res:
            return {
                'ret_code':
                    1,
                'task': task_id,
                'exception':
                    ValueError(f'No result is received for task {task_id}')
            }

        try:
            if res['ret_code'] != 0 or env.verbosity >= 3:
                _show_err_and_out(task_id, res)
        except Exception as e:
            # if ret_code does not exist...
            return {'ret_code': 1, 'output': {}, 'exception': e}
        return res


class RemoteHost:
    '''A remote host class that manages how to communicate with remote host'''

    def __init__(self, config: Dict[str, Union[str, int, List[str]]]) -> None:
        self.config = config
        self.cm_opts = self._get_control_master_options()
        self.alias = self.config['alias']
        #
        self.address = self.config['address']
        self.port = self.config.get('port', 22)
        self.shared_dirs = self._get_shared_dirs()
        self.path_map = self._get_path_map()

    def _get_shared_dirs(self) -> List[Any]:
        value = self.config.get('shared', [])
        if isinstance(value, str):
            return [value]
        if isinstance(value, Sequence):
            return value
        raise ValueError(
            'Option shared can only be a string or a list of strings')

    def _get_path_map(self) -> Dict[str, str]:
        res: Dict = {}
        # if user-specified path_map, it overrides CONFIG
        path_map = self.config.get('path_map', [])
        #
        if not path_map:
            return res
        if isinstance(path_map, str):
            path_map = [path_map]
        if isinstance(path_map, Sequence):
            for v in path_map:
                if ' -> ' not in v:
                    raise ValueError(
                        f'Path map should be separated as from -> to, {v} specified'
                    )
                elif v.count(' -> ') > 1:
                    raise ValueError(
                        f'Path map should be separated as from -> to, {v} specified'
                    )
                res[v.split(' -> ')[0]] = v.split(' -> ')[1]
        elif isinstance(path_map, dict):
            for k, v in path_map.items():
                res[k] = v
        else:
            raise ValueError(
                f'Unacceptable path_mapue for configuration path_map: {path_map}'
            )
        return res

    def _get_control_master_options(self):
        master_dir = os.path.join(
            os.path.expanduser('~'), '.ssh', 'controlmasters')
        if not os.path.isdir(master_dir):
            try:
                os.makedirs(master_dir, exist_ok=True)
            except Exception as e:
                env.logger.debug(
                    f'Failed to create ssh control master directory {master_dir}: {e}'
                )
                return ''
        return f'-o "ControlMaster=auto" -o "ControlPath={master_dir}/%r@%h:%p" -o "ControlPersist=10m"'

    def _get_send_cmd(self, rename=False):
        if rename:
            return 'ssh ' + self.cm_opts + ''' -q {host} -p {port} "mkdir -p {dest:dpq}" && ''' + \
                '''rsync -a --no-g -e 'ssh ''' + self.cm_opts + ''' -p {port}' {source:aep} "{host}:{dest:dep}" && ''' + \
                '''ssh ''' + self.cm_opts + ''' -q {host} -p {port} "mv {dest:dep}/{source:b} {dest:ep}" '''
        return 'ssh ' + self.cm_opts + ''' -q {host} -p {port} "mkdir -p {dest:dpq}" && rsync -a --no-g -e 'ssh -p {port}' {source:aep} "{host}:{dest:dep}"'''

    def _get_receive_cmd(self, rename=False):
        if rename:
            return '''rsync -a --no-g -e 'ssh ''' + self.cm_opts + ''' -p {port}' {host}:{source:e} "{dest:adep}" && ''' + \
                '''mv "{dest:adep}/{source:b}" "{dest:aep}"'''
        return '''rsync -a --no-g -e 'ssh ''' + self.cm_opts + ''' -p {port}' {host}:{source:e} "{dest:adep}"'''

    def _get_execute_cmd(self, under_workdir=True) -> str:
        return self.config.get(
            'execute_cmd', 'ssh ' + self.cm_opts +
            """ -q {host} -p {port} "bash --login -c '""" +
            ('[ -d {workdir} ] || mkdir -p {workdir}; cd {workdir} && '
             if under_workdir else '') + ''' {cmd}'" ''')

    def _get_query_cmd(self):
        return self.config.get(
            'query_cmd', '''ssh ''' + self.cm_opts +
            ''' -q {host} -p {port} "bash --login -c 'sos status {task} -v 0'" '''
        )

    def is_shared(self, path):
        fullpath = os.path.abspath(os.path.expanduser(path))
        for sdir in self.shared_dirs:
            if fullpath.startswith(sdir):
                # issue 710, if a directory is both under path_map and shared, then it is not considered to be shared.
                if not any(
                        fullpath.startswith(mdir)
                        for mdir in self.path_map.keys()):
                    return True
        return False

    def _map_path(self, source):
        result = {}
        cwd = os.getcwd()
        if isinstance(source, (str, path)):
            dest = os.path.abspath(os.path.expanduser(source))
            # we use samefile to avoid problems with case-insensitive file system #522
            # we also use the "cwd" name to avoid wrong case for cwd. For example,
            # if the cwd = '/Users/user/Project'
            # then, dest = '/USERS/USER/PROJECT/a.txt'
            # would be converted to '/Users/user/Project/a.txt' before path mapping
            if os.path.exists(dest[:len(cwd)]) and os.path.samefile(
                    dest[:len(cwd)], cwd):
                dest = cwd + dest[len(cwd):]
            matched = [
                k for k in self.path_map.keys()
                if os.path.exists(dest[:len(k)]) and
                os.path.samefile(dest[:len(k)], k)
            ]
            if matched:
                # pick the longest key that matches
                k = max(matched, key=len)
                dest = self.path_map[k] + dest[len(k):]
            else:
                env.logger.warning(
                    f'Path {source} is not under any specified paths of localhost and is mapped to {dest} on remote host.'
                )
            result[source] = dest.replace('\\', '/')
        elif isinstance(source, (Sequence, set, sos_targets)):
            for src in source:
                result.update(self._map_path(src))
        else:
            env.logger.debug(f'Ignore unmappable source {source}')
            return {source: source}
        return result

    #
    # Interface functions
    #
    def _map_var(self, source):
        cwd = os.getcwd()
        if isinstance(source, path):
            source = str(source)
        if isinstance(source, str):
            dest = os.path.abspath(os.path.expanduser(source))
            # we use samefile to avoid problems with case-insensitive file system #522
            # we also use the "cwd" name to avoid wrong case for cwd. For example,
            # if the cwd = '/Users/user/Project'
            # then, dest = '/USERS/USER/PROJECT/a.txt'
            # would be converted to '/Users/user/Project/a.txt' before path mapping
            if os.path.exists(dest[:len(cwd)]) and os.path.samefile(
                    dest[:len(cwd)], cwd):
                dest = cwd + dest[len(cwd):]
            matched = [
                k for k in self.path_map.keys()
                if os.path.exists(dest[:len(k)]) and
                os.path.samefile(dest[:len(k)], k)
            ]
            if matched:
                # pick the longest key that matches
                k = max(matched, key=len)
                dest = self.path_map[k] + dest[len(k):]
            else:
                env.logger.debug(
                    f'Path {source} is not under any specified paths of localhost and is mapped to {dest} on remote host.'
                )
            return dest.replace('\\', '/')
        elif isinstance(source, (Sequence, set, sos_targets)):
            ret = [self._map_var(x) for x in source]
            return [x for x in ret if x is not None]
        else:
            env.logger.debug(f'Ignore unmappable source {source}')
            return source

    def _reverse_map_var(self, dest):
        if isinstance(dest, path):
            dest = str(dest)
        if isinstance(dest, str):
            matched = [
                l for l, r in self.path_map.items()
                if dest.startswith(r) and (len(r) == len(dest) or r.endswith(
                    '/') or r.endswith('\\') or dest[len(r)] in ('/', '\\'))
            ]
            if matched:
                # pick the longest key that matches
                k = max(matched, key=len)
                dest = k + dest[len(self.path_map[k]):]
            else:
                env.logger.debug(
                    f'Path {dest} is not under any specified paths of localhost and is mapped to {dest} on local host.'
                )
            return dest.replace('\\', '/')
        elif isinstance(dest, (Sequence, set, sos_targets)):
            ret = [self._reverse_map_var(x) for x in dest]
            return [x for x in ret if x is not None]
        else:
            env.logger.debug(f'Ignore unmappable source {dest}')
            return dest

    def _remote_abs(self, path):
        # return an absolute path relative to remote host
        path = str(path)
        if os.path.isabs(path):
            return path
        return os.path.join(self._map_var(os.getcwd()), path)

    def send_to_host(self, items):
        # we only copy files and directories, not other types of targets
        if isinstance(items, str):
            items = [items]
        elif isinstance(items, path):
            items = [str(items)]
        elif isinstance(items, Sequence):
            ignored = [x for x in items if not isinstance(x, (str, path))]
            if ignored:
                env.logger.info(f'``Ignore`` {ignored}')
            items = [x for x in items if isinstance(x, (str, path))]
        elif isinstance(items, dict):
            items = items
        else:
            env.logger.warning(
                f'Unrecognized items to be sent to host: {items}')
            return {}

        if isinstance(items, Sequence):
            from .utils import find_symbolic_links
            new_items = []
            for item in items:
                links = find_symbolic_links(item)
                for link, realpath in links.items():
                    env.logger.info(
                        f'Adding {realpath} for symbolic link {link}')
                new_items.extend(links.values())
            items.extend(new_items)

            sending = self._map_path(items)
        else:
            sending = items
        sent = {}
        for source in sorted(sending.keys()):
            dest = self._remote_abs(sending[source])
            if self.is_shared(source):
                if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                    env.log_to_file(
                        'TASK', f'Skip sending {source} on shared file system')
            else:
                if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                    env.log_to_file(
                        'TASK', f'Sending ``{source}`` to {self.alias}:{dest}')
                cmd = cfg_interpolate(
                    self._get_send_cmd(
                        rename=os.path.basename(source) != os.path.basename(
                            dest)), {
                                'source': sos_targets(str(source).rstrip('/')),
                                'dest': sos_targets(dest),
                                'host': self.address,
                                'port': self.port
                            })
                if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                    env.log_to_file('TASK', cmd)
                ret = subprocess.call(
                    cmd,
                    shell=True,
                    stderr=subprocess.DEVNULL,
                    stdout=subprocess.DEVNULL)
                if (ret != 0):
                    raise RuntimeError(
                        f'Failed to copy {source} to {self.alias} using command "{cmd}". The remote host might be unavailable.'
                    )
            sent[source] = dest
        return sent

    def receive_from_host(self, items):
        if isinstance(items, dict):
            # specify as local:remote
            # needs remote:local
            receiving = {self._remote_abs(y): str(x) for x, y in items.items()}
        else:
            # y could be path
            receiving = {
                self._remote_abs(y): str(x)
                for x, y in self._map_path(items).items()
            }
        #
        received = {}
        for source in sorted(receiving.keys()):
            dest = receiving[source]
            dest_dir = os.path.dirname(dest)
            if dest_dir and not os.path.isdir(dest_dir):
                try:
                    os.makedirs(dest_dir)
                except Exception as e:
                    env.logger.error(
                        f'Failed to create destination directory {dest_dir}: {e}'
                    )
            if self.is_shared(dest) and os.path.basename(
                    source) == os.path.basename(dest):
                env.logger.debug(
                    f'Skip retrieving ``{dest}`` from shared file system')
                received[dest] = source
            else:
                cmd = cfg_interpolate(
                    self._get_receive_cmd(
                        rename=os.path.basename(source) != os.path.basename(
                            dest)), {
                                'source': sos_targets(str(source).rstrip('/')),
                                'dest': sos_targets(dest),
                                'host': self.address,
                                'port': self.port
                            })
                if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
                    env.log_to_file('TASK', cmd)
                try:
                    ret = subprocess.call(
                        cmd,
                        shell=True,
                        stderr=subprocess.DEVNULL,
                        stdout=subprocess.DEVNULL)
                    if (ret != 0):
                        raise RuntimeError(f'command return {ret}')
                    received[dest] = source
                except Exception as e:
                    raise RuntimeError(
                        f'Failed to copy {source} from {self.alias} using command "{cmd}": {e}'
                    )
        return received

    #
    # Interface
    #
    def prepare_task(self, task_id):
        try:
            self._prepare_task(task_id)
            return True
        except Exception as e:
            env.logger.error(e)
            return False

    def _prepare_task(self, task_id):
        task_file = os.path.join(
            os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
        if not os.path.isfile(task_file):
            raise ValueError(f'Missing task definition {task_file}')
        tf = TaskFile(task_id)
        params, old_runtime = tf.get_params_and_runtime()
        task_vars = params.sos_dict
        runtime = {
            '_runtime': {
                x: task_vars['_runtime'][x]
                for x in ('verbosity', 'sig_mode', 'run_mode', 'walltime',
                          'cores', 'mem')
                if x in task_vars['_runtime']
            },
            task_id: {}
        }

        if self.config.get('max_mem', None) is not None and task_vars['_runtime'].get('mem', None) is not None \
                and self.config['max_mem'] < task_vars['_runtime']['mem']:
            raise ValueError(
                f'Task {task_id} requested more mem ({task_vars["_runtime"]["mem"]}) than allowed max_mem ({self.config["max_mem"]})'
            )
        if self.config.get('max_cores', None) is not None and task_vars['_runtime'].get('cores', None) is not None \
                and self.config['max_cores'] < task_vars['_runtime']['cores']:
            raise ValueError(
                f"Task {task_id} requested more cores ({task_vars['_runtime']['cores']}) than allowed max_cores ({self.config['max_cores']})"
            )
        if self.config.get('max_walltime', None) is not None and task_vars['_runtime'].get('walltime', None) is not None \
                and expand_time(self.config['max_walltime']) < expand_time(task_vars['_runtime']['walltime']):
            raise ValueError(
                f'Task {task_id} requested more walltime ({task_vars["_runtime"]["walltime"]}) than allowed max_walltime ({self.config["max_walltime"]})'
            )

        if task_vars['_input'] and not isinstance(task_vars['_input'],
                                                  Undetermined):
            sent = self.send_to_host(task_vars['_input'])
            if sent:
                env.logger.info(
                    f'{task_id} ``sent`` {short_repr(sent.keys())} to {self.alias}'
                )
        if task_vars['_depends'] and not isinstance(task_vars['_depends'],
                                                    Undetermined):
            sent = self.send_to_host(task_vars['_depends'])
            if sent:
                env.logger.info(
                    f'{task_id} ``sent`` {short_repr(sent.keys())} to {self.alias}'
                )
        if 'to_host' in task_vars['_runtime']:
            sent = self.send_to_host(task_vars['_runtime']['to_host'])
            if sent:
                env.logger.info(
                    f'{task_id} ``sent`` {short_repr(sent.keys())} to {self.alias}'
                )

        # map variables
        if 'workdir' in task_vars['_runtime']:
            runtime['_runtime']['workdir'] = self._map_var(
                task_vars['_runtime']['workdir'])
        else:
            runtime['_runtime']['workdir'] = self._map_var(os.getcwd())

        mapped_vars = {
            '_input', '_output', '_depends', 'input', 'output', 'depends'
        }
        if 'mapp_vars' in task_vars['_runtime']:
            if isinstance(task_vars['_runtime']['mapped_vars_vars'], str):
                mapped_vars.add(task_vars['_runtime']['mapped_vars_vars'])
            elif isinstance(task_vars['_runtime']['mapped_vars_vars'],
                            (set, Sequence)):
                mapped_vars |= set(task_vars['_runtime']['mapped_vars_vars'])
            else:
                raise ValueError(
                    f'Unacceptable value for runtime option mapped_vars_vars: {task_vars["_runtime"]["mapped_vars_vars"]}'
                )

        for var in mapped_vars:
            if var not in task_vars:
                # input, output, depends might not exist
                continue
            if not task_vars[var]:
                continue
            elif isinstance(task_vars[var], str):
                runtime[task_id][var] = self._map_var(task_vars[var])
                env.log_to_file(
                    'TASK',
                    f'On {self.alias}: ``{var}`` = {short_repr(task_vars[var])}'
                )
            elif isinstance(task_vars[var], (Sequence, set)):
                runtime[task_id][var] = type(task_vars[var])(
                    self._map_var(task_vars[var]))
                env.log_to_file(
                    'TASK',
                    f'On {self.alias}: ``{var}`` = {short_repr(task_vars[var])}'
                )
            else:
                env.logger.warning(
                    f'Failed to map {var} of type {task_vars[var].__class__.__name__}'
                )

        # master task??
        if hasattr(params, 'task_stack'):
            for tid, tdef in params.task_stack:
                runtime[tid] = {}
                for var in mapped_vars:
                    if var not in tdef.sos_dict:
                        # input, output, depends might not exist
                        continue
                    if not tdef.sos_dict[var]:
                        continue
                    elif isinstance(tdef.sos_dict[var], str):
                        runtime[tid][var] = self._map_var(tdef.sos_dict[var])
                    elif isinstance(tdef.sos_dict[var], (Sequence, set)):
                        runtime[tid][var] = type(tdef.sos_dict[var])(
                            self._map_var(tdef.sos_dict[var]))
                    else:
                        env.logger.warning(
                            f'Failed to map {var} of type {tdef.sos_dict[var].__class__.__name__}'
                        )

        # server restrictions #488
        for key in ('max_mem', 'max_cores', 'max_walltime'):
            if key in self.config:
                runtime['_runtime'][key] = format_HHMMSS(
                    self.config[key]
                ) if key == 'max_walltime' else self.config[key]

        # only update task file if there are runtime information
        if len(runtime) > 1 or runtime['_runtime'] or runtime != old_runtime:
            tf.runtime = runtime

        tf.status = 'pending'
        self.send_task_file(task_file)

    def send_task_file(self, task_file):
        send_cmd = cfg_interpolate(
            'ssh ' + self.cm_opts +
            ' -q {address} -p {port} "[ -d ~/.sos/tasks ] || mkdir -p ~/.sos/tasks" && '
            + 'rsync --ignore-existing -a --no-g -e "ssh ' + self.cm_opts +
            ' -q -p {port}" {job_file:ap} {address}:.sos/tasks/', {
                'job_file': sos_targets(task_file),
                'address': self.address,
                'port': self.port
            })
        # use scp for this simple case
        try:
            subprocess.check_call(send_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f'Failed to copy job {task_file} to {self.alias} using command {send_cmd}: {e}'
            )

    def check_output(self, cmd: object, under_workdir=False) -> object:
        if isinstance(cmd, list):
            cmd = subprocess.list2cmdline(cmd)
        try:
            cmd = cfg_interpolate(
                self._get_execute_cmd(under_workdir=under_workdir), {
                    'host': self.address,
                    'port': self.port,
                    'cmd': cmd,
                    'workdir': self._map_var(os.getcwd())
                })
        except Exception as e:
            raise ValueError(
                f'Failed to run command {cmd}: {e} ({env.sos_dict["CONFIG"]})')
        if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
            env.log_to_file('TASK', f'Executing command ``{cmd}``')
        try:
            return subprocess.check_output(cmd, shell=True).decode()
        except Exception as e:
            env.logger.debug(f'Check output of {cmd} failed: {e}')
            raise

    def check_call(self, cmd, under_workdir=False, **kwargs):
        if isinstance(cmd, list):
            cmd = subprocess.list2cmdline(cmd)
        try:
            cmd = cfg_interpolate(
                self._get_execute_cmd(under_workdir=under_workdir), {
                    'host': self.address,
                    'port': self.port,
                    'cmd': cmd,
                    'workdir': self._map_var(os.getcwd())
                })
        except Exception as e:
            raise ValueError(f'Failed to run command {cmd}: {e}')
        if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
            env.log_to_file('TASK', f'Executing command ``{cmd}``')
        try:
            return subprocess.check_call(cmd, shell=True, **kwargs)
        except Exception as e:
            env.logger.debug(f'Check output of {cmd} failed: {e}')
            raise

    def run_command(self, cmd, wait_for_task, realtime=False, **kwargs):
        if isinstance(cmd, list):
            cmd = subprocess.list2cmdline(cmd)
        try:
            cmd = cfg_interpolate(
                self._get_execute_cmd(under_workdir=False), {
                    'host': self.address,
                    'port': self.port,
                    'cmd': cmd,
                    'workdir': self._map_var(os.getcwd())
                })
        except Exception as e:
            raise ValueError(f'Failed to run command {cmd}: {e}')
        if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
            env.log_to_file('TASK', f'Executing command ``{cmd}``')
        if realtime:
            from .utils import pexpect_run
            return pexpect_run(cmd)
        elif wait_for_task or sys.platform == 'win32':
            # keep proc persistent to avoid a subprocess is still running warning.
            return subprocess.Popen(cmd, shell=True, **kwargs)
        else:
            p = DaemonizedProcess(cmd, **kwargs)
            p.start()
            p.join()

    def receive_result(self, task_id: str) -> Dict[str, int]:
        # for filetype in ('res', 'status', 'out', 'err'):
        sys_task_dir = os.path.join(os.path.expanduser('~'), '.sos', 'tasks')
        # use -p to preserve modification times so that we can keep the job status locally.
        receive_cmd = cfg_interpolate(
            "scp -P {port} -p -q {address}:.sos/tasks/{task}.* {sys_task_dir}",
            {
                'port': self.port,
                'address': self.address,
                'task': task_id,
                'sys_task_dir': sys_task_dir
            })
        # it is possible that local files are readonly (e.g. a pluse file) so we first need to
        # make sure the files are readable and remove them. Also, we do not want any file that is
        # obsolete to appear as new after copying
        for lfile in glob.glob(
                os.path.join(
                    os.path.expanduser('~'), '.sos', 'tasks', task_id + '.*')):
            if not os.access(lfile, os.W_OK):
                os.chmod(lfile, stat.S_IREAD | stat.S_IWRITE)
            os.remove(lfile)
        if 'TASK' in env.config['SOS_DEBUG'] or 'ALL' in env.config['SOS_DEBUG']:
            env.log_to_file('TASK', receive_cmd)
        ret = subprocess.call(
            receive_cmd,
            shell=True,
            stderr=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL)
        if (ret != 0):
            # this time try to get error message
            ret = subprocess.call(receive_cmd, shell=True)
            if (ret != 0):
                raise RuntimeError(
                    'Failed to retrieve result of job {} from {} with cmd\n{}'
                    .format(task_id, self.alias, receive_cmd))

        tf = TaskFile(task_id)
        params = tf.params
        res = tf.result

        if not res:
            env.logger.debug(
                f'Result for {task_id} is not received (no result)')
            return {
                'ret_code': 1,
                'task': task_id,
                'exception': RuntimeError(f'Task {task_id} failed or aborted'),
                'output': sos_targets()
            }

        if ('ret_code' in res and res['ret_code'] != 0) or ('succ' in res and
                                                            res['succ'] != 0):
            _show_err_and_out(task_id, res)
            env.logger.info(f'Ignore remote results for failed job {task_id}')
            return res

        if env.verbosity >= 3:
            _show_err_and_out(task_id, res)
        # do we need to copy files? We need to consult original task file
        # not the converted one
        job_dict = params.sos_dict
        if '_output' in job_dict and job_dict['_output'] and not isinstance(
                job_dict['_output'],
                Undetermined) and env.config['run_mode'] != 'dryrun':
            received = self.receive_from_host(
                [x for x in job_dict['_output'] if isinstance(x, (str, path))])
            if received:
                env.logger.info(
                    f'{task_id} ``received`` {short_repr(received.keys())} from {self.alias}'
                )
        if 'from_host' in job_dict[
                '_runtime'] and env.config['run_mode'] != 'dryrun':
            if isinstance(job_dict['_runtime']['from_host'], (dict, str)):
                received = self.receive_from_host(
                    job_dict['_runtime']['from_host'])
                if received:
                    env.logger.info(
                        f'{task_id} ``received`` {short_repr(received.keys())} from {self.alias}'
                    )
            else:
                env.logger.warning(
                    f"Expecting a dictionary from from_host: {job_dict['_runtime']['from_host']} received"
                )
        # we need to translate result from remote path to local
        if 'output' in res:
            if '_output' not in job_dict:
                # this should exist, but let us check it for robustness
                env.logger.warning('Missing _output in task dict')
                res['output'] = sos_targets()
            elif job_dict['_output'].undetermined():
                res['output'] = sos_targets(
                    self._reverse_map_var(res['output']))
            else:
                res['output'] = job_dict['_output']
        if 'subtasks' in res:
            for tid, subparams in params.task_stack:
                if tid in res['subtasks'] and 'output' in res['subtasks'][tid]:
                    if '_output' not in subparams.sos_dict:
                        # this should not happen
                        env.logger.warning('Missing _output in subparams')
                        res['subtasks'][tid]['output'] = sos_targets()
                    elif subparams.sos_dict['_output'].undetermined():
                        res['subtasks'][tid]['output'] = sos_targets(
                            self._reverse_map_var(
                                res['subtasks'][tid]['output']))
                    else:
                        res['subtasks'][tid]['output'] = subparams.sos_dict[
                            '_output']
        return res


#
# host instances are shared by all tasks so there should be only one
# instance for each host.
#


class Host:
    host_instances: Dict = {}

    def __init__(self, alias: Optional[str] = '',
                 start_engine: bool = True) -> None:
        # a host started from Jupyter notebook might not have proper stdout
        # (a StringIO) and cannot be used for DaemonizedFork).
        self._get_config(alias)
        self._get_host_agent(start_engine)

    # for test purpose
    @classmethod
    def reset(cls) -> None:
        for host in cls.host_instances.values():
            del host._task_engine
        cls.host_instances = {}

    def _get_local_host(self) -> str:
        if 'CONFIG' not in env.sos_dict or 'hosts' not in env.sos_dict['CONFIG']:
            from .utils import load_config_files
            load_config_files()
        # look for an entry with gethost
        if 'hosts' not in env.sos_dict['CONFIG']:
            return 'localhost'
        # try host name
        hostname = socket.gethostname().lower()
        for host, host_info in env.sos_dict['CONFIG']['hosts'].items():
            # find by key hostname
            if 'hostname' in host_info and host_info['hostname'].lower(
            ) == hostname:
                return host
        for host, host_info in env.sos_dict['CONFIG']['hosts'].items():
            # find by key hostname
            if 'hostname' in host_info and (
                    host_info['hostname'].lower().split('.')[0] == hostname or
                    host_info['hostname'].lower() == hostname.split('.')[0]):
                return host
        for host, host_info in env.sos_dict['CONFIG']['hosts'].items():
            # find by alias
            if host.lower() == hostname:
                return host
        for host, host_info in env.sos_dict['CONFIG']['hosts'].items():
            # find by address
            if 'address' in host_info and (
                    host_info['address'].split('@')[-1].lower() == hostname or
                    host_info['address'].split(
                        '.', 1)[0].split('@')[-1].lower() == hostname):
                return host
        # try IP Address
        hostname = socket.gethostname()
        ips = socket.gethostbyname_ex(hostname)[2]
        ips = [ip for ip in ips if not ip.startswith("127.")]
        # if the IP matches any address?
        for host, host_info in env.sos_dict['CONFIG']['hosts'].items():
            # find by key hostname
            if 'address' in host_info and any(
                    ip == host_info['address'].split('@')[-1] for ip in ips):
                return host
        # now check if a key localhost is defined
        if 'localhost' in env.sos_dict['CONFIG']:
            if env.sos_dict['CONFIG']['localhost'] not in env.sos_dict[
                    'CONFIG']['hosts']:
                raise ValueError(
                    f"Undefined localhost {env.sos_dict['CONFIG']['localhost']}"
                )
            return env.sos_dict['CONFIG']['localhost']
        raise ValueError(
            "No localhost could be identified from hostname, ip address, or a localhost key in config file"
        )

    def _get_remote_host(self, alias: Optional[str]) -> str:
        # get a remote host specified by Alias
        if 'CONFIG' not in env.sos_dict or 'hosts' not in env.sos_dict['CONFIG']:
            from .utils import load_config_files
            load_config_files()
        if not alias or alias == 'localhost':
            return self._get_local_host()
        if not isinstance(alias, str):
            raise ValueError(f'A string is expected for host {alias}')
        if 'hosts' not in env.sos_dict['CONFIG'] or alias not in env.sos_dict[
                'CONFIG']['hosts']:
            raise ValueError(f'Undefined remote host {alias}')
        return alias

    def _get_config(self, alias: Optional[str]) -> None:
        LOCAL = self._get_local_host()
        REMOTE = self._get_remote_host(alias)
        self.alias = REMOTE

        # now we need to find definition for local and remote host
        if LOCAL == 'localhost' and REMOTE == 'localhost':
            self.config = {
                'address': 'localhost',
                'alias': 'localhost',
            }
        elif 'hosts' in env.sos_dict['CONFIG']:
            if LOCAL not in env.sos_dict['CONFIG']['hosts']:
                raise ValueError(f'No hosts definition for local host {LOCAL}')
            if REMOTE not in env.sos_dict['CONFIG']['hosts']:
                raise ValueError(
                    f'No hosts definition for remote host {REMOTE}')

            # now we have definition for local and remote hosts
            cfg = env.sos_dict['CONFIG']['hosts']
            self.config = {
                x: y
                for x, y in cfg[self.alias].items()
                if x not in ('paths', 'shared')
            }

            # if local and remote hosts are the same
            if LOCAL == REMOTE or (
                    'address' in env.sos_dict['CONFIG']['hosts'][REMOTE] and
                    env.sos_dict['CONFIG']['hosts'][REMOTE]['address'] ==
                    'localhost'):
                # there would be no path map
                self.config['path_map'] = []
                self.config['shared'] = ['/']
                # override address setting to use localhost
                self.config['address'] = 'localhost'
            else:
                if 'address' not in env.sos_dict['CONFIG']['hosts'][REMOTE]:
                    raise ValueError(
                        f'No address defined for remote host {REMOTE}')
                self.config['path_map'] = []

                def normalize_value(x):
                    x = cfg_interpolate(x)
                    return x if x.endswith(os.sep) else (x + os.sep)

                if 'shared' in cfg[LOCAL] and 'shared' in cfg[REMOTE]:
                    common = set(cfg[LOCAL]['shared'].keys()) & set(
                        cfg[REMOTE]['shared'].keys())
                    if common:
                        self.config['shared'] = [
                            normalize_value(cfg[LOCAL]['shared'][x])
                            for x in common
                        ]
                        self.config['path_map'] = [
                            f'{normalize_value(cfg[LOCAL]["shared"][x])} -> {normalize_value(cfg[REMOTE]["shared"][x])}'
                            for x in common
                        ]
                # if paths are defined for both local and remote host, define path_map
                if ('paths' in cfg[LOCAL] and
                        cfg[LOCAL]['paths']) and ('paths' in cfg[REMOTE] and
                                                  cfg[REMOTE]['paths']):
                    if any(k not in cfg[REMOTE]['paths']
                           for k in cfg[LOCAL]['paths'].keys()):
                        env.logger.debug(
                            f'One or more local paths {", ".join(cfg[LOCAL]["paths"].keys())} cannot be mapped to remote host {REMOTE} with paths {",".join(cfg[REMOTE]["paths"].keys())}'
                        )
                    #
                    self.config['path_map'].extend([
                        f'{normalize_value(cfg[LOCAL]["paths"][x])} -> {normalize_value(cfg[REMOTE]["paths"][x])}'
                        for x in cfg[LOCAL]['paths'].keys()
                        if x in cfg[REMOTE]['paths']
                    ])
        elif LOCAL == REMOTE:
            # now we have checked local and remote are not defined, but they are the same, so
            # it is safe to assume that they are both local hosts
            self.config = {
                'address': 'localhost',
                'alias': LOCAL,
            }
        else:
            raise ValueError(
                f'Undefined local and remote hosts {LOCAL} and {REMOTE}.')
        #
        self.config['alias'] = self.alias
        self.description = self.config.get('description', '')

        # standardize parameters max_walltime, max_cores, and max_mem for the host
        if 'max_walltime' in self.config:
            self.config['max_walltime'] = format_HHMMSS(
                self.config['max_walltime'])
        if 'max_cores' in self.config:
            if not isinstance(self.config['max_cores'], int):
                raise ValueError('An integer is expected for max_cores')
        if 'max_mem' in self.config:
            self.config['max_mem'] = expand_size(self.config['max_mem'])

    def _get_host_agent(self, start_engine: bool) -> None:
        if 'queue_type' not in self.config:
            self._task_engine_type = 'process'
        else:
            self._task_engine_type = self.config['queue_type'].strip()
        if self.alias not in self.host_instances:
            if self.config['address'] == 'localhost':
                self.host_instances[self.alias] = LocalHost(self.config)
            else:
                self.host_instances[self.alias] = RemoteHost(self.config)

            if self._task_engine_type == 'process':
                task_engine = BackgroundProcess_TaskEngine(
                    self.host_instances[self.alias])
            else:
                task_engine = None

                available_engines = []
                for entrypoint in pkg_resources.iter_entry_points(
                        group='sos_taskengines'):
                    try:
                        if entrypoint.name == self._task_engine_type:
                            task_engine = entrypoint.load()(
                                self.host_instances[self.alias])
                            break
                        available_engines.append(entrypoint.name)
                    except Exception as e:
                        raise RuntimeError(
                            f'Failed to load task engine {self._task_engine_type}: {e}'
                        )

                if task_engine is None:
                    raise RuntimeError(
                        f'Failed to locate task engine type {self._task_engine_type}. Available engine types are {", ".join(available_engines)}'
                    )

            self.host_instances[self.alias]._task_engine = task_engine
            # the task engine is a thread and will run continously
            if start_engine:
                self.host_instances[self.alias]._task_engine.start()

        self._host_agent = self.host_instances[self.alias]
        # for convenience
        self._task_engine = self._host_agent._task_engine

    # public interface
    #
    # based on Host definition
    #
    def send_to_host(self, items):
        return self._host_agent.send_to_host(items)

    def receive_from_host(self, items):
        return self._host_agent.receive_from_host(items)

    def map_var(self, rvars):
        return self._host_agent._map_var(rvars)

    def submit_task(self, task_id: str) -> str:
        return self._task_engine.submit_task(task_id)

    def check_status(self, tasks: List[str]) -> List[str]:
        # find the task engine
        return [self._task_engine.check_task_status(task) for task in tasks]

    def retrieve_results(self, tasks: List[str]) -> Dict[str, Union[
            Dict[str, Union[int, str, Dict[int, Dict[Any, Any]], float]],
            Dict[str, int], Dict[str, Union[int, str, Dict[file_target, str],
                                            Dict[int, Dict[Any, Any]], float]],
            Dict[str, Union[int, str, Dict[int, Dict[str, int]], float]]]]:
        return {task: self._host_agent.receive_result(task) for task in tasks}


def list_queues(cfg, hosts=[], verbosity=1):
    env.verbosity = 1
    all_hosts = cfg.get('hosts', [])
    if not all_hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    for host in hosts:
        if host not in all_hosts:
            env.logger.warning(f'Undefined host {host}')
    host_description = [['Alias', 'Address', 'Queue Type', 'Description'],
                        ['-----', '-------', '----------', '-----------']]
    for host in sorted([x for x in hosts
                        if x in all_hosts] if hosts else all_hosts):
        try:
            h = Host(host, start_engine=False)
        except Exception as e:
            if verbosity == 0:
                print(f'{host} ({e})')
            elif verbosity in (1, 2):
                host_description.append([host, '?', '?', str(e)])
            else:
                print(f'Queue:       {host}')
                print(f'Error:       {str(e)}')
                if isinstance(cfg['hosts'][host], dict):
                    print('Configuration:')
                    for key in cfg['hosts'][host].keys():
                        print(
                            f'  {(key + ":").ljust(24)} {cfg["hosts"][host][key]}'
                        )
                print()
            continue
        if verbosity == 0:
            print(h.alias)
        elif verbosity in (1, 2):
            host_description.append([
                h.alias, h._host_agent.address, h._task_engine_type,
                h.description
            ])
        else:
            print(f'Queue:       {h.alias}')
            print(f'Address:     {h._host_agent.address}')
            print(f'Queue Type:  {h._task_engine_type}')
            print(f'Description: {h.description}')
            print('Configuration:')
            keys = sorted(h.config.keys())
            for key in keys:
                print(f'  {(key + ":").ljust(24)} {h.config[key]}')
            print()
    if verbosity in (1, 2):
        width = [(len(x) for x in row) for row in host_description]
        max_width = [max(x) for x in zip(*width)]
        print('\n'.join(' '.join([t.ljust(w)
                                  for t, w in zip(row, max_width)])
                        for row in host_description))


def status_of_queues(cfg, hosts=[], verbosity=1):
    env.verbosity = 1
    all_hosts = cfg.get('hosts', [])
    if not all_hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    for host in hosts:
        if host not in all_hosts:
            env.logger.warning(f'Undefined host {host}')
    host_description = [[
        'Alias', 'Address', 'Queue Type', 'Running', 'Pending', 'Completed'
    ], ['-----', '-------', '----------', '-------', '-------', '---------']]
    for host in sorted([x for x in hosts
                        if x in all_hosts] if hosts else all_hosts):
        try:
            h = Host(host, start_engine=True)
            status = h._task_engine.query_tasks(tasks=[], verbosity=0)
            if not status:
                raise ValueError("Failed to query status")
        except Exception as e:
            if verbosity == 0:
                print(f'{host} ({e})')
            elif verbosity in (1, 2):
                host_description.append([host, '?', '?', '?', '?', '?'])
            else:
                print(f'Queue:       {host}')
                print(f'Error:       {str(e)}')
                if isinstance(cfg['hosts'][host], dict):
                    print('Configuration:')
                    for key in cfg['hosts'][host].keys():
                        print(
                            f'  {(key + ":").ljust(24)} {cfg["hosts"][host][key]}'
                        )
                print()
            continue
        status = [x.strip() for x in status.splitlines() if x.strip()]
        running = str(status.count('running'))
        pending = str(status.count('pending'))
        completed = str(status.count('completed'))

        if verbosity == 0:
            print(f'{h.alias} {running} {pending} {completed}')
        elif verbosity in (1, 2):
            host_description.append([
                h.alias, h._host_agent.address, h._task_engine_type, running,
                pending, completed
            ])
        else:
            print(f'Queue:       {h.alias}')
            print(f'Address:     {h._host_agent.address}')
            print(f'Queue Type:  {h._task_engine_type}')
            print(f'Description: {h.description}')
            for k in set(status):
                if k not in ('running', 'pending', 'completed'):
                    print(f'{(k+":").ljust(26)} {status.count(k)}')
            print('Configuration:')
            keys = sorted(h.config.keys())
            for key in keys:
                print(f'  {(key + ":").ljust(24)} {h.config[key]}')
            print()
    if verbosity in (1, 2):
        width = [(len(x) for x in row) for row in host_description]
        max_width = [max(x) for x in zip(*width)]
        print('\n'.join(' '.join([t.ljust(w)
                                  for t, w in zip(row, max_width)])
                        for row in host_description))


def test_ssh(host):
    if host.address == 'localhost':
        return 'OK'
    address, port = host.address, host.port
    try:
        cmd = cfg_interpolate('ssh -q {host} -p {port} true', {
            'host': address,
            'port': port
        })
        env.logger.info(cmd)
        p = pexpect.spawn(cmd)
        # could be prompted for Password or password, so use assword
        i = p.expect([
            "(?i)are you sure you want to continue connecting", "[pP]assword:",
            pexpect.EOF
        ],
                     timeout=5)
        if i == 0:
            p.sendline('yes')
            p.expect([
                "(?i)are you sure you want to continue connecting",
                "[pP]assword:", pexpect.EOF
            ],
                     timeout=5)
        if i == 1:
            p.close(force=True)
            stty_sane()
            return f'ssh connection to {address} was prompted for password. Please set up public key authentication to the remote host before continue.'
        if i == 2:
            if not p.before:
                return "OK"
            p.close(force=True)
            return p.before.decode()
    except pexpect.TIMEOUT:
        return f'ssh connection to {address} time out with prompt: {str(p.before)}'
    except Exception as e:
        return f'Failed to check remote connection {address}:{port}: {e}'
    return "OK"


def test_scp(host):
    if host.address == 'localhost':
        return 'OK'
    # test send task file
    import random
    tID = random.randint(1, 100000)
    task_filename = os.path.join(
        os.path.expanduser('~'), '.sos', f'test_{tID}.tmp')
    with open(task_filename, 'w') as test_task:
        test_task.write('test task')
    # test scp
    try:
        host.send_task_file(task_filename)
    except Exception as e:
        return str(e)
    # test remove file using ssh
    try:
        host.check_call(f'rm -f ~/.sos/test_{tID}.tmp')
    except Exception as e:
        return str(e)
    # test rsync
    return 'OK'


def test_sos(host):
    # test the execution of sos commands
    try:
        ret = host.check_call(
            'sos -h', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if ret == 0:
            return 'OK'
        else:
            return 'sos not installed or not accessible on host.'
    except Exception as e:
        return str(e)


def test_paths(host):
    if host.address == 'localhost':
        return 'OK'
    # shared means, if localhost creates a file, it should be
    # instantly available on the remote host
    if not host.path_map:
        return 'No path_map between local and remote host.'
    import random
    tID = random.randint(1, 100000)
    for local, remote in host.path_map.items():
        if local in host.shared_dirs:
            # will be tested by 'shared'
            continue
        # now, let us see if two directory has the same files?
        if not os.path.isdir(local):
            return f'Mapped directory {local} does not exist.'
        # remote?
        try:
            host.check_output(f'ls -a {path(remote):q}')
        except:
            return f'Failed to access shared directory {remote} on remote host.'

        # test if local directory is writable
        try:
            with open(os.path.join(local, f".sos_test_{tID}.txt"),
                      'w') as tFile:
                tFile.write(f'{tID}')
        except:
            return f'Failed to write to mapped directory {local}'
        # test if file can be sent
        try:
            host.send_to_host(os.path.join(local, f".sos_test_{tID}.txt"))
        except Exception as e:
            return f'Failed to send files under {local} to remote host under {remote}: {e}'
        # the file should be available on remote host
        try:
            remote_content = host.check_output(
                f'cat {remote}/.sos_test_{tID}.txt')
        except Exception as e:
            return f'Failed to send files under {local} to remote host under {remote}: {e}'
        if remote_content != str(tID):
            return f'Content of file sent does not match: {tID} sent, {remote_content} received'
        # test retrieving files
        # remove local file
        os.remove(os.path.join(local, f".sos_test_{tID}.txt"))
        # copy file back
        try:
            host.receive_from_host(os.path.join(local, f".sos_test_{tID}.txt"))
        except Exception as e:
            return f'Failed to receive file from remote host {remote}: {e}'
        #
        if not os.path.isfile(os.path.join(local, f".sos_test_{tID}.txt")):
            return f'Failed to receive file from remote host {remote}: file does not exist'
        # check file content?
        with open(os.path.join(local, f".sos_test_{tID}.txt"), 'r') as tFile:
            remote_content = tFile.read()
        if remote_content != str(tID):
            return f'Content of received file does not match: {tID} expected, {remote_content} received.'
        # if everything ok, remove local and remote test files
        os.remove(os.path.join(local, f".sos_test_{tID}.txt"))
        #
        try:
            remote_content = host.check_output(
                f'rm {remote}/.sos_test_{tID}.txt')
        except Exception as e:
            return f'Failed to remove test file on remote host: {e}'
    return 'OK'


def test_shared(host):
    if host.address == 'localhost':
        return 'OK (localhost)'
    # shared means, if localhost creates a file, it should be
    # instantly available on the remote host
    for local in host.shared_dirs:
        if local not in host.path_map:
            return f'shared directory {local} not in path_map'
        # now, let us see if two directory has the same files?
        if not os.path.isdir(local):
            return f'shared directory {local} does not exist.'
        local_files = os.listdir(local)
        # remote?
        remote = host.path_map[local]
        try:
            remote_files = host.check_output(f'ls -a {path(remote):q}')
        except:
            return f'Failed to access shared directory {remote} on remote host.'
        remote_files = [
            x for x in remote_files.splitlines() if x not in ('.', '..')
        ]
        #
        if sorted(local_files) != sorted(remote_files):
            return f'shared directory {local} has different content on remote host under {remote}'

    return f'OK (shared {" ".join(host.shared_dirs)})'


def stty_sane():
    try:
        subprocess.check_call('stty sane', shell=True)
    except:
        pass


def test_queue(host):
    try:
        h = Host(host, start_engine=False)
    except Exception:
        return [host, '?', '?', '-', '-', '-', '-', '-']
    ssh_res = test_ssh(h._host_agent)
    return [
        h.alias, h._host_agent.address, h._task_engine_type, ssh_res,
        test_scp(h._host_agent) if ssh_res.startswith('OK') else '-',
        test_sos(h._host_agent) if ssh_res.startswith('OK') else '-',
        test_paths(h._host_agent) if ssh_res.startswith('OK') else '-',
        test_shared(h._host_agent) if ssh_res.startswith('OK') else '-'
    ]


def test_queues(cfg, hosts=[], verbosity=1):
    env.verbosity = 1
    all_hosts = cfg.get('hosts', [])
    if not all_hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    host_description = [[
        'Alias', 'Address', 'Queue Type', 'ssh', 'scp', 'sos', 'paths', 'shared'
    ], [
        '-----', '-------', '----------', '---', '---', '---', '-----', '------'
    ]]
    for host in hosts:
        if host not in all_hosts:
            env.logger.warning(f'Undefined host {host}')
    hosts = [x for x in hosts if x in all_hosts] if hosts else all_hosts
    if not hosts:
        return
    from multiprocessing import Pool
    pool = Pool(min(len(hosts), 10))
    host_description.extend(pool.map(test_queue, hosts))
    if verbosity == 0:
        # just print succ or self
        for hd in host_description:
            print(f'{hd[0]} {"OK" if all(x=="OK" for x in hd[3:]) else "FAIL"}')
    elif verbosity in (1, 2):
        shortened = host_description[:2]
        for row in host_description[2:]:
            shortened.append(row[:3] + [
                'OK' if x.startswith('OK') else ('-' if x == '-' else 'FAIL')
                for x in row[3:]
            ])
        width = [(len(x) for x in row) for row in shortened]
        max_width = [max(x) for x in zip(*width)]
        print('\n'.join(' '.join([t.ljust(w)
                                  for t, w in zip(row, max_width)])
                        for row in shortened))
        if any('FAILED' in row for row in shortened):
            print(
                '\nUse command "sos remote --test host -v3" to check details of hosts with failed tests.'
            )
    else:
        for row in host_description[2:]:
            print(f'Alias:       {row[0]}')
            print(f'Address:     {row[1]}')
            print(f'Queue Type:  {row[2]}')
            print(f'ssh:         {row[3]}')
            print(f'scp:         {row[4]}')
            print(f'sos:         {row[5]}')
            print(f'paths:       {row[6]}')
            print(f'shared:      {row[7]}')
            print()


def copy_public_key(host, agent, password):
    try:
        if password is None:
            import getpass
            password = getpass.getpass(
                f'Please enter password for {agent.address}: ')
        cmd = f"scp -P {agent.port if agent.port else 22} {os.path.expanduser('~')}/.ssh/id_rsa.pub {agent.address}:id_rsa.pub.{host}"
        env.logger.info(cmd)
        p = pexpect.spawn(cmd, echo=False)
        i = p.expect([
            "(?i)are you sure you want to continue connecting", "[pP]assword:",
            pexpect.EOF
        ])
        if i == 0:
            p.sendline('yes')
            p.expect([
                "(?i)are you sure you want to continue connecting",
                "[pP]assword:", pexpect.EOF
            ],
                     timeout=5)

        if i == 1:
            p.sendline(password)
            i = p.expect(['assword:', pexpect.EOF])
            if i == 0:
                p.close(force=True)
                return f'Incorrect password specified (you can try to specify it with command line option --password)'
        if i == 2:
            p.close()
            return f'Failed to copy public key to {agent.address}'
    except Exception as e:
        p.close()
        return f'Failed to copy public key to {host}: {e}'
    #
    # ssh
    try:
        cmd = f"ssh {agent.address} -p {agent.port} '[ -d .ssh ] || mkdir .ssh && chmod 700 .ssh; cat id_rsa.pub.{host} >> .ssh/authorized_keys; rm -f id_rsa.pub.{host}'"
        env.logger.info(cmd)
        p = pexpect.spawn(cmd, echo=False)
        i = p.expect([
            "(?i)are you sure you want to continue connecting", 'assword:',
            pexpect.EOF
        ])
        if i == 0:
            p.sendline('yes')
            p.expect([
                "(?i)are you sure you want to continue connecting",
                "[pP]assword:", pexpect.EOF
            ],
                     timeout=5)

        if i == 1:
            p.sendline(password)
            i = p.expect(['assword:', pexpect.EOF])
            if i == 0:
                p.close(force=True)
                return f'Incorrect password specified (you can try to specify it with command line option --password)'
        elif i != 1:
            p.close()
            return f'Failed to append public key to .ssh/authorized_keys'
    except Exception as e:
        p.close()
        return f'Failed to append public key to .ssh/authorized_keys: {e}'
    return 'OK'


def setup_remote_access(cfg, hosts=[], password='', verbosity=1):
    env.verbosity = verbosity
    all_hosts = cfg.get('hosts', [])
    if not all_hosts and not hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    for host in hosts:
        if host not in all_hosts:
            env.logger.warning(
                f'Treating undefined host {host} as address of a remote host.')
    # public_key
    public_key = os.path.join(os.path.expanduser('~'), '.ssh', 'id_rsa.pub')

    for host in sorted(hosts if hosts else all_hosts):
        try:
            if host in all_hosts:
                h = Host(host, start_engine=False)
                host_agent = h._host_agent
            else:
                from argparse import Namespace
                host_agent = Namespace(address=host, port=22)
        except Exception as e:
            env.logger.error(f'Failed to start task engine {host}: {e}')
            continue

        if os.path.isfile(public_key):
            env.logger.info('Using existing public key .ssh/id_rsa.pub')
        else:
            env.logger.info(f'Public key {public_key} is found. Creating one.')
            try:
                subprocess.check_call(
                    'echo | ssh-keygen -t rsa',
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)
                if not os.path.isfile(public_key):
                    raise RuntimeError('public key not found after ssh-keygen')
            except Exception as e:
                raise RuntimeError(f'Failed to create public key: {e}')
        #
        # can ssh?
        response = test_ssh(host_agent)
        if response.startswith('OK'):
            env.logger.info(
                f'Public key access is already enabled for host ``{host}`` with address ``{host_agent.address}``'
            )
            continue
        elif 'Could not resolve hostname' in response:
            env.logger.error(response)
            sys.exit(1)
        #
        response = copy_public_key(host, host_agent, password)
        stty_sane()
        if not response.startswith('OK'):
            env.logger.error(response)
            sys.exit(1)
        # file copied, check ssh again.
        response = test_ssh(host_agent)
        if response.startswith('OK'):
            env.logger.info(
                f'Public key access is successfully set up for host ``{host}`` with address ``{host_agent.address}``'
            )
            continue
        else:
            env.logger.error(
                f'Failed to connect to {host} after passing public key. Possible problems include permission of .ssh and home directories.'
            )


def login_host(cfg, host):
    try:
        h = Host(host, start_engine=False)
    except Exception as e:
        raise ValueError(f'Unrecognized or invalid host {host}: {e}')

    address, port = h._host_agent.address, h._host_agent.port
    try:
        env.logger.info(f'Running ``ssh {address} -p {port}``')
        os.execvp('ssh', ['ssh', address, '-p', str(port)])
    except Exception as e:
        raise RuntimeError(f'Failed to log in to {host}: {e}')


def run_command_on_hosts(cfg, hosts, cmd, verbosity):
    env.verbosity = verbosity
    if not hosts:
        hosts = cfg.get('hosts', [])
    if not hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    for host in hosts:
        # runing command on all hosts
        try:
            env.logger.info(f'Running ``{" ".join(cmd)}`` on ``{host}``')
            h = Host(host, start_engine=False)
            print(h._host_agent.check_output(cmd, under_workdir=True))
        except Exception as e:
            from .utils import get_traceback
            if verbosity and verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error(e)


def push_to_hosts(cfg, hosts, items, verbosity):
    env.verbosity = verbosity
    if not hosts:
        hosts = cfg.get('hosts', [])
    if not hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    for host in hosts:
        try:
            env.logger.info(f'Pushing ``{" ".join(items)}`` to ``{host}``')

            h = Host(host, start_engine=False)
            #
            sent = h.send_to_host(items)
            #
            print('{} item{} sent:\n{}'.format(
                len(sent), ' is' if len(sent) <= 1 else 's are', '\n'.join([
                    '{} => {}'.format(x, sent[x]) for x in sorted(sent.keys())
                ])))
        except Exception as e:
            from .utils import get_traceback
            if verbosity and verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error(e)
            sys.exit(1)


def pull_from_host(cfg, hosts, items, verbosity):
    env.verbosity = verbosity
    if not hosts:
        hosts = cfg.get('hosts', [])
    if not hosts:
        env.logger.warning(
            "No remote host or task queue is defined in ~/.sos/hosts.yml.")
        return
    if len(hosts) > 1:
        raise ValueError('Can only pull from a single remote host.')
    try:
        env.logger.info(f'Pulling ``{" ".join(items)}`` from ``{hosts[0]}``')

        host = Host(hosts[0], start_engine=False)
        #
        received = host.receive_from_host(items)
        #
        print('{} item{} received:\n{}'.format(
            len(received), ' is' if len(received) <= 1 else 's are', '\n'.join([
                '{} <= {}'.format(x, received[x])
                for x in sorted(received.keys())
            ])))
    except Exception as e:
        from .utils import get_traceback
        if verbosity and verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)