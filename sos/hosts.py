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
import os
import sys
import stat

import subprocess
import multiprocessing as mp
import pickle
import shutil
import glob
import pkg_resources
from collections.abc import Sequence

from .utils import env, pickleable, short_repr, load_config_files, expand_size, format_HHMMSS, expand_time
from .sos_eval import interpolate, Undetermined
from .sos_task import BackgroundProcess_TaskEngine, TaskParams, loadTask

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
            env.logger.error('_Fork #1 failed: {0}\n'.format(err))
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
            env.logger.error('_Fork #2 failed: {0}\n'.format(err))
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
        except:
            # #493
            pass

        # fork a new process
        subprocess.Popen(self.cmd, shell=True, close_fds=True)
        return

class LocalHost:
    '''For local host, no path map, send and receive ...'''

    def __init__(self, alias='localhost'):
        self.alias = alias
        self.address = 'localhost'
        # we checkk local jobs more aggressively
        self.config = {'alias': 'localhost', 'status_check_interval': 2}
        self._procs = []

    def prepare_task(self, task_id):
        def_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.def')
        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
        # add server restriction on task file
        params = loadTask(def_file)
        task_vars = params.sos_dict

        if 'max_mem' not in self.config and 'max_cores' not in self.config and 'max_walltime' not in self.config:
            shutil.copyfile(def_file, task_file)
        else:

            task_vars['_runtime']['max_mem'] = self.config.get('max_mem', None)
            task_vars['_runtime']['max_cores'] = self.config.get('max_cores', None)
            task_vars['_runtime']['max_walltime'] = self.config.get('max_walltime', None)
            if task_vars['_runtime']['max_walltime'] is not None:
                task_vars['_runtime']['max_walltime'] = format_HHMMSS(task_vars['_runtime']['max_walltime'])

            if self.config.get('max_mem', None) is not None and task_vars['_runtime'].get('mem', None) is not None \
                    and self.config['max_mem'] < task_vars['_runtime']['mem']:
                raise ValueError('Task {} requested more mem ({}) than allowed max_mem ({})'.format(
                    task_id, task_vars['_runtime']['mem'], self.config['max_mem']))
            if self.config.get('max_cores', None) is not None and task_vars['_runtime'].get('cores', None) is not None \
                    and self.config['max_cores'] < task_vars['_runtime']['cores']:
                raise ValueError('Task {} requested more cores ({}) than allowed max_cores ({})'.format(
                    task_id, task_vars['_runtime']['cores'], self.config['max_cores']))
            if self.config.get('max_walltime', None) is not None and task_vars['_runtime'].get('walltime', None) is not None \
                    and expand_time(self.config['max_walltime']) < expand_time(task_vars['_runtime']['walltime']):
                raise ValueError('Task {} requested more walltime ({}) than allowed max_walltime ({})'.format(
                    task_id, task_vars['_runtime']['walltime'], self.config['max_walltime']))

            new_param = TaskParams(
                name = params.name,
                task = params.task,
                sos_dict = params.sos_dict,
                sigil = params.sigil
            )
            new_param.save(task_file)
        #
        if 'to_host' in task_vars['_runtime'] and isinstance(task_vars['_runtime']['to_host'], dict):
            for l, r in task_vars['_runtime']['to_host'].items():
                if l != r:
                    shutil.copy(l, r)
        return True

    def send_task_file(self, task_file):
        # on the same file system, no action is needed.
        pass

    def check_output(self, cmd):
        # get the output of command
        try:
            return subprocess.check_output(cmd, shell=True).decode()
        except Exception as e:
            env.logger.warning('Check output of {} failed: {}'.format(cmd, e))
            raise

    def run_command(self, cmd, wait_for_task):
        # run command but does not wait for result.
        if wait_for_task or sys.platform == 'win32':
            self._procs.append(subprocess.Popen(cmd, shell=True))
        else:
            p = DaemonizedProcess(cmd)
            p.start()
            p.join()

    def receive_result(self, task_id):
        res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.res')
        try:
            with open(res_file, 'rb') as result:
                res = pickle.load(result)
        except Exception as e:
            env.logger.warning('Failed to receive result for task {}: {}'.format(task_id, e))
            return {'ret_code': 1, 'output': {}}

        sys_task_dir = os.path.join(os.path.expanduser('~'), '.sos', 'tasks')
        task_file = os.path.join(sys_task_dir, task_id + '.def')
        params = loadTask(task_file)
        job_dict = params.sos_dict

        if 'from_host' in job_dict['_runtime'] and isinstance(job_dict['_runtime']['from_host'], dict):
            for l, r in job_dict['_runtime']['from_host'].items():
                if l != r:
                    shutil.copy(r, l)
        return res


class RemoteHost:
    '''A remote host class that manages how to communicate with remote host'''
    def __init__(self, config):
        self.config = config
        self.alias = self.config['alias']
        #
        self.address = self.config['address']
        self.port = self.config.get('port', 22)
        self.shared_dirs = self._get_shared_dirs()
        self.path_map = self._get_path_map()
        self.execute_cmd = self._get_execute_cmd()
        self._procs = []

        self.task_dir = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', self.alias)
        if not os.path.isdir(self.task_dir):
            os.mkdir(self.task_dir)

    def _get_shared_dirs(self):
        value = self.config.get('shared', [])
        if isinstance(value, str):
            return [value]
        elif isinstance(value, Sequence):
            return value
        else:
            raise ValueError('Option shared can only be a string or a list of strings')

    def _get_path_map(self):
        res = {}
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
                    raise ValueError('Path map should be separated as from -> to, {} specified'.format(v))
                elif v.count(' -> ') > 1:
                    raise ValueError('Path map should be separated as from -> to, {} specified'.format(v))
                res[v.split(' -> ')[0]] = v.split(' -> ')[1]
        elif isinstance(path_map, dict):
            for k,v in path_map.items():
                res[k] = v
        else:
            raise ValueError('Unacceptable path_mapue for configuration path_map: {}'.format(path_map))
        return res

    def _get_send_cmd(self, rename=False):
        if rename:
            return '''ssh -q ${host} -p ${port} "mkdir -p ${dest!dpq}" && ''' + \
                   '''rsync -av -e 'ssh -p ${port}' ${source!aep} "${host}:${dest!dep}" && ''' + \
                   '''ssh -q ${host} -p ${port} "mv ${dest!dep}/${source!b} ${dest!ep}" '''
        else:
            return '''ssh -q ${host} -p ${port} "mkdir -p ${dest!dpq}" && rsync -av -e 'ssh -p ${port}' ${source!aep} "${host}:${dest!dep}"'''

    def _get_receive_cmd(self, rename=False):
        if rename:
            return '''rsync -av -e 'ssh -p ${port}' ${host}:${source!e} "${dest!adep}" && ''' + \
                    '''mv "${dest!adep}/${source!b}" "${dest!aep}"'''
        else:
            return '''rsync -av -e 'ssh -p ${port}' ${host}:${source!e} "${dest!adep}"'''

    def _get_execute_cmd(self):
        return self.config.get('execute_cmd',
            '''ssh -q ${host} -p ${port} "bash --login -c '${cmd}'" ''')

    def _get_query_cmd(self):
        return self.config.get('query_cmd',
            '''ssh -q ${host} -p ${port} "bash --login -c 'sos status ${task} -v 0'" ''')

    def is_shared(self, path):
        fullpath = os.path.abspath(os.path.expanduser(path))
        for dir in self.shared_dirs:
            if fullpath.startswith(dir):
                return True
        return False

    def _map_path(self, source):
        result = {}
        cwd = os.getcwd()
        if isinstance(source, str):
            dest = os.path.abspath(os.path.expanduser(source))
            # we use samefile to avoid problems with case-insensitive file system #522
            # we also use the "cwd" name to avoid wrong case for cwd. For example,
            # if the cwd = '/Users/user/Project'
            # then, dest = '/USERS/USER/PROJECT/a.txt'
            # would be converted to '/Users/user/Project/a.txt' before path mapping
            if os.path.exists(dest[:len(cwd)]) and os.path.samefile(dest[:len(cwd)], cwd):
                dest = cwd + dest[len(cwd):]
            matched = [k for k in self.path_map.keys() if os.path.exists(dest[:len(k)]) and os.path.samefile(dest[:len(k)], k)]
            if matched:
                # pick the longest key that matches
                k = max(matched, key=len)
                dest = self.path_map[k] + dest[len(k):]
            else:
                env.logger.warning('Path {} is not under any specified paths of localhost and is mapped to {} on remote host.'.format(source, dest))
            result[source] = dest.replace('\\', '/')
        elif isinstance(source, Sequence):
            for src in source:
                result.update(self._map_path(src))
        else:
            env.logger.debug('Ignore unmappable source {}'.format(source))
            return {source: source}
        return result

    #
    # Interface functions
    #
    def _map_var(self, source):
        cwd = os.getcwd()
        if isinstance(source, str):
            dest = os.path.abspath(os.path.expanduser(source))
            # we use samefile to avoid problems with case-insensitive file system #522
            # we also use the "cwd" name to avoid wrong case for cwd. For example,
            # if the cwd = '/Users/user/Project'
            # then, dest = '/USERS/USER/PROJECT/a.txt'
            # would be converted to '/Users/user/Project/a.txt' before path mapping
            if os.path.exists(dest[:len(cwd)]) and os.path.samefile(dest[:len(cwd)], cwd):
                dest = cwd + dest[len(cwd):]
            matched = [k for k in self.path_map.keys() if os.path.exists(dest[:len(k)]) and os.path.samefile(dest[:len(k)], k)]
            if matched:
                # pick the longest key that matches
                k = max(matched, key=len)
                dest = self.path_map[k] + dest[len(k):]
            else:
                env.logger.warning('Path {} is not under any specified paths of localhost and is mapped to {} on remote host.'.format(source, dest))
            return dest.replace('\\', '/')
        elif isinstance(source, Sequence):
            ret = [self._map_var(x) for x in source]
            return [x for x in ret if x is not None]
        else:
            env.logger.debug('Ignore unmappable source {}'.format(source))
            return source

    def _send_to_host(self, items):
        # we only copy files and directories, not other types of targets
        if isinstance(items, str):
            items = [items]
        elif isinstance(items, Sequence):
            items = [x for x in items if isinstance(x, str)]
        elif isinstance(items, dict):
            for x,y in items.items():
                if not isinstance(x, str):
                    env.logger.warning('Unrecognized item to be sent to host: {}'.format(x))
                if not isinstance(y, str):
                    env.logger.warning('Unrecognized item to be sent to host: {}'.format(y))
            items = {x:y for x,y in items.items() if isinstance(x, str) and isinstance(y, str)}
        else:
            env.logger.warning('Unrecognized items to be sent to host: {}'.format(items))
            return

        if isinstance(items, Sequence):
            from .utils import find_symbolic_links
            new_items = []
            for item in items:
                links = find_symbolic_links(item)
                for link, realpath in links.items():
                    env.logger.info('Adding {} for symbolic link {}'.format(realpath, link))
                new_items.extend(links.values())
            items.extend(new_items)

            sending = self._map_path(items)
        else:
            sending = items

        for source in sorted(sending.keys()):
            if self.is_shared(source):
                env.logger.debug('Skip sending {} on shared file system'.format(source))
            else:
                dest = sending[source]
                env.logger.info('Sending ``{}`` to {}:{}'.format(source, self.alias, dest))
                cmd = interpolate(self._get_send_cmd(rename=os.path.basename(source) != os.path.basename(dest)),
                        '${ }', {'source': source.rstrip('/'), 'dest': dest, 'host': self.address, 'port': self.port})
                env.logger.debug(cmd)
                ret = subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                if (ret != 0):
                    raise RuntimeError('Failed to copy {} to {} using command "{}". The remote host might be unavailable.'.format(source, self.alias, cmd))

    def _receive_from_host(self, items):
        if isinstance(items, dict):
            # specify as local:remote
            # needs remote:local
            receiving = {y:x for x,y in items.items()}
        else:
            receiving = {y:x for x,y in self._map_path(items).items()}
        #
        received = []
        for source in sorted(receiving.keys()):
            dest = receiving[source]
            dest_dir = os.path.dirname(dest)
            if dest_dir and not os.path.isdir(dest_dir):
                try:
                    os.path.makedirs(dest_dir)
                except Exception as e:
                    env.logger.error('Failed to create destination directory {}'.format(dest_dir))
            if self.is_shared(dest) and os.path.basename(source) == os.path.basename(dest):
                env.logger.debug('Skip retrieving ``{}`` from shared file system'.format(dest))
            else:
                cmd = interpolate(self._get_receive_cmd(rename=os.path.basename(source) != os.path.basename(dest)),
                    '${ }', {'source': source.rstrip('/'), 'dest': dest, 'host': self.address, 'port': self.port})
                env.logger.debug(cmd)
                try:
                    ret = subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                    if (ret != 0):
                        raise RuntimeError('command return {}'.format(ret))
                    received.append(dest)
                except Exception as e:
                    raise  RuntimeError('Failed to copy {} from {} using command "{}": {}'.format(source, self.alias, cmd, e))
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
        def_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.def')
        params = loadTask(def_file)
        task_vars = params.sos_dict

        if self.config.get('max_mem', None) is not None and task_vars['_runtime'].get('mem', None) is not None \
                and self.config['max_mem'] < task_vars['_runtime']['mem']:
            raise ValueError('Task {} requested more mem ({}) than allowed max_mem ({})'.format(
                task_id, task_vars['_runtime']['mem'], self.config['max_mem']))
        if self.config.get('max_cores', None) is not None and task_vars['_runtime'].get('cores', None) is not None \
                and self.config['max_cores'] < task_vars['_runtime']['cores']:
            raise ValueError('Task {} requested more cores ({}) than allowed max_cores ({})'.format(
                task_id, task_vars['_runtime']['cores'], self.config['max_cores']))
        if self.config.get('max_walltime', None) is not None and task_vars['_runtime'].get('walltime', None) is not None \
                and expand_time(self.config['max_walltime']) < expand_time(task_vars['_runtime']['walltime']):
            raise ValueError('Task {} requested more walltime ({}) than allowed max_walltime ({})'.format(
                task_id, task_vars['_runtime']['walltime'], self.config['max_walltime']))

        if task_vars['_input'] and not isinstance(task_vars['_input'], Undetermined):
            self._send_to_host(task_vars['_input'])
            env.logger.info('{} ``send`` {}'.format(task_id, short_repr(task_vars['_input'])))
        if task_vars['_depends'] and not isinstance(task_vars['_depends'], Undetermined):
            self._send_to_host(task_vars['_depends'])
            env.logger.info('{} ``send`` {}'.format(task_id, short_repr(task_vars['_depends'])))
        if 'to_host' in task_vars['_runtime']:
            if isinstance(task_vars['_runtime']['to_host'], dict):
                th = {}
                for x,y in task_vars['_runtime']['to_host'].items():
                    if y.startswith('/'):
                        th[x] = y
                    elif y.startswith('~'):
                        th[x] = self._map_var(task_vars['_runtime']['home_dir']) + y[1:]
                    else:
                        th[x] = self._map_var(task_vars['_runtime']['cur_dir']) + '/' + y
                self._send_to_host(th)
            else:
                self._send_to_host(task_vars['_runtime']['to_host'])
            env.logger.info('{} ``send`` {}'.format(task_id, short_repr(task_vars['_runtime']['to_host'])))

        # map variables
        mvars = ['_input', '_output', '_depends', 'input', 'output', 'depends', '_runtime',
            '_local_input_{}'.format(task_vars['_index']),
            '_local_output_{}'.format(task_vars['_index'])] + list(task_vars['__signature_vars__'])
        seen = set()
        # avoid variable to be handled twice
        vars = [x for x in mvars if not (x in seen or seen.add(x))]
        preserved = set()
        if 'preserved_vars' in task_vars['_runtime']:
            if isinstance(task_vars['_runtime']['preserved_vars'], str):
                preserved.add(task_vars['_runtime']['preserved_vars'])
            elif isinstance(task_vars['_runtime']['preserved_vars'], (set, Sequence)):
                preserved |= set(task_vars['_runtime']['preserved_vars'])
            else:
                raise ValueError('Unacceptable value for runtime option preserved_vars: {}'.format(task_vars['_runtime']['preserved_vars']))
        env.logger.debug('Translating {}'.format(vars))
        for var in vars:
            if var in preserved:
                env.logger.debug('Value of variable {} is preserved'.format(var))
            elif var == '_runtime':
                task_vars[var]['cur_dir'] = self._map_var(task_vars[var]['cur_dir'])
                task_vars[var]['home_dir'] = self._map_var(task_vars[var]['home_dir'])
                if 'workdir' in task_vars[var]:
                    task_vars[var]['workdir'] = self._map_var(task_vars[var]['workdir'])
            elif var in task_vars:
                if isinstance(task_vars[var], (type(None), int)) or not pickleable(task_vars[var], var):
                    continue
                try:
                    old_var = task_vars[var]
                    task_vars[var] = self._map_var(task_vars[var])
                    if not task_vars[var]:
                        continue
                    # looks a bit suspicious
                    if isinstance(old_var, str) and old_var != task_vars[var] and not os.path.exists(os.path.expanduser(old_var)) \
                            and os.sep not in old_var:
                        env.logger.warning('On {}: ``{}`` = {}'.format(self.alias, var, short_repr(task_vars[var])))
                    else:
                        env.logger.info('On {}: ``{}`` = {}'.format(self.alias, var, short_repr(task_vars[var])))
                except Exception as e:
                    env.logger.debug('Failed to map variable {}: {}'.format(var, e))
            else:
                env.logger.debug('Variable {} not in env.'.format(var))

        # server restrictions #488
        task_vars['_runtime']['max_mem'] = self.config.get('max_mem', None)
        task_vars['_runtime']['max_cores'] = self.config.get('max_cores', None)
        task_vars['_runtime']['max_walltime'] = self.config.get('max_walltime', None)
        if task_vars['_runtime']['max_walltime'] is not None:
            task_vars['_runtime']['max_walltime'] = format_HHMMSS(task_vars['_runtime']['max_walltime'])

        new_param = TaskParams(
            name = params.name,
            task = params.task,
            sos_dict = params.sos_dict,
            sigil = params.sigil
        )
        task_file = os.path.join(self.task_dir, task_id + '.task')
        new_param.save(task_file)
        self.send_task_file(task_id + '.task')

    def send_task_file(self, task_file):
        job_file = os.path.join(self.task_dir, task_file)
        send_cmd = interpolate('ssh -q ${address} -p ${port} "[ -d ~/.sos/tasks ] || mkdir -p ~/.sos/tasks" && scp -q -P ${port} ${job_file!ap} ${address}:.sos/tasks/',
                '${ }', {'job_file': job_file, 'address': self.address, 'port': self.port})
        # use scp for this simple case
        try:
            subprocess.check_output(send_cmd, shell=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError('Failed to copy job {} to {} using command {}: {}'.format(task_file, self.alias, send_cmd, e))

    def check_output(self, cmd):
        try:
            cmd = interpolate(self.execute_cmd, '${ }', {
                'host': self.address, 'port': self.port,
                'cmd': cmd})
        except Exception as e:
            raise ValueError('Failed to run command {}: {}'.format(cmd, e))
        env.logger.debug('Executing command ``{}``'.format(cmd))
        try:
            return subprocess.check_output(cmd, shell=True).decode()
        except Exception as e:
            env.logger.debug('Check output of {} failed: {}'.format(cmd, e))
            raise

    def run_command(self, cmd, wait_for_task):
        try:
            cmd = interpolate(self.execute_cmd, '${ }', {
                'host': self.address, 'port': self.port,
                'cmd': cmd})
        except Exception as e:
            raise ValueError('Failed to run command {}: {}'.format(cmd, e))
        env.logger.debug('Executing command ``{}``'.format(cmd))

        if wait_for_task or sys.platform == 'win32':
            # keep proc persistent to avoid a subprocess is still running warning.
            self._procs.append(subprocess.Popen(cmd, shell=True))
        else:
            p = DaemonizedProcess(cmd)
            p.start()
            p.join()

    def receive_result(self, task_id):
        # for filetype in ('res', 'status', 'out', 'err'):
        sys_task_dir = os.path.join(os.path.expanduser('~'), '.sos', 'tasks')
        # use -p to preserve modification times so that we can keep the job status locally.
        receive_cmd = interpolate("scp -P ${port} -p -q ${address}:.sos/tasks/${task}.* ${sys_task_dir}",
                '${ }', {'port': self.port, 'address': self.address, 'task': task_id, 'sys_task_dir': sys_task_dir})
        # it is possible that local files are readonly (e.g. a pluse file) so we first need to
        # make sure the files are readable and remove them. Also, we do not want any file that is
        # obsolete to appear as new after copying
        for lfile in glob.glob(os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.*')):
            if os.path.splitext(lfile)[-1] == '.def':
                continue
            if not os.access(lfile, os.W_OK):
                os.chmod(lfile, stat.S_IREAD | stat.S_IWRITE)
            os.remove(lfile)
        env.logger.debug(receive_cmd)
        ret = subprocess.call(receive_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if (ret != 0):
            # this time try to get error message
            ret = subprocess.call(receive_cmd, shell=True)
            if (ret != 0):
                raise RuntimeError('Failed to retrieve result of job {} from {} with cmd\n{}'.format(task_id, self.alias, receive_cmd))
        # show results? Not sure if this is a good idea but helps debugging at this point
        if env.verbosity >= 2:
            out_file = os.path.join(sys_task_dir, task_id + '.out')
            err_file = os.path.join(sys_task_dir, task_id + '.err')
            if os.path.isfile(out_file):
                env.logger.info('{}.out:'.format(task_id))
                with open(out_file) as out:
                    print(out.read())
            if os.path.isfile(err_file):
                env.logger.info('{}.err:'.format(task_id))
                with open(err_file) as err:
                    print(err.read())

        res_file = os.path.join(sys_task_dir, task_id + '.res')
        if not os.path.isfile(res_file):
            env.logger.debug('Result for {} is not received'.format(task_id))
            return {'ret_code': 1, 'output': {}}

        with open(res_file, 'rb') as result:
            res = pickle.load(result)

        if ('ret_code' in res and res['ret_code'] != 0) or ('succ' in res and res['succ'] != 0):
            env.logger.info('Ignore remote results for failed job {}'.format(task_id))
        else:
            # do we need to copy files? We need to consult original task file
            # not the converted one
            task_file = os.path.join(sys_task_dir, task_id + '.def')
            params = loadTask(task_file)
            job_dict = params.sos_dict
            #
            if job_dict['_output'] and not isinstance(job_dict['_output'], Undetermined):
                received = self._receive_from_host([x for x in job_dict['_output'] if isinstance(x, str)])
                if received:
                    env.logger.info('{} ``received`` {}'.format(task_id, short_repr(received)))
            if 'from_host' in job_dict['_runtime']:
                if isinstance(job_dict['_runtime']['from_host'], dict):
                    fh = {}
                    for x,y in job_dict['_runtime']['from_host'].items():
                        if y.startswith('/'):
                            fh[x] = y
                        elif y.startswith('~'):
                            fh[x] = self._map_var(job_dict['_runtime']['home_dir']) + y[1:]
                        else:
                            fh[x] = self._map_var(job_dict['_runtime']['cur_dir']) + '/' + y
                    received = self._receive_from_host(fh)
                else:
                    received = self._receive_from_host(job_dict['_runtime']['from_host'])
                if received:
                    env.logger.info('{} ``received`` {}'.format(task_id, short_repr(received)))
        return res


#
# host instances are shared by all tasks so there should be only one
# instance for each host.
#

class Host:
    host_instances = {}

    def __init__(self, alias='', start_engine=True):
        # a host started from Jupyter notebook might not have proper stdout
        # (a StringIO) and cannot be used for DaemonizedFork).
        self._get_config(alias)
        self._get_host_agent(start_engine)

    # for test purpose
    @classmethod
    def reset(cls):
        for host in cls.host_instances.values():
            del host._task_engine
        cls.host_instances = {}


    @classmethod
    def not_wait_for_tasks(cls):
        return all(host._task_engine.wait_for_task is False for host in cls.host_instances.values())

    def _get_config(self, alias):
        if 'CONFIG' not in env.sos_dict:
            from .utils import load_config_files
            env.sos_dict.set('CONFIG', load_config_files())
        if not alias or alias == 'localhost':
            # if no alias is specified, we are using localhost  ->  localhost
            if 'localhost' not in env.sos_dict['CONFIG']:
                # true default ... we are running localhost  ->  localhost without definition
                self.alias = 'localhost'
                LOCAL = 'localhost'
            else:
                # use local host ... it is possible that localhost is a queue system
                self.alias = env.sos_dict['CONFIG']['localhost']
                LOCAL = self.alias
        elif not isinstance(alias, str):
            raise ValueError('An alias or host address is expected. {} provided.'.format(self.alias))
        else:
            # specified "remote" host.
            # but then we would require a definition for localhost
            if 'localhost' not in env.sos_dict['CONFIG']:
                raise ValueError('localhost undefined in sos configuration file when a remote host {} is specified.'.format(alias))
            self.alias = alias
            LOCAL = env.sos_dict['CONFIG']['localhost']

        # just to make it clear that alias refers to remote_host
        REMOTE = self.alias
        # now we need to find definition for local and remote host
        if 'hosts' not in env.sos_dict['CONFIG']:
            if LOCAL == 'localhost' and REMOTE == 'localhost':
                self.config = {
                        'address': 'localhost',
                        'alias': 'localhost',
                }
            else:
                raise ValueError('No hosts definitions for local and remote hosts {} and {}'.format(LOCAL, REMOTE))
        else:
            if LOCAL not in env.sos_dict['CONFIG']['hosts']:
                raise ValueError('No hosts definition for local host {}'.format(LOCAL))
            if REMOTE not in env.sos_dict['CONFIG']['hosts']:
                raise ValueError('No hosts definition for remote host {}'.format(REMOTE))

            # now we have definition for local and remote hosts
            cfg = env.sos_dict['CONFIG']['hosts']
            self.config = {x:y for x,y in cfg[self.alias].items() if x not in ('paths', 'shared')}
            # if local and remote hosts are the same
            if LOCAL == REMOTE:
                # there would be no path map
                self.config['path_map'] = []
                self.config['shared'] = ['/']
                # override address setting to use localhost
                self.config['address'] = 'localhost'
            else:
                if 'address' not in env.sos_dict['CONFIG']['hosts'][REMOTE]:
                    raise ValueError('No address defined for remote host {}'.format(REMOTE))
                self.config['path_map'] = []
                def append_slash(x):
                    return x if x.endswith(os.sep) else (x + os.sep)
                if 'shared' in cfg[LOCAL] and 'shared' in cfg[REMOTE]:
                    common = set(cfg[LOCAL]['shared'].keys()) & set(cfg[REMOTE]['shared'].keys())
                    if common:
                        self.config['shared'] = [append_slash(cfg[LOCAL]['shared'][x]) for x in common]
                        self.config['path_map'] = ['{} -> {}'.format(append_slash(cfg[LOCAL]['shared'][x]), append_slash(cfg[REMOTE]['shared'][x])) \
                            for x in common]
                # if paths are defined for both local and remote host, define path_map
                if ('paths' in cfg[LOCAL] and cfg[LOCAL]['paths']) and ('paths' in cfg[REMOTE] and cfg[REMOTE]['paths']):
                    if any(k not in cfg[REMOTE]['paths'] for k in cfg[LOCAL]['paths'].keys()):
                        raise ValueError('One or more local paths {} cannot be mapped to remote host {} with paths {}'.format(
                            ','.join(cfg[LOCAL]['paths'].keys()), REMOTE, ','.join(cfg[REMOTE]['paths'].keys())))
                    #
                    self.config['path_map'].extend(['{} -> {}'.format(append_slash(cfg[LOCAL]['paths'][x]), append_slash(cfg[REMOTE]['paths'][x])) \
                        for x in cfg[LOCAL]['paths'].keys()])
        #
        self.config['alias'] = self.alias
        self.description = self.config.get('description', '')

        # standardize parameters max_walltime, max_cores, and max_mem for the host
        if 'max_walltime' in self.config:
            self.config['max_walltime'] = format_HHMMSS(self.config['max_walltime'])
        if 'max_cores' in self.config:
            if not isinstance(self.config['max_cores'], int):
                raise ValueError('An integer is expected for max_cores')
        if 'max_mem' in self.config:
            self.config['max_mem'] = expand_size(self.config['max_mem'])

    def _get_host_agent(self, start_engine):
        if self.alias not in self.host_instances:
            if self.config['address'] == 'localhost':
                self.host_instances[self.alias] = LocalHost()
            else:
                self.host_instances[self.alias] = RemoteHost(self.config)

            if 'queue_type' not in self.config:
                self._task_engine_type = 'process'
                task_engine = BackgroundProcess_TaskEngine(self.host_instances[self.alias])
            else:
                self._task_engine_type = self.config['queue_type'].strip()
                task_engine = None

                available_engines = []
                for entrypoint in pkg_resources.iter_entry_points(group='sos_taskengines'):
                    try:
                        available_engines.append(entrypoint.name)
                        if entrypoint.name == self._task_engine_type:
                            task_engine = entrypoint.load()(self.host_instances[self.alias])
                            break
                    except Exception as e:
                        env.logger.warning('Failed to load task engine {}: {}'.format(self._task_engine_type, e))

                if task_engine is None:
                    raise RuntimeError('Failed to locate task engine type {}. Available engine types are {}'.format(
                        self._task_engine_type, ', '.join(available_engines)))

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

    def map_var(self, vars):
        return self._host_agent._map_var(vars)

    def submit_task(self, task_id):
        return self._task_engine.submit_task(task_id)

    def check_status(self, tasks):
        # find the task engine
        return [self._task_engine.check_task_status(task) for task in tasks]

    def retrieve_results(self, tasks):
        return {task: self._host_agent.receive_result(task) for task in tasks}


def list_queues(config_file, verbosity = 1):
    cfg = load_config_files(config_file)
    env.sos_dict.set('CONFIG', cfg)
    hosts = cfg.get('hosts', [])
    if not hosts:
        sys.exit("No hosts is defined.")
    host_description = [['Alias', 'Address', 'Queue Type', 'Description'],
                        ['-----', '-------', '----------', '-----------']]
    for host in sorted(hosts):
        try:
            h = Host(host, start_engine=False)
        except Exception as e:
            env.logger.warning('Invalid remote host {} from localhost: {}'.format(host, e))
            continue
        if verbosity == 0:
            print(h.alias)
        elif verbosity in (1, 2):
            host_description.append([
                h.alias, h._host_agent.address, h._task_engine_type, h.description])
        else:
            print('Queue:       {}'.format(h.alias))
            print('Address:     {}'.format(h._host_agent.address))
            print('Queue Type:  {}'.format(h._task_engine_type))
            print('Description: {}'.format(h.description))
            print('Configuration:')
            keys = sorted(h.config.keys())
            for key in keys:
                print('  {} {}'.format((key + ':').ljust(24), h.config[key]))

            if verbosity == 4:
                if 'template_file' in keys:
                    template_file = h.config['template_file']
                    if not os.path.isfile(os.path.expanduser(template_file)):
                        env.warning('Missing template_file {}'.format(template_file))
                    else:
                        print('------ begin of {} -------------'.format(template_file))
                        with open(os.path.expanduser(template_file)) as tfile:
                            print(tfile.read())
                        print('------ end of file ---------------')
            print()
    if verbosity in (1, 2):
        width = [(len(x) for x in row) for row in host_description]
        max_width = [max(x) for x in zip(*width)]
        print('\n'.join(' '.join([t.ljust(w) for t,w in zip(row, max_width)]) for row in host_description))
