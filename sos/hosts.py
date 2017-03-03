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

import subprocess
import pickle
import pkg_resources
from collections.abc import Sequence

from .utils import env, pickleable, short_repr
from .sos_eval import interpolate, Undetermined
from .sos_task import BackgroundProcess_TaskEngine, TaskParams

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

class LocalHost:
    '''For local host, no path map, send and receive ...'''

    def __init__(self, alias='localhost'):
        self.alias = alias
        # we checkk local jobs more aggressively
        self.config = {'alias': 'localhost', 'status_check_interval': 2}

    def prepare_task(self, task_id):
        return task_id

    def send_task(self, task):
        # on the same file system, no action is needed.
        pass

    def check_output(self, cmd):
        # get the output of command
        try:
            return subprocess.check_output(cmd, shell=True).decode()
        except Exception as e:
            env.logger.warning('Check output of {} failed: {}'.format(cmd, e))
            return ''

    def run_command(self, cmd):
        # run command but does not wait for result.
        p = subprocess.Popen(cmd, shell=True)
        #ret = p.wait()
        #if (ret != 0):
        #    raise RuntimeError('Failed to execute {}'.format(cmd))
        #return ret
        return p

    def receive_result(self, task_id):
        res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.res')
        with open(res_file, 'rb') as result:
            res = pickle.load(result)
        return res


class RemoteHost:
    '''A remote host class that manages how to communicate with remote host'''
    def __init__(self, config):
        self.config = config
        self.alias = self.config['alias']
        #
        self.address = self.config['address']
        self.shared_dirs = self._get_shared_dirs()
        self.path_map = self._get_path_map()
        self.send_cmd = self._get_send_cmd()
        self.receive_cmd = self._get_receive_cmd()
        self.execute_cmd = self._get_execute_cmd()

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
        # use ordered map so that users can control the order
        # in which substitution happens.
        from collections import OrderedDict
        res = OrderedDict()
        # if user-specified path_map, it overrides CONFIG
        path_map = self.config.get('path_map', [])
        #
        if not path_map:
            return res
        if isinstance(path_map, str):
            path_map = [path_map]
        if isinstance(path_map, Sequence):
            for v in path_map:
                if ':' not in v:
                    raise ValueError('Path map should be separated as from:to, {} specified'.format(v))
                elif v.count(':') > 1:
                    raise ValueError('Path map should be separated as from:to, {} specified'.format(v))
                res[v.split(':')[0]] = v.split(':')[1]
        elif isinstance(path_map, dict):
            for k,v in path_map.items():
                res[k] = v
        else:
            raise ValueError('Unacceptable path_mapue for configuration path_map: {}'.format(path_map))
        return res

    def _get_send_cmd(self):
        return self.config.get('send_cmd',
            '''ssh -q ${host} "mkdir -p ${dest!dq}"; rsync -av ${source!ae} "${host}:${dest!de}"''')

    def _get_receive_cmd(self):
        return self.config.get('receive_cmd',
            'mkdir -p ${dest!dq}; rsync -av ${host}:${source!ae} "${dest!de}"')

    def _get_execute_cmd(self):
        return self.config.get('execute_cmd',
            '''ssh -q ${host} "bash --login -c '${cmd}'" ''')

    def _get_query_cmd(self):
        return self.config.get('query_cmd',
            '''ssh -q ${host} "bash --login -c 'sos status ${task} -v 0'" ''')

    def is_shared(self, path):
        fullpath = os.path.abspath(os.path.expanduser(path))
        for dir in self.shared_dirs:
            if fullpath.startswith(dir):
                return True
        return False

    def _map_path(self, source):
        result = {}
        if isinstance(source, str):
            dest = os.path.abspath(os.path.expanduser(source))
            for k,v in self.path_map.items():
                if dest.startswith(k):
                    dest = v + dest[len(k):]
            result[source] = dest
        elif isinstance(source, Sequence):
            for src in source:
                result.update(self._map_path(src))
        else:
            raise ValueError('Unacceptable parameter {} to option to_host'.format(source))
        return result

    #
    # Interface functions
    #
    def _map_var(self, source):
        if isinstance(source, str):
            dest = os.path.abspath(os.path.expanduser(source))
            for k,v in self.path_map.items():
                if dest.startswith(k):
                    dest = v + dest[len(k):]
            return dest
        elif isinstance(source, Sequence):
            return [self.map_var(x) for x in source]
        else:
            raise ValueError('Cannot map variables {} of type {}'.format(source, type(source).__name__))

    def _send_to_host(self, items):
        sending = self._map_path(items)
        for source in sorted(sending.keys()):
            if self.is_shared(source):
                env.logger.debug('Skip sending {} on shared file system'.format(source))
            else:
                dest = sending[source]
                env.logger.info('Sending ``{}`` to {}:{}'.format(source, self.alias, dest))
                cmd = interpolate(self.send_cmd, '${ }', {'source': source, 'dest': dest, 'host': self.address})
                env.logger.debug(cmd)
                ret = subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                if (ret != 0):
                    raise RuntimeError('Failed to copy {} to {}'.format(source, self.alias))

    def _receive_from_host(self, items):
        receiving = {y:x for x,y in self._map_path(items).items()}
        #
        for source in sorted(receiving.keys()):
            dest = receiving[source]
            if self.is_shared(dest):
                env.logger.info('Skip retrieving ``{}`` from shared file system'.format(dest))
            else:
                env.logger.info('Receiving ``{}`` from {}:{}'.format(dest, self.alias, source))
                cmd = interpolate(self.receive_cmd, '${ }', {'source': source, 'dest': dest, 'host': self.address})
                try:
                    ret = subprocess.call(cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
                    if (ret != 0):
                        raise RuntimeError('Failed to copy {} from {}'.format(source, self.alias))
                except Exception as e:
                    raise  RuntimeError('Failed to copy {} from {}: {}'.format(source, self.alias, e))

    #
    # Interface
    #
    #

    def prepare_task(self, task_id):
        task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
        with open(task_file, 'rb') as task:
            params = pickle.load(task)
            task_vars = params.data[1]

        if task_vars['_input'] and not isinstance(task_vars['_input'], Undetermined):
            self._send_to_host(task_vars['_input'])
        if task_vars['_depends'] and not isinstance(task_vars['_depends'], Undetermined):
            self._send_to_host(task_vars['_depends'])
        if 'to_host' in task_vars['_runtime']:
            self._send_to_host(task_vars['_runtime']['to_host'])

        # map variables
        vars = ['_input', '_output', '_depends', 'input', 'output', 'depends', '__report_output__', '_runtime',
            '_local_input_{}'.format(task_vars['_index']),
            '_local_output_{}'.format(task_vars['_index'])] + list(task_vars['__signature_vars__'])
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
                if 'workdir' in task_vars[var]:
                    task_vars[var]['workdir'] = self._map_var(task_vars[var]['workdir'])
            elif var in task_vars and pickleable(task_vars[var], var):
                try:
                    task_vars[var] = self._map_var(task_vars[var])
                    if task_vars[var] != task_vars[var]:
                        env.logger.info('On {}: ``{}`` = {}'.format(self.host.alias, var, short_repr(task_vars[var])))
                    else:
                        env.logger.debug('{}: {} is kept'.format(var, short_repr(task_vars[var])))
                except Exception as e:
                    env.logger.debug(e)
            else:
                env.logger.debug('Variable {} not in env.'.format(var))

        new_param = TaskParams(
            name = params.name,
            data = (
                params.data[0],
                task_vars,
                params.data[2],
            )
        )
        job_file = os.path.join(self.task_dir, task_id + '.task')
        with open(job_file, 'wb') as jf:
            try:
                pickle.dump(new_param, jf)
            except Exception as e:
                env.logger.warning(e)
                raise


    def send_task(self, task_id):
        job_file = os.path.join(self.task_dir, task_id + '.task')
        send_cmd = 'ssh -q {1} "[ -d ~/.sos/tasks ] || mkdir -p ~/.sos/tasks ]"; scp -q {0} {1}:.sos/tasks/'.format(job_file, self.address)
        # use scp for this simple case
        ret = subprocess.call(send_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if (ret != 0):
            raise RuntimeError('Failed to copy job {} to {} using command {}'.format(task_id, self.alias, send_cmd))

    def check_output(self, cmd):
        try:
            cmd = interpolate(self.execute_cmd, '${ }', {
                'host': self.address,
                'cmd': cmd})
        except Exception as e:
            raise ValueError('Failed to run command {}: {}'.format(cmd, e))
        env.logger.debug('Executing command ``{}``'.format(cmd))
        try:
            return subprocess.check_output(cmd, shell=True).decode()
        except Exception as e:
            env.logger.warning('Check output of {} failed: {}'.format(cmd, e))
            return ''

    def run_command(self, cmd):
        try:
            cmd = interpolate(self.execute_cmd, '${ }', {
                'host': self.address,
                'cmd': cmd})
        except Exception as e:
            raise ValueError('Failed to run command {}: {}'.format(cmd, e))
        env.logger.debug('Executing command ``{}``'.format(cmd))
        p = subprocess.Popen(cmd, shell=True)
        return p

    def receive_result(self, task_id):
        receive_cmd = 'scp -q {}:.sos/tasks/{}.res {}'.format(self.address, task_id, self.task_dir)
        ret = subprocess.call(receive_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if (ret != 0):
            raise RuntimeError('Failed to retrieve result of job {} from {}'.format(task_id, self.alias))
        res_file = os.path.join(self.task_dir, task_id + '.res')
        with open(res_file, 'rb') as result:
            res = pickle.load(result)

        if res['succ'] != 0:
            env.logger.error('Remote job failed.')
        else:
            # do we need to copy files? We need to consult original task file
            # not the converted one
            task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
            with open(task_file, 'rb') as task:
                params = pickle.load(task)
                job_dict = params.data[1]
            #
            if job_dict['_output'] and not isinstance(job_dict['_output'], Undetermined):
                self._receive_from_host(job_dict['_output'])
            if 'from_host' in job_dict['_runtime']:
                self._receive_from_host(job_dict['_runtime']['from_host'])
        return res


#
# host instances are shared by all tasks so there should be only one
# instance for each host.
#

class Host:
    host_instances = {}

    def __init__(self, alias=''):
        self._get_config(alias)
        self._get_host_agent()

    def _get_config(self, alias):
        if not alias:
            self.alias = 'localhost'
            self.config = {'alias': 'localhost'}
        elif isinstance(alias, str):
            # read from configuration file
            self.alias = alias
            if 'hosts' in env.sos_dict['CONFIG'] and \
                alias in env.sos_dict['CONFIG']['hosts']:
                self.config = env.sos_dict['CONFIG']['hosts'][alias]
                if 'address' not in self.config:
                    self.config['address'] = alias
                self.config['alias'] = alias
            else:
                self.config = { 'address': alias, 'alias': alias }
        elif isinstance(alias, dict):
            if 'address' not in alias:
                raise ValueError('Please define at least "address" for host specification')
            self.config = alias
            self.alias = self.config.get('alias', self.config['address'])

    def _get_host_agent(self):
        if self.alias not in self.host_instances:
            if self.alias == 'localhost':
                self.host_instances[self.alias] = LocalHost()
            else:
                self.host_instances[self.alias] = RemoteHost(self.config)

            if 'task_engine' not in self.config:
                task_engine_type = 'background_execution'
                task_engine = BackgroundProcess_TaskEngine(self.host_instances[self.alias])
            else:
                task_engine_type = self.config['task_engine']
                task_engine = None

                available_engines = []
                for entrypoint in pkg_resources.iter_entry_points(group='sos_taskengines'):
                    try:
                        available_engines.append(entrypoint.name)
                        if entrypoint.name == task_engine_type:
                            task_engine = entrypoint.load()(self.host_instances[self.alias])
                    except Exception as e:
                        env.logger.debug('Failed to load task engine {}: {}'.format(task_engine_type, e))

                if task_engine is None:
                    raise RuntimeError('Failed to locate task engine type {}. Available engine types are {}'.format(
                        task_engine, ', '.join(available_engines)))

            self.host_instances[self.alias]._task_engine = task_engine
            # the task engine is a thread and will run continously
            self.host_instances[self.alias]._task_engine.start()

        self._host_agent = self.host_instances[self.alias]
        # for convenience
        self._task_engine = self._host_agent._task_engine
        #
        # task engine can be unavailable because of time need to start it
        #
        #if not self._task_engine.is_alive():
        #    # wait a bit
        #    time.sleep(1)
        #    if not self._task_engine.is_alive():
        #        env.logger.warning('Restart non-working task engine')
        #        self._task_engine.start()

    # public interface
    #
    # based on Host definition
    #
    def send_to_host(self, items):
        return self._host_agent.send_to_host(items)

    def receive_from_host(self, items):
        return self._host_agent.receive_from_host(items)

    def map_var(self, vars):
        return self._host_agent.map_var(vars)

    def submit_task(self, task_id):
        #
        # Here we need to prepare the task for execution.
        #
        self._host_agent.prepare_task(task_id)
        #
        self._host_agent.send_task(task_id)
        env.logger.info('{} ``queued``'.format(task_id))
        return self._task_engine.submit_task(task_id)

    def kill_task(self, task_id):
        return self._task_engine.kill_task(task_id)

    def check_status(self, tasks):
        # find the task engine
        return [self._task_engine.check_task_status(task) for task in tasks]

    def retrieve_results(self, tasks):
        return {task: self._host_agent.receive_result(task) for task in tasks}


