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
import time
import math
import pickle
import pkg_resources
from collections.abc import Sequence

from .utils import env
from .sos_eval import interpolate
from .sos_task import BackgroundProcess_TaskEngine

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
        self.config = {'alias': 'localhost'}

    def send_to_host(self, items):
        pass

    def receive_from_host(self, items):
        pass

    def map_var(self, vars):
        return vars

    def send_task(self, task):
        # on the same file system, no action is needed.
        pass

    def run_command(self, cmd):
        # run command but does not wait for result.
        p = subprocess.Popen(cmd, shell=True)
        #ret = p.wait()
        #if (ret != 0):
        #    raise RuntimeError('Failed to execute {}'.format(cmd))
        #return ret
        return p

    def check_output(self, cmd):
        # get the output of command
        try:
            return subprocess.check_output(cmd, shell=True)
        except Exception as e:
            env.logger.warning('Check output of {} failed: {}'.format(cmd, e))
            return ''

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
    def map_var(self, source):
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

    def send_to_host(self, items):
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

    def receive_from_host(self, items):
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

    def send_task(self, task_id):
        job_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
        send_cmd = 'ssh -q {1} "[ -d ~/.sos/tasks ] || mkdir -p ~/.sos/tasks ]"; scp -q {0} {1}:.sos/tasks/'.format(job_file, self.address)
        # use scp for this simple case
        ret = subprocess.call(send_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if (ret != 0):
            raise RuntimeError('Failed to copy job {} to {} using command {}'.format(task_id, self.alias, send_cmd))

    def receive_result(self, task_id):
        receive_cmd = 'scp -q {}:.sos/tasks/{}.res {}/.sos/tasks'.format(self.address, task_id, os.path.expanduser('~'))
        ret = subprocess.call(receive_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if (ret != 0):
            raise RuntimeError('Failed to retrieve result of job {} from {}'.format(task_id, self.alias))
        res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.res')
        with open(res_file, 'rb') as result:
            res = pickle.load(result)
        return res

    def check_output(self, cmd):
        try:
            cmd = interpolate(self.execute_cmd, '${ }', {
                'host': self.address,
                'cmd': cmd})
        except Exception as e:
            raise ValueError('Failed to run command {}: {}'.format(cmd, e))
        env.logger.debug('Executing command ``{}``'.format(cmd))
        try:
            return subprocess.check_output(cmd, shell=True)
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

#
# host instances are shared by all tasks so there should be only one
# instance for each host.
#

class Host:
    host_instances = {}

    @classmethod
    def pending_tasks(cls):
        tasks = []
        for host in cls.host_instances:
            tasks.extend(host._task_engine.pending_tasks)
        return tasks

    @classmethod
    def running_tasks(cls):
        tasks = []
        for host in cls.host_instances:
            tasks.extend(host._task_engine.running_tasks)
        return tasks

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
        if not self._task_engine.is_alive():
            # wait a bit
            time.sleep(1)
            if not self._task_engine.is_alive():
                env.logger.warning('Restart non-working task engine')
                self._task_engine.start()

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
        self._host_agent.send_task(task_id)
        env.logger.info('{} ``queued``'.format(task_id))
        return self._task_engine.submit_task(task_id)
        
    def kill_task(self, task_id):
        return self._task_engine.kill_task(task_id)

    def wait_task(self, task_id):
        st = time.time()
        while True:
            status = self._task_engine.check_task_status(task_id)
            if status not in ('pending', 'running', 'completed-old', 'failed-old', 'failed-missing-output', 'failed-old-missing-output'):
                break
            elapsed = time.time() - st
            # the longer the wait, the less frequent the check
            time.sleep(max(2, math.log(elapsed, 1.3)))
        if status == 'completed':
            env.logger.info('{} ``completed``'.format(task_id))
            return self._host_agent.receive_result(task_id)
        raise RuntimeError('Job returned with status {}'.format(status))


