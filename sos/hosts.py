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

class LocalHost:
    '''For local host, no path map, send and receive ...'''
    def __init__(self):
        self.alias = 'localhost'

    def send_to_host(self, items):
        pass

    def receive_from_host(self, items):
        pass

    def map_vars(self, vars):
        return vars

    def send_task(self, task):
        pass

    def run_command(self, cmd):
        p = subprocess.Popen(cmd, shell=True)
        #ret = p.wait()
        #if (ret != 0):
        #    raise RuntimeError('Failed to execute {}'.format(cmd))
        #return ret
        return p

    def check_output(self, cmd):
        return subprocess.check_output(cmd, shell=True)

    def receive_result(self, task_id):
        res_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.res')
        with open(res_file, 'rb') as result:
            res = pickle.load(result)
        return res 

class RemoteHost:
    '''A remote host class that manages how to communicate with remote host'''
    def __init__(self, config):
        self.config = config
        #
        self.address = self.config['address']
        self.shared_dirs = self._get_shared_dirs()
        self.path_map = self._get_path_map()
        self.send_cmd = self._get_send_cmd()
        self.receive_cmd = self._get_receive_cmd()
        self.execute_cmd = self._get_execute_cmd()
        self.query_cmd = self._get_query_cmd()

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
            '''ssh ${host} "mkdir -p ${dest!dq}"; rsync -av ${source!ae} "${host}:${dest!de}"''')

    def _get_receive_cmd(self):
        return self.config.get('receive_cmd',
            'mkdir -p ${dest!dq}; rsync -av ${host}:${source!ae} "${dest!de}"')

    def _get_execute_cmd(self):
        return self.config.get('execute_cmd',
            '''ssh ${host} "nohup bash --login -c '${cmd}' > ~/.sos/tasks/${task}.out 2> ~/.sos/tasks/${task}.err" & ''')

    def _get_query_cmd(self):
        return self.config.get('query_cmd',
            '''ssh ${host} "bash --login -c 'sos status ${task} -v 0'" ''')

    def is_shared(self, path):
        fullpath = os.path.abspath(os.path.expanduser(path))
        for dir in self.shared_dirs:
            if fullpath.startswith(dir):
                return True
        return False

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
        send_cmd = 'scp {} {}:.sos/tasks'.format(job_file, self.address)
        # use scp for this simple case
        ret = subprocess.call(send_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        if (ret != 0):
            raise RuntimeError('Failed to copy job {} to {}'.format(task_id, self.alias))

    def receive_result(self, task_id):
        receive_cmd = 'scp {}:.sos/tasks/{}.res {}/.sos/tasks'.format(self.address, task_id, os.path.expanduser('~'))
        ret = subprocess.checcall(receive_cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
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
        env.logger.info('Executing command ``{}``'.format(cmd))
        env.logger.debug(cmd)
        return subprocess.check_output(cmd, shell=True)

    def run_command(self, cmd):
        try:
            cmd = interpolate(self.execute_cmd, '${ }', {
                'host': self.address,
                'cmd': cmd})
        except Exception as e:
            raise ValueError('Failed to run command {}: {}'.format(cmd, e))
        env.logger.info('Executing command ``{}``'.format(cmd))
        env.logger.debug(cmd)
        p = subprocess.Popen(cmd, shell=True)
        ret = p.wait()
        if (ret != 0):
            raise RuntimeError('Failed to execute {}'.format(cmd))
        return ret


class Host:
    def __init__(self, host = ''):
        if not host:
            self.alias = 'losthost'
            self.config = {'alias': 'localhost'}
        elif isinstance(host, str):
            # read from configuration file
            self.alias = host
            if 'hosts' in env.sos_dict['CONFIG'] and \
                self.alias in env.sos_dict['CONFIG']['hosts']:
                self.config = env.sos_dict['CONFIG']['hosts'][host]
                if 'address' not in self.config:
                    self.config['address'] = host
            else:
                self.config = { 'address': host }
        elif isinstance(host, dict):
            if 'address' not in host:
                raise ValueError('Please define at least "address" for host specification')
            self.config = host
            self.alias = self.config.get('alias', self.config['address'])

        if self.alias == 'losthost':
            self._host_agent = LocalHost()
        else:
            self._host_agene = RemoteHost(self.config)

        if 'task_engine' not in self.config:
            self._task_engine_type = 'background_execution'
            self._task_engine = BackgroundProcess_TaskEngine(self._host_agent)
        else:
            self._task_engine_type = self.config['task_engine']
            self._task_engine = None

            available_engines = []
            for entrypoint in pkg_resources.iter_entry_points(group='sos_taskengines'):
                try:
                    available_engines.append(entrypoint.name)
                    if entrypoint.name == self._task_engine_type:
                        self.task_engine = entrypoint.load()(self._host_agent)
                except Exception as e:
                    env.logger.debug('Failed to load task engine {}: {}'.format(self._task_engine_type, e))

            if self._task_engine is None:
                raise RuntimeError('Failed to locate task engine type {}. Available engine types are {}'.format(
                    self.task_engine, ', '.join(available_engines)))

    def send_to_host(self, items):
        return self._host_agent.send_to_host(items)

    def receive_from_host(self, items):
        return self._host_agent.receive_from_host(items)

    def map_vars(self, vars):
        return self._host_agent.map_vars(vars)

    def submit_task(self, task_id):
        self._host_agent.send_task(task_id)
        env.logger.info('{} ``submitted``'.format(task_id))
        return self._task_engine.submit_task(task_id)
        
    def query_task(self, task_id):
        return self._task_engine.query_task(task_id)

    def kill_task(self, task_id):
        return self._task_engine.kill_task(task_id)

    def wait_task(self, task_id):
        st = time.time()
        while True:
            status = self.query_task(task_id).decode().strip()
            if status not in ('pending', 'running', 'completed-old', 'failed-old', 'failed-missing-output', 'failed-old-missing-output'):
                break
            elapsed = time.time() - st
            # the longer the wait, the less frequent the check
            time.sleep(max(2, math.log(elapsed, 1.3)))
        if status == 'completed':
            env.logger.info('{} ``completed``'.format(task_id))
            return self._host_agent.receive_result(task_id)
        raise RuntimeError('Job returned with status {}'.format(status))


