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
import pickle

from sos.utils import env
from sos.sos_eval import SoS_exec

class TaskParams(object):
    '''A parameter object that encaptulates parameters sending to
    task executors. This would makes the output of workers, especially
    in the web interface much cleaner (issue #259)'''
    def __init__(self, name, data):
        self.name = name
        self.data = data

    def __repr__(self):
        return self.name


def execute_task(task_file, verbosity=None, sigmode=None):
    '''A function that execute specified task within a local dictionary
    (from SoS env.sos_dict). This function should be self-contained in that
    it can be handled by a task manager, be executed locally in a separate
    process or remotely on a different machine.'''
    with open(task_file, 'rb') as task:
        params = pickle.load(task)

    task, global_def, global_sigil, sos_dict, sigil = params.data
    if verbosity is not None:
        env.verbosity = verbosity
    if sigmode is not None:
        env.sigmode = sigmode
    env.register_process(os.getpid(), 'spawned_job with {} {}'
        .format(sos_dict['_input'], sos_dict['_output']))
    try:
        # set current directory if specified
        orig_dir = os.getcwd()
        if '_runtime' in sos_dict and 'workdir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['workdir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['workdir']))
                except Exception as e:
                    raise RuntimeError('Failed to create workdir {}'.format(sos_dict['_runtime']['workdir']))
            os.chdir(os.path.expanduser(sos_dict['_runtime']['workdir']))
        # set environ ...
        # we join PATH because the task might be executed on a different machine
        if '_runtime' in sos_dict and 'env' in sos_dict['_runtime']:
            for key, value in sos_dict['_runtime']['env'].items():
                if 'PATH' in key and key in os.environ:
                    os.environ[key] = value + os.pathsep + os.environ[key]
                else:
                    os.environ[key] = value

        env.sos_dict.quick_update(sos_dict)
        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        # define functions to_host, which uses CONFIG and _runtime and cannot be defined in runtime
        SoS_exec('''
from collections import Sequence, OrderedDict
def to_host(source):
    global _runtime
    global CONFIG
    path_map = OrderedDict()
    host = _runtime['on_host']

    if not host:
        return source

    if 'path_map' in _runtime:
        val = _runtime['path_map']
    elif 'path_map' in CONFIG['hosts'][host]:
        val = CONFIG['hosts'][host]['path_map']
    else:
        val = None
    if val is not None:
        if isinstance(val, str):
            val = [val]
        if isinstance(val, Sequence):
            for v in val:
                if ':' not in v or v.count(':') > 1:
                    raise ValueError('Path map should be separated as from:to, {} specified'.format(v))
                path_map[v.split(':')[0]] = v.split(':')[1]
        elif isinstance(val, dict):
            for k,v in val.items():
                path_map[k] = v
        else:
            raise ValueError('Unacceptable value for configuration path_map: {}'.format(val))

    def map_path(source):
        if os.path.isabs(source):
            dest = source
        else:
            dest = os.path.join(_runtime['cur_dir'], source)
        for k,v in path_map.items():
            if dest.startswith(k):
                dest = v + dest[len(k):]
        return dest

    if isinstance(source, str):
        return map_path(source)
    elif isinstance(source, Sequence):
        return [map_path(x) for x in source]
    else:
        raise ValueError('Unacceptable parameter {} to function to_host'.format(source))
''', None)
        # re-execute global definition because some of the definitions in the
        # global section might not be pickaleable (e.g. functions) and cannot
        # be passed to this separate process.
        if global_def:
            SoS_exec(global_def, global_sigil)
        # step process
        #if signature is None:
        #    env.sos_dict.set('__step_sig__', None)
        #else:
        #    env.sos_dict.set('__step_sig__', os.path.basename(signature.proc_info).split('.')[0])
        SoS_exec(task, sigil)
        os.chdir(orig_dir)
    except Exception as e:
        return {'succ': 1, 'exception': e, 'path': os.environ['PATH']}
    except KeyboardInterrupt:
        raise RuntimeError('KeyboardInterrupt from {}'.format(os.getpid()))
    finally:
        env.sos_dict.set('__step_sig__', None)

    #if signature is not None:
    #    signature.write(env.sos_dict['_local_input_{}'.format(env.sos_dict['_index'])],
    #        env.sos_dict['_local_output_{}'.format(env.sos_dict['_index'])])
    #    signature.release()
    env.deregister_process(os.getpid())
    return {'succ': 0, 'output': env.sos_dict['_output'], 'path': os.environ['PATH']}
