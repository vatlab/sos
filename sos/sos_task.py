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
from io import StringIO
from tokenize import generate_tokens

from sos.utils import env
from sos.sos_eval import SoS_exec

from .target import textMD5, RuntimeInfo
from .monitor import ProcessMonitor

from collections import OrderedDict

class TaskParams(object):
    '''A parameter object that encaptulates parameters sending to
    task executors. This would makes the output of workers, especially
    in the web interface much cleaner (issue #259)'''
    def __init__(self, name, data):
        self.name = name
        self.data = data

    def __repr__(self):
        return self.name

def execute_task(task_id, verbosity=None, sigmode=None, monitor_interval=5):
    '''A function that execute specified task within a local dictionary
    (from SoS env.sos_dict). This function should be self-contained in that
    it can be handled by a task manager, be executed locally in a separate
    process or remotely on a different machine.'''
    env.logger.info('Executing task {}'.format(task_id))
    # start a monitoring file, which would be killed after the job
    # is done (killed etc)
    m = ProcessMonitor(task_id, interval=monitor_interval)
    m.start()

    task_file = os.path.join(os.path.expanduser('~'), '.sos', 'tasks', task_id + '.task')
    with open(task_file, 'rb') as task:
        params = pickle.load(task)

    task, sos_dict, sigil = params.data
    if verbosity is not None:
        env.verbosity = verbosity
    if sigmode is not None:
        env.sig_mode = sigmode
    env.register_process(os.getpid(), 'spawned_job with {} {}'
        .format(sos_dict['_input'], sos_dict['_output']))

    env.sos_dict.quick_update(sos_dict)

    skipped = False
    if env.sig_mode == 'ignore' or env.sos_dict['_output'] is None:
        sig = None
    else:
        tokens = [x[1] for x in generate_tokens(StringIO(task).readline)]
        # try to add #task so that the signature can be different from the step
        # if everything else is the same
        sig = RuntimeInfo(textMD5('#task\n' + ' '.join(tokens)), task,
            env.sos_dict['_input'], env.sos_dict['_output'], env.sos_dict['_depends'], env.sos_dict['__signature_vars__'])
        sig.lock()

        idx = env.sos_dict['_index']
        if env.sig_mode == 'default':
            matched = sig.validate()
            if isinstance(matched, dict):
                # in this case, an Undetermined output can get real output files
                # from a signature
                env.sos_dict.set('_input', matched['input'])
                env.sos_dict.set('_depends', matched['depends'])
                env.sos_dict.set('_output', matched['output'])
                env.sos_dict.set('_local_input', matched['local_output'])
                env.sos_dict.set('_local_output', matched['local_output'])
                env.sos_dict['local_input'].extend(env.sos_dict['_local_input'])
                env.sos_dict['local_output'].extend(env.sos_dict['_local_output'])
                env.sos_dict.update(matched['vars'])
                env.logger.info('Task ``{}`` (index={}) is ``ignored`` due to saved signature'.format(env.sos_dict['step_name'], idx))
                skipped = True
        elif env.sig_mode == 'assert':
            matched = sig.validate()
            if isinstance(matched, str):
                raise RuntimeError('Signature mismatch: {}'.format(matched))
            else:
                env.sos_dict.set('_input', matched['input'])
                env.sos_dict.set('_depends', matched['depends'])
                env.sos_dict.set('_output', matched['output'])
                env.sos_dict.set('_local_input', matched['local_output'])
                env.sos_dict.set('_local_output', matched['local_output'])
                env.sos_dict['local_input'].extend(env.sos_dict['_local_input'])
                env.sos_dict['local_output'].extend(env.sos_dict['_local_output'])
                env.sos_dict.update(matched['vars'])
                env.logger.info('Step ``{}`` (index={}) is ``ignored`` with matching signature'.format(env.sos_dict['step_name'], idx))
                skipped = True
        elif env.sig_mode == 'build':
            # build signature require existence of files
            if sig.write(
                env.sos_dict['_local_input_{}'.format(idx)],
                env.sos_dict['_local_output_{}'.format(idx)],
                rebuild=True):
                env.logger.info('Task ``{}`` (index={}) is ``ignored`` with signature constructed'.format(env.sos_dict['step_name'], idx))
                skipped = True
        elif env.sig_mode == 'force':
            skipped = False
        else:
            raise RuntimeError('Unrecognized signature mode {}'.format(env.sig_mode))

    if skipped:
        return {'succ': 0, 'output': env.sos_dict['_output'], 'path': os.environ['PATH']}

    try:
        # go to 'cur_dir'
        orig_dir = os.getcwd()
        if '_runtime' in sos_dict and 'cur_dir' in sos_dict['_runtime']:
            if not os.path.isdir(os.path.expanduser(sos_dict['_runtime']['cur_dir'])):
                try:
                    os.makedirs(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
                except Exception as e:
                    raise RuntimeError('Failed to create cur_dir {}'.format(sos_dict['_runtime']['cur_dir']))
            os.chdir(os.path.expanduser(sos_dict['_runtime']['cur_dir']))
        # go to user specified workdir
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
                    new_path = OrderedDict()
                    for p in value.split(os.pathsep):
                        new_path[p] = 1
                    for p in value.split(os.environ[key]):
                        new_path[p] = 1
                    os.environ[key] = os.pathsep.join(new_path.keys())
                else:
                    os.environ[key] = value

        SoS_exec('import os, sys, glob', None)
        SoS_exec('from sos.runtime import *', None)
        # step process
        SoS_exec(task, sigil)
        os.chdir(orig_dir)
    except Exception as e:
        env.logger.error('Task {} terminated with error: {}'.format(task_id, e))
        return {'succ': 1, 'exception': e, 'path': os.environ['PATH']}
    except KeyboardInterrupt:
        raise RuntimeError('KeyboardInterrupt from {}'.format(os.getpid()))
    finally:
        env.sos_dict.set('__step_sig__', None)

    if sig:
        sig.write(env.sos_dict['_local_input_{}'.format(env.sos_dict['_index'])],
            env.sos_dict['_local_output_{}'.format(env.sos_dict['_index'])])
        sig.release()
    env.deregister_process(os.getpid())
    return {'succ': 0, 'output': env.sos_dict['_output'], 'path': os.environ['PATH']}
