#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
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
import re
import subprocess
import tempfile
import pipes
from .utils import env

#
# A decoration function that allows SoS to replace all SoS actions
# with a null action.
#
def SoS_Action(run_mode='run'):
    run_mode = [run_mode] if isinstance(run_mode, str) else run_mode
    def runtime_decorator(func):
        def action_wrapper(*args, **kwargs):
            if env.run_mode not in run_mode:
                return 0
            else:
                return func(*args, **kwargs)
        return action_wrapper
    return runtime_decorator


class SoS_ExecuteScript:
    def __init__(self, script, interpreter, suffix):
        self.script = script
        self.interpreter = interpreter
        self.script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=suffix, delete=False).name
        with open(self.script_file, 'w') as script_file:
            script_file.write(self.script)

    def run(self):
        if '{}' in self.interpreter:
            cmd = self.interpreter.replace('{}', pipes.quote(self.script_file))
        else:
            cmd = self.interpreter + ' ' + pipes.quote(self.script_file)
        try:
            p = subprocess.Popen(cmd, shell=True)
            env.register_process(p.pid, 'Runing {}'.format(self.script_file))
            ret = p.wait()
        finally:
            env.deregister_process(p.pid)
        if ret != 0:
            raise RuntimeError('Failed to execute script')



class SoS_ExecuteRScript(SoS_ExecuteScript):
    '''SoS_Execute in-line R script using Rscript as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='Rscript',
            suffix='.R')

class SoS_ExecuteShellScript(SoS_ExecuteScript):
    '''SoS_Execute in-line shell script using bash as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='bash',
            suffix='.sh')

class SoS_ExecutePythonScript(SoS_ExecuteScript):
    '''SoS_Execute in-line python script using python as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='python',
            suffix='.py')


class SoS_ExecutePython3Script(SoS_ExecuteScript):
    '''SoS_Execute in-line python script using python as interpreter. Please
    check action SoS_ExecuteScript for more details.
    '''
    def __init__(self, script=''):
        SoS_ExecuteScript.__init__(self, script=script, interpreter='python3',
            suffix='.py')

@SoS_Action(run_mode='run')
def run(script):
    SoS_ExecuteShellScript(script).run()

@SoS_Action(run_mode='run')
def python(script):
    SoS_ExecutePythonScript(script).run()

@SoS_Action(run_mode='run')
def python3(script):
    SoS_ExecutePython3Script(script).run()

@SoS_Action(run_mode='run')
def R(script):
    SoS_ExecuteRScript(script).run()


try:
    from shutil import which
except ImportError:
    # this function is only define in Python 3.3 +
    def which(cmd, mode=os.F_OK | os.X_OK, path=None):
        """Given a command, mode, and a PATH string, return the path which
        conforms to the given mode on the PATH, or None if there is no such
        file.

        `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
        of os.environ.get("PATH"), or can be overridden with a custom search
        path.

        """
        # Check that a given file can be accessed with the correct mode.
        # Additionally check that `file` is not a directory, as on Windows
        # directories pass the os.access check.
        def _access_check(fn, mode):
            return (os.path.exists(fn) and os.access(fn, mode)
                    and not os.path.isdir(fn))

        # Short circuit. If we're given a full path which matches the mode
        # and it exists, we're done here.
        if _access_check(cmd, mode):
            return cmd

        path = (path or os.environ.get("PATH", os.defpath)).split(os.pathsep)
        files = [cmd]

        seen = set()
        for dir in path:
            dir = os.path.normcase(dir)
            if not dir in seen:
                seen.add(dir)
                for thefile in files:
                    name = os.path.join(dir, thefile)
                    if _access_check(name, mode):
                        return name
        return None


@SoS_Action(run_mode=['dryrun', 'run'])
def check_command(cmds):
    '''Check the existence of command `cmd` and raise an error if
    command does not exist. `cmd` can be one command or a list of
    commands.'''
    cmds = [cmds] if isinstance(cmds, str) else cmds
    #
    for cmd in cmds:
        name = which(cmd)
        if not name:
            raise RuntimeError('Command ``{}`` not found!'.format(cmd))
        env.logger.info('Command ``{}`` is located as ``{}``.'.format(cmd, name))
    return 0


@SoS_Action(run_mode=['dryrun', 'run'])
def fail_if(expr, msg=''):
    '''Raise an exception with `msg` if condition `expr` is False'''
    if expr:
        raise RuntimeError(msg)


@SoS_Action(run_mode=['dryrun', 'run'])
def warn_if(expr, msg=''):
    '''Yield an warning message `msg` if `expr` is False '''
    if expr:
        env.logger.warning(msg)

@SoS_Action(run_mode=['dryrun', 'run'])
def search_output(cmd, pattern):
    '''Raise an exception if output of `cmd` does not match specified `pattern`.
    Multiple patterns can be specified as a list of patterns.'''
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode()
    env.logger.trace('Output of command ``{}`` is ``{}``'.format(cmd, output))
    #
    pattern = [pattern] if isinstance(pattern, str) else pattern
    if all([re.search(x, output, re.MULTILINE) is None for x in pattern]):
        raise RuntimeError('Output of command ``{}`` does not match specified regular expression ``{}``.'
            .format(cmd, ' or '.join(pattern)))
