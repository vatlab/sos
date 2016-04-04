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
import shlex
from shutil import which
from .utils import env, interpolate, glob_wildcards

__all__ = ['SoS_Action', 'SoS_ExecuteScript',
    'check_command', 'fail_if', 'warn_if',
    'run', 'bash', 'csh', 'tcsh', 'zsh', 'sh',
    'python', 'python3',
    'perl', 'ruby', 'node', 'JavaScript',
    'R', 'check_R_library',
    ]

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
        action_wrapper.run_mode = run_mode
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

@SoS_Action(run_mode=['run'])
def execute_script(script, interpreter, suffix):
    return SoS_ExecuteScript(script, interpreter, suffix).run()

@SoS_Action(run_mode=['dryrun', 'run'])
def check_command(cmd, pattern = None):
    '''Raise an exception if output of `cmd` does not match specified `pattern`.
    Multiple patterns can be specified as a list of patterns.
    When pattern is None, check the existence of command `cmd`
    and raise an error if command does not exist.'''
    ret_val = 0
    if pattern is None and len(shlex.split(cmd)) == 1:
        name = which(cmd)
        if not name:
            raise RuntimeError('Command ``{}`` not found!'.format(cmd))
        env.logger.info('Command ``{}`` is located as ``{}``.'.format(cmd, name))
    else:
        try:
            output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode()
        except subprocess.CalledProcessError as e:
            ret_val = e.returncode
            output = e.output
            env.logger.warning(e)
        #
        env.logger.trace('Output of command ``{}`` is ``{}``'.format(cmd, output))
        #
        if pattern:
            pattern = [pattern] if isinstance(pattern, str) else pattern
            if all([re.search(x, output, re.MULTILINE) is None for x in pattern]):
                raise RuntimeError('Output of command ``{}`` does not match specified regular expression ``{}``.'
                    .format(cmd, ' or '.join(pattern)))
    return ret_val

@SoS_Action(run_mode=['dryrun', 'run'])
def fail_if(expr, msg=''):
    '''Raise an exception with `msg` if condition `expr` is False'''
    if expr:
        raise RuntimeError(msg)
    return 0

@SoS_Action(run_mode=['dryrun', 'run'])
def warn_if(expr, msg=''):
    '''Yield an warning message `msg` if `expr` is False '''
    if expr:
        env.logger.warning(msg)
    return 0

@SoS_Action(run_mode='run')
def run(script):
    return SoS_ExecuteScript(script, 'bash', '.sh').run()

@SoS_Action(run_mode='run')
def bash(script):
    return SoS_ExecuteScript(script, 'bash', '.sh').run()

@SoS_Action(run_mode='run')
def csh(script):
    return SoS_ExecuteScript(script, 'csh', '.csh').run()

@SoS_Action(run_mode='run')
def tcsh(script):
    return SoS_ExecuteScript(script, 'tcsh', '.sh').run()

@SoS_Action(run_mode='run')
def zsh(script):
    return SoS_ExecuteScript(script, 'zsh', '.zsh').run()

@SoS_Action(run_mode='run')
def sh(script):
    return SoS_ExecuteScript(script, 'sh', '.sh').run()

@SoS_Action(run_mode='run')
def python(script):
    return SoS_ExecuteScript(script, 'python', '.py').run()

@SoS_Action(run_mode='run')
def python3(script):
    return SoS_ExecuteScript(script, 'python3', '.py').run()

@SoS_Action(run_mode='run')
def perl(script):
    return SoS_ExecuteScript(script, 'perl', '.pl').run()

@SoS_Action(run_mode='run')
def ruby(script):
    return SoS_ExecuteScript(script, 'ruby', '.rb').run()

@SoS_Action(run_mode='run')
def node(script):
    return SoS_ExecuteScript(script, 'node', '.js').run()

@SoS_Action(run_mode='run')
def JavaScript(script):
    return SoS_ExecuteScript(script, 'node', '.js').run()

@SoS_Action(run_mode='run')
def R(script):
    return SoS_ExecuteScript(script, 'Rscript --default-packages=methods,utils,stats', '.R').run()

@SoS_Action(run_mode=['dryrun', 'run'])
def check_R_library(name, version = None):
    '''Check existence and version match of R library.
    cran and bioc packages are unique yet might overlap with github.
    Therefore if the input name is {repo}/{pkg} the package will be
    installed from github if not available, else from cran or bioc
    '''
    output_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False).name
    if len(glob_wildcards('{repo}/{pkg}', [name])['repo']):
        # package is from github
        check_R_library('devtools')
        install_script = interpolate('''
        options(warn=-1)
        package_repo <- ${name!r}
        package <- basename(package_repo)
        if (require(package, character.only=TRUE, quietly=TRUE)) {
            write(paste(package, packageVersion(package), "AVAILABLE"), file="${output_file}")
        } else {
            devtools::install_github(package_repo)
            # if it still does not exist, write the package name to output
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "INSTALLED"), file="${output_file}")
            } else {
                write(paste(package, "NA", "MISSING"), file="${output_file}")
                quit("no")
            }
        }
        cur_version <- packageVersion(package)
        ''', '${ }', locals())
    else:
        # package is from cran or bioc
        install_script = interpolate('''
        options(warn=-1)
        package <- ${name!r}
        if (require(package, character.only=TRUE, quietly=TRUE)) {
            write(paste(package, packageVersion(package), "AVAILABLE"), file="${output_file}")
        } else {
            install.packages(package, repos="http://cran.us.r-project.org",
                quiet=FALSE)
            # if the package still does not exist
            if (!require(package, character.only=TRUE, quietly=TRUE)) {
                source("http://bioconductor.org/biocLite.R")
                biocLite(package, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
            }
            # if it still does not exist, write the package name to output
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "INSTALLED"), file="${output_file}")
            } else {
                write(paste(package, "NA", "MISSING"), file="${output_file}")
                quit("no")
            }
        }
        cur_version <- packageVersion(package)
        ''', '${ }', locals())
    version_script = ''
    if version is not None:
        version = [version] if isinstance(version, str) else version
        operators = []
        for idx, value in enumerate(version):
            value = str(value)
            if value.endswith('+'):
                operators.append('>=')
                version[idx] = value[:-1]
            elif value.endswith('-'):
                operators.append('<')
                version[idx] = value[:-1]
            else:
                operators.append('==')
        # check version and mark version mismatch
        # if current version satisfies any of the
        # requirement the check program quits
        for x, y in zip(version, operators):
            version_script += '''
            if (cur_version {1} {0}) {{
              quit("no")
            }}
            '''.format(repr(x), y)
        version_script += 'write(paste(package, cur_version, "VERSION_MISMATCH"), file = {})'.\
          format(repr(output_file))
    SoS_ExecuteScript(install_script + version_script, 'Rscript --default-packages=methods,utils,stats', '.R').run()
    ret_val = 0
    with open(output_file) as tmp:
        for line in tmp:
            lib, version, status = line.split()
            if status.strip() == "MISSING":
                raise RuntimeError('R Library {} is not available and cannot be installed.'.format(lib))
            elif status.strip() == 'AVAILABLE':
                env.logger.info('R library {} ({}) is available'.format(lib, version))
            elif status.strip() == 'INSTALLED':
                env.logger.info('R library {} ({}) has been installed'.format(lib, version))
            elif status.strip() == 'VERSION_MISMATCH':
                env.logger.warning('R library {} ({}) does not satisfy version requirement!'.format(lib, version))
                ret_val = 1
            else:
                raise RuntimeError('This should not happen: {}'.format(line))
    try:
        os.remove(self.output_file)
    except:
        pass
    return ret_val
