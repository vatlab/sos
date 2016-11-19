#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)

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

import pkg_resources
from .sos_script import SoS_Script
from .utils import logger, get_output, sos_handle_parameter_
from .actions import SoS_Action, execute_script, sos_run, \
    fail_if, warn_if, stop_if, download, run, bash, csh, tcsh, zsh, sh, \
    python, python3, perl, ruby, node, JavaScript, \
    docker_build, docker_commit, report, pandoc
from .sos_eval import interpolate, sos_namespace_
from .pattern import expand_pattern
from .target import dynamic, executable, sos_variable, env_variable
from .main import runfile

# silent pyflakes
SoS_Script
logger, get_output, sos_handle_parameter_
SoS_Action, execute_script, sos_run
fail_if, warn_if, stop_if, download, run, bash, csh, tcsh, zsh, sh
python, python3, perl, ruby, node, JavaScript
docker_build, docker_commit, report, pandoc
interpolate, sos_namespace_
expand_pattern
dynamic, executable, sos_variable, env_variable
runfile


for entrypoint in pkg_resources.iter_entry_points(group='sos_targets'):
    # Grab the function that is the actual plugin.
    name = entrypoint.name
    try:
        plugin = entrypoint.load()
        globals()[name] = plugin
    except Exception as e:
        print('Failed to load target {}: {}'.format(entrypoint.name, e))

for entrypoint in pkg_resources.iter_entry_points(group='sos_actions'):
    name = entrypoint.name
    try:
        plugin = entrypoint.load()
        globals()[name] = plugin
    except Exception as e:
        print('Failed to load action {}: {}'.format(entrypoint.name, e))
