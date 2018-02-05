#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
from .utils import logger, get_output, sos_handle_parameter_
from .eval import interpolate, sos_namespace_
from .pattern import expand_pattern
from .targets import path, paths

# silent pyflakes
logger, get_output, sos_handle_parameter_
interpolate, sos_namespace_
expand_pattern, path, paths


sos_symbols_ = {
    'logger', 'get_output', 'sos_handle_parameter_'
    'interpolate', 'sos_namespace_',
    'expand_pattern', 'runfile'
}

def _load_group(group):
    global sos_symbols_
    for _entrypoint in pkg_resources.iter_entry_points(group=group):
        # import all targets and actions from entry_points
        # Grab the function that is the actual plugin.
        _name = _entrypoint.name
        sos_symbols_.add(_name)
        try:
            _plugin = _entrypoint.load()
            globals()[_name] = _plugin
        except Exception as e:
            # look for sos version requirement
            if 'Requirement.parse' in str(e):
                import re
                from ._version import __version__
                from pkg_resources import parse_version
                m = re.search("Requirement.parse\('sos>=([^)]*)'\)", str(e))
                if m:
                    if parse_version(__version__) < parse_version(m.group(1)):
                        logger.warning(
                            f'Failed to load target {_entrypoint.name}: please upgrade your version of sos from {__version__} to at least version {m.group(1)}')
                        continue
            if _name == 'run':
                # this is critical so we print the warning
                logger.warning(f'Failed to load target {_entrypoint.name}: {e}')
            else:
                logger.trace(f'Failed to load target {_entrypoint.name}: {e}')

_load_group('sos_targets')
_load_group('sos_actions')
_load_group('sos_functions')

