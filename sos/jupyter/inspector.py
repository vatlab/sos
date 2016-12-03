#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
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
import glob
from sos.utils import env
from sos.sos_syntax import SOS_USAGES

class SoS_VariableInspector(object):
    def __init__(self, kernel):
        self.kernel = kernel

    def inspect(self, name, line, pos):
        try:
            format_dict, _ = self.kernel.preview_var(name)
            return format_dict
        except Exception as e:
            return {}

def log(obj):
    with open(os.path.expanduser('~/a.txt'), 'a') as a:
        a.write('{}\n'.format(obj))

class SoS_SyntaxInspector(object):
    def __init__(self, kernel):
        self.kernel = kernel

    def inspect(self, name, line, pos):
        log(name)
        log(line)
        log(pos)
        log(self.kernel.ALL_MAGICS)
        if line.startswith('%') and name in self.kernel.ALL_MAGICS and pos <= len(name) + 1:
            if hasattr(self.kernel, 'get_{}_parser'.format(name)):
                log('had attr')
                parser = getattr(self.kernel, 'get_{}_parser'.format(name))()
                return {'text/plain': parser.format_help()}
            else:
                return {'text/plain': 'Magic %{}'.format(name) }
        elif line.startswith(name + ':') and pos <= len(name):
            # input: etc
            if name in SOS_USAGES:
                return {'text/plain': SOS_USAGES[name]}
            elif name in env.sos_dict and hasattr(env.sos_dict[name], '__doc__'):
                # action?
                return {'text/plain': env.sos_dict[name].__doc__}
            else:
                return {}
        else:
            return {}

class SoS_Inspector(object):
    def __init__(self, kernel):
        self.inspectors = [
            SoS_SyntaxInspector(kernel),
            SoS_VariableInspector(kernel),
        ]

    def inspect(self, name, line, pos):
        for c in self.inspectors:
            try:
                data = c.inspect(name, line, pos)
                if data:
                    return data
            except Exception as e:
                log(e)
                continue
        # No match
        return {}

