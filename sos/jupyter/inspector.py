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
import pydoc
from types import ModuleType
from sos.utils import env, short_repr
from sos.sos_eval import SoS_eval, get_default_global_sigil

class SoS_FileInspector(object):
    def __init__(self, kernel):
        self.kernel = kernel

    def inspect(self, name, line, pos):
        l = 0
        r = 10000000
        #
        for char in (' ', '\t', '"', "'", '=', '('):
            try:
                idx = line[pos:].index(char)
                if idx < r:
                    r = idx
            except:
                pass
            try:
                idx = line[:pos].rindex(char)
                if idx > l:
                    l = idx
            except:
                pass
        filename = line[l+1:pos + r]
        if os.path.isfile(os.path.expanduser(filename)):
        else:
            return {}

class SoS_VariableInspector(object):
    def __init__(self, kernel):
        self.kernel = kernel

    def inspect(self, name, line, pos):
        try:
            if name in env.sos_dict:
                obj = env.sos_dict[name]
            else:
                obj = SoS_eval(name, sigil=get_default_global_sigil())
            if callable(obj) or isinstance(obj, ModuleType):
                return {'text/plain': pydoc.getdoc(obj)}
            else:
                format_dict, md_dict = self.kernel.format_obj(obj)
                return format_dict
        except Exception as e:
            log(e)
            return {}
    
def log(obj):
    with open(os.path.expanduser('~/a.txt'), 'a') as a:
        a.write('{}\n'.format(obj))

class SoS_Inspector(object):
    def __init__(self, kernel):
        self.inspectors = [
            SoS_VariableInspector(kernel),
            SoS_FileInspector(kernel),
        ]

    def inspect(self, name, line, pos):
        for c in self.inspectors:
            try:
                data = c.inspect(name, line, pos)
                log(name)
                if data:
                    log(data)
                    return data
            except Exception as e:
                raise
        # No match
        return {}

