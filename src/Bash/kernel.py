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

from collections.abc import Sequence
from sos.utils import short_repr, env

#
#  support for %get
#
#  Converting a Python object to a R expression that will be executed
#  by the R kernel.
#
#
def _Bash_repr(obj):
    if isinstance(obj, bool):
        return 'TRUE' if obj else 'FALSE'
    elif isinstance(obj, (int, float)):
        return repr(obj)
    elif isinstance(obj, str):
        return obj
    elif isinstance(obj, Sequence):
        if len(obj) == 0:
            return ''
        return ' '.join(_Bash_repr(x) for x in obj)
    elif obj is None:
        return ''
    elif isinstance(obj, dict):
        return ' '.join(_Bash_repr(x) for x in obj.keys())
    elif isinstance(obj, set):
        return ' '.join(_Bash_repr(x) for x in obj)
    else:
        return repr('Unsupported datatype {}'.format(short_repr(obj)))

class sos_Bash:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'bash'
        self.background_color = '#E6EEFF'
        self.init_statements = ''

    def get_vars(self, names):
        for name in names:
            stmt = 'export {}={!r}'.format(name, _Bash_repr(env.sos_dict[name]))
            self.sos_kernel.run_cell(stmt, True, False, on_error='Failed to get variable {}'.format(name))

    def put_vars(self, items, to_kernel=None):
        # first let us get all variables with names starting with sos
        response = self.sos_kernel.get_response('env', ('stream'))
        if self.sos_kernel._debug_mode:
            self.sos_kernel.warn('Response {}'.format(response))
        response = [x[1]['text'].split('=', 1) for x in response]
        all_vars = {x:y.strip() for x,y in response if x.startswith('sos') or x in items}

        for item in items:
            if item not in all_vars:
                self.sos_kernel.warn('Variable not exist: {}'.format(item))

        return all_vars

        

