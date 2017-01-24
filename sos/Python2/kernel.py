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

import pickle
from sos.utils import env

class sos_Python2:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'python2'
        self.background_color = '#F6FAEA'
        self.init_statements = ''

    def sos_to_lan(self, name, obj):
        # try to dump data in a python 2 compatible fashion so that python 2 can load it
        if self.sos_kernel._debug_mode:
            self.sos_kernel.warn("RUN import pickle\nglobals().update(pickle.loads({!r}))".format(pickle.dumps({name:obj}, protocol=2, fix_imports=True)))
        return name, "import pickle\nglobals().update(pickle.loads({!r}))".format(pickle.dumps({name:obj}, protocol=2, fix_imports=True))

    def lan_to_sos(self, items):
        # python 2 uses protocol 2, which python 3 should be able to pick up
        stmt = 'import pickle\n__vars__={{ {} }}\n__vars__.update({{x:y for x,y in locals().items() if x.startswith("sos")}})\npickle.dumps(__vars__)'.format(','.join('"{0}":{0}'.format(x) for x in items))
        response = self.sos_kernel.get_response(stmt, ['execute_result'])
        try:
            ret = pickle.loads(eval(response['data']['text/plain']).encode('utf-8'))
            if self.sos_kernel._debug_mode:
                self.sos_kernel.warn('Get: {}'.format(ret))
            return ret
        except Exception as e:
            self.sos_kernel.warn('Failed to import variables {}: {}'.format(items, e))
            return {}
