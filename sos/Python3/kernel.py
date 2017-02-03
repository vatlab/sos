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

import pickle

class sos_Python3:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'python3'
        self.background_color = '#EAFAF1'
        self.init_statements = ''

    def sos_to_lan(self, name, obj):
        return name, "import pickle\nglobals().update(pickle.loads({!r}))".format(pickle.dumps({name:obj}))

    def lan_to_sos(self, items):
        stmt = 'import pickle\n__vars__={{ {} }}\n__vars__.update({{x:y for x,y in locals().items() if x.startswith("sos")}})\npickle.dumps(__vars__)'.format(','.join('"{0}":{0}'.format(x) for x in items))
        response = self.sos_kernel.get_response(stmt, ['execute_result'])[0][1]
        try:
            ret = pickle.loads(eval(response['data']['text/plain']))
            return ret
        except Exception as e:
            self.sos_kernel.warn('Failed to import variables {}: {}'.format(items, e))
            return {}
