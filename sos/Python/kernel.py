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

class sos_Python:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'python3'
        self.init_statements = ''

    def sos_to_lan(self, name, obj):
        return name, "import pickle\nglobals().update(pickle.loads({!r}))".format(pickle.dumps({name:obj}))

    def lan_to_sos(self, items):
        default_items = [x for x in env.sos_dict.keys() if x.startswith('sos') and x not in self.sos_kernel.original_keys]
        all_items = default_items + items
        if not all_items:
            return {}
        response = self.sos_kernel.get_response('import pickle\npickle.dumps({{ {} }})'.format(','.join('"{0}":{0}'.format(x) for x in all_items)),
            ['execute_result'])
        return pickle.loads(eval(response['data']['text/plain']))
