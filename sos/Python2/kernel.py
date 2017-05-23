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
import sys
from collections import OrderedDict

__init_statement__ = '''
def __version_info__(module):
    # return the version of Python module
    try:
        code = ("import %s; version=str(%s.__version__)" %
                (module, module))
        ns_g = ns_l = {}
        exec(compile(code, "<string>", "exec"), ns_g, ns_l)
        return ns_l["version"]
    except Exception as e:
        import pkg_resources
        try:
            return pkg_resources.require(module)[0].version
        except Exception as e:
            return 'na'

def __loaded_modules__():
    from types import ModuleType
    res = {}
    for key,value in globals().items():
        if isinstance(value, ModuleType):
            res[value.__name__] = __version_info__(value.__name__)
    return {x:y for x,y in res.items() if y != 'na'}
'''

class sos_Python2:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'python2'
        self.background_color = '#F6FAEA'
        self.init_statements = __init_statement__

    def sos_to_lan(self, name, obj):
        # try to dump data in a python 2 compatible fashion so that python 2 can load it
        if self.sos_kernel._debug_mode:
            self.sos_kernel.warn("RUN import pickle\nglobals().update(pickle.loads({!r}))".format(pickle.dumps({name:obj}, protocol=2, fix_imports=True)))
        return name, "import pickle\nglobals().update(pickle.loads({!r}))".format(pickle.dumps({name:obj}, protocol=2, fix_imports=True))

    def lan_to_sos(self, items):
        # python 2 uses protocol 2, which python 3 should be able to pick up
        stmt = 'import pickle\n__vars__={{ {} }}\n__vars__.update({{x:y for x,y in locals().items() if x.startswith("sos")}})\npickle.dumps(__vars__)'.format(','.join('"{0}":{0}'.format(x) for x in items))
        response = self.sos_kernel.get_response(stmt, ['execute_result'])[0][1]
        try:
            ret = pickle.loads(eval(response['data']['text/plain']).encode('utf-8'))
            if self.sos_kernel._debug_mode:
                self.sos_kernel.warn('Get: {}'.format(ret))
            return ret
        except Exception as e:
            self.sos_kernel.warn('Failed to import variables {}: {}'.format(items, e))
            return {}

    def sessioninfo(self):
        res = OrderedDict()
        res['Version'] = sys.version
        modules = self.sos_kernel.get_response('__loaded_modules__()', ['execute_result'])[0][1]
        res.update(eval(modules['data']['text/plain']))
        return res
