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
from sos.utils import env, short_repr

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
    res = []
    for key,value in globals().items():
        if isinstance(value, ModuleType):
            res.append([value.__name__, __version_info__(value.__name__)])
    return [(x,y) for x,y in res if y != 'na']
'''


class sos_Python2:

    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'python2'
        self.background_color = '#F6FAEA'
        self.init_statements = __init_statement__

    def get_vars(self, names):
        self.sos_kernel.run_cell("import pickle", True, False)
        for name in names:
            stmt = "globals().update(pickle.loads({!r}))\n".format(pickle.dumps({name: env.sos_dict[name]}, protocol=2, fix_imports=True))
            self.sos_kernel.run_cell(stmt, True, False, on_error='Failed to get variable {} from SoS to Python2'.format(name))

    def load_pickled(self, item):
        if isinstance(item, bytes):
            return pickle.loads(item)
        elif isinstance(item, str):
            return pickle.loads(item.encode('utf-8'))
        else:
            self.warn('Cannot restore from result of pickle.dumps: {}'.format(short_repr(item)))
            return {}

    def put_vars(self, items, to_kernel=None):
        # python 2 uses protocol 2, which python 3 should be able to pick up
        stmt = 'import pickle\n__vars__={{ {} }}\n__vars__.update({{x:y for x,y in locals().items() if x.startswith("sos")}})\npickle.dumps(__vars__)'.format(','.join('"{0}":{0}'.format(x) for x in items))
        response = self.sos_kernel.get_response(stmt, ['execute_result'])[0][1]
        if to_kernel == 'Python2':
            # to self, this should allow all variables to be passed
            return 'import pickle\nglobals().update(pickle.loads({}))'.format(response['data']['text/plain'])
        else:
            try:
                ret = self.load_pickled(eval(response['data']['text/plain']))
                if self.sos_kernel._debug_mode:
                    self.sos_kernel.warn('Get: {}'.format(ret))
                return ret
            except Exception as e:
                self.sos_kernel.warn('Failed to import variables {}: {}'.format(items, e))
                return {}

    def sessioninfo(self):
        modules = self.sos_kernel.get_response('import pickle;import sys;res=[("Version", sys.version)];res.extend(__loaded_modules__());pickle.dumps(res)', ['execute_result'])[0][1]
        return self.load_pickled(eval(modules['data']['text/plain']))
