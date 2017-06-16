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

from sos.utils import short_repr, env
import json

JS_init_statement = '''
__get_sos_vars = function() {
    vars = []
    for(v in this) { 
        if (v.startsWith('sos') )
        { 
            vars.push(v);
        }
    }
    return vars;
}

'''

#
#  support for %get
#
#  Converting a Python object to a JSON format to be loaded by JavaScript
#
def _JS_repr(obj):
    try:
        # for JSON serlializable type, we simply dump it
        return json.dumps(obj)
    except:
        import numpy
        import pandas
        if isinstance(obj, (numpy.intc, numpy.intp, numpy.int8, numpy.int16, numpy.int32, numpy.int64,\
                numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.float16, numpy.float32, \
                numpy.float64, numpy.matrixlib.defmatrix.matrix, numpy.ndarray)):
            return json.dumps(obj.tolist())
        elif isinstance(obj, pandas.DataFrame):
            return obj.to_json(orient='index')
        elif isinstance(obj, set):
            return json.dumps(list(obj))
        else:
            return 'Unsupported seralizable data {} with type {}'.format(short_repr(obj), obj.__class__.__name__)


class sos_JavaScript:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'javascript'
        self.background_color = '#00ff80'
        self.init_statements = JS_init_statement

    def get_vars(self, names):
        for name in names:
            self.sos_kernel.run_cell('{} = {}'.format(name, _JS_repr(env.sos_dict[name])), True, False)

    def put_vars(self, items, to_kernel=None):
        # first let us get all variables with names starting with sos
        response = self.sos_kernel.get_response('__get_sos_vars()', ('execute_result'))[0][1]
        expr = response['data']['text/plain']
        items += eval(expr)

        if not items:
            return {}

        py_repr = 'JSON.stringify({{ {} }})'.format(','.join('"{0}":{0}'.format(x) for x in items))
        response = self.sos_kernel.get_response(py_repr, ('execute_result'))[0][1]
        expr = response['data']['text/plain']
        try:
            return json.loads(eval(expr))
        except Exception as e:
            self.sos_kernel.warn('Failed to convert {} to Python object: {}'.format(expr, e))
            return {}
