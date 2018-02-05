#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
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
import ast
from .utils import env,  text_repr

def interpolate(text, global_dict=None, local_dict=None):
    '''Evaluate expressions in `text` '''
    # step 1, make it a f-string (add quotation marks and f
    # step 2, evaluate as a string
    try:
        return eval('f' + text_repr(text), global_dict, local_dict)
    except Exception as e:
        raise ValueError(f'Failed to interpolate {text}: {e}')

def cfg_interpolate(text, local_dict={}):
    # handle nested interpolate ...
    while True:
        res = interpolate(text, local_dict, env.sos_dict.get('CONFIG', {}))
        if res == text:
            break
        else:
            text = res
    return res

def accessed_vars(statement, filename='<string>', mode='exec'):
    '''Parse a Python statement and analyze the symbols used. The result
    will be used to determine what variables a step depends upon.'''
    try:
        return {node.id for node in ast.walk(ast.parse(statement, filename, mode)) if isinstance(node, ast.Name)}
    except:
        # try to treat them as parameters
        try:
            return {node.id for node in ast.walk(ast.parse('__NULLFUNC__(' + statement + ')', filename, mode)) if isinstance(node, ast.Name)}
        except:
            raise RuntimeError(f'Failed to parse statement: {statement}')

def SoS_eval(expr):
    '''Evaluate an expression with sos dict.'''
    return eval(expr, env.sos_dict._dict)

def _is_expr(expr):
    try:
        compile(expr, '<string>', 'eval')
        return True
    except Exception:
        return False

class StatementHash(object):
    stmt_hash = {}
    def __init__(self):
        pass

    def hash(self, script):
        h = hash(script)
        StatementHash.stmt_hash[h] = script
        return f'script_{h}'

    def script(self, hash):
        return StatementHash.stmt_hash[int(hash[7:])]

stmtHash = StatementHash()

def SoS_exec(script, _dict=None, return_result=True):
    '''Execute a statement.'''
    if _dict is None:
        _dict = env.sos_dict._dict

    if not return_result:
        exec(compile(script, filename=stmtHash.hash(script), mode='exec'), _dict)
        return None

    try:
        stmts = list(ast.iter_child_nodes(ast.parse(script)))
        if not stmts:
            return
        if isinstance(stmts[-1], ast.Expr):
            # the last one is an expression and we will try to return the results
            # so we first execute the previous statements
            if len(stmts) > 1:
                exec(compile(ast.Module(body=stmts[:-1]), filename=stmtHash.hash(script), mode="exec"), _dict)
            # then we eval the last one
            res = eval(compile(ast.Expression(body=stmts[-1].value), filename=stmtHash.hash(script), mode="eval"), _dict)
        else:
            # otherwise we just execute the entire code
            exec(compile(script, filename=stmtHash.hash(script), mode='exec'), _dict)
            res = None
    except SyntaxError as e:
        raise SyntaxError(f"Invalid code {script}: {e}")

    #if check_readonly:
    #    env.sos_dict.check_readonly_vars()
    return res

#
# dynamic expression that cannot be resolved during parsing
# at prepare mode etc, and has to be resolved at run time.
#
class Undetermined(object):
    def __init__(self, expr=''):
        if not isinstance(expr, str):
            raise RuntimeError(f'Undetermined expression has to be a string: "{expr}" passed')
        self.expr = expr.strip()

    def value(self):
        return SoS_eval(self.expr)

    def __repr__(self):
        return f'Undetermined({self.expr!r})'

    def __hash__(self):
        raise RuntimeError('Undetermined expression should be evaluated before used. '
            'This is certainly a bug so please report this to SoS developer.')

    def targets(self):
        return self

class sos_namespace_(object):
    '''A namespace that is created by evaluating statements
    and use the results as attributes of the object.'''
    def __init__(self, stmts):
        # we need to define functions defined by sos ...
        exec('from sos.runtime import *', self.__dict__)
        # the results of the statments will be saved as
        # attribute of this object.
        SoS_exec(stmts, _dict=self.__dict__)

class on_demand_options(object):
    '''Expression that will be evaluated upon request.'''
    def __init__(self, items):
        self._expressions = {}
        if items:
            self._expressions.update(items)

    def set(self, key, value):
        self._expressions[key] = repr(value)

    def __contains__(self, key):
        return key in self._expressions

    def __setitem__(self, key, value):
        self._expressions[key] = value

    def __getitem__(self, key):
        # first check if the value if cached
        if key not in self._expressions:
            raise KeyError(key)
        try:
            return SoS_eval(self._expressions[key])
        except Exception as e:
            if key == 'skip':
                raise ValueError(f'Failed to evaluate option {key} with value {self._expressions[key]}: Only constant values are allowed for section option skip')
            else:
                raise ValueError(f'Failed to evaluate option {key} with value {self._expressions[key]}: {e}')

    def __repr__(self):
        return repr(self._expressions)
