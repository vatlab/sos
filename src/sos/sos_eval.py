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
import os
import sys
import re
import collections
import ast

from .utils import env, Error, short_repr, DelayedAction

def interpolate(text, local_dict=None, global_dict=None):
    '''Evaluate expressions in `text` '''
    return text

def cfg_interpolate(text, local_dict={}):
    return interpolate(text, '${ }', local_dict, env.sos_dict.get('CONFIG', {}))

def accessed_vars(statement):
    '''Parse a Python statement and analyze the symbols used. The result
    will be used to determine what variables a step depends upon.'''
    return {node.id for node in ast.walk(ast.parse(statement)) if isinstance(node, ast.Name)}

def SoS_eval(expr):
    '''Evaluate an expression with sos dict.'''
    return eval(expr, env.sos_dict._dict)

def _is_expr(expr):
    try:
        compile(expr, '<string>', 'eval')
        return True
    except Exception:
        return False

def SoS_exec(stmts, _dict=None):
    '''Execute a statement after modifying (convert ' ' string to raw string,
    interpolate expressions) strings.'''
    # the trouble here is that we have to execute the statements line by line
    # because the variables defined. The trouble is in cases such as class
    # definition
    #
    # class A:
    #     def __init__(self):
    #         pass
    # # this is already correct syntax but the remaining piece is not.
    #     def another_one(self):
    #         pass
    #
    # We therefore has to be a bit more clever on this.
    #
    # we first group all lines by their own group
    # we then try to find syntaxly valid groups
    code_group = [x for x in stmts.split('\n')]
    idx = 0
    if _dict is None:
        _dict = env.sos_dict._dict
    while True:
        try:
            # test current group
            compile(code_group[idx], filename = '<string>', mode='exec')
            # if it is ok, go next
            idx += 1
            if idx == len(code_group):
                break
        except Exception:
            # error happens merge the next line
            if idx < len(code_group) - 1:
                code_group[idx] += '\n' + code_group[idx + 1]
                code_group.pop(idx + 1)
            else:
                # if no next group, expand previously correct one
                if idx == 0:
                    raise RuntimeError('Failed to find syntax correct group: {}'.format(stmts))
                # break myself again
                code_group = code_group[: idx] + code_group[idx].split('\n') + code_group[idx+1:]
                # go back
                idx -= 1
                code_group[idx] += '\n' + code_group[idx + 1]
                code_group.pop(idx+1)
    #
    # execute statements one by one
    res = None
    executed = ''
    code_group = [x for x in code_group if x.strip()]
    for idx,code in enumerate(code_group):
        stmts = code
        if not stmts.strip():
            continue
        try:
            if idx + 1 == len(code_group) and _is_expr(stmts):
                res = eval(stmts, _dict)
            else:
                exec(stmts, _dict)
        except Exception as e:
            env.logger.debug('Failed to execute {}'.format(short_repr(stmts)))
            raise
        #finally:
        #    del act
        executed += stmts + '\n'
    env.sos_dict.check_readonly_vars()
    return res

#
# dynamic expression that cannot be resolved during parsing
# at prepare mode etc, and has to be resolved at run time.
#
class Undetermined(object):
    def __init__(self, expr=''):
        if not isinstance(expr, str):
            raise RuntimeError('Undetermined expression has to be a string: "{}" passed'.format(expr))
        self.expr = expr.strip()

    def value(self):
        return SoS_eval(self.expr)

    def __repr__(self):
        return 'Undetermined({!r})'.format(self.expr)

    def __hash__(self):
        raise RuntimeError('Undetermined expression should be evaluated before used. '
            'This is certainly a bug so please report this to SoS developer.')

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
            raise ValueError('Failed to evaluate option {} with value {}: {}'
                .format(key, self._expressions[key], e))

    def __repr__(self):
        return repr(self._expressions)
