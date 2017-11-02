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
from io import StringIO
from shlex import quote
from subprocess import list2cmdline
from tokenize import generate_tokens, untokenize, NAME, STRING, INDENT

from .utils import env, Error, short_repr, DelayedAction
from .sos_syntax import FORMAT_SPECIFIER, SIMPLE_SUB

# function interpolate is needed because it is required by the SoS
# script (not seen but translated to have this function)
__all__ = ['interpolate']


class InterpolationError(Error):
    '''Exception raised for interpolation related errors'''
    def __init__(self, text, msg):
        if '\n' in text:
            # if there is newline in original text, output it in separate
            # lines
            msg = '{}:\n{}\n'.format(msg, text)
        else:
            msg = '{}: {}'.format(msg, text)
        Error.__init__(self, msg)
        self.args = (text, msg)

class UnresolvableObject(Error):
    def __init__(self, obj):
        super(UnresolvableObject, self).__init__('Unresolvable object ({}) during string interpolation, use converter R if needed.'.format(obj))

def sos_compile(expr, *args, **kwargs):
    ''' Compiling an statement but ignores tab error (mixed tab and space)'''
    try:
        compile(expr, *args, **kwargs)
    except TabError:
        # if tab error, try to fix it by replacing \t with 4 spaces
        result = []
        for toknum, tokval, _, _, _  in generate_tokens(StringIO(expr).readline):
            if toknum == INDENT and '\t' in tokval:
                tokval = tokval.replace('\t', '    ')
            # the resusting string is put back to the expression (or statement)
            result.append((toknum, tokval))
        # other compiling errors are still raised
        compile(untokenize(result), *args, **kwargs)

def interpolate(text, sigil, local_dict=None, global_dict=None):
    '''Evaluate expressions in `text` marked by specified `sigil` using provided
    global and local dictionaries, and replace the expressions with their formatted strings.'''
    return text

def cfg_interpolate(text, local_dict={}):
    return interpolate(text, '${ }', local_dict, env.sos_dict.get('CONFIG', {}))

#def cfg_get(key, local_dict={}, default=None):
#    '''Get config from CONFIG dictionary'''
#    cfg = env.sos_dict.get('CONFIG', {})
#    if isinstance(key, str):
#        if key in cfg:
#            val = cfg[key]
#        else:
#            env.logger.warning('Failed to get {} from CONFIG: {} returned'.format(key, default))
#            return default
#    else:
#        val = cfg
#        for k in key:
#            if k in val:
#                val = val[k]
#            else:
#                env.logger.warning('Failed to get {} from CONFIG: {} returned'.format(key, default))
#                return default
#    if isinstance(val, str):
#        return cfg_interpolate(val, local_dict)
#    else:
#        env.logger.warning('Failed to interpolate {} from CONFIG: {} returned'.format(key, val))
#        return val

default_global_sigil = '${ }'

def set_default_global_sigil(val):
    global default_global_sigil
    if val is not None and val.count(' ') != 1:
        raise ValueError('A sigil should be specified as None or two strings separated by a space')
    default_global_sigil = val

def get_default_global_sigil():
    global default_global_sigil
    return default_global_sigil

accessed_vars_cache = {}
def accessed_vars(statement, sigil):
    '''Parse a Python statement and analyze the symbols used. The result
    will be used to determine what variables a step depends upon.'''
    global accessed_vars_cache
    if statement in accessed_vars_cache:
        return accessed_vars_cache[statement]

    if sigil is None:
        left_sigil = None
    else:
        left_sigil = sigil.split(' ')[0]
    result = set()
    prev_tok = None
    for toknum, tokval, _, _, _  in generate_tokens(StringIO(statement).readline):
        if toknum == NAME and prev_tok != '.':
            if tokval not in ('None', 'True', 'False'):
                result.add(tokval)
        if toknum == STRING and left_sigil is not None and left_sigil in tokval:
            #FIXME: get accessed_vars without string interpolation
            result = set() # |= ss.accessed_vars
        prev_tok = tokval
    return result

def SoS_eval(expr, sigil):
    '''Evaluate an expression with sos dict.'''
    return eval(expr, env.sos_dict._dict)

def _is_expr(expr):
    try:
        sos_compile(expr, '<string>', 'eval')
        return True
    except Exception:
        return False

def SoS_exec(stmts, sigil, _dict=None):
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
            sos_compile(code_group[idx], filename = '<string>', mode='exec')
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

    def value(self, sigil):
        return SoS_eval(self.expr, sigil)

    def __repr__(self):
        return 'Undetermined({!r})'.format(self.expr)

    def __hash__(self):
        raise RuntimeError('Undetermined expression should be evaluated before used. '
            'This is certainly a bug so please report this to SoS developer.')

class sos_namespace_(object):
    '''A namespace that is created by evaluating statements
    and use the results as attributes of the object.'''
    def __init__(self, stmts, sigil):
        # we need to define functions defined by sos ...
        exec('from sos.runtime import *', self.__dict__)
        # the results of the statments will be saved as
        # attribute of this object.
        SoS_exec(stmts, _dict=self.__dict__, sigil=sigil)

class on_demand_options(object):
    '''Expression that will be evaluated upon request.'''
    def __init__(self, items, sigil):
        self._expressions = {}
        if items:
            self._expressions.update(items)
        self._sigil = sigil

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
            return SoS_eval(self._expressions[key], self._sigil)
        except Exception as e:
            raise ValueError('Failed to evaluate option {} with value {}: {}'
                .format(key, self._expressions[key], e))

    def __repr__(self):
        return repr(self._expressions)
