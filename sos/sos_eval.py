#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
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
import os
import re
import collections
from io import StringIO
from shlex import quote
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

def sos_compile(expr, *args, **kwargs):
    ''' Compiling an statement but ignores tab error (mixed tab and space'''
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

#
# String intepolation
#
class SoS_String:
    '''This class implements SoS string that support interpolation using
    a local and a global dictionary. It is similar to format string of
    Python 3.6 with configurable sigil (default to ${ }) and allows the
    resulting object of types other than simple int, None etc. It also supports
    a special conversion flag !q that properly quote a filename so that
    it can be used safely within a shell script. '''

    LEFT_PATTERNS = {
        # if not preceded by a backslash
        '${': re.compile(r'(?<!\\)\$\{')
        }

    CONVERTERS = {
        'e': os.path.expanduser,
        'a': lambda x: os.path.abspath(os.path.expanduser(x)),
        'd': os.path.dirname,
        'b': os.path.basename,
        'n': lambda x: os.path.splitext(x)[0],
        'q': quote,
        'r': repr,
        's': str,
        # these are handled elsewhere
        ',': lambda x: x,
        '!': lambda x: x,
        }

    def __init__(self, sigil, local_dict={}, trace_vars=False):
        # do not check sigil here because the function will be called quite frequently
        # the sigil will be checked when it is entered in SoS script.
        self.l, self.r = sigil.split(' ')
        self.default_sigil = sigil == '${ }'
        if self.l not in self.LEFT_PATTERNS:
            self.LEFT_PATTERNS[self.l] = re.compile(r'(?<!\\){}'.format(re.escape(self.l)))
        self.left_pattern = self.LEFT_PATTERNS[self.l]
        self.error_count = 0
        self.local_dict = local_dict
        self.my_eval = eval
        if trace_vars:
            self.accessed_vars = set()

    def interpolate(self, text):
        '''Intepolate string with local and global dictionary'''
        #
        # We could potentially parse the text and find all interpolation text,
        # but we cannot really do it because of possible nested interpolation
        #
        # 'in' test is 10 times faster than split so we do this test first.
        # here we do not consider the case of \${
        if self.l not in text:
            return text
        #
        # this function uses a direct substitution method to handle simple
        # cases ${var}. Performance test shows that it can cut string interpolation
        # time roughtly in half.
        if self.default_sigil:
            text = self.direct_interpolate(text)
            # 'in' test is 10 times faster than split so we do this test first.
            if self.l not in text:
                return text
        #
        # split by left sigil
        #
        # '${a} part1 ${ expr2 ${ nested }} and another ${expr2 {}} and done'
        #
        # '' 'a} part1 ' ' expr2 ' 'nested }} and another' 'expr2 {}} and done'
        #
        pieces = self.left_pattern.split(text, 1)
        if len(pieces) == 1:
            return pieces[0].replace('\\' + self.l, self.l)
        # the first piece must be before sigil and be completed text
        #env.logger.trace('"{}" interpolated to "{}"'.format(text, res))
        return (pieces[0] + self._interpolate(pieces[1])).replace('\\' + self.l, self.l)

    def direct_interpolate(self, text):
        pieces = SIMPLE_SUB.split(text)
        # replace pieces 1, 3, 5, ... etc with their values
        for i in range(1, len(pieces), 2):
            if hasattr(self, 'accessed_vars'):
                self.accessed_vars |= accessed_vars(pieces[i], '${ }')
                pieces[i] = ''
            else:
                pieces[i] = self._repr(eval(pieces[i], env.sos_dict._dict, self.local_dict))
        return ''.join(pieces)

    def _interpolate(self, text, start_nested=0):
        '''Intepolate an expression with unknown ending location. We cannot split
        the string because of potentially nested expression. start_nested indicates
        that the location of nested expression can only happen after this point.
        '''
        # no matching }, must be wrong
        if self.r not in text:
            raise InterpolationError(text[:20], "Missing ending sigil {}".format(self.r))
        #
        # location of first ending sigil
        i = text.index(self.r)
        if i == 0:
            # if no expression is found
            return self.interpolate(text[len(self.r):])
        #
        # substr contains ${
        if self.l in text[start_nested:]:
            k = text.index(self.l)
        else:
            # k is very far away
            k = len(text) + 100
        #
        # nested
        if k < i:
            #                     i
            # something  ${ nested} }
            #            k
            #
            try:
                # expolate string recursively.
                return self._interpolate(text[:k] + self._interpolate(text[k + len(self.l):]))
            except Exception as e:
                self.error_count += 1
                if self.error_count > 10:
                    raise
                # This is for the case where inner sigil is actually part of the syntax. For example, if
                # sigil = []
                #
                #     [[x*2 for x in [a, b]]]
                #
                # will first try to evaluate
                #
                #     x*2 for x in [a, b]
                #
                # without success. This part will then try to evaluate
                #
                #     [x*2 x for x in [a, b]]
                #
                # namely keeping [] as python expression
                #
                return self._interpolate(text, start_nested=k + len(self.l))
        else:
            # non-nested case, proceed from left to right
            #
            #            i
            # something {} } ${ another }
            #                k
            j = i
            while True:
                try:
                    # is syntax correct?
                    mo = FORMAT_SPECIFIER.match(text[:j])
                    if mo:
                        expr, fmt, conversion = mo.group('expr', 'format_spec', 'conversion')
                    else:
                        expr = text[:j]
                        fmt = None
                        conversion = None
                    # if the syntax is correct
                    sos_compile(expr, '<string>', 'eval')
                    try:
                        if hasattr(self, 'accessed_vars'):
                            self.accessed_vars |= accessed_vars(expr, self.l + ' ' + self.r)
                            return self.interpolate(text[j+len(self.r):])
                        else:
                            result = eval(expr, env.sos_dict._dict, self.local_dict)
                            return self._repr(result, fmt, conversion) + self.interpolate(text[j+len(self.r):])
                    except Exception as e:
                        raise InterpolationError(expr, e)
                    # evaluate the expression and interpolate the next expression
                except Exception as e:
                    self.error_count += 1
                    if self.error_count > 10:
                        raise
                    if self.r not in text[j+1:]:
                        if hasattr(self, 'accessed_vars'):
                            return self.interpolate(text[len(self.r):])
                        else:
                            raise InterpolationError(text[:j], e)
                    j = text.index(self.r, j+1)
                    #                           j
                    # something {} } ${ another }
                    #                k
                    if j > k:
                        return self._interpolate(text[:k] + self._interpolate(text[k+len(self.l):]))

    def _format(self, obj, fmt, conversion):
        '''Format an object in basic type (not list etc)
        '''
        # handling special !q conversion flag
        if conversion:
            for c in conversion:
                obj = self.CONVERTERS[c](obj)
        return ('{' + (fmt if fmt else '') + '}').format(obj)

    def _repr(self, obj, fmt=None, conversion=None):
        '''Format an object. fmt will be applied to all elements if obj is not
        in a basic type. Callable object cannot be outputed (an InterpolationError
        will be raised).
        '''
        if isinstance(obj, str):
            return obj if fmt is None and conversion is None else self._format(obj, fmt, conversion)
        elif isinstance(obj, collections.Iterable):
            # the object might be nested...
            sep = ', ' if conversion and ',' in conversion else ' '
            return sep.join([self._repr(x, fmt, conversion) for x in obj])
        elif isinstance(obj, (collections.Callable, Undetermined)):
            raise InterpolationError(repr(obj), 'Cannot interpolate callable object.')
        else:
            return repr(obj) if fmt is None and conversion is None else self._format(obj, fmt, conversion)

def interpolate(text, sigil, local_dict={}):
    '''Evaluate expressions in `text` marked by specified `sigil` using provided
    global and local dictionaries, and replace the expressions with their formatted strings.'''
    return SoS_String(sigil, local_dict).interpolate(text)


default_global_sigil = '${ }'

def set_default_global_sigil(val):
    global default_global_sigil
    if val is not None and val.count(' ') != 1:
        raise ValueError('A sigil should be specified as None or two strings separated by a space')
    default_global_sigil = val

def get_default_global_sigil():
    global default_global_sigil
    return default_global_sigil
    
def ConvertString(s, sigil):
    '''Convert a unicode string to a raw string and interpolate expressions
    within it by parsing the python expression and statement BEFORE they are
    evaluated (or executed).

    FIXME: the expression might have a dynamic option which should prevent
    string interpolation. Not sure how to handle this option right now.
    '''
    # if no sigil is specified, only check space/tab...
    indent_space = False
    indent_tab = False
    result = []
    if sigil is None:
        if '\t' not in s:
            return s
        for toknum, tokval, _, _, _  in generate_tokens(StringIO(s).readline):
            if toknum == INDENT:
                if '\t' in tokval:
                    tokval = tokval.replace('\t', '    ')
                    indent_tab = True
                elif ' ' in tokval:
                    indent_space = True
            # the resusting string is put back to the expression (or statement)
            result.append((toknum, tokval))
        if indent_space and indent_tab:
            env.logger.warning('Tabs converted to 4 spaces due to mixed use of tab and space in statement {}'.format(short_repr(s)))
    else:
        left_sigil = sigil.split(' ')[0]
        if left_sigil not in s and not '\t' in s:
            return s
        # tokenize the input syntax.
        for toknum, tokval, _, _, _  in generate_tokens(StringIO(s).readline):
            if toknum == STRING:
                # if this item is a string that uses triple single quote
                # if tokval.startswith("'''"):
                #     # we convert it to a raw string
                #     tokval = u'r' + tokval
                # we then perform interpolation on the string and put it back to expression
                if (tokval.startswith('"') or tokval.startswith('r"') or tokval.startswith('u"')) and left_sigil in tokval:
                    tokval = 'interpolate(' + tokval + ", \'" + sigil + "', locals())"
            if toknum == INDENT:
                if '\t' in tokval:
                    tokval = tokval.replace('\t', '    ')
                    indent_tab = True
                elif ' ' in tokval:
                    indent_space = True
            # the resusting string is put back to the expression (or statement)
            result.append((toknum, tokval))
    if indent_space and indent_tab:
        env.logger.warning('Tabs converted to 4 spaces due to mixed use of tab and space in statement {}'.format(short_repr(s)))
    return untokenize(result)

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
            result.add(tokval)
        if toknum == STRING and left_sigil is not None and left_sigil in tokval:
            # if it is a string, check if variables used during
            # string interpolation
            ss = SoS_String(sigil, {}, True)
            ss.interpolate(eval(tokval))
            result |= ss.accessed_vars
        prev_tok = tokval
    return result

def param_of(name, text):
    '''return parameters of parameter of name
    for example:
        name='input'
        text='func(input=1, output=2)'
    returns 1 as the parameter of input. '''
    params = re.split(r'({}\s*=\s*)'.format(name), text)
    exprs = []
    for param in params[2::2]:
        expr = ''
        try:
            for _, tokval, _, _, _ in generate_tokens(StringIO(param).readline):
                try:
                    expr += tokval
                    compile(expr, '<string>', 'eval')
                    exprs.append(expr)
                    break
                except:
                    continue
        except:
            continue
    return exprs

def SoS_eval(expr, sigil):
    '''Evaluate an expression after modifying (convert ' ' string to raw string,
    interpolate expressions) strings.'''
    expr = ConvertString(expr, sigil)
    return eval(expr, env.sos_dict._dict)

def _is_expr(expr):
    try:
        sos_compile(expr, '<string>', 'eval')
        return True
    except:
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
        except:
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
        stmts = ConvertString(code, sigil)
        if not stmts.strip():
            continue
        if env.run_mode == 'prepare':
            env.logger.trace('Preparing statement:\n{}'.format(stmts))
        #else:
        #    env.logger.trace('Executing statement:\n{}'.format(stmts))
        #
        try:
            if env.run_mode == 'interactive':
                act = DelayedAction(env.logger.info, 'Running {}'.format(short_repr(code)))
            else:
                act = None
            if idx + 1 == len(code_group) and _is_expr(stmts):
                res = eval(stmts, _dict)
            else:
                exec(stmts, _dict)
        finally:
            del act
        executed += stmts + '\n'
        # check if the statement has altered any readonly variables
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
