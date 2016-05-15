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
import sys
import re
import signal
import collections
import token
from io import StringIO
from shlex import quote
from tokenize import generate_tokens, untokenize
from contextlib import contextmanager

from .utils import env, Error, shortRepr

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

    # Format specifier that can be used at the end of the string to convert
    # result (!s|r|q) and control output format (:.2f etc). The pattern
    # is constructed according to Python format mini language.
    _FORMAT_SPECIFIER_TMPL = r'''
        ^                                   # start of expression
        (?P<expr>.*?)                       # any expression
        (?P<conversion>!\s*                 # conversion starting with !
        [srqabe,]+                          # conversion, q, a, b, and , are added by SoS
        )?
        (?P<format_spec>:\s*                # format_spec starting with :
        (?P<fill>.?[<>=^])?                 # optional fill|align
        (?P<sign>[-+ ])?                    # optional sign
        \#?                                 #
        0?                                  #
        (?P<width>\d+)?                     # optional width
        (?P<precision>\.\d+)?               # optional precision
        (?P<type>[bcdeEfFgGnosxX%])?        # optional type
        )?                                  # optional format_spec
        \s*$                                # end of tring
        '''

    # DOTALL makes . matchs also to newline so this supports multi-line expression
    FORMAT_SPECIFIER = re.compile(_FORMAT_SPECIFIER_TMPL, re.VERBOSE | re.DOTALL)

    def __init__(self, sigil = '${ }', local_dict={}):
        # do not check sigil here because the function will be called quite frequently
        # the sigil will be checked when it is entered in SoS script.
        self.l, self.r = sigil.split(' ')
        self.error_count = 0
        self.local_dict = local_dict

    def interpolate(self, text):
        '''Intepolate string with local and global dictionary'''
        #
        # We could potentially parse the text and find all interpolation text,
        # but we cannot really do it because of possible nested interpolation
        #
        # split by left sigil
        #
        # '${a} part1 ${ expr2 ${ nested }} and another ${expr2 {}} and done'
        #
        # '' 'a} part1 ' ' expr2 ' 'nested }} and another' 'expr2 {}} and done'
        #
        # 'in' test is 10 times faster than split so we do this test first.
        if self.l not in text:
            return text
        pieces = text.split(self.l, 1)
        # the first piece must be before sigil and be completed text
        #env.logger.trace('"{}" interpolated to "{}"'.format(text, res))
        return pieces[0] + self._interpolate(pieces[1])

    def _interpolate(self, text, start_nested=0):
        '''Intepolate an expression with unknown ending location. We cannot split
        the string because of potentially nested expression. start_nested indicates
        that the location of nested expression can only happen after this point.
        '''
        # no matching }, must be wrong
        if self.r not in text:
            raise InterpolationError(text[:20], "Missing {}".format(self.r))
        #
        # location of first ending sigil
        i = text.index(self.r)
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
                    mo = self.FORMAT_SPECIFIER.match(text[:j])
                    if mo:
                        expr, fmt, conversion = mo.group('expr', 'format_spec', 'conversion')
                    else:
                        expr = text[:j]
                        fmt = None
                        conversion = None
                    # if the syntax is correct
                    compile(expr, '<string>', 'eval')
                    try:
                        result = eval(expr, env.sos_dict._dict, self.local_dict)
                    except Exception as e:
                        raise InterpolationError(expr, e)
                    # evaluate the expression and interpolate the next expression
                    return self._repr(result, fmt, conversion) + self.interpolate(text[j+len(self.r):])
                except Exception as e:
                    self.error_count += 1
                    if self.r not in text[j+1:]:
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
            if isinstance(obj, str):
                if 'e' in conversion:
                    obj = os.path.expanduser(obj)
                if 'a' in conversion:
                    obj = os.path.abspath(os.path.expanduser(obj))
                if 'b' in conversion:
                    obj = os.path.basename(obj)
                if 'q' in conversion:
                    # special SoS conversion for shell quotation.
                    obj = quote(obj)
            if 'r' in conversion:
                obj = repr(obj)
            if 's' in conversion:
                obj = str(obj)
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
            sep = ',' if conversion and ',' in conversion else ' '
            return sep.join([self._repr(x, fmt, conversion) for x in obj])
        elif isinstance(obj, collections.Callable):
            raise InterpolationError(repr(obj), 'Cannot interpolate callable object.')
        else:
            return repr(obj) if fmt is None and conversion is None else self._format(obj, fmt, conversion)

def interpolate(text, sigil='${ }', local_dict={}):
    '''Evaluate expressions in `text` marked by specified `sigil` using provided
    global and local dictionaries, and replace the expressions with their formatted strings.'''
    return SoS_String(sigil, local_dict).interpolate(text)


def ConvertString(s, sigil):
    '''Convert a unicode string to a raw string and interpolate expressions
    within it by parsing the python expression and statement BEFORE they are
    evaluated (or executed).

    FIXME: the expression might have a dynamic option which should prevent
    string interpolation. Not sure how to handle this option right now.
    '''
    result = []
    left_sigil = sigil.split(' ')[0]
    # tokenize the input syntax.
    if sys.version_info.major == 2:
        g = generate_tokens(StringIO(s.decode()).readline)
    else:
        # python 3 string is already unicode string
        g = generate_tokens(StringIO(s).readline)
    for toknum, tokval, _, _, _  in g:
        if toknum == token.STRING:
            # if this item is a string that uses triple single quote
            if tokval.startswith("'''"):
                # we convert it to a raw string
                tokval = u'r' + tokval
            # we then perform interpolation on the string and put it back to expression
            if left_sigil in tokval:
                tokval = 'interpolate(' + tokval + ", \'" + sigil + "', locals())"
        # the resusting string is put back to the expression (or statement)
        result.append((toknum, tokval))
    return untokenize(result)


class TimeoutException(Exception):
    def __init__(self, msg=''):
        self.msg = msg

@contextmanager
def time_limit(seconds, msg=''):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out for option {}".format(msg))
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def SoS_eval(expr, sigil='${ }'):
    '''Evaluate an expression after modifying (convert ' ' string to raw string,
    interpolate expressions) strings.'''
    expr = ConvertString(expr, sigil)
    try:
        if env.run_mode == 'inspect':
            # make sure that the expression can be completed in 5 seconds
            with time_limit(env.sos_dict['CONFIG'].get('sos_inspect_timeout', 5), expr):
                return eval(expr, env.sos_dict._dict)
        else:
            return eval(expr, env.sos_dict._dict)
    except Exception as e:
        if env.run_mode != 'run':
            env.sos_dict['__execute_errors__'].append(expr, e)
            return None
        else:		
            raise

def _is_expr(expr):
    try:
        compile(expr, '<string>', 'eval')
        return True
    except:
        return False

def SoS_exec(stmts, sigil='${ }'):
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
    while True:
        try:
            # test current group
            compile(code_group[idx], filename = '<string>', mode='exec')
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
        if env.run_mode == 'inspect':
            env.logger.trace('Checking\n{}'.format(stmts))
        elif env.run_mode == 'prepare':
            env.logger.trace('Preparing\n{}'.format(stmts))
        else:
            env.logger.trace('Executing\n{}'.format(stmts))
        #
        try:
            if env.run_mode == 'inspect':
                # make sure that the expression can be completed in 5 seconds
                with time_limit(env.sos_dict['CONFIG'].get('sos_inspect_timeout', 5), stmts):
                    if idx + 1 == len(code_group) and _is_expr(stmts):
                        res = eval(stmts, env.sos_dict._dict)
                    else:
                        exec(stmts, env.sos_dict._dict)
            elif idx + 1 == len(code_group) and _is_expr(stmts):
                res = eval(stmts, env.sos_dict._dict)
            else:
                exec(stmts, env.sos_dict._dict)
        except Exception as e:
            if env.run_mode != 'run':
                if isinstance(e, InterpolationError):
                    if env.run_mode == 'inspect':
                        # this should not matter in inspect mode because many variables do not yet
                        # exist...
                        env.logger.debug('Failed to interpolate {}: {}'.format(shortRepr(stmts), e))
                else:
                    env.sos_dict['__execute_errors__'].append(stmts, e)
            else:
                raise
        executed += stmts + '\n'
        # check if the statement has altered any readonly variables
        env.sos_dict.check_readonly_vars()
    return res

#
# dynamic expression that cannot be resolved during parsing
# at inspect mode etc, and has to be resolved at run time.
#
class Undetermined(object):
    def __init__(self, expr):
        if not isinstance(expr, str):
            raise RuntimeError('Undetermined expression has to be a string: "{}" passed'.format(expr))
        self.expr = expr.strip()

    def value(self, sigil='${ }'):
        return SoS_eval(self.expr, sigil)

    def __repr__(self):
        return 'Undetermined({!r})'.format(self.expr)

    def __hash__(self):
        raise RuntimeError('Undetermined expression should be evaluated before used. '
            'This is certainly a bug so please report this to SoS developer.')


