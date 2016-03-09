#!/usr/bin/env python
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
import re
import logging
import collections

import token
from tokenize import generate_tokens, untokenize
from io import StringIO

try:
    # python 3.3
    from shlex import quote
except ImportError:
    # python 2.7
    from pipes import quote

class ColoredFormatter(logging.Formatter):
    ''' A logging formatter that uses color to differntiate logging messages
    and emphasize texts. Texts that would be empahsized are quoted with
    double backslashes (`` ``).
    '''
    def __init__(self, msg):
        logging.Formatter.__init__(self, msg)
        #
        # color for different logging levels. The current terminal color
        # is used for INFO
        self.LEVEL_COLOR = {
            'TRACE': 'DARK_CYAN',
            'DEBUG': 'BLUE',
            'WARNING': 'PURPLE',
            'ERROR': 'RED',
            'CRITICAL': 'RED_BG',
        }
        self.COLOR_CODE={
            'ENDC':0,  # RESET COLOR
            'BOLD':1,
            'UNDERLINE':4,
            'BLINK':5,
            'INVERT':7,
            'CONCEALD':8,
            'STRIKE':9,
            'GREY30':90,
            'GREY40':2,
            'GREY65':37,
            'GREY70':97,
            'GREY20_BG':40,
            'GREY33_BG':100,
            'GREY80_BG':47,
            'GREY93_BG':107,
            'DARK_RED':31,
            'RED':91,
            'RED_BG':41,
            'LIGHT_RED_BG':101,
            'DARK_YELLOW':33,
            'YELLOW':93,
            'YELLOW_BG':43,
            'LIGHT_YELLOW_BG':103,
            'DARK_BLUE':34,
            'BLUE':94,
            'BLUE_BG':44,
            'LIGHT_BLUE_BG':104,
            'DARK_MAGENTA':35,
            'PURPLE':95,
            'MAGENTA_BG':45,
            'LIGHT_PURPLE_BG':105,
            'DARK_CYAN':36,
            'AUQA':96,
            'CYAN_BG':46,
            'LIGHT_AUQA_BG':106,
            'DARK_GREEN':32,
            'GREEN':92,
            'GREEN_BG':42,
            'LIGHT_GREEN_BG':102,
            'BLACK':30,
        }

    def colorstr(self, astr, color):
        return '\033[{}m{}\033[{}m'.format(color, astr,
            self.COLOR_CODE['ENDC'])

    def emphasize(self, msg, level_color=0):
        # display text within `` and `` in green
        return re.sub(r'``([^`]*)``', '\033[32m\\1\033[{}m'.format(level_color), str(msg))

    def format(self, record):
        level_name = record.levelname
        if level_name in self.LEVEL_COLOR:
            level_color = self.COLOR_CODE[self.LEVEL_COLOR[level_name]]
            record.color_levelname = self.colorstr(level_name, level_color)
            record.color_name = self.colorstr(record.name, self.COLOR_CODE['BOLD'])
            record.color_msg = self.colorstr(self.emphasize(record.msg, level_color), level_color)
        else:
            # for INFO, use default color
            record.color_levelname = record.levelname
            record.color_msg = self.emphasize(record.msg)
        return logging.Formatter.format(self, record)

class RuntimeEnvironments(object):
    '''A singleton object that provides runtime environment for SoS.
    Atributes of this object include:

    logger: 
        a logging object 

    verbosity:
        a verbosity level object that sets the verbosity level of the logger
    
    logfile:
        name of logfile for the logger. default to no logfile.

    '''
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(RuntimeEnvironments, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        # logger
        self._logger = None
        self._verbosity = '1'
        self._logfile = None
        self._set_logger()

    #
    # attribute logger
    #
    def _set_logger(self, unused=None):
        if not hasattr(logging, 'TRACE'):
            logging.TRACE = 5
            logging.addLevelName(logging.TRACE, "TRACE")
        # create a logger, but shutdown the previous one
        if self._logger is not None:
            self._logger.handlers = []
        self._logger = logging.getLogger()
        self._logger.setLevel(logging.DEBUG)
        # output to standard output
        cout = logging.StreamHandler()
        levels = {
            '0': logging.WARNING,
            '1': logging.INFO,
            '2': logging.DEBUG,
            '3': logging.TRACE,
            None: logging.INFO
        }
        #
        cout.setLevel(levels[self._verbosity])
        cout.setFormatter(ColoredFormatter('%(color_levelname)s: %(color_msg)s'))
        self._logger.addHandler(cout)
        self._logger.trace = lambda msg, *args: self._logger._log(logging.TRACE, msg, args)
        # output to a log file
        if self._logfile is not None:
            ch = logging.FileHandler(self._logfile, mode = 'a')
            # debug informaiton and time is always written to the log file
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
            self._logger.addHandler(ch)
    #
    # atribute logger
    #
    @property
    def logger(self):
        return self._logger
    #
    # attribute verbosity
    #
    def _set_verbosity(self, v):
        if v in ['0', '1', '2', '3']:
            self._verbosity = v
            # reset logger to appropriate logging level
            self._set_logger()
    #
    verbosity = property(lambda self: self._verbosity, _set_verbosity)
    #
    # attribute logfile
    #
    def _set_logfile(self, f):
        self._logfile = f
        # reset logger to include log file
        self._set_logger()
    #
    logfile = property(lambda self: self._logfile, _set_logfile)


# set up environment variable and a default logger
env = RuntimeEnvironments()

# exception classes
class Error(Exception):
    '''Base class for SoS_ScriptParser exceptions.'''

    def _get_message(self):
        '''Getter for 'message'; needed only to override deprecation in
        BaseException.'''
        return self.__message

    def _set_message(self, value):
        '''Setter for 'message'; needed only to override deprecation in
        BaseException.'''
        self.__message = value

    # BaseException.message has been deprecated since Python 2.6.  To prevent
    # DeprecationWarning from popping up over this pre-existing attribute, use
    # a new property that takes lookup precedence.
    message = property(_get_message, _set_message)

    def __init__(self, msg=''):
        self.message = msg
        Exception.__init__(self, msg)

    def __repr__(self):
        return self.message

    __str__ = __repr__


class InterpolationError(Error):
    """Base class for interpolation-related exceptions."""

    def __init__(self, text, msg):
        #if len(text) > 20:
        #    msg = '{}: "{}..."'.format(msg, text[:20])
        #else:
        msg = '{}: "{}"'.format(msg, text)
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
        (?P<specifier>
        (?P<conversion>!\s*                 # conversion starting with !
        [s|r|q]                             # conversion, q is added by SoS                       
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
        )?                                  # optional specifier
        \s*$                                # end of tring
        '''
    
    # DOTALL makes . matchs also to newline so this supports multi-line expression
    FORMAT_SPECIFIER = re.compile(_FORMAT_SPECIFIER_TMPL, re.VERBOSE | re.DOTALL)

    def __init__(self, globals, locals, sigil = '${ }'):
        if sigil.count(' ') != 1 or sigil.startswith(' ') or sigil.endswith(' '):
            raise ValueError('Incorrect sigil "{}"'.format(sigl))
        self.sigil = sigil.split(' ')
        if self.sigil[0] == self.sigil[1]:
            raise ValueError('Incorrect sigl "{}"'.format(sigil))
        #
        self.globals = globals
        self.locals = locals

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
        pieces = text.split(self.sigil[0], 1)
        if len(pieces) == 1:
            # nothing to split, so we are done
            return text
        else:
            # the first piece must be before sigil and be completed text
            return pieces[0] + self._interpolate(pieces[1])

    def _interpolate(self, text, start_nested=0):
        '''Intepolate an expression with unknown ending location. We cannot split
        the string because of potentially nested expression. start_nested indicates
        that the location of nested expression can only happen after this point.
        '''
        # no matching }, must be wrong
        if self.sigil[1] not in text:
            raise InterpolationError(text[:20], "Missing {}".format(self.sigil[1]))
        # 
        # location of first ending sigil
        i = text.index(self.sigil[1])
        #
        # substr contains ${
        if self.sigil[0] in text[start_nested:]:
            k = text.index(self.sigil[0])
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
                return self._interpolate(text[:k] + self._interpolate(text[k + len(self.sigil[0]):]))
            except Exception as e:
                # This is for the case where inner sigil is actually part of the syntax. For example, if
                # sigil = []
                #
                # [[x*2 for x in [a, b]]]
                #
                # will first try to evaluate 
                #
                # x*2 for x in [a, b]
                #
                # without success. This part will then try to evaluate
                #
                # [x*2 x for x in [a, b]]
                #
                # namely keeping [] as python expression
                #
                return self._interpolate(text, start_nested=k + len(self.sigil[0]))
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
                        expr = mo.group('expr')
                        fmt = mo.group('specifier')
                    else:
                        expr = text[:j]
                        fmt = None
                    # if the syntax is correct
                    compile(expr, '<string>', 'eval')
                    try:
                        result = eval(expr, self.globals, self.locals)
                    except Exception as e:
                        raise InterpolationError(expr, e)
                    # evaluate the expression and interpolate the next expression
                    return self._repr(result, fmt) + self.interpolate(text[j+len(self.sigil[1]):])
                except Exception as e:
                    if self.sigil[1] not in text[j+1:]:
                        raise InterpolationError(text[:j], e)
                    j = text.index(self.sigil[1], j+1)
                    #                           j
                    # something {} } ${ another }
                    #                k
                    if j > k:
                        return self._interpolate(text[:k] + self._interpolate(text[k+len(self.sigil[0]):]))

    def _format(self, obj, fmt):
        '''Format an object in basic type (not list etc)
        '''
        # handling special !q conversion flag
        if fmt.startswith('!q'):
            # special SoS conversion for shell quotation.
            return self._format(quote(obj), fmt[2:])
        else:
            # use 
            return ('{' + fmt + '}').format(obj)
        
    def _repr(self, obj, fmt=None):
        '''Format an object. fmt will be applied to all elements if obj is not
        in a basic type. Callable object cannot be outputed (an InterpolationError
        will be raised).
        '''
        if isinstance(obj, basestring):
            return obj if fmt is None else self._format(obj, fmt)
        elif isinstance(obj, collections.Iterable):
            # the object might be nested...
            return ' '.join([self._repr(x, fmt) for x in obj])
        elif isinstance(obj, collections.Callable):
            raise InterpolationError(repr(obj), 'Cannot interpolate callable object.')
        else:
            return repr(obj) if fmt is None else self._format(obj, fmt)

def interpolate(text, globals={}, locals={}, sigil='${ }'):
    # all string is raw string....
    return SoS_String(globals, locals, sigil).interpolate(text)


class _WorkflowDict(dict):
    """A dictionary object that
    1. restrict data type to string, list of string, or dictionary of string.
    2. Generate logging message for debugging purposes.
    3. Allow access of values from attributes.
    """
    def __init__(self):
        dict.__init__(self)

    def __setitem__(self, key, value):
        if not (isinstance(value, basestring) or  \
            (isinstance(value, list) and all(isinstance(x, basestring) for x in value)) or \
            (isinstance(value, dict) and all(isinstance(x, basestring) for x in value.keys()) and 
                all(isinstance(x, basestring) for x in value.values())) ):
            raise ValueError('SoS variable can only be string or list or dictionary of strings.')
        dict.__setitem__(self, key, value)
        if isinstance(value, str) or len(value) <= 2 or len(str(value)) < 50:
            env.logger.debug('Workflow variable ``{}`` is set to ``{}``'.format(key, str(value)))
        else: # should be a list or tuple
            val = str(value).split(' ')[0] + ' ...] ({} items)'.format(len(value))
            env.logger.debug('Workflow variable ``{}`` is set to ``{}``'.format(key, val))

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def __getattr__(self, key):
        return self.__dict__[key]


def ConvertString(s, globals, locals, sigil):
    '''Convert a unicode string to a raw string and interpolate expressions
    within it by parsing the python expression and statement BEFORE they are
    evaluated (or executed).
    '''
    result = []
    # tokenize the input syntax.
    g = generate_tokens(StringIO(s.decode()).readline)
    for toknum, tokval, _, _, _  in g:
        if toknum == token.STRING:
            # if this item is a string that uses single or triple single quote
            if tokval.startswith("'"):
                # we convert it to a raw string
                tokval = u'r' + tokval 
            # we then perform interpolation on the string and put it back to expression
            tokval = repr(interpolate(eval(tokval), globals, locals, sigil))
        # the resusting string is put back to the expression (or statement)
        result.append((toknum, tokval))
    return untokenize(result)

def SoS_eval(expr, globals, locals, sigil='${ }'):
    '''Evaluate an expression after modifying (convert ' ' string to raw string,
    interpolate expressions) strings.'''
    expr = ConvertString(expr, globals, locals, sigil)
    return eval(expr, globals, locals)

def SoS_exec(stmts, globals, locals, sigil='${ }'):
    '''Execute a statement after modifying (convert ' ' string to raw string,
    interpolate expressions) strings.'''
    stmts = ConvertString(stmts, globals, locals, sigil)
    exec(stmts, globals, locals)
