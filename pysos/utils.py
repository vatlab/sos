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
import copy
import time
import types
import logging
import signal
import glob
import collections
import traceback
import hashlib
import gzip
import zipfile
import tarfile
import pickle
import yaml
import token
import psutil
import urllib
import pycurl
import blessings
import subprocess
from io import StringIO
from shlex import quote
from tokenize import generate_tokens, untokenize
from contextlib import contextmanager
from html.parser import HTMLParser
from itertools import chain

# function interpolate is needed because it is required by the SoS
# script (not seen but translated to have this function)
__all__ = ['logger', 'interpolate', 'get_output', 'expand_pattern']


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


def shortRepr(obj, noneAsNA=False):
    '''Return a short representation of obj for clarity.'''
    if obj is None:
        return 'unspecified' if noneAsNA else 'None'
    elif isinstance(obj, (str, int, float, bool)) or (isinstance(obj, collections.Sequence) \
        and len(obj) <= 2) or len(str(obj)) < 50:
        return repr(obj)
    elif isinstance(obj, collections.Sequence): # should be a list or tuple
        return repr(obj).split(' ')[0] + ' ...] ({} items)'.format(len(obj))
    elif isinstance(obj, dict):
        first_key = obj.keys()[0]
        return '{{{!r}:{!r}, ...}} ({} items)'.format(first_key, obj[first_key], len(obj))
    else:
        return '{}...'.format(repr(obj)[:40])

#
# SoS Workflow dictionary
#
class WorkflowDict(object):
    """A dictionary object that
    1. Generate logging message for debugging purposes.
    2. Generate warning message if ALLCAP variables are changed.

    IMPORTANT:

    Python does not allow the passing of a derived class of dict as globals
    to eval or exec. Doing so will result in strange behavior such as __builtins__
    not found. We then have to embed a real dictionary in WorkflowDict instead of
    deriving a dict from it.
    """
    def __init__(self, *args, **kwargs):
        self._dict = dict(*args, **kwargs)
        self._readonly_vars = {}

    def set(self, key, value):
        '''A short cut to set value to key without triggering any logging
        or warning message.'''
        self._check_readonly(key, value)
        self._dict[key] = value

    def quick_update(self, obj):
        '''Update without readonly check etc. For fast internal update'''
        self._dict.update(obj)

    def update(self, obj):
        '''Redefine update to trigger logging message'''
        for k,v in obj.items():
            self._check_readonly(k, v)
        #
        self._dict.update(obj)
        for k, v in obj.items():
            if env.verbosity > 2:
                self._log(k, v)

    def __contains__(self, key):
        return key in self._dict

    def __getattr__(self, attr):
        # for attributes that cannot be found, default to dictionary attribute
        # (e.g. keys, pop, get...)
        return getattr(self._dict, attr)

    def __getitem__(self, key):
        return self._dict[key]

    def __setitem__(self, key, value):
        '''Set value to key, trigger logging and warning messages if needed'''
        if env.verbosity > 2:
            self._log(key, value)
        if env.run_mode == 'dryrun':
            self._warn(key, value)
        if key in ('input', 'output', 'depends', '_input', '_output', '_depends', '_runtime'):
            raise ValueError('Variable {} can only be set by SoS'.format(key))
        self.set(key, value)

    def check_readonly_vars(self):
        for key in env.readonly_vars:
            if key in self._readonly_vars:
                if self._dict[key] != self._readonly_vars[key]:
                    raise RuntimeError('Variable {} is readonly and cannot be changed from {} to {}.'
                        .format(key, self._dict[key], self._readonly_vars[key]))
            elif key in self._dict:
                self._readonly_vars[key] = self._dict[key]

    def _check_readonly(self, key, value):
        if key in env.readonly_vars:
            if key not in self._readonly_vars:
                self._readonly_vars[key] = value
            # if the key already exists
            if key in self._dict:
                if self._dict[key] != self._readonly_vars[key]:
                    raise RuntimeError('Variable {} is readonly and cannot be changed from {} to {}.'
                    .format(key, self._dict[key], self._readonly_vars[key]))
                if value != self._dict[key]:
                    raise RuntimeError('Variable {} is readonly and cannot be changed from {} to {}.'
                        .format(key, self._dict[key], value))

    def _log(self, key, value):
        env.logger.debug('``{}`` = ``{}``'.format(key, shortRepr(value)))

    def _warn(self, key, value):
        if key.isupper() and key in self._dict and self._dict[key] != value:
            env.logger.warning('Changing readonly variable {} from {} to {}'
                .format(key, self._dict[key], value))
        if key.startswith('_') and not key.startswith('__') and key not in ('_input', '_output', '_step', '_index', '_depends', '_runtime'):
            env.logger.warning('{}: Variables with leading underscore is reserved for SoS temporary variables.'.format(key))

    def clone_pickleable(self):
        '''Return a copy of the existing dictionary but keep only the ones that are pickleable'''
        return {x:copy.deepcopy(y) for x,y in self._dict.items() if pickleable(y)}
#
# Runtime environment
#
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
        self.reset()

    def reset(self):
        # logger
        self._logger = None
        self._verbosity = 2
        self._logfile = None
        self._set_logger()
        #
        # run mode, this mode controls how SoS actions behave
        #
        self.run_mode = 'run'
        #
        # signature mode can be
        #
        # default              (save signature, skip if signature match)
        # ignore_signature     (ignore existing signature but still saves signature)
        # assert_signature     (should
        self.sig_mode = 'default'
        #
        # global dictionaries used by SoS during the
        # execution of SoS workflows
        self.sos_dict = WorkflowDict()
        # variables that are defined in global and parameters sections and are readonly
        self.readonly_vars = set()
        #
        # a list of variables that will be sent back from subprocess
        # in addition to aliased stepinfo. This is designed for testing
        # purposes only
        self.shared_vars = set()
        # maximum number of concurrent jobs
        self.max_jobs = 1
        self.running_jobs = 0
        # this directory will be used by a lot of processes
        self.exec_dir = os.getcwd()
        if not os.path.isdir('.sos'):
            os.mkdir('.sos')

    def register_process(self, pid, msg=''):
        '''Register a process used by this SoS instance. It will also be
        used to check resource used.'''
        proc_file = os.path.join(self.exec_dir, '.sos/proc_{}'.format(pid))
        self.logger.trace('Register {} {}'.format(pid, msg))
        with open(proc_file, 'w') as p:
            p.write(msg)

    def deregister_process(self, pid):
        proc_file = os.path.join(self.exec_dir, '.sos/proc_{}'.format(pid))
        self.logger.trace('Deregister {} at {}'.format(pid, proc_file))
        if os.path.isfile(proc_file):
            os.remove(proc_file)

    def cleanup(self):
        '''Clean up all running processes'''
        for p in glob.glob(os.path.join(self.exec_dir, '.sos/proc_*')):
            pid = int(os.path.basename(p)[5:])
            try:
                env.logger.trace('Killing {} and all its children'.format(pid))
                # psutil might not exist if SoS is not properly installed
                # but we are not acting like the end of world here
                parent = psutil.Process(pid)
                for child in parent.children(recursive=True):
                    child.kill()
                parent.kill()
            except Exception as e:
                env.logger.debug(e)
            os.remove(p)

    #
    # attribute logger
    #
    def _set_logger(self, unused=None):
        if not hasattr(logging, 'TRACE'):
            logging.TRACE = 5
            logging.addLevelName(logging.TRACE, "TRACE")
        # create a logger, we current use the regular logger but we should
        # switch to multiprocessing.get_logger if we notice trouble in, for example,
        # logging from multiple processes.
        self._logger = logging.getLogger()
        # clear previous handler
        for handler in self._logger.handlers:
            self._logger.removeHandler(handler)
        self._logger.setLevel(logging.DEBUG)
        # output to standard output
        cout = logging.StreamHandler()
        levels = {
            0: logging.ERROR,
            1: logging.WARNING,
            2: logging.INFO,
            3: logging.DEBUG,
            4: logging.TRACE,
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
        if v in [0, 1, 2, 3, 4]:
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
logger = env.logger


#
# String formatting
#
def getTermWidth():
    '''Get the width of current user terminal to properly wrap SoS
    output when well-formatted output is required.
    '''
    width = blessings.Terminal().width
    return 75 if width is None else width

class _DeHTMLParser(HTMLParser):
    '''This parser analyzes input text, removes HTML tags such as
    <p>, <br>, <ul>, <li> etc and returns properly formatted texts.
    '''
    def __init__(self):
        HTMLParser.__init__(self)
        self.__text = []

    def handle_data(self, data):
        text = data.strip()
        if len(text) > 0:
            text = re.sub('[ \t\r\n]+', ' ', text)
            self.__text.append(text + ' ')

    def handle_starttag(self, tag, attrs):
        if tag == 'p':
            self.__text.append('\n\n\n\n')
        elif tag == 'br':
            self.__text.append('\n\n')
        elif tag == 'ul':
            self.__text.append('')
        elif tag == 'li':
            self.__text.append('\n\n  * ')

    def handle_endtag(self, tag):
        if tag == 'ul':
            self.__text.append('\n\n')
        if tag == 'li':
            self.__text.append('\n\n')

    def handle_startendtag(self, tag, attrs):
        if tag == 'br':
            self.__text.append('\n\n')

    def text(self):
        return ''.join(self.__text).strip()

def dehtml(text):
    '''Remove HTML tag in input text and format the texts
    accordingly. '''
    try:
        parser = _DeHTMLParser()
        parser.feed(text)
        parser.close()
        return parser.text()
    except Exception as e:
        env.logger.warning('Failed to dehtml text: {}'.format(e))
        return text

#
# String interpolation and expression evluation
#

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
        (?P<specifier>
        (?P<conversion>!\s*                 # conversion starting with !
        [s|r|q]?                            # conversion, q is added by SoS
        (?P<sep>,)?                         # optional character to join sequences
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
                        expr, fmt, sep = mo.group('expr', 'specifier', 'sep')
                        if sep:
                            if fmt.startswith('!,'):
                                fmt = fmt[2:]
                            else:
                                fmt = fmt[:2] + fmt[3:]
                        else:
                            sep = ' '
                    else:
                        expr = text[:j]
                        fmt = None
                        sep = None
                    # if the syntax is correct
                    compile(expr, '<string>', 'eval')
                    try:
                        result = eval(expr, env.sos_dict._dict, self.local_dict)
                    except Exception as e:
                        raise InterpolationError(expr, e)
                    # evaluate the expression and interpolate the next expression
                    return self._repr(result, fmt, sep) + self.interpolate(text[j+len(self.r):])
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

    def _repr(self, obj, fmt=None, sep=' '):
        '''Format an object. fmt will be applied to all elements if obj is not
        in a basic type. Callable object cannot be outputed (an InterpolationError
        will be raised).
        '''
        if isinstance(obj, str):
            return obj if fmt is None else self._format(obj, fmt)
        elif isinstance(obj, collections.Iterable):
            # the object might be nested...
            return sep.join([self._repr(x, fmt, sep) for x in obj])
        elif isinstance(obj, collections.Callable):
            raise InterpolationError(repr(obj), 'Cannot interpolate callable object.')
        else:
            return repr(obj) if fmt is None else self._format(obj, fmt)

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
        if env.run_mode == 'dryrun':
            # make sure that the expression can be completed in 5 seconds
            with time_limit(env.sos_dict['CONFIG'].get('sos_dryrun_timeout', 5), expr):
                return eval(expr, env.sos_dict._dict)
        else:
            return eval(expr, env.sos_dict._dict)
    except Exception as e:
        if env.run_mode != 'run':
            env.sos_dict['__execute_errors__'].append(expr, e)
            return None
        else:
            raise

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
    executed = ''
    for code in code_group:
        stmts = ConvertString(code, sigil)
        if not stmts.strip():
            continue
        env.logger.trace('Executing\n{}'.format(stmts))
        try:
            if env.run_mode == 'dryrun':
                # make sure that the expression can be completed in 5 seconds
                with time_limit(env.sos_dict['CONFIG'].get('sos_dryrun_timeout', 5), stmts):
                    exec(stmts, env.sos_dict._dict)
            else:
                exec(stmts, env.sos_dict._dict)
        except Exception as e:
            if env.run_mode != 'run':
                env.sos_dict['__execute_errors__'].append(stmts, e)
            else:
                raise
        executed += stmts + '\n'
        # check if the statement has altered any readonly variables
        env.sos_dict.check_readonly_vars()
    #
    return executed



#
# Runtime signature
#

def textMD5(text):
    '''Get md5 of a piece of text'''
    m = hashlib.md5()
    m.update(text.encode())
    return m.hexdigest()

def fileMD5(filename, partial=True):
    '''Calculate partial MD5, basically the first and last 32M
    of the file for large files. This should signicicantly reduce
    the time spent on the creation and comparison of file signature
    when dealing with large bioinformat ics datasets. '''
    filesize = os.path.getsize(filename)
    # calculate md5 for specified file
    md5 = hashlib.md5()
    block_size = 2**20  # buffer of 1M
    try:
        if (not partial) or filesize < 2**26:
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    if not data:
                        break
                    md5.update(data)
        else:
            count = 64
            # otherwise, use the first and last 32M
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    count -= 1
                    if count == 32:
                        f.seek(-2**25, 2)
                    if not data or count == 0:
                        break
                    md5.update(data)
    except IOError as e:
        sys.exit('Failed to read {}: {}'.format(filename, e))
    return md5.hexdigest()


class FileSignature:
    '''Record file MD5 information to sign downloaded files, also add
    decompressed files in case the file is decompressed.'''
    def __init__(self, filename, workdir='.'):
        #
        sig_name = os.path.realpath(os.path.expanduser(filename))
        self.filenames = [sig_name]
        #
        # If the output path is outside of the current working directory
        rel_path = os.path.relpath(sig_name, os.path.realpath(workdir))
        # if this file is not relative to cache, use global signature file
        if rel_path.startswith('../'):
            self.sig_file = os.path.join(os.path.expanduser('~/.sos/.runtime'), sig_name.lstrip(os.sep) + '.file_info')
        else:
            # if this file is relative to cache, use local directory
            self.sig_file = os.path.join('.sos/.runtime', rel_path + '.file_info')
        # path to file
        sig_path = os.path.split(self.sig_file)[0]
        if not os.path.isdir(sig_path):
            try:
                os.makedirs(sig_path)
            except Exception as e:
                raise RuntimeError('Failed to create runtime directory {}: {}'.format(sig_path, e))

    def add(self, filename):
        '''add related files to the same signature'''
        self.filenames.append(os.path.abspath(os.path.expanduser(filename)))

    def write(self):
        '''Write .file_info file with signature'''
        with open(self.sig_file, 'w') as md5:
            for filename in self.filenames:
                md5.write('{}\t{}\n'.format(filename, fileMD5(filename)))

    def validate(self):
        '''Check if file matches its signature'''
        if not os.path.isfile(self.sig_file):
            return False
        with open(self.sig_file) as md5:
            for line in md5:
                f, m = line.rsplit('\t', 1)
                if not os.path.isfile(f):
                    return False
                if fileMD5(f) != m.strip():
                    env.logger.debug('MD5 mismatch {}'.format(f))
                    return False
        return True

class RuntimeInfo:
    '''Record run time information related to a number of output files. Right now only the
    .exe_info files are used.
    '''
    def __init__(self, script, input_files=[], output_files=[], dependent_files = [], index=None, pid=None, workdir='.'):
        '''Runtime information for specified output files
        index:
            in case of partial output, output files can be the same form (dynamic) so we need index to differntiate

        output_files:
            intended output file

        pid:
            process id.

        workdir:
            Current working directory.,
        '''
        self.script = script if isinstance(script, str) else ''.join(script)
        self.input_files = [input_files] if isinstance(input_files, str) else input_files
        self.output_files = [output_files] if isinstance(output_files, str) else output_files
        self.dependent_files = [dependent_files] if isinstance(dependent_files, str) else dependent_files
        if self.input_files is None:
            raise RuntimeError('Cannot create runtime signature for unknown input')
        if self.output_files is None:
            raise RuntimeError('Cannot create runtime signature for unknown output')
        #
        if self.output_files and not isinstance(self.output_files[0], Undetermined):
            sig_name = os.path.realpath(os.path.expanduser(self.output_files[0])) + textMD5('{} {} {} {}'.format(script, input_files, output_files, dependent_files))
        else:
            sig_name = textMD5('{} {} {} {} {}'.format(script, input_files, output_files, dependent_files, index))
        #
        # If the output path is outside of the current working directory
        rel_path = os.path.relpath(sig_name, os.path.realpath(workdir))
        # if this file is not relative to cache, use global signature file
        if rel_path.startswith('../'):
            info_file = os.path.join(os.path.expanduser('~/.sos/.runtime'), sig_name.lstrip(os.sep))
        else:
            # if this file is relative to cache, use local directory
            info_file = os.path.join('.sos/.runtime', rel_path)
        # path to file
        sig_path = os.path.split(info_file)[0]
        if not os.path.isdir(sig_path):
            try:
                os.makedirs(sig_path)
            except Exception as e:
                raise RuntimeError('Failed to create runtime directory {}: {}'.format(sig_path, e))
        env.logger.trace('Using signature file {} for output {} and index {}'.format(info_file, output_files, index))
        if pid is None:
            self.pid = os.getpid()
        else:
            self.pid = pid
        self.proc_info = '{}.exe_info'.format(info_file)

    def set(self, files, file_type):
        # add signature file if input and output files are dynamic
        env.logger.trace('Set {} of signature to {}'.format(file_type, files))
        if file_type == 'input':
            self.input_files = files
        elif file_type == 'output':
            self.output_files = files
        elif file_type == 'depends':
            self.depends_files = files
        else:
            raise RuntimeError('Invalid signature file type')

    def write(self):
        '''Write .exe_info file with signature of script, input, output and dependent files.'''
        if not self.proc_info:
            return
        env.logger.trace('Write signature {}'.format(self.proc_info))
        with open(self.proc_info, 'w') as md5:
            md5.write('{}\n'.format(textMD5(self.script)))
            md5.write('# input\n')
            for f in self.input_files:
                if isinstance(f, Undetermined):
                    raise ValueError('Cannot write signature for undetermined input file')
                md5.write('{}\t{}\n'.format(f, fileMD5(os.path.realpath(os.path.expanduser(f)))))
            md5.write('# output\n')
            for f in self.output_files:
                if isinstance(f, Undetermined):
                    raise ValueError('Cannot write signature for undetermined output file')
                md5.write('{}\t{}\n'.format(f, fileMD5(os.path.realpath(os.path.expanduser(f)))))
            md5.write('# dependent\n')
            for f in self.dependent_files:
                if isinstance(f, Undetermined):
                    raise ValueError('Cannot write signature for undetermined dependent file')
                md5.write('{}\t{}\n'.format(f, fileMD5(os.path.realpath(os.path.expanduser(f)))))

    def validate(self):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        if not self.proc_info or not os.path.isfile(self.proc_info):
            env.logger.trace('Fail because of no signature file {}'.format(self.proc_info))
            return False
        env.logger.trace('Validating {}'.format(self.proc_info))
        #
        # file not exist?
        self.sig_files = self.input_files + self.output_files + self.dependent_files
        if not all(isinstance(x, Undetermined) or os.path.isfile(x) for x in self.sig_files):
            env.logger.trace('Fail because of missing one of the files {}'.format(', '.join(self.sig_files)))
            return False
        #
        files_checked = {os.path.realpath(os.path.expanduser(x)):False for x in self.sig_files if not isinstance(x, Undetermined)}
        res = {'input': [], 'output': [], 'depends': []}
        cur_type = 'input'
        with open(self.proc_info) as md5:
            cmdMD5 = md5.readline().strip()   # command
            if textMD5(self.script) != cmdMD5:
                env.logger.trace('Fail because of command change')
                return False
            for line in md5:
                if line.startswith('#'):
                    if line == '# input\n':
                        cur_type = 'input'
                    elif line == '# output\n':
                        cur_type = 'output'
                    elif line == '# dependent\n':
                        cur_type = 'depends'
                    else:
                        env.logger.error('Unrecognized line in sig file {}'.format(line))
                    continue
                try:
                    f, m = line.rsplit('\t', 1)
                    freal = os.path.realpath(os.path.expanduser(f))
                    if not os.path.isfile(freal):
                        env.logger.debug('File {} not exist'.format(f))
                        return False
                    res[cur_type].append(f)
                except Exception as e:
                    env.logger.debug('Wrong md5 line {} in {}: {}'.format(line, self.proc_info, e))
                    continue
                if fileMD5(freal) != m.strip():
                    env.logger.debug('MD5 mismatch {}'.format(f))
                    return False
                # for dynamic files, they are in sig file but not in self.sig_files
                if freal in files_checked:
                    files_checked[freal] = True
        #
        if not all(files_checked.values()):
            env.logger.warning('No MD5 signature for {}'.format(', '.join(x for x,y in files_checked.items() if not y)))
            return False
        env.logger.trace('Signature matches and returns {}'.format(res))
        return res


#
# The following is adapted from snakemake so that the syntax of the pattern parameter
# would be the same as snakemake.
#
# https://bitbucket.org/snakemake/snakemake/src/22eff35627401d3bb243e068c5dce97107e7090b/snakemake/io.py?at=master&fileviewer=file-view-default
#
# I will re-implement it if there is any license issue with the code
#

_wildcard_regex = re.compile(
    "\{\s*(?P<name>\w+?)(\s*,\s*(?P<constraint>([^\{\}]+|\{\d+(,\d+)?\})*))?\s*\}")


def regex(filepattern):
    f = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        f.append(re.escape(filepattern[last:match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "If multiple wildcards of the same name "
                    "appear in a string, eventual constraints have to be defined "
                    "at the first occurence and will be inherited by the others.")
            f.append("(?P={})".format(wildcard))
        else:
            wildcards.add(wildcard)
            f.append("(?P<{}>{})".format(wildcard, match.group("constraint") if
                                         match.group("constraint") else ".+"))
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)


def glob_wildcards(pattern, files=None):
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = os.path.dirname(pattern[:first_wildcard.start(
    )]) if first_wildcard else os.path.dirname(pattern)
    if not dirname:
        dirname = "."

    names = [match.group('name')
             for match in _wildcard_regex.finditer(pattern)]
    res = {x: [] for x in names}
    pattern = re.compile(regex(pattern))

    if files is None:
        files = ((os.path.join(dirpath, f) if dirpath != "." else f)
                 for dirpath, dirnames, filenames in os.walk(dirname)
                 for f in chain(filenames, dirnames))

    for f in files:
        match = re.match(pattern, f)
        if match:
            for name, value in match.groupdict().items():
                res[name].append(value)
    return res

def apply_wildcards(pattern,
                    wildcards,
                    fill_missing=False,
                    fail_dynamic=False,
                    dynamic_fill=None,
                    keep_dynamic=False):
    def format_match(match):
        name = match.group("name")
        try:
            value = wildcards[name]
            if fail_dynamic and value == dynamic_fill:
                raise RuntimeError(name)
            return str(value)  # convert anything into a str
        except KeyError as ex:
            if keep_dynamic:
                return "{{{}}}".format(name)
            elif fill_missing:
                return dynamic_fill
            else:
                raise RuntimeError('Wildcard apply error: {} ({})'.format(ex, wildcards))

    return re.sub(_wildcard_regex, format_match, pattern)

def extract_pattern(pattern, ifiles):
    '''This function match pattern to a list of input files, extract and return
    pieces of filenames as a list of variables with keys defined by pattern.'''
    res = glob_wildcards(pattern, [])
    for ifile in ifiles:
        matched = glob_wildcards(pattern, [ifile])
        for key in matched.keys():
            if not matched[key]:
                env.logger.warning('Filename {} does not match pattern {}. None returned.'.format(ifile, pattern))
                res[key].append(None)
            else:
                res[key].extend(matched[key])
    return res

def expand_pattern(pattern):
    '''This function expand patterns against the current namespace
    and return a list of filenames'''
    ofiles = []
    sz = None
    res = glob_wildcards(pattern, [])
    sz = None
    wildcard = [{}]
    for key in res.keys():
        if key not in env.sos_dict:
            raise ValueError('Undefined variable {} in pattern {}'.format(key, pattern))
        if not isinstance(env.sos_dict[key], str) and isinstance(env.sos_dict[key], collections.Sequence):
            if sz is None:
                sz = len(env.sos_dict[key])
                wildcard = [copy.deepcopy(wildcard[0]) for x in range(sz)]
            elif sz != len(env.sos_dict[key]):
                raise ValueError('Variables in output pattern should have the same length (other={}, len({})={})'
                    .format(sz, key, len(env.sos_dict[key])))
            for idx, value in enumerate(env.sos_dict[key]):
                wildcard[idx][key] = value
        else:
            for v in wildcard:
                v[key] = env.sos_dict[key]
    #
    for card in wildcard:
        ofiles.append(apply_wildcards(pattern, card, fill_missing=False,
           fail_dynamic=False, dynamic_fill=None, keep_dynamic=False))
    return ofiles

def get_traceback():
    output = StringIO()
    exc_type, exc_value, exc_traceback = sys.exc_info()
    #print "*** print_tb:"
    traceback.print_tb(exc_traceback, limit=1, file=output)
    #print "*** print_exception:"
    traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=5, file=output)
    result = output.getvalue()
    output.close()
    return result
    #print "*** print_exc:"
    #traceback.print_exc()
    #print "*** format_exc, first and last line:"
    #formatted_lines = traceback.format_exc().splitlines()
    #print formatted_lines[0]
    #print formatted_lines[-1]
    #print "*** format_exception:"
    #print repr(traceback.format_exception(exc_type, exc_value,
    #                                      exc_traceback))
    #print "*** extract_tb:"
    #print repr(traceback.extract_tb(exc_traceback))
    #print "*** format_tb:"
    #print repr(traceback.format_tb(exc_traceback))
    #print "*** tb_lineno:", exc_traceback.tb_lineno


def pickleable(obj):
    if isinstance(obj, (str, bool, int, float, complex, bytes)):
        return True
    if isinstance(obj, (types.ModuleType, WorkflowDict)):
        return False
    try:
        pickle.dumps(obj)
        return True
    except:
        return False

class ProgressBar:
    '''A text-based progress bar, it differs from regular progress bar in that
    1. it can start from the middle with init count
    '''
    def __init__(self, message, totalCount = None, disp=True, index=None):
        if not disp:
            self.update = self.empty
            self.curlUpdate = self.empty
            self.progress = self.empty
            self.outputProgress = self.empty
            self.done = self.empty
            self.main = ''
            self.finished = 0
            return
        self.index = index
        if self.index is not None:
            self.term = blessings.Terminal(stream=sys.stderr)
        self.main = message
        self.main_start_time = time.time()
        self.message = self.main
        # get terminal width
        self.term_width = getTermWidth()
        #
        # total count
        self.count = 0
        # total initial count
        self.init_count = self.count
        #
        self.finished = 0
        self.reset('', totalCount)

    def reset(self, msg='', totalCount = None):
        if msg:
            self.message = '{} - {}'.format(self.main, msg)
        self.finished += self.count
        self.count = 0
        self.totalCount = totalCount
        self.min_progress_count = None if self.totalCount is None else self.totalCount / 1000
        self.last_progress_count = 0
        self.start_time = None
        self.last_time = None
        self.outputProgress()

    def empty(self, *args, **kwargs):
        return

    def update(self, count):
        '''completed count jobs'''
        # do not update if the diferent is less than 0.1% of the total count.
        # this is to avoid excess of calling the time() function
        if self.totalCount is not None and (count - self.count) < self.min_progress_count:
            return
        self.count = count
        self.outputProgress()

    def curlUpdate(self, total, existing, upload_t, upload_d):
        '''Update called from pycurl'''
        self.count = existing
        self.totalCount = total
        self.outputProgress()

    def progress(self, count):
        self.last_progress_count += count
        if self.last_progress_count > self.min_progress_count:
            self.count += self.last_progress_count
            self.outputProgress()
            self.last_progress_count = 0

    def outputProgress(self):
        '''Output progress'''
        if not self.start_time:
            self.start_time = time.time()
            self.last_time = self.start_time
        cur_time = time.time()
        # stop update progress bar more than once per second.
        if self.count > 0 and self.count > self.init_count and \
            self.count != self.totalCount and cur_time - self.last_time < 1:
            return
        msg = ['', '', '', '', '', '', '']
        # message
        msg[0] = self.message + ':'
        self.last_time = cur_time
        second_elapsed = cur_time - self.start_time
        if second_elapsed < 0.0001 or self.count == 0:
            msg[4] = ''
        else:
            cps = (self.count - self.init_count) / second_elapsed
            # speed
            if cps > 1000000:
                msg[4] = ' {:.1f}M/s'.format(cps/1000000)
            elif cps > 1000:
                msg[4] = ' {:.1f}K/s'.format(cps/1000)
            elif cps > 0.05:
                msg[4] = ' {:.1f}/s'.format(cps)
            elif cps > 1e-6:
                msg[4] = ' {:.1f}s each'.format(1. / cps)
            else:
                msg[4] = ' 0.0/s'
        # estimated time left
        if self.totalCount:
            perc = min(1, float(self.count) / self.totalCount)
            init_perc = min(1, float(self.init_count) / self.totalCount)
            time_left = (second_elapsed / (perc - init_perc) * (1 - perc)) if perc > init_perc else 0
            msg[5] += ' in {}{}'.format('' if time_left < 86400 else '{} day{} '
                .format(int(time_left/86400), 's' if time_left > 172800 else ''),
                time.strftime('%H:%M:%S', time.gmtime(time_left)))
        # percentage / progress
        if self.count > 0:
            msg[3] = ' {:,}'.format(int(self.count))
            m3Len = len(msg[3])
        else:
            msg[3] = ' '
            m3Len = 1
        if self.totalCount:
            # percentage
            perc = min(1, float(self.count) / self.totalCount)
            msg[1] = ' {:5.1f}%'.format(perc * 100)
            width = self.term_width - len(msg[0]) - len(msg[1]) - m3Len - len(msg[4]) - len(msg[5])
            if width > 5:
                front = int(perc * (width - 4))
                back = width - 4 - front
                msg[2] = ' [{}>{}]'.format('=' * front, ' ' * back)
        else:
            width = self.term_width - len(msg[0]) - len(msg[1]) - m3Len - len(msg[4])
            msg[6] = ' '*width
        # use stderr to avoid messing up process output
        if self.index is None:
            sys.stderr.write('\r' + ''.join(msg))
        else:
            with self.term.location( 0, self.term.height - self.index - 1):
                sys.stderr.write('\r' + ''.join(msg))

    def done(self, done_msg=''):
        '''Finish, output a new line'''
        # if an message is given, display and quit.
        if self.totalCount and not done_msg:
            self.count = self.totalCount
        #
        msg = ['', '', '', '', '', '']
        # message
        msg[0] = self.main + ':'
        second_elapsed = time.time() - self.main_start_time
        cps = 0 if second_elapsed < 0.0001 else (self.finished + self.count) / second_elapsed
        # speed
        if cps > 1000000:
            msg[4] = ' {:.1f}M/s'.format(cps/1000000)
        elif cps > 1000:
            msg[4] = ' {:.1f}K/s'.format(cps/1000)
        elif cps > 0.05:
            msg[4] = ' {:.1f}/s'.format(cps)
        elif cps > 1e-6:
            msg[4] = ' {:.1f}s each'.format(1. / cps)
        else:
            msg[4] = ' 0.0/s'
        #
        msg[3] = ' {:,}'.format(self.finished + self.count)
        m3Len = len(msg[3])
        msg[5] += ' in {}{}'.format('' if second_elapsed < 86400 else '{} day{} '
            .format(int(second_elapsed/86400), 's' if second_elapsed > 172800 else ''),
                time.strftime('%H:%M:%S', time.gmtime(second_elapsed)))
        # percentage / progress
        if self.totalCount:
            # percentage
            msg[1] = ' 100%'
            width = self.term_width - len(msg[0]) - len(msg[1]) - m3Len - len(msg[4]) - len(msg[5])
            if width > 4:
                front = int(width - 3)
                msg[2] = ' [{}]'.format('=' * front)
        #
        if done_msg:
            msg[0] = done_msg
            msg[1] = ''
            msg[2] = ''
            msg[3] = ''
            msg[4] = ''
            msg[5] = ''
        if self.index is None:
            sys.stderr.write('\r' + ''.join(msg) + '\n')
            sys.stderr.flush()
        else:
            with self.term.location(0, self.term.height - self.index - 1):
                sys.stderr.write('\r' + ''.join(msg))
                sys.stderr.flush()

#
# download file with progress bar
#

def downloadURL(URL, dest, decompress=False, index=None):
    # use libcurl
    dest = os.path.abspath(os.path.expanduser(dest))
    dest_dir, filename = os.path.split(dest)
    #
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    if not os.path.isdir(dest_dir):
        raise RuntimeError('Failed to create destination directory to download {}'.format(URL))
    #
    message = filename
    if len(message) > 30:
        message = message[:10] + '...' + message[-16:]
    #
    dest_tmp = dest + '.tmp_{}'.format(os.getpid())
    term_width = getTermWidth()
    try:
        env.logger.debug('Download {} to {}'.format(URL, dest))
        prog = ProgressBar(message, disp=env.verbosity > 1, index=index)
        sig = FileSignature(dest)
        if os.path.isfile(dest) and sig.validate():
            prog.done(message + ': \033[32m use existing {}\033[0m'.format(' '*(term_width - len(message) - 15)))
            return True
        #
        with open(dest_tmp, 'wb') as f:
            c = pycurl.Curl()
            c.setopt(pycurl.URL, str(URL))
            c.setopt(pycurl.WRITEFUNCTION, f.write)
            c.setopt(pycurl.SSL_VERIFYPEER, False)
            c.setopt(pycurl.NOPROGRESS, False)
            c.setopt(pycurl.PROGRESSFUNCTION, prog.curlUpdate)
            c.perform()
        if c.getinfo(pycurl.HTTP_CODE) == 404:
            prog.done(message + ':\033[91m 404 Error {}\033[0m'.format(' '*(term_width - len(message) - 12)))
            try:
                os.remove(dest_tmp)
            except OSError:
                pass
            return False
        os.rename(dest_tmp, dest)
        decompressed = 0
        if decompress:
            if zipfile.is_zipfile(dest):
                zip = zipfile.ZipFile(dest)
                zip.extractall(dest_dir)
                names = zip.namelist()
                for name in names:
                    if not os.path.isfile(os.path.join(dest_dir, name)):
                        return False
                    else:
                        sig.add(os.path.join(dest_dir, name))
                        decompressed += 1
            elif tarfile.is_tarfile(dest):
                with tarfile.open(dest, 'r:*') as tar:
                    tar.extractall(dest_dir)
                    # only extract files
                    files = [x.name for x in tar.getmembers() if x.isfile()]
                    for name in files:
                        if not os.path.isfile(os.path.join(dest_dir, name)):
                            return False
                        else:
                            sig.add(os.path.join(dest_dir, name))
                            decompressed += 1
            elif dest.endswith('.gz'):
                decomp = dest[:-3]
                with gzip.open(dest, 'rb') as fin, open(decomp, 'wb') as fout:
                    buffer = fin.read(100000)
                    while buffer:
                        fout.write(buffer)
                        buffer = fin.read(100000)
                sig.add(decomp)
                decompressed += 1
        decompress_msg = '' if not decompressed else ' ({} file{} decompressed)'.format(
            decompressed, '' if decompressed <= 1 else 's')
        prog.done(message + ':\033[32m downloaded{} {}\033[0m'.format(decompress_msg,
            ' '*(term_width - len(message) - 13 - len(decompress_msg))))
    except Exception as e:
        env.logger.error(e)
        return False
    finally:
        # if there is something wrong still remove temporary file
        if os.path.isfile(dest_tmp):
            os.remove(dest_tmp)
    sig.write()
    return os.path.isfile(dest)

class frozendict(dict):
    '''A fronzen dictionary that disallow changing of its elements
    Copied from http://code.activestate.com/recipes/414283/
    '''
    def _blocked_attribute(obj):
        raise RuntimeError("Cannot modify a readonly dictionary.")
    _blocked_attribute = property(_blocked_attribute)

    __delitem__ = __setitem__ = clear = _blocked_attribute
    pop = popitem = setdefault = update = _blocked_attribute

    def __new__(cls, *args):
        new = dict.__new__(cls)
        dict.__init__(new, *args)
        return new

    def __init__(self, *args):
        pass

    def __hash__(self):
        try:
            return self._cached_hash
        except AttributeError:
            h = self._cached_hash = hash(tuple(sorted(self.items())))
            return h

    def __getattr__(self, key):
        return dict.__getitem__(self, key)

    def __setattr__(self, key, value):
        raise RuntimeError('Cannot modify a readonly dictionary')

    def __repr__(self):
        return "frozendict(%s)" % dict.__repr__(self)

#
# dynamic expression that cannot be resolved during parsing
# at dryrun mode etc, and has to be resolved at run time.
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

#
# A utility function that returns output of a command
def get_output(cmd):
    try:
        output = subprocess.check_output(cmd, shell=True).decode()
    except subprocess.CalledProcessError as e:
        if e.output.decode():
            env.logger.error(e.output.decode())
        raise RuntimeError(e)
    return output

#
# search a path and locate script and other files
#
def locate_script(filename, start=''):
    #
    attemp = os.path.expanduser(filename)
    if os.path.isfile(attemp):
        return ('', attemp)
    #
    token = urllib.parse.urlparse(filename)
    # if no scheme or netloc, the URL is not acceptable
    if all([getattr(token, qualifying_attr) for qualifying_attr in  ('scheme', 'netloc')]):
        try:
            local_filename, headers = urllib.request.urlretrieve(filename)
            with open(local_filename) as script:
                content = script.read()
            #
            return (content, filename)
        except Exception as e:
            env.logger.error(e)
            raise ValueError('Failed to open {}'.format(filename))
    #
    # a search path
    pathes = [start]
    sos_config_file = os.path.expanduser('~/.sos/config.json')
    if os.path.isfile(sos_config_file):
        try:
            with open(sos_config_file) as config:
                cfg = yaml.safe_load(config)
        except Exception as e:
            raise RuntimeError('Failed to parse global sos config file {}, is it in JSON format?'.format(sos_config_file))
        #
        pathes.extend(cfg.get('sos_path', []))
    #
    sos_config_file = '.sos/config.json'
    if os.path.isfile(sos_config_file):
        try:
            with open(sos_config_file) as config:
                cfg = yaml.safe_load(config)
        except Exception as e:
            raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format?'.format(sos_config_file))
        #
        pathes.extend(cfg.get('sos_path', []))
    #
    for path in pathes:
        if not path:
            continue
        attemp = os.path.join(os.path.expanduser(path), os.path.expanduser(filename))
        if os.path.isfile(attemp):
            return ('', attemp)
        # is it an URL?
        token = urllib.parse.urlparse(path)
        # if no scheme or netloc, the URL is not acceptable
        if all([getattr(token, qualifying_attr) for qualifying_attr in  ('scheme', 'netloc')]):
            url = path + ('' if path.endswith('/') else '/') + filename
            try:
                local_filename, headers = urllib.request.urlretrieve(url)
                with open(local_filename) as script:
                    content = script.read()
                return content, url
            except Exception as e:
                pass
    #
    raise ValueError('Failed to locate {}'.format(filename))
