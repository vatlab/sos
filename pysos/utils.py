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
import os
import sys
import re
import copy
import types
import logging
import collections
import hashlib
import shutil
import token
from tokenize import generate_tokens, untokenize
from io import StringIO

try:
    # python 3.3
    from shlex import quote
except ImportError:
    # python 2.7
    from pipes import quote

try:
    # python 3.x
    from html.parser import HTMLParser
except ImportError:
    # python 2.7
    from HTMLParser import HTMLParser

#
# Logging
#

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


def shortRepr(obj):
    '''Return a short representation of obj for clarity.'''
    if isinstance(obj, (str, int, float, bool)) or (isinstance(obj, collections.Sequence) \
        and len(obj) <= 2) or len(str(obj)) < 50:
        return str(obj)
    elif isinstance(obj, collections.Sequence): # should be a list or tuple
        return str(obj).split(' ')[0] + ' ...] ({} items)'.format(len(obj))
    elif isinstance(obj, dict):
        first_key = obj.keys()[0]
        return '{{{}:{}, ...}} ({} items)'.format(first_key, obj[first_key], len(obj))
    else:
        return '{}...'.format(repr(obj)[:40])

#
# SoS Workflow dictionary
#
class WorkflowDict(dict):
    """A dictionary object that
    1. Generate logging message for debugging purposes.
    2. Generate warning message if ALLCAP variables are changed.
    """
    _protect_vars_assigned = False
     
    class protect_vars_assigned:
        '''A environment that change the env.protected_vars_assigned to true
        when statement is executed in this mode. This is friendier with
        excetion because _protect_vars_assigned would be turned off as soon
        as the statement finishes, or if an exception raises.
        '''
        def __init__(self, wf_dict):
            self.wf_dict = wf_dict

        def __enter__(self):
            self.wf_dict._protect_vars_assigned = True

        def __exit__(self,  etype, value, traceback):
            self.wf_dict._protect_vars_assigned = False

    def readonly_assignment(self):
        return self.protect_vars_assigned(self)
    
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self._readonly = set()

    def set(self, key, value):
        '''A short cut to set value to key without triggering any logging
        or warning message.'''
        self._check_readonly(key, value)
        dict.__setitem__(self, key, value)
        if self._protect_vars_assigned:
            self._readonly.add(key)

    def update(self, obj):
        '''Redefine update to trigger logging message'''
        for k,v in obj.items():
            self._check_readonly(k, v)
        #
        dict.update(self, obj)
        for k, v in obj.items():
            if env.verbosity > 2:
                self._log(k, v)
            if self._protect_vars_assigned:
                self._readonly.add(k)

    def clear(self):
        dict.clear(self)
        self._readonly = set()

    def __setitem__(self, key, value):
        '''Set value to key, trigger logging and warning messages if needed'''
        if env.verbosity > 2:
            self._log(key, value)
        if env.run_mode == 'dryrun':
            self._warn(key, value)
        self.set(key, value)

    def _check_readonly(self, key, value):
        if key in self._readonly and value != dict.__getitem__(self, key):
            raise RuntimeError('Variable {} is readonly and cannot be changed from {} to {}.'
                .format(key, dict.__getitem__(self, key), value))

    def _log(self, key, value):
        env.logger.debug('Workflow variable ``{}`` is set to ``{}``'.format(key, shortRepr(value)))

    def _warn(self, key, value):
        if key.isupper() and dict.__contains__(self, key) and dict.__getitem__(self, key) != value:
            env.logger.warning('Changing readonly variable {} from {} to {}'
                .format(key, dict.__getitem__(self, key), value))
        if key.startswith('_') and key not in ('_input', '_output', '_step', '_index', '_depends'):
            env.logger.warning('{}: Variables with leading underscore is reserved for SoS temporary variables.'.format(key))

    def clone_pickleable(self):
        '''Return a copy of the existing dictionary but keep only the ones that are pickleable'''
        # FIXME: not well tested
        return {x:copy.deepcopy(y) for x,y in self.items() if not callable(y) and not isinstance(y, types.ModuleType)}
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
        # local and global dictionaries used by SoS during the
        # execution of SoS workflows
        self.context_stack = []
        self.locals = WorkflowDict()
        self.globals = globals()
        # 
        # maximum number of concurrent jobs
        self.max_jobs = 1
        self.running_jobs = 0

    class ContextStack:
        '''A context stack and pushes existing workflow dict (env.locals)
        to a stack and make the new workflow the current dict. The 
        context will be poped as soon as the with statement ends and/or
        an exception is raised.
        '''
        def __init__(self, environ, new_dict):
            #
            self.environ = environ
            self.new_dict = new_dict

        def __enter__(self):
            # archive old items
            self.environ.context_stack.append(self.environ.locals)
            self.environ.locals = self.new_dict

        def __exit__(self, etype, value, traceback):
            self.environ.locals = self.environ.context_stack.pop()

    def push_context(self, wf_dict):
        return self.ContextStack(self, wf_dict)

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


#
# String formatting
#


from array import array
try:
    from fcntl import ioctl
    import termios
except ImportError:
    pass

def getTermWidth():
    '''Get the width of current user terminal to properly wrap SoS
    output when well-formatted output is required.
    '''
    try:
        h, w = array('h', ioctl(sys.stderr, termios.TIOCGWINSZ, '\0' * 8))[:2]
        return w
    except:
        return 78

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

    def __init__(self, sigil = '${ }'):
        # do not check sigil here because the function will be called quite frequently
        # the sigil will be checked when it is entered in SoS script.
        self.l, self.r = sigil.split(' ')
        self.error_count = 0

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
                        result = eval(expr, env.globals, env.locals)
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
        if isinstance(obj, basestring):
            return obj if fmt is None else self._format(obj, fmt)
        elif isinstance(obj, collections.Iterable):
            # the object might be nested...
            return sep.join([self._repr(x, fmt, sep) for x in obj])
        elif isinstance(obj, collections.Callable):
            raise InterpolationError(repr(obj), 'Cannot interpolate callable object.')
        else:
            return repr(obj) if fmt is None else self._format(obj, fmt)

def interpolate(text, sigil='${ }'):
    '''Evaluate expressions in `text` marked by specified `sigil` using provided
    global and local dictionaries, and replace the expressions with their formatted strings.'''
    return SoS_String(sigil).interpolate(text)


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
                tokval = 'interpolate(' + tokval + ", \'" + sigil + "')"
        # the resusting string is put back to the expression (or statement)
        result.append((toknum, tokval))
    return untokenize(result)

def SoS_eval(expr, sigil='${ }'):
    '''Evaluate an expression after modifying (convert ' ' string to raw string,
    interpolate expressions) strings.'''
    expr = ConvertString(expr, sigil)
    return eval(expr, env.globals, env.locals)

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
                    raise RuntimeError('Failed to find syntax correct group')
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
        exec(stmts, env.globals, env.locals)
        executed += stmts + '\n'
    #
    env.logger.trace('Executed\n{}'.format(executed))
    return executed



#
# Runtime signature
#

def textMD5(text):
    '''Get md5 of a piece of text'''
    m = hashlib.md5()
    m.update(text.encode())
    return m.hexdigest()

def partialMD5(filename):
    '''Calculate partial MD5, basically the first and last 32M
    of the file for large files. This should signicicantly reduce
    the time spent on the creation and comparison of file signature
    when dealing with large bioinformat ics datasets. '''
    filesize = os.path.getsize(filename)
    # calculate md5 for specified file
    md5 = hashlib.md5()
    block_size = 2**20  # buffer of 1M
    try:
        if filesize < 2**26:
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


class RuntimeInfo:
    '''Record run time information related to a number of output files. Right now only the
    .exe_info files are used.
    '''
    def __init__(self, script, input_files=[], output_files=[], dependent_files = [], pid=None, workdir='.'):
        '''Runtime information for specified output files
        output_files:
            intended output file

        pid:
            process id.

        workdir:
            Current working directory.,
        '''
        self.script = script if isinstance(script, basestring) else ''.join(script)
        self.input_files = [input_files] if isinstance(input_files, basestring) else input_files
        self.output_files = [output_files] if isinstance(output_files, basestring) else output_files
        self.dependent_files = dependent_files if isinstance(dependent_files, basestring) else dependent_files
        #
        sig_name = os.path.realpath(os.path.expanduser(self.output_files[0])) + textMD5('{} {} {} {}'.format(script, input_files, output_files, dependent_files))
        #
        if not output_files:
            self.sig_file = None
            self.proc_out = None
            self.proc_err = None
            self.proc_lck = None
            self.proc_info = None
            self.proc_cmd = None
            self.proc_prog = None
            self.proc_done = None
            self.manifest = None
        else:
            #
            # If the output path is outside of the current working directory
            rel_path = os.path.relpath(sig_name, os.path.realpath(workdir))
            # if this file is not relative to cache, use global signature file
            if rel_path.startswith('../'):
                self.sig_file = os.path.join(os.path.expanduser('~/.sos/.runtime'), sig_name.lstrip(os.sep))
            else:
                # if this file is relative to cache, use local directory
                self.sig_file = os.path.join('.sos/.runtime', rel_path)
            # path to file
            sig_path = os.path.split(self.sig_file)[0]
            if not os.path.isdir(sig_path):
                try:
                    os.makedirs(sig_path)
                except Exception as e:
                    raise RuntimeError('Failed to create runtime directory {}: {}'.format(sig_path, e))
            env.logger.trace('Using signature file {} for output {}'.format(self.sig_file, output_files))
            if pid is None:
                self.pid = os.getpid()
            else:
                self.pid = pid
            self.proc_out = '{}.out_{}'.format(self.sig_file, self.pid)
            self.proc_err = '{}.err_{}'.format(self.sig_file, self.pid)
            self.proc_lck = '{}.lck'.format(self.sig_file)
            self.proc_info = '{}.exe_info'.format(self.sig_file)
            self.proc_cmd = '{}.cmd'.format(self.sig_file)
            self.proc_done = '{}.done_{}'.format(self.sig_file, self.pid)
            self.proc_prog = '{}.working_{}'.format(self.sig_file, self.pid)
            self.manifest = '{}.manifest'.format(self.sig_file)

    def write(self):
        '''Write .exe_info file with signature of script, input, output and dependent files.'''
        if not self.proc_info:
            return
        with open(self.proc_info, 'w') as md5:
            md5.write('{}\n'.format(textMD5(self.script)))
            for f in self.input_files + self.output_files + self.dependent_files:
                f = os.path.realpath(os.path.expanduser(f))
                md5.write('{}\t{}\n'.format(f, partialMD5(f)))

    def validate(self):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        if not self.proc_info or not os.path.isfile(self.proc_info):
            return False
        #
        # duplicated files are only tested once.
        all_files = [os.path.realpath(os.path.expanduser(x)) for x in set(self.input_files + self.output_files + self.dependent_files)]
        # file not exist?
        if not all(os.path.isfile(x) for x in all_files):
            return False
        #
        files_checked = {x:False for x in all_files}
        with open(self.proc_info) as md5:
            cmdMD5 = md5.readline().strip()   # command
            if textMD5(self.script) != cmdMD5:
                return False
            for line in md5:
                try:
                    f, m = line.rsplit('\t', 1)
                    f = os.path.realpath(os.path.expanduser(f))
                except Exception as e:
                    env.logger.error('Wrong md5 line {} in {}: {}'.format(line, md5file, e))
                    continue
                if f not in files_checked:
                    env.logger.waring('{} not need to be checked'.format(f))
                    continue
                if partialMD5(f) != m.strip():
                    env.logger.error('MD5 mismatch {}'.format(f))
                    return False
                files_checked[f] = True
        #
        if not all(files_checked.values()):
            env.logger.warning('No MD5 signature for {}'.format(', '.join(x for x,y in files_checked.items() if not y)))
            return False
        return True

    def clear(self, types=['out', 'err', 'done']):
        if self.sig_file is None:
            return
        for filename in sum([glob.glob(self.sig_file + '.{}_*'.format(x)) for x in types], []):
            try:
                os.remove(filename)
            except Exception as e:
                env.logger.warning('Fail to remove {}: {}'.format(filename, e))


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
                raise WildcardError(name)
            return str(value)  # convert anything into a str
        except KeyError as ex:
            if keep_dynamic:
                return "{{{}}}".format(name)
            elif fill_missing:
                return dynamic_fill
            else:
                raise RuntimeError('Wildcard apply error {}'.format(ex))

    return re.sub(_wildcard_regex, format_match, pattern)

