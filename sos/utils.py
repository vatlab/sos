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
import time
import os
import sys
import re
import copy
import types
import logging
import math
import collections
import traceback
import threading
import base64
import pickle
import yaml
import urllib
import urllib.parse
import urllib.request
import argparse
from tqdm import tqdm as ProgressBar
from collections.abc import Sequence
from io import StringIO, FileIO

from html.parser import HTMLParser
import fasteners

__all__ = ['logger', 'get_output']


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
        if sys.platform == 'win32':
            return astr
        else:
            return '\033[{}m{}\033[0m'.format(color, astr)

    def emphasize(self, msg, level_color=0):
        # display text within `` and `` in green
        if sys.platform == 'win32':
            return str(msg).replace('``', '')
        elif level_color == 0:
            return re.sub(r'``([^`]*)``', '\033[32m\\1\033[0m', str(msg))
        else:
            return re.sub(r'``([^`]*)``', '\033[0m\033[32m\\1\033[0m\033[{}m'.format(level_color), str(msg))

    def format(self, record):
        level_name = record.levelname
        if level_name in self.LEVEL_COLOR:
            level_color = self.COLOR_CODE[self.LEVEL_COLOR[level_name]]
            record.color_levelname = self.colorstr(level_name, level_color)
            record.color_name = self.colorstr(record.name, level_color)
            record.color_msg = self.colorstr(self.emphasize(record.msg, level_color), level_color)
        else:
            # for INFO, use default color
            record.color_levelname = record.levelname
            record.color_msg = self.emphasize(record.msg)
        return logging.Formatter.format(self, record)


def short_repr(obj, noneAsNA=False):
    '''Return a short representation of obj for clarity.'''
    if obj is None:
        return 'unspecified' if noneAsNA else 'None'
    elif isinstance(obj, str) and len(obj) > 50:
        return '{}...{}'.format(obj[:30].replace('\n', '\\n'), obj[-10:].replace('\n', '\\n'))
    elif isinstance(obj, (str, int, float, bool)) or (isinstance(obj, collections.Sequence) \
        and len(obj) <= 2) or len(str(obj)) < 50:
        return repr(obj)
    elif isinstance(obj, collections.Sequence): # should be a list or tuple
        return '[{}, ...] ({} items)'.format(short_repr(obj[0]), len(obj))
    elif isinstance(obj, dict):
        if obj:
            first_key = list(obj.keys())[0]
            return '{{{!r}:{!r}, ...}} ({} items)'.format(first_key, short_repr(obj[first_key]), len(obj))
        else:
            return '{}'
    else:
        return '{}...'.format(repr(obj)[:40])

#
# SoS Workflow dictionary
#
class WorkflowDict(object):
    """A dictionary object that keeps all SoS workflow objects.

    IMPORTANT:

    Python does not allow the passing of a derived class of dict as globals
    to eval or exec. Doing so will result in strange behavior such as __builtins__
    not found. We then have to embed a real dictionary in WorkflowDict instead of
    deriving a dict from it.
    """
    def __init__(self, *args, **kwargs):
        self._dict = dict(*args, **kwargs)

    def set(self, key, value):
        '''A short cut to set value to key without triggering any logging
        or warning message.'''
        self._dict[key] = value

    def quick_update(self, obj):
        '''Update without sanity check etc. For fast internal update'''
        self._dict.update(obj)

    def update(self, obj):
        '''Redefine update to trigger logging message'''
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
        if env.config['run_mode'] == 'prepare':
            self._warn(key, value)
        if key in ('input', 'output', 'depends', '_input', '_output', '_depends', '_runtime'):
            raise ValueError('Variable {} can only be set by SoS'.format(key))
        self.set(key, value)

    def _log(self, key, value):
        env.logger.debug('Set ``{}`` = ``{}``'.format(key, short_repr(value)))

    def _warn(self, key, value):
        if key.startswith('_') and not key.startswith('__') and key not in ('_input', '_output', '_step', '_index', '_depends', '_runtime'):
            env.logger.warning('{}: Variables with leading underscore is reserved for SoS temporary variables.'.format(key))

    def clone_selected_vars(self, selected=[]):
        return {x:copy.deepcopy(y) for x,y in self._dict.items() if (not selected or x in selected) and pickleable(y, x)}

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

    _exec_dir = None
    def reset(self):
        # logger
        self._logger = None
        self._verbosity = 2
        self._logfile = None
        self._set_logger()
        #
        # run mode, this mode controls how SoS actions behave
        #
        self.config = {
            'config_file': None,
            'output_dag': None,
            'wait_for_task': None,
            'resume_mode': False,
            'default_queue': '',
            'max_procs': 4,
            'max_running_jobs': None,
            'sig_mode': 'default',
            'run_mode': 'run',
        }

        #
        # global dictionaries used by SoS during the
        # execution of SoS workflows
        self.sos_dict = WorkflowDict()
        # parameters of the workflow, which will be handled differently
        self.parameter_vars = set()
        #
        # maximum number of concurrent jobs
        self.running_jobs = 0
        # this directory will be used by a lot of processes
        self.exec_dir = os.getcwd()
        #
        os.makedirs(os.path.join(os.path.expanduser('~'), '.sos', 'tasks'), exist_ok=True)
        os.makedirs(os.path.join(os.path.expanduser('~'), '.sos', '.runtime'), exist_ok=True)
        os.makedirs(os.path.join('.sos', '.runtime'), exist_ok=True)

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
        while self._logger.hasHandlers():
            self._logger.removeHandler(self._logger.handlers[0])
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
    # attribute exec_dir
    def _assure_runtime_dir(self, dir):
        os.makedirs(os.path.join(dir, '.sos'), exist_ok=True)

    def _set_exec_dir(self, dir):
        if not os.path.isdir(dir):
            raise RuntimeError('Exec dir {} does not exist.'.format(dir))
        self._assure_runtime_dir(dir)
        self._exec_dir = dir

    def _get_exec_dir(self):
        if self._exec_dir is None:
            raise RuntimeError('Exec dir is not set')
        self._assure_runtime_dir(self._exec_dir)
        return self._exec_dir
   
    exec_dir = property(_get_exec_dir, _set_exec_dir)
        
    #
    # attribute logger
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


class AbortExecution(Error):
    '''Abort a step and continue'''
    def __init__(self, msg):
        Error.__init__(self, msg)

def get_traceback():
    output = StringIO()
    exc_type, exc_value, exc_traceback = sys.exc_info()
    #print "*** print_tb:"
    traceback.print_tb(exc_traceback, limit=1, file=output)
    #print "*** print_exception:"
    try:
        traceback.print_exception(exc_type, exc_value, exc_traceback,
                              limit=5, file=output)
    except:
        # the above function call can fail under Python 3.4 for some
        # exception but we do not really care if that happens
        pass
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


def pickleable(obj, name):
    if isinstance(obj, (str, bool, int, float, complex, bytes)):
        return True
    if isinstance(obj, (types.ModuleType, WorkflowDict)):
        return False
    try:
        pickle.dumps(obj)
        return True
    except Exception as e:
        env.logger.debug('Object {} with value {} is not passed because it is not pickleable: {}'.format(name, short_repr(obj), e))
        return False

class ProgressFileObj(FileIO):
    '''A wrapper of a file object that update a progress bar
    during file read.
    '''
    def __init__(self, prog, *args, **kwargs):
        FileIO.__init__(self, *args, **kwargs)
        self.prog = prog

    def read(self, n, *args):
        self.prog.update(n)
        return FileIO.read(self, n, *args)

def stable_repr(obj):
    if isinstance(obj, str):
        return repr(obj)
    elif isinstance(obj, collections.abc.Mapping):
        items = [stable_repr(k) + ':' + stable_repr(obj[k]) for k in obj.keys()]
        return '{' + ', '.join(sorted(items)) + '}'
    elif isinstance(obj, collections.abc.Set):
        items = [stable_repr(x) for x in obj]
        return '{' + ', '.join(sorted(items)) + '}'
    elif isinstance(obj, collections.abc.Sequence):
        return '[' + ', '.join(stable_repr(k) for k in obj) + ']'
    else:
        return repr(obj)

# frozen dict makes it extremely difficult to transfer CONFIG from
# one process to another...
class frozendict(dict):
    '''A fronzen dictionary that disallow changing of its elements
    Copied from http://code.activestate.com/recipes/414283/
    '''
#    def _blocked_attribute(obj):
#        raise RuntimeError("Cannot modify a readonly dictionary.")
#    _blocked_attribute = property(_blocked_attribute)
#
#    __delitem__ = __setitem__ = clear = _blocked_attribute
#    pop = popitem = setdefault = update = _blocked_attribute
#
#    def __new__(cls, *args):
#        new = dict.__new__(cls)
#        dict.__init__(new, *args)
#        return new
#
    def __init__(self, *args, **kwargs):
        super(frozendict, self).__init__(*args, **kwargs)
#
#    def __hash__(self):
#        try:
#            return self._cached_hash
#        except AttributeError:
#            h = self._cached_hash = hash(tuple(sorted(self.items())))
#            return h
#
    def __getattr__(self, key):
        return dict.__getitem__(self, key)

    def __setattr__(self, key, value):
        raise RuntimeError('Cannot modify a readonly dictionary')

    def __repr__(self):
        return stable_repr(self)

    def __getstate__(self):
        return {x:y for x,y in self.items()}

    def __setstate__(self, other):
        self.update(other)

    def __deepcopy__(self, memo):
        return frozendict(self)

#
# A utility function that returns output of a command
def get_output(cmd, show_command=False, prompt='$ '):
    import subprocess
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, shell=True).decode()
    except subprocess.CalledProcessError as e:
        if e.output.decode():
            env.logger.error(e.output.decode())
        raise RuntimeError(e)
    if show_command:
        return '{}{}\n{}'.format(prompt, cmd, output)
    else:
        return output

#
# search a path and locate script and other files
#
def locate_script(filename, start=''):
    #
    attemp = os.path.abspath(os.path.expanduser(filename))
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
    sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yml')
    if os.path.isfile(sos_config_file):
        try:
            with open(sos_config_file) as config:
                cfg = yaml.safe_load(config)
        except Exception as e:
            raise RuntimeError('Failed to parse global sos config file {}, is it in JSON format?'.format(sos_config_file))
        #
        pathes.extend(cfg.get('sos_path', []))
    #
    sos_config_file = 'config.yml'
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

def text_repr(text):
    """return a valid string representation of text, but requires that
    it is double quoted so that sos can interpolate text in script style
    """
    # in the simple case, we can just use r""" """
    if '"""' not in text and not text.endswith('"'):
        return 'r"""' + text + '"""'
    # if things need to be quoted, let us first use repr
    # to quote them
    r = repr(text)
    # if the result happens to be double quoted, good
    # although it appears to me that Python 3 only use single quote.
    if r.startswith('"'):
        return r
    # otherwise we have to manually change single quote to double quote
    #
    # The problem here is that the representation can have a bunch of \'
    # and I have to hope that \' will be correctly interpreted in " "
    # strings
    return '"' + r.replace('"', r'\"')[1:-1] + '"'

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ int(c) if c.isdigit() else c for c in re.split('(\d+)', text) ]

def transcribe(text, action=None):
    if action is not None:
        text = '{}:\n{}'.format(action, '    ' + text.replace('\n', '\n    ') + '\n')
    with open(os.path.join(env.exec_dir, '.sos', 'transcript.txt'), 'a') as trans:
        trans.write(text)

def dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: None
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], dict)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]

def PrettyRelativeTime(time_diff_secs):
    secs = int(time_diff_secs)
    rec = [
        (secs // (3600*24), 'day'), 
        (secs % (3600*24) // 3600, 'hr'),
        (secs % 3600 // 60, 'min'),
        (secs % 60, 'sec')]
    txt = ' '.join(['{} {}'.format(x,y) for x,y in rec if x > 0])
    return txt if txt else '0 sec'

# display file size in K, M, G etc automatically. Code copied from
# http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
# for its compact size
def pretty_size(n,pow=0,b=1024,u='B',pre=['']+[p+'i'for p in'KMGTPEZY']):
    pow,n=min(int(math.log(max(n*b**pow,1),b)),len(pre)-1),n*b**pow
    return "%%.%if %%s%%s"%abs(pow%(-pow-1))%(n/b**float(pow),pre[pow],u)

def expand_size(size):
    m = re.match(r'\s*([+-]?)([\.\d]*)\s*(\S+)\s*', size)
    if not m:
        raise ValueError('Invalid size specified: {}'.format(size))
    sign, num, unit = m.groups()
    sign = -1 if sign == '-' else 1
    if not unit:
        return sign * int(num)
    if not num:
        num = 1
    s = {x + 'I' :1024**(idx+1) for idx,x in enumerate('KMGTPEZY')}
    s.update({x: 1000**(idx+1) for idx,x in enumerate('KMGTPEZY')})
    unit = unit[:-1].upper() if unit[-1].upper().endswith('B') else unit.upper()
    if unit not in s:
        raise ValueError('Invalid size specified: {}'.format(size))
    return sign * int(float(num) * s[unit])

def find_symbolic_links(item):
    item = os.path.expanduser(item)
    if os.path.isfile(item):
        return {}
    elif os.path.islink(item):
        if not os.path.exists(item):
            env.logger.warning('Non-existent symbolic link {}'.format(item))
        return {item: os.path.realpath(item)}
    else:
        result = {}
        for x in os.listdir(item):
            result.update(find_symbolic_links(os.path.join(item, x)))
        return result

class ActivityNotifier(threading.Thread):
    def __init__(self, msg, delay=5):
        threading.Thread.__init__(self)
        self.msg = msg
        self.start_time = time.time()
        self.delay = delay
        self.event = threading.Event()
        self.start()

    def run(self):
        prog = None
        while True:
            self.event.wait(self.delay)
            if self.event.is_set():
                if prog:
                    prog.close()
                break
            if not prog:
                prog = ProgressBar(desc=self.msg, position=0, bar_format='{desc}', total=100000000)
            second_elapsed = time.time() - self.start_time
            prog.set_description(self.msg + ' ({}{})'.format(
                    '' if second_elapsed < 86400 else '{} day{} '
                    .format(int(second_elapsed/86400), 's' if second_elapsed > 172800 else ''),
                    time.strftime('%H:%M:%S', time.gmtime(second_elapsed)) ))
            prog.update(1)

    def stop(self):
        self.event.set()

class DelayedAction:
    '''Call the passed function with param after a few seconds. It is most often
    used to display certain message only if an action takes a long time.

        action = delayedAction(env.logger.info, 'This might take a while', 5)
        some_action_that_might_take_a_while
        del action

    if the action finishes very quick, the message will not be displayed.
    '''
    def __init__(self, func, param, delay=5):
        self.timer = threading.Timer(delay, func, (param,))
        self.timer.start()

    def __del__(self):
        self.timer.cancel()

class ArgumentError(Error):
    """Raised when an invalid argument is passed."""
    def __init__(self, msg):
        Error.__init__(self, msg)
        self.args = (msg, )

def _parse_error(msg):
    '''This function will replace error() function in argparse module so that SoS
    can hijack errors raised from it.'''
    raise ArgumentError(msg)

def sos_handle_parameter_(key, defvalue):
    '''Parse command line arguments and set values to parameters section.
    NOTE: parmeters will not be handled if it is already defined in
    the environment. This makes the parameters variable.
    '''
    if key in env.sos_dict['sos_symbols_']:
        env.logger.warning('Parameter {} overrides a SoS function.'.format(key))

    env.parameter_vars.add(key)
    if not env.sos_dict['__args__']:
        if isinstance(defvalue, type):
            raise ArgumentError('Argument {} of type {} is required'.format(key, defvalue.__name__))
        return defvalue
    # if the parameter is passed from action sos_run
    if isinstance(env.sos_dict['__args__'], dict):
        if key in env.sos_dict['__args__']:
            return env.sos_dict['__args__'][key]
    #
    parser = argparse.ArgumentParser()
    # thre is a possibility that users specify --cut-off instead of --cut_off for parameter
    # cut_off. It owuld be nice to allow both.
    #
    # Argparse would produce cut_off for both definition of --cut-off and --cut_off, however
    # you can only use the matching input...

    if isinstance(defvalue, type) or defvalue is None:
        if defvalue == bool:
            feature_parser = parser.add_mutually_exclusive_group(required=True)
            feature_parser.add_argument('--{}'.format(key), dest=key, action='store_true')
            feature_parser.add_argument('--no-{}'.format(key), dest=key, action='store_false')
            if '_' in key:
                feature_parser.add_argument('--{}'.format(key.replace('_', '-')), dest=key, action='store_true')
                feature_parser.add_argument('--no-{}'.format(key.replace('_', '-')), dest=key, action='store_false')
        else:
            if defvalue is None:
                defvalue = str
            # if only a type is specified, it is a required document of required type
            if '_' in key:
                feature_parser = parser.add_mutually_exclusive_group(required=True)
                feature_parser.add_argument('--{}'.format(key), dest=key, type=str if hasattr(defvalue, '__iter__') else defvalue,
                    help='', nargs='+' if defvalue != str and hasattr(defvalue, '__iter__') else '?')
                feature_parser.add_argument('--{}'.format(key.replace('_', '-')), dest=key, type=str if hasattr(defvalue, '__iter__') else defvalue,
                    help='', nargs='+' if defvalue != str and hasattr(defvalue, '__iter__') else '?')
            else:
                parser.add_argument('--{}'.format(key), dest=key, type=str if hasattr(defvalue, '__iter__') else defvalue,
                    help='', required=True, nargs='+' if defvalue != str and hasattr(defvalue, '__iter__') else '?')
    else:
        if isinstance(defvalue, bool):
            feature_parser = parser.add_mutually_exclusive_group(required=False)
            feature_parser.add_argument('--{}'.format(key), dest=key, action='store_true')
            feature_parser.add_argument('--no-{}'.format(key), dest=key, action='store_false')
            if '_' in key:
                feature_parser.add_argument('--{}'.format(key.replace('_', '-')), dest=key, action='store_true')
                feature_parser.add_argument('--no-{}'.format(key.replace('_', '-')), dest=key, action='store_false')
            feature_parser.set_defaults(key=defvalue)
        else:
            if isinstance(defvalue, str):
                deftype = str
            elif isinstance(defvalue, Sequence):
                if len(defvalue) > 0:
                    deftype = type(defvalue[0])
                else:
                    deftype = str
            else:
                deftype = type(defvalue)
            if '_' in key:
                feature_parser = parser.add_mutually_exclusive_group(required=False)
                feature_parser.add_argument('--{}'.format(key), dest=key, type=deftype,
                    nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, str) else '?',
                    default=defvalue)
                feature_parser.add_argument('--{}'.format(key.replace('_', '-')), dest=key, type=deftype,
                    nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, str) else '?',
                    default=defvalue)
            else:
                parser.add_argument('--{}'.format(key), dest=key, type=deftype,
                    nargs='*' if isinstance(defvalue, Sequence) and not isinstance(defvalue, str) else '?',
                    default=defvalue)
    #
    parser.error = _parse_error
    parsed, unknown = parser.parse_known_args(env.sos_dict['__args__']['__args__'] if isinstance(env.sos_dict['__args__'], dict) else env.sos_dict['__args__'])
    return vars(parsed)[key]


def load_config_files(filename=None):
    cfg = {}
    # global site file
    sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'hosts.yml')
    if os.path.isfile(sos_config_file):
        with fasteners.InterProcessLock('/tmp/sos_config_'):
            try:
                with open(sos_config_file) as config:
                    cfg = yaml.safe_load(config)
            except Exception as e:
                raise RuntimeError('Failed to parse global sos hosts file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
    else:
        env.logger.debug("{} does not exist".format(sos_config_file))
    # global config file
    sos_config_file = os.path.join(os.path.expanduser('~'), '.sos', 'config.yml')
    if os.path.isfile(sos_config_file):
        with fasteners.InterProcessLock('/tmp/sos_config_'):
            try:
                with open(sos_config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse global sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
    else:
        env.logger.debug("{} does not exist".format(sos_config_file))
    # local config file
    sos_config_file = 'config.yml'
    if os.path.isfile(sos_config_file):
        with fasteners.InterProcessLock('/tmp/sos_config_'):
            try:
                with open(sos_config_file) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse local sos config file {}, is it in YAML/JSON format? ({})'.format(sos_config_file, e))
    else:
        env.logger.debug("{} does not exist".format(sos_config_file))
    # user-specified configuration file.
    if filename is not None:
        if not os.path.isfile(filename):
            raise RuntimeError('Config file {} not found'.format(filename))
        with fasteners.InterProcessLock('/tmp/sos_config_'):
            try:
                with open(filename) as config:
                    dict_merge(cfg, yaml.safe_load(config))
            except Exception as e:
                raise RuntimeError('Failed to parse config file {}, is it in YAML/JSON format? ({})'.format(filename, e))
    return cfg
 
def format_HHMMSS(v):
    if isinstance(v, int):
        h, m, s = v // 3600, v %3600//60, v % 60
    elif isinstance(v, str):
        # the time can be spacified as age. 
        try:
            return format_HHMMSS(expand_time(v))
        except Exception as e:
            raise ValueError('walltime should be specified as a integer with unit s (default), h, m, d or string in the format of HH:MM:SS. "{}" specified ({})'.format(v, e))
    else:
        raise ValueError('walltime should be specified as a integer with unit s (default), h, m, d or string in the format of HH:MM:SS. "{}" of type {} specified'.format(v, v.__class__.__name__))
    return '{:02d}:{:02d}:{:02d}'.format(h,m,s)


def expand_time(v, default_unit='s'):
    # expand walltime from '00:00:00' format to second
    if isinstance(v, str):
        try:
            sign = {'+': 1, '-': -1}[v[1]]
            v = v[1:]
        except:
            # if there is no sign, assume +
            sign = 1

        # format HHMMSS
        if ':' in v:
            try:
                h, m, s = map(int, v.split(':'))
                return sign * (h * 60 * 60 + m * 60 + s)
            except Exception as e:
                raise ValueError('Input of option walltime should be an integer with unit s (default), h, m, d or a string in the format of HH:MM:SS. {} specified: {}'.format(v, e))
        #
        try:
            unit = {'s': 1, 'm': 60, 'h': 3600, 'd': 3600*24}[v[-1]]
            v = v[:-1]
        except:
            unit =  {'s': 1, 'm': 60, 'h': 3600, 'd': 3600*24}[default_unit]
        #
        try:
            return sign * unit * int(v)
        except:
            raise ValueError('Unacceptable time for parameter age, expecting [+/-] num [s|m|h|d] or HH:MM:SS (e.g. +5h): {} provided'.format(v))
    elif isinstance(v, int):
        return v
    else:
         raise ValueError('Input of option walltime should be an integer with unit s (default), h, m, d or a string in the format of HH:MM:SS. {} specified.'.format(v))



def tail_of_file(filename, n, offset=None, ansi2html=False):
    """Reads a n lines from f with an offset of offset lines.  The return
    value is a tuple in the form ``(lines, has_more)`` where `has_more` is
    an indicator that is `True` if there are more lines in the file.
    """
    avg_line_length = 74
    to_read = n + (offset or 0)

    with open(filename) as f:
        while 1:
            try:
                f.seek(-(avg_line_length * to_read), 2)
            except IOError:
                # woops.  apparently file is smaller than what we want
                # to step back, go to the beginning instead
                f.seek(0)
            pos = f.tell()
            lines = f.read().splitlines()
            if len(lines) >= to_read or pos == 0:
                if ansi2html:
                    return convertAnsi2html('\n'.join(lines[-to_read:offset and -offset or None]))
                else:
                    return '\n'.join(lines[-to_read:offset and -offset or None])
            avg_line_length *= 1.3

def sample_of_file(filename, n):
    '''Draw a sample of n lines from filename, largely evenly.'''
    with open(filename) as f:
        lines = f.readlines()
        if len(lines) <= n:
            return ''.join(lines)
        else:
            m = len(lines)
            return ''.join([lines[x*m//n + m//(2*n)] for x in range(n)])


def linecount_of_file(filename):
    f = open(filename, 'rb')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.raw.read

    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)

    return lines

def isPrimitive(obj):
    # test if object is of primitive types (string, number, sequence etc
    # http://stackoverflow.com/questions/6391694/check-if-a-variables-type-is-primitive
    return not hasattr(obj, '__dict__')


def save_var(name, var):
    if isPrimitive(var):
        return '{}={}\n'.format(name, stable_repr(var))
    else:
         # for more complex type, we use pickle + base64
         return '{}:={}\n'.format(name, base64.b64encode(pickle.dumps(var)))

def load_var(line):
    from .target import remote
    assert remote
    key, value = line.split('=', 1)
    if key.endswith(':'):
        return key[:-1], pickle.loads(base64.b64decode(eval(value.strip())))
    else:
        return key, eval(value.strip())


def version_info(module):
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

def loaded_modules(namespace={}):
    from types import ModuleType
    res = {}
    for key,value in namespace.items():
        if isinstance(value, ModuleType):
            res[value.__name__] = version_info(value.__name__)
    return [(x,y) for x,y in res.items() if y != 'na']

def convertAnsi2html(txt):
    # 94 is blue, debug
    # 32 is darkgreen, emphasize
    # 95 is purple, warning
    # 91 is red, error
    # 36 is cray, trace
    return txt.replace('\033[94m', '<font color="">'). \
        replace('\033[32m', '<font color="DarkGreen">'). \
        replace('\033[36m', '<font color="cyan">'). \
        replace('\033[95m', '<font color="purple">'). \
        replace('\033[91m', '<font color="red">'). \
        replace('\033[0m', '</font>'). \
        replace('\n', '<br>')
