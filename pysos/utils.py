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

import logging
import re

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
    'A singleton object that provides runtime environment for SoS'
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
    logger = property(lambda self: self._logger, _set_logger)
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
env.logger = None


## import os
## import sys
## import glob
## import logging
## import subprocess
## import urllib
## import getpass
## import time
## import tempfile
## import tokenize
## import gzip
## import copy
## import threading
## import re
## import shlex
## import stat
## import signal
## import random
## import shutil
## import hashlib
## import ConfigParser
## from HTMLParser import HTMLParser
## import tarfile
## import binascii
## from collections import namedtuple
## from itertools import chain
## import site_options
## 
## class RuntimeEnvironments(object):
##     # the following make RuntimeEnvironments a singleton class
##     _instance = None
##     def __new__(cls, *args, **kwargs):
##         if not cls._instance:
##             # *args, **kwargs are not passed to avoid
##             # DeprecationWarning: object.__new__() takes no parameters
##             # cls._instance = super(Singleton, cls).__new__(cls, *args, **kwargs) 
##             cls._instance = super(RuntimeEnvironments, cls).__new__(cls) #, *args, **kwargs)
##         return cls._instance
## 
##     def __init__(self):
##         self._search_path = self.persistent_options['search_path'][0]
##         # logger
##         self._logger = None
##     # 
##     # attribute check_update
##     #
##     #def _set_check_update(self, v):
##     #    if v in ['1', True, 'T', 'True', 'Y', 'Yes']:
##     #        self._check_update = True
##     #    else:
##     #        self._check_update = False
##     #
##     #check_update = property(lambda self: self._check_update, _set_check_update)
##     #
##     # attribute term_width
##     def _set_term_width(self, v):
##         try:
##             self._term_width = int(v)
##         except:
##             self._term_width = None
##     #
##     term_width = property(lambda self: self._term_width, _set_term_width)
##     #
##     # attribute logfile_verbosity
##     #
##     def _set_logfile_verbosity(self, v):
##         if v in ['0', '1', '2']:
##             self._logfile_verbosity = v
##     #
##     logfile_verbosity = property(lambda self: self._logfile_verbosity, _set_logfile_verbosity)
##     #
##     #
##     # attribute verbosity
##     #
##     def _set_verbosity(self, v):
##         if v in ['0', '1', '2', '3']:
##             self._verbosity = v
##     #
##     verbosity = property(lambda self: self._verbosity, _set_verbosity)
##        def _set_temp_dir(self, path=None):
##         # user can explicitly set a path ('None' could be saved by a previous version of vtools)
##         if path not in [None, 'None', '']:
##             path = os.path.expanduser(path)
##             if not os.path.isdir(path):
##                 raise ValueError('Temp directory {} does not exist'.format(path))
##             if os.path.isdir(path) and (
##                     (not os.access(path, os.R_OK)) or (not os.access(path, os.W_OK)) or
##                     (os.stat(path).st_mode & stat.S_ISVTX == 512)):
##                 raise ValueError('Cannot set temporary directory to directory {} because '.format(path) + \
##                     'it is not empty or is not writable or deletable. Please clear this directory or use '
##                     'command "vtools admin --set_runtime_option temp_dir=DIR" to set it to another path, '
##                     'or a random path (empty DIR).')
##             self._temp_dir = path
##             # create a random subdirectory in this directory
##             while True:
##                 subdir = os.path.join(path, '_tmp_{}'.format(random.randint(1, 1000000)))
##                 if not os.path.isdir(subdir):
##                     if self._proj_temp_dir is not None and os.path.isdir(self._proj_temp_dir):
##                         try:
##                             shutil.rmtree(env._proj_temp_dir)
##                         except:
##                             pass
##                     self._proj_temp_dir = subdir
##                     os.makedirs(subdir)
##                     break
##         else:
##             # the usual case
##             if self._temp_dir is None:
##                 self._proj_temp_dir = tempfile.mkdtemp() 
##             try:
##                 if not os.path.isdir(os.path.expanduser(self._proj_temp_dir)):
##                     os.makedirs(os.path.expanduser(self._proj_temp_dir))
##                 while True:
##                     subdir = os.path.join(self._proj_temp_dir, '_tmp_{}'.format(random.randint(1, 1000000)))
##                     if not os.path.isdir(subdir):
##                         if self._proj_temp_dir is not None and os.path.isdir(self._proj_temp_dir):
##                             try:
##                                 shutil.rmtree(env._proj_temp_dir)
##                             except:
##                                 pass
##                         self._proj_temp_dir = subdir
##                         os.makedirs(subdir)
##                         break
##             except:
##                 sys.stderr.write('Failed to create a temporary directory {}.\n'.format(self._proj_temp_dir))
##                 self._proj_temp_dir = tempfile.mkdtemp()
##     #
##     def _get_temp_dir(self):
##         if self._proj_temp_dir is None:
##             self._set_temp_dir()
##         return os.path.expanduser(self._proj_temp_dir)
##     #
##     temp_dir = property(_get_temp_dir, _set_temp_dir)
##         # attribute search_path
##     def _set_search_path(self, val):
##         if val not in ['None', None]:
##             self._search_path = val
##     #
##     search_path = property(lambda self: self._search_path, _set_search_path)
##     #
##     # user stash
##     def _set_user_stash(self, val):
##         if val not in ['None', None]:
##             self._user_stash = val
##     #
##     user_stash = property(lambda self: self._user_stash, _set_user_stash)
##     #
##     #
##     # attribute logger
##     class ColoredFormatter(logging.Formatter):
##         # A variant of code found at http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored
##         def __init__(self, msg):
##             logging.Formatter.__init__(self, msg)
##             self.LEVEL_COLOR = {
##                 'TRACE': 'DARK_CYAN',
##                 'DEBUG': 'BLUE',
##                 'WARNING': 'PURPLE',
##                 'ERROR': 'RED',
##                 'CRITICAL': 'RED_BG',
##                 }
##             self.COLOR_CODE={
##                 'ENDC':0,  # RESET COLOR
##                 'BOLD':1,
##                 'UNDERLINE':4,
##                 'BLINK':5,
##                 'INVERT':7,
##                 'CONCEALD':8,
##                 'STRIKE':9,
##                 'GREY30':90,
##                 'GREY40':2,
##                 'GREY65':37,
##                 'GREY70':97,
##                 'GREY20_BG':40,
##                 'GREY33_BG':100,
##                 'GREY80_BG':47,
##                 'GREY93_BG':107,
##                 'DARK_RED':31,
##                 'RED':91,
##                 'RED_BG':41,
##                 'LIGHT_RED_BG':101,
##                 'DARK_YELLOW':33,
##                 'YELLOW':93,
##                 'YELLOW_BG':43,
##                 'LIGHT_YELLOW_BG':103,
##                 'DARK_BLUE':34,
##                 'BLUE':94,
##                 'BLUE_BG':44,
##                 'LIGHT_BLUE_BG':104,
##                 'DARK_MAGENTA':35,
##                 'PURPLE':95,
##                 'MAGENTA_BG':45,
##                 'LIGHT_PURPLE_BG':105,
##                 'DARK_CYAN':36,
##                 'AUQA':96,
##                 'CYAN_BG':46,
##                 'LIGHT_AUQA_BG':106,
##                 'DARK_GREEN':32,
##                 'GREEN':92,
##                 'GREEN_BG':42,
##                 'LIGHT_GREEN_BG':102,
##                 'BLACK':30,
##             }
## 
##         def colorstr(self, astr, color):
##             return '\033[{}m{}\033[{}m'.format(self.COLOR_CODE[color], astr,
##                 self.COLOR_CODE['ENDC'])
## 
##         def emphasize(self, msg, in_color):
##             # display text within `` and `` in green
##             # This is done for levelname not in self.LEVEL_COLOR, e.g.
##             # for info that uses native color. The text will not be 
##             # visible if someone is using a green background
##             if in_color == 0:
##                 return re.sub(r'``([^`]*)``', '\033[32m\\1\033[0m', str(msg))
##             else:
##                 return re.sub(r'``([^`]*)``', '\033[32m\\1\033[{}m'.format(self.COLOR_CODE[in_color]), str(msg))
## 
##         def format(self, record):
##             record = copy.copy(record)
##             levelname = record.levelname
##             if levelname in self.LEVEL_COLOR:
##                 record.levelname = self.colorstr(levelname, self.LEVEL_COLOR[levelname])
##                 record.name = self.colorstr(record.name, 'BOLD')
##                 record.msg = self.colorstr(self.emphasize(record.msg,
##                     self.LEVEL_COLOR[levelname]), self.LEVEL_COLOR[levelname])
##             else:
##                 record.msg = self.emphasize(record.msg, 0)
##             return logging.Formatter.format(self, record)
## 
##     def _set_logger(self, logfile=None):
##         # create a logger, but shutdown the previous one
##         if not hasattr(logging, 'TRACE'):
##             logging.TRACE = 5
##             logging.addLevelName(logging.TRACE, "TRACE")
##         #
##         if self._logger is not None:
##             self._logger.handlers = []
##         self._logger = logging.getLogger()
##         self._logger.setLevel(logging.DEBUG)
##         # output to standard output
##         cout = logging.StreamHandler()
##         levels = {
##             '0': logging.WARNING,
##             '1': logging.INFO,
##             '2': logging.DEBUG,
##             '3': logging.TRACE,
##             None: logging.INFO
##         }
##         #
##         cout.setLevel(levels[self._verbosity])
##         cout.setFormatter(self.ColoredFormatter('%(levelname)s: %(message)s'))
##         self._logger.addHandler(cout)
##         self._logger.trace = lambda msg, *args: self._logger._log(logging.TRACE, msg, args)
##         # output to a log file
##         if logfile is not None:
##             ch = logging.FileHandler(logfile.lstrip('>'), mode = ('a' if logfile.startswith('>>') else 'w'))
##             # NOTE: debug informaiton is always written to the log file
##             ch.setLevel(levels[self._logfile_verbosity])
##             ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
##             self._logger.addHandler(ch)
##     #
##     logger = property(lambda self: self._logger, _set_logger)
## 
## # the singleton object of RuntimeEnvironments
## env = RuntimeEnvironments()
## # create a default logger without logging to file, this makes sure a logger
## # will be usable even when a project is failed to create
## env.logger = None
## 
## # the following is copied from shutils.which from Python 3.3
## def which(cmd, mode=os.F_OK | os.X_OK, path=None):
##     """Given a command, mode, and a PATH string, return the path which
##     conforms to the given mode on the PATH, or None if there is no such
##     file.
## 
##     `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
##     of os.environ.get("PATH"), or can be overridden with a custom search
##     path.
## 
##     """
##     # Check that a given file can be accessed with the correct mode.
##     # Additionally check that `file` is not a directory, as on Windows
##     # directories pass the os.access check.
##     def _access_check(fn, mode):
##         return (os.path.exists(fn) and os.access(fn, mode)
##                 and not os.path.isdir(fn))
## 
##     # Short circuit. If we're given a full path which matches the mode
##     # and it exists, we're done here.
##     if _access_check(cmd, mode):
##         return cmd
## 
##     path = (path or os.environ.get("PATH", os.defpath)).split(os.pathsep)
##     files = [cmd]
## 
##     seen = set()
##     for dir in path:
##         dir = os.path.normcase(dir)
##         if not dir in seen:
##             seen.add(dir)
##             for thefile in files:
##                 name = os.path.join(dir, thefile)
##                 if _access_check(name, mode):
##                     return name
##     return None
## 
## #
## # Well, it is not easy to do reliable download
## # 
## def downloadURL(URL, dest, quiet, message=None):
##     # use libcurl? Recommended but not always available
##     if 'VTOOLS_ENV' in os.environ and 'NOWEB' in os.environ['VTOOLS_ENV']:
##          raise RuntimeError('Failed to download from {}: no internet connection (set by NOWEB in VTOOLS_ENV environment variable)'.format(URL))
##     env.logger.trace('Download {}'.format(URL))
##     filename = os.path.split(urlparse.urlsplit(URL).path)[-1]
##     # message during downloading
##     if message is None:
##         message = filename
##     if len(message) > 30:
##         message = message[:10] + '...' + message[-16:]
##     if os.path.isdir(dest):
##         dest = os.path.join(dest, filename)
##     #
##     try:
##         import pycurl
##         if not quiet:
##             prog = ProgressBar(message)
##         dest_tmp = TEMP(dest)
##         with open(dest_tmp, 'wb') as f:
##             c = pycurl.Curl()
##             c.setopt(pycurl.URL, str(URL))
##             c.setopt(pycurl.WRITEFUNCTION, f.write)
##             if not quiet:
##                 c.setopt(pycurl.NOPROGRESS, False)
##                 c.setopt(pycurl.PROGRESSFUNCTION, prog.curlUpdate)
##             c.perform()
##         if not quiet:
##             prog.done()
##         if c.getinfo(pycurl.HTTP_CODE) == 404:
##             try:
##                 os.remove(dest_tmp)
##             except OSError:
##                 pass
##             raise RuntimeError('ERROR 404: Not Found.')
##         os.rename(dest_tmp, dest)
##         if os.path.isfile(dest):
##             return dest
##         else:
##             raise RuntimeError('Failed to download {} using pycurl'.format(URL))
##     except Exception:
##         # no pycurl module, or when download fail
##         pass
##     # use wget? Almost universally available under linux
##     if which('wget'):
##         # for some strange reason, passing wget without shell=True can fail silently.
##         dest_tmp = TEMP(dest)
##         p = subprocess.Popen('wget {} -O {} {}'.format('-q' if quiet else '',
##             dest_tmp, URL), shell=True)
##         ret = p.wait()
##         os.rename(dest_tmp, dest)
##         if ret == 0 and os.path.isfile(dest):
##             return dest
##         else:
##             try:
##                 os.remove(dest_tmp)
##             except OSError:
##                 pass
##             raise RuntimeError('Failed to download {} using wget'.format(URL))
##     # use python urllib?
##     if not quiet:
##         prog = ProgressBar(message)
##     try:
##         urllib.URLopener().open(URL)
##     except IOError as error_code:
##         if error_code[1] == 404:
##             raise RuntimeError('ERROR 404: Not Found.')
##         else:
##             raise RuntimeError('Unknown error has happend: {}'.format(error_code[1]))
##     else:
##         dest_tmp = TEMP(dest)
##         urllib.urlretrieve(URL, dest_tmp, reporthook=None if quiet else prog.urllibUpdate)
##         os.rename(dest_tmp, dest)
##     if not quiet:
##         prog.done()
##     # all methods tried
##     if os.path.isfile(dest):
##         return dest
##     # if all failed
##     raise RuntimeError('Failed to download {}'.format(fileToGet))
## 
## 
## def downloadFile(fileToGet, dest_dir = None, quiet = False, checkUpdate = False,
##     message=None):
##     '''Download file from URL to filename.'''
##     # two special cases. Move files around to avoid re-download these files.
##     if fileToGet == 'reference/hg18.crr':
##         if os.path.isfile(os.path.join(env.local_resource, 'ftp.completegenomics.com/ReferenceFiles/build36.crr')) and \
##             not os.path.isfile(os.path.join(env.local_resource, 'reference/hg18.crr')):
##             shutil.move(os.path.join(env.local_resource, 'ftp.completegenomics.com/ReferenceFiles/build36.crr'),
##                 os.path.join(env.local_resource, 'reference/hg18.crr'))
##     elif fileToGet == 'reference/hg19.crr':
##         if os.path.isfile(os.path.join(env.local_resource, 'ftp.completegenomics.com/ReferenceFiles/build37.crr')) and \
##             not os.path.isfile(os.path.join(env.local_resource, 'reference/hg19.crr')):
##             shutil.move(os.path.join(env.local_resource, 'ftp.completegenomics.com/ReferenceFiles/build37.crr'),
##                 os.path.join(env.local_resource, 'reference/hg19.crr'))
##     #
##     # if a complete URL is given, DO NOT download from variant tools repository
##     # 
##     # Downloaded file will look similar to
##     #
##     # ~/.variant_tools/ftp.completegenomics.com/refgenome/build36.crr
##     # 
##     # unless a specific dest_dir is given. NO md5 check is possible.
##     #
##     # for backward compatibility, remove http://vtools.houstonbioinformatics.org and
##     # use new server
##     if fileToGet.startswith('http://vtools.houstonbioinformatics.org/'):
##         fileToGet = fileToGet[len('http://vtools.houstonbioinformatics.org/'):]
##     #
##     if '://' in fileToGet:
##         filename = os.path.split(urlparse.urlsplit(fileToGet).path)[-1]
##         # get filename from URL
##         local_fileToGet = fileToGet.split('://', 1)[1]
##         # use root local_resource directory if dest_dir is None
##         if dest_dir is not None:
##             dest = os.path.join(dest_dir, filename)
##             if (not checkUpdate) and os.path.isfile(dest):
##                 env.logger.trace('Using existing file {}'.format(dest))
##                 return dest
##         else:
##             # look for the file in local resource directory
##             dest_dir = os.path.join(env.shared_resource, os.path.split(local_fileToGet)[0])
##             dest = os.path.join(env.shared_resource, local_fileToGet)
##             # if the file is there, return it directly
##             if (not checkUpdate) and os.path.isfile(dest):
##                 env.logger.trace('Using existing file {}'.format(dest))
##                 return dest
##             # if the shared resource is not writable, write to local_resource
##             if not os.access(env.shared_resource, os.W_OK):
##                 dest_dir = os.path.join(env.local_resource, os.path.split(local_fileToGet)[0])
##                 dest = os.path.join(env.local_resource, local_fileToGet)
##                 # if exists in local user-specific .variant_tools, return it
##                 if (not checkUpdate) and os.path.isfile(dest):
##                     env.logger.trace('Using existing file {}'.format(dest))
##                     return dest
##         # start to write to it
##         if not os.path.isdir(dest_dir):
##             os.makedirs(dest_dir)
##         #
##         try:
##             env.logger.trace('Downloading {} to {}'.format(fileToGet, dest))
##             return downloadURL(fileToGet, dest, quiet, message)
##         except Exception as e:
##             raise ValueError('Failed to download URL {}: {}'.format(fileToGet, e))
##     # 
##     # otherwise, download from variant tools repository, but first let us check
##     # if the file is in the repository
##     #
##     filename = os.path.split(fileToGet)[-1]
##     local_fileToGet = fileToGet
##     resource = ResourceManager()
##     resource.getLocalManifest()
##     if fileToGet not in resource.manifest:
##         # update to the latest manifest and see if we still 
##         # cannot find the file
##         resource.getRemoteManifest()
##         if fileToGet not in resource.manifest:
##             # look in user stash if avail
##             if env.user_stash is not None:
##                 for us in env.user_stash.split(';'):
##                     if not os.path.isdir(os.path.expanduser(us)):
##                         env.logger.warning('Stash directory ({}) does not exist. Check ~/.variant_tools/user_options.py for details.'
##                             .format(us))
##                     # we usually download file in directories such as annoDB/mydb
##                     # so we allow both structure and unstructured stash directory
##                     for usf in  [os.path.expanduser(os.path.join(us, x)) for x in (fileToGet, os.path.basename(fileToGet))]:
##                         if os.path.isfile(usf):
##                             return usf
##             raise RuntimeError('Failed to download {} because it is not in the variant tools online repository or local stash directories.'.format(fileToGet))
##     #
##     fileSig = resource.manifest[fileToGet]
##     if dest_dir is not None:
##         dest = os.path.join(dest_dir, os.path.split(filename)[-1])
##         # if exists in local user-specific .variant_tools, return it
##         if (not checkUpdate) and os.path.isfile(dest):
##             env.logger.trace('Using existing file {}'.format(dest))
##             if calculateMD5(dest, partial=True) != fileSig[1]:
##                 env.logger.warning('MD5 signature mismatch: {} (signature {}, calculated {})'
##                     .format(fileToGet, fileSig[1], calculateMD5(dest, partial=True)))
##             return dest
##     else:
##         # look for the file in shared resource directory
##         dest_dir = os.path.join(env.shared_resource, os.path.split(local_fileToGet)[0])
##         dest = os.path.join(env.shared_resource, local_fileToGet)
##         # if the file is there, return it directly
##         if (not checkUpdate) and os.path.isfile(dest):
##             env.logger.trace('Using existing file {}'.format(dest))
##             if calculateMD5(dest, partial=True) != fileSig[1]:
##                 env.logger.warning('MD5 signature mismatch: {} (signature {}, calculated {})'
##                     .format(fileToGet, fileSig[1], calculateMD5(dest, partial=True)))
##             return dest
##         # if the share resource is not writable, write to ~/.variant_tools
##         if not os.access(env.shared_resource, os.W_OK):
##             dest_dir = os.path.join(env.local_resource, os.path.split(local_fileToGet)[0])
##             dest = os.path.join(env.local_resource, local_fileToGet)
##             # if exists in local user-specific .variant_tools, return it
##             if (not checkUpdate) and os.path.isfile(dest):
##                 env.logger.trace('Using existing file {}'.format(dest))
##                 if calculateMD5(dest, partial=True) != fileSig[1]:
##                     env.logger.warning('MD5 signature mismatch: {} (signature {}, calculated {})'
##                         .format(fileToGet, fileSig[1], calculateMD5(dest, partial=True)))
##                 return dest
##     #
##     if not os.path.isdir(dest_dir):
##         os.makedirs(dest_dir)
##     # 
##     # if the file is in the repository, try to find a mirror
##     servers = [fileSig[4][2*i] for i in range(len(fileSig[4])//2)]
##     weights = [fileSig[4][2*i+1] for i in range(len(fileSig[4])//2)]
##     #
##     # if there is a local server, use it regardless of weight
##     for server in servers:
##         # if the path is a local file repository, do not check for mirrors
##         if server.startswith('file://') or '://' not in server:
##             source_file = '{}/{}'.format(server, local_fileToGet)
##             if source_file.startswith('file://'):
##                 source_file = source_file[len('file://'):]
##             #
##             if os.path.isfile(source_file):
##                 env.logger.trace('Copying {} to {}'.format(source_file, dest))
##                 shutil.copyfile(source_file, dest)
##                 if calculateMD5(dest, partial=True) != fileSig[1]:
##                     env.logger.warning('MD5 signature mismatch: {}'
##                         .format(dest))
##                 return dest
##             else:
##                 env.logger.warning('Cannot locate {} from a local file server {}.'.format(fileToGet, server))
##     #
##     # no local file server
##     while servers:
##         if len(servers) == 1:
##             idx = 0
##         else:
##             r = random.random() * sum(weights)
##             s = 0
##             for idx in range(len(servers)):
##                 s += weights[idx]
##                 if s >= r:
##                     break
##         #
##         try:
##             env.logger.trace('Download {} from {}'.format(fileToGet, servers[idx]))
##             downloaded = downloadURL('{}/{}'.format(servers[idx], fileToGet), dest, quiet, message)
##             if calculateMD5(downloaded, partial=True) != fileSig[1]:
##                 env.logger.warning('Downloaded file {} is different from remote copy. You might '
##                     'want to remove local file and try again to update your local copy.'
##                     .format(downloaded))
##             return downloaded
##         except Exception as e:
##             env.logger.warning('Failed to download {} from {}: {}'.format(fileToGet, servers[idx], e))
##             servers.pop(idx)
##             weights.pop(idx)
##     # failed to get file
##     raise Exception('Failed to download file {}'.format(fileToGet))
## 
## 
## class FileInfo:
##     def __init__(self, filename):
##         self.filename = filename
##         self._size = None
##         self._md5 = None
##         self._mtime = None
##         self._ctime = None
##         self._first_line = None
## 
##     def save(self):
##         '''Create a .file_info file with information about the original
##         file.'''
##         with open(self.filename + '.file_info', 'w') as info:
##             info.write('{}\n{}\n{}\n{}\n{}\n'.format(
##                 os.path.getsize(self.filename),
##                 calculateMD5(self.filename, partial=True),
##                 os.path.getctime(self.filename),
##                 os.path.getmtime(self.filename),
##                 self._getFirstLine(self.filename)))
## 
##     def load(self):
##         try:
##             with open(self.filename + '.file_info') as info:
##                 self._size = int(info.readline().strip())
##                 self._md5 = info.readline().strip()
##                 self._ctime = float(info.readline().strip())
##                 self._mtime = float(info.readline().strip())
##                 try:
##                     self._first_line = info.readline().strip()
##                 except:
##                     # early version of file_info does not have first line
##                     pass
##         except Exception as e:
##             raise ValueError('Corrupted file info file {}.file_info: {}'.format(self.filename, e))
##     
##     def _getFirstLine(self, filename):
##         try:
##             with openFile(filename) as input:
##                 return input.readline().decode()
##         except:
##             return ''
## 
##     def md5(self):
##         if self._md5 is None:
##             if os.path.isfile(self.filename):
##                 self._md5 = calculateMD5(self.filename, partial=True)
##             else:
##                 self.load()
##         return self._md5
## 
##     def size(self):
##         if self._size is None:
##             if os.path.isfile(self.filename):
##                 self._size = os.path.getsize(self.filename)
##             else:
##                 self.load()
##         return self._size
## 
##     def mtime(self):
##         if self._mtime is None:
##             if os.path.isfile(self.filename):
##                 self._mtime = os.path.getmtime(self.filename)
##             else:
##                 self.load()
##         return self._mtime
## 
##     def firstline(self):
##         if self._first_line is None:
##             if os.path.isfile(self.filename):
##                 self._first_line = self._getFirstLine(self.filename)
##             else:
##                 self.load()
##         return self._first_line
## 
## def existAndNewerThan(ofiles, ifiles, md5file=None, pipeline=None):
##     '''Check if ofiles is newer than ifiles. The oldest timestamp
##     of ofiles and newest timestam of ifiles will be used if 
##     ofiles or ifiles is a list. If a md5file is specified,
##     timestamp will be ignored if md5 signature of all ofiles
##     and ifiles match.'''
##     # if there is no input or output file, ofiles cannot be newer than ifiles.
##     if not ofiles or ifiles == ofiles:
##         return False
##     _ifiles = [ifiles] if not isinstance(ifiles, list) else ifiles
##     _ifiles = [x for x in _ifiles if x != env.null_input]
##     _ofiles = [ofiles] if not isinstance(ofiles, list) else ofiles
##     #
##     # file exist?
##     for ifile in _ifiles:
##         if not (os.path.isfile(ifile) or os.path.isfile(ifile + '.file_info')):
##             raise RuntimeError('Input file {} is not found.'.format(ifile))
##     # out file does not exist
##     if not all([os.path.isfile(x) or os.path.isfile(x + '.file_info') for x in _ofiles]):
##         return False
##     #
##     # compare timestamp of input and output files
##     ifiles_checked = {os.path.realpath(x):False for x in _ifiles}
##     md5matched = []
##     if md5file:
##         nFiles = [0]
##         with open(md5file) as md5:
##             md5.readline()   # command
##             line = md5.readline()
##             if not line.startswith('#Start:'):
##                 env.logger.warning('Invalid exe_info file {}'.format(md5file))
##                 return False
##             for line in md5:
##                 if line.startswith('#'):
##                     if not line.startswith('#End:'):
##                         env.logger.warning('Invalid exe_info file {}'.format(md5file))
##                         return False
##                     nFiles.append(0)
##                     continue
##                 # stdout and stderr are separated from md5 by newlines
##                 if not line.strip():
##                     break
##                 try:
##                     f_raw, s, m = line.split('\t')
##                     f = substituteVars(f_raw, pipeline.VARS, pipeline.GLOBALS)
##                     nFiles[-1] += 1
##                     s = int(s)
##                 except Exception as e:
##                     env.logger.error('Wrong md5 line {} in {}: {}'.format(line, md5file, e))
##                     continue
##                 # we do not check if f is one of _ifiles or _ofiles because presentation
##                 # of files might differ
##                 if not any([os.path.realpath(f) == x for x in ifiles_checked.keys()]):
##                     if not any([os.path.realpath(f) == os.path.realpath(x) for x in _ofiles]):
##                         env.logger.warning('{} in exe_info is not an required input or putput file.'.format(f))
##                 else:
##                     ifiles_checked[os.path.realpath(f)] = True
##                 #
##                 if not (os.path.isfile(f) or os.path.isfile(f + '.file_info')):
##                     env.logger.warning('{} in {} does not exist.'.format(f, md5file))
##                     return False
##                 try:
##                     f_info = FileInfo(f)
##                     if f_info.size() != s:
##                         env.logger.warning(
##                             'Size of existing file differ from recorded file: {}'
##                             .format(f))
##                         return False
##                     if f_info.md5() != m.strip():
##                         env.logger.warning(
##                             'md5 of existing file differ from recorded file: {}'
##                             .format(f))
##                         return False
##                 except Exception as e:
##                     env.logger.warning(e)
##                     return False
##                 md5matched.append(f)
##             #
##             if not all(ifiles_checked.values()):
##                 env.logger.error('Input or dependent file {} is not recorded in exe_info file.'.format(', '.join([x for x,y in ifiles_checked.items() if not y])))
##                 return False
##         if len(nFiles) != 2 or nFiles[1] == 0:
##             env.logger.warning('Corrupted exe_info file {}'.format(md5file))
##             return False    
##     #
##     def samefile(x,y):
##         if x == y:
##             return True
##         if os.path.isfile(x):
##             if os.path.isfile(y):
##                 return os.path.samefile(x, y)
##             elif os.path.isfile(y + '.file_info'):
##                 return True
##             return False
##         else:
##             if os.path.isfile(y):
##                 return True
##             else:
##                 return False
##     # check if all files have matching signature, do not check timestamp
##     if all([any([samefile(x, y) for y in md5matched]) for x in _ifiles]) \
##         and all([any([samefile(x, y) for y in md5matched]) for x in _ofiles]):
##         return True
##     if not _ifiles:
##         return True
##     # md5 not available 
##     output_timestamp = min([FileInfo(x).mtime() for x in _ofiles])
##     input_timestamp = max([FileInfo(x).mtime() for x in _ifiles])
##     if output_timestamp < input_timestamp:
##         env.logger.debug('Ignoring older existing output file {}.'
##             .format(', '.join(_ofiles)))
##         return False
##     else:
##         return True
## 
## def physicalMemory():
##     '''Get the amount of physical memory in the system'''
##     # MacOSX?
##     if platform.platform().startswith('Darwin'):
##         # FIXME
##         return None
##     elif platform.platform().startswith('Linux'):
##         try:
##             res = subprocess.check_output('free').decode().split('\n')
##             return int(res[1].split()[1])
##         except Exception as e:
##             return None
## 
## class VariableSubstitutor:
##     def __init__(self, text, asString):
##         self.text = text
##         self.asString = asString
## 
##     def var_expr(self, var):
##         if type(var) == str:
##             # tries to be clever and quote filenames with space
##             if os.path.isfile(var) and ' ' in var:
##                 return "'{}'".format(var)
##             else:
##                 return var
##         elif type(var) == list:
##             if self.asString:
##                 return ' '.join([self.var_expr(x) for x in var])
##             else:
##                 return [self.var_expr(x) for x in var]
##         else:
##             env.logger.debug('Return value of pipeline variable is not string or list of strings: {}'.format(var))
##             return str(var)
## 
##     def _substitute(self, text, PipelineVars, PipelineGlobals):
##         if float(PipelineVars['pipeline_format']) <= 1.0 is None:
##             # for the first version of pipeline specification file, newlines are
##             # replaced with ' '. The newer version (1.0+) keeps newline to faciliate the
##             # inclusion of multi-line scripts etc.
##             text =  ' '.join(text.split())
##         # now, find ${}, excluding simple {}, {{}} etc
##         pieces = re.split('(\${(?:[^{}]|{[^{}]*{[^}]*}[^{}]*}|{[^}]*})*})', text)
##         for idx, piece in enumerate(pieces):
##             if piece.startswith('${') and piece.endswith('}'):
##                 KEY = piece[2:-1].lower()
##                 # if the KEY is in the format of ${VAR}
##                 if re.match('^\s*[\w\d_]+\s*$', KEY):
##                     if KEY.strip() in PipelineVars:
##                         pieces[idx] = self.var_expr(PipelineVars[KEY.strip()])
##                     else:
##                         env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
##                             .format(piece, KEY))
##                     continue
##                 #
##                 # if the KEY is in the format of ${VAR[0]} or ${VAR[2:]}
##                 match = re.match('^([\w\d_]+)\s*((\[[\s\d:-]+\])+)$', KEY)
##                 if match:
##                     KEY_name = match.group(1)
##                     KEY_index = match.group(2)
##                     if KEY_name in PipelineVars:
##                         VAL = PipelineVars[KEY_name]
##                         # handle index
##                         #
##                         # split index by [][]
##                         for sub_index in re.split('\]\s*\[', KEY_index):
##                             try:
##                                 sub_index = sub_index.strip().lstrip('[').rstrip(']')
##                                 if sub_index.count(':') == 0:
##                                     VAL = VAL[int(sub_index)]
##                                 elif sub_index.count(':') == 1:
##                                     idx1, idx2 = sub_index.split(':')
##                                     idx1 = int(idx1) if idx1.strip() else None
##                                     idx2 = int(idx2) if idx2.strip() else None
##                                     VAL = VAL[idx1:idx2]
##                                 elif sub_index.count(':') == 2:
##                                     idx1, idx2, idx3 = sub_index.split(':')
##                                     idx1 = int(idx1) if idx1.strip() else None
##                                     idx2 = int(idx2) if idx2.strip() else None
##                                     idx3 = int(idx3) if idx3.strip() else None
##                                     VAL = VAL[idx1:idx2:idx3]
##                                 else:
##                                     raise ValueError('Invalid index string {}'.format(KEY_index))
##                             except Exception as e:
##                                 env.logger.warning("Failed to interpret {} as a pipeline varialbe: {}"
##                                     .format(piece, e))
##                         pieces[idx] = self.var_expr(VAL)
##                     else:
##                         env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
##                             .format(piece, KEY))
##                     continue
##                 #
##                 # now, the lambda function form
##                 #
##                 if ':' in KEY:
##                     # a lambda function?
##                     try:
##                         FUNC = eval('lambda {}'.format(piece[2:-1].replace('\\\n', '').replace('\n', ' ')))
##                     except Exception as e:
##                         env.logger.warning('Failed to interpret {} as a pipeline variable: {}'
##                             .format(piece, e))
##                         continue
##                     KEY = KEY.split(':', 1)[0].strip()
##                     try:
##                         # allow the use of user defined functions in expressions
##                         globals().update(PipelineGlobals)
##                         if not KEY:
##                             # if there is no KEY, this is a lamba function without parameter
##                             pieces[idx] = self.var_expr(FUNC())
##                         elif ',' not in KEY:
##                             # single varialbe
##                             if KEY in PipelineVars:
##                                 VAL = PipelineVars[KEY]
##                             else:
##                                 env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
##                                     .format(piece, KEY))
##                                 continue
##                             pieces[idx] = self.var_expr(FUNC(VAL))
##                         else:
##                             # several parameters
##                             KEYS = KEY.split(',')
##                             VAL = []
##                             for KEY in KEYS:
##                                 # single varialbe
##                                 if KEY in PipelineVars:
##                                     VAL.append(PipelineVars[KEY])
##                                 else:
##                                     env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
##                                         .format(piece, KEY))
##                                     continue
##                             pieces[idx] = self.var_expr(FUNC(*VAL))
##                     except Exception as e:
##                         env.logger.warning('Failed to interpret {} as a pipeline variable: {}'
##                             .format(piece, e))
##                         continue
##                 else:
##                     env.logger.warning('Failed to interpret {} as a pipeline variable'
##                             .format(piece))
##                     continue
##         #
##         if self.asString:
##             if float(PipelineVars['pipeline_format']) <= 1.0:
##                 # now, join the pieces together, but remove all newlines
##                 return ' '.join(''.join(pieces).split())
##             else:
##                 return ''.join(pieces)
##         else:
##             pieces = [x for x in pieces if x]
##             if not pieces:
##                 return ''
##             elif len(pieces) == 1:
##                 return pieces[0]
##             else:
##                 if all([isinstance(x, str) or len(x) <= 1 for x in pieces]):
##                     return ''.join([x if isinstance(x, str) else (x[0] if len(x) == 1 else '') for x in pieces])
##                 else:
##                     raise ValueError('Variables must be string type (or list of length 1) in variable assignment if text and expressions are mixed: {} evalulated as {}'
##                         .format(text, pieces))
## 
## 
##     def substituteWith(self, PipelineVars, PipelineGlobals):
##         count = 1
##         while count < 10:
##             if isinstance(self.text, str):
##                 new_text = self._substitute(self.text, PipelineVars, PipelineGlobals)
##             else:
##                 new_text = [self._substitute(x, PipelineVars, PipelineGlobals) for x in self.text]
##             if new_text == self.text:
##                 return new_text
##             else:
##                 self.text = new_text
##             count += 1
##         raise ValueError('Failed to evaluate pipeline varialbe {}. Perhpas the variable is nested.'.format(self.text))
##     
## def substituteVars(text, PipelineVars, PipelineGlobals, asString=True):
##     # if asString is to, the return value is forced to be string
##     # Otherwise, the evaluate values are returned.
##     return VariableSubstitutor(text, asString).substituteWith(PipelineVars, PipelineGlobals)
## 
## 
## class RuntimeFiles: 
##     def __init__(self, output_files=[], pid=None):
##         if not output_files:
##             self.sig_file = None
##             self.proc_out = None
##             self.proc_err = None
##             self.proc_lck = None
##             self.proc_info = None
##             self.proc_cmd = None
##             self.proc_prog = None
##             self.proc_done = None
##             self.manifest = None
##         else:
##             if isinstance(output_files, list):
##                 output_file = output_files[0]
##             elif isinstance(output_files, str):
##                 output_file = output_files
##             else:
##                 raise ValueError('Invalid output file specification: {}'.format(output_files))
##             #
##             # what is the relative 
##             # The parental directory of cache?
##             cache_parent = os.path.dirname(env.cache_dir.rstrip(os.sep))
##             #
##             # is the file relative to this cache_parent?
##             rel_path = os.path.relpath(os.path.realpath(os.path.expanduser(output_file)), cache_parent)
##             # if this file is not relative to cache, use global signature file
##             if rel_path.startswith('../'):
##                 self.sig_file = os.path.join(env.local_resource, '.runtime', os.path.realpath(os.path.expanduser(output_file)).lstrip(os.sep))
##             else:
##                 # if this file is relative to cache, use cache to store signature file
##                 self.sig_file = os.path.join(env.cache_dir, '.runtime', rel_path)
##             # path to file
##             sig_path = os.path.split(self.sig_file)[0]
##             if not os.path.isdir(sig_path):
##                 try:
##                     os.makedirs(sig_path)
##                 except Exception as e:
##                     raise RuntimeError('Failed to create runtime directory {}: {}'.format(sig_path, e))
##             env.logger.trace('Using signature file {} for output {}'.format(self.sig_file, output_file))
##             if pid is None:
##                 self.pid = os.getpid()
##             else:
##                 self.pid = pid
##             self.proc_out = '{}.out_{}'.format(self.sig_file, self.pid)
##             self.proc_err = '{}.err_{}'.format(self.sig_file, self.pid)
##             self.proc_lck = '{}.lck'.format(self.sig_file)
##             self.proc_info = '{}.exe_info'.format(self.sig_file)
##             self.proc_cmd = '{}.cmd'.format(self.sig_file)
##             self.proc_done = '{}.done_{}'.format(self.sig_file, self.pid)
##             self.proc_prog = '{}.working_{}'.format(self.sig_file, self.pid)
##             self.manifest = '{}.manifest'.format(self.sig_file)
##             #
##             # now if there is an old signature file, let us move it to the new location
##             if os.path.isfile('{}.exe_info'.format(output_file)):
##                 env.logger.info('Moving {}.exe_info from older version of variant tools to local cache'.format(output_file))
##                 shutil.move('{}.exe_info'.format(output_file), self.proc_info)
##             
##     def clear(self, types=['out', 'err', 'done']):
##         if self.sig_file is None:
##             return
##         for filename in sum([glob.glob(self.sig_file + '.{}_*'.format(x)) for x in types], []):
##             try:
##                 os.remove(filename)
##             except Exception as e:
##                 env.logger.warning('Fail to remove {}: {}'.format(filename, e))
## 
## 
