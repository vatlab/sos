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
import shlex
import shutil
import fasteners
from copy import deepcopy
import pkg_resources
from shlex import quote
import subprocess
from pathlib import Path, WindowsPath, PosixPath
try:
    from xxhash import xxh64 as hash_md5
except ImportError:
    from hashlib import md5 as hash_md5

from collections.abc import Sequence, Iterable

from .utils import env, Error, short_repr, stable_repr, save_var, load_var, isPrimitive, TimeoutInterProcessLock
from .eval import Undetermined

__all__ = ['dynamic', 'executable', 'env_variable', 'sos_variable']

class UnknownTarget(Error):
    def __init__(self, target):
        Error.__init__(self, 'Target unavailable: %s' % target)
        self.target = target

class RemovedTarget(Error):
    def __init__(self, target):
        Error.__init__(self, 'Target removed: %s' % target)
        self.target = target

class UnavailableLock(Error):
    """Raised when there are errors in prepare mode. Such errors are not raised
    immediately, but will be collected and raised at the end """

    def __init__(self, signature):
        Error.__init__(self, f'Failed to obtain lock {signature[2]} for input {short_repr(signature[0])} and output {short_repr(signature[1])}. It is likely ' +
            'that these files are protected by another SoS process or concurrant task that is generating the same set of files. Please manually remove the lockfile ' +
            'if you are certain that no other process is using the lock.')
        self.input = signature[0]
        self.output = signature[1]
        self.sig_file = signature[2]

#
# Runtime signature
#
def textMD5(text):
    '''Get md5 of a piece of text'''
    m = hash_md5()
    if isinstance(text, str):
        m.update(text.encode())
    else:
        m.update(text)
    return m.hexdigest()

def fileMD5(filename, partial=True):
    '''Calculate partial MD5, basically the first and last 8M
    of the file for large files. This should signicicantly reduce
    the time spent on the creation and comparison of file signature
    when dealing with large bioinformat ics datasets. '''
    filesize = os.path.getsize(filename)
    # calculate md5 for specified file
    md5 = hash_md5()
    block_size = 2**20  # buffer of 1M
    try:
        # 2**24 = 16M
        if (not partial) or filesize < 2**24:
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    if not data:
                        break
                    md5.update(data)
        else:
            count = 16
            # otherwise, use the first and last 32M
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    count -= 1
                    if count == 8:
                        # 2**23 = 8M
                        f.seek(-2**23, 2)
                    if not data or count == 0:
                        break
                    md5.update(data)
    except IOError as e:
        sys.exit(f'Failed to read {filename}: {e}')
    return md5.hexdigest()


class BaseTarget:
    '''A base class for all targets (e.g. a file)'''
    def __init__(self, *args):
        self._sigfile = None

    def target_exists(self, mode='any'):
        # mode should be 'any', 'target', or 'signature'
        raise RuntimeError('Undefined base function')

    def target_name(self):
        # name of the target, which should be able to differentiate
        # this object with other targets of the same type.
        raise RuntimeError('Undefined base function')

    def target_signature(self, mode='any'):
        # signature of the content of the target, which should be
        # able to detect changes of the content of target
        #
        # if mode == 'target', the target has to exist and the signature
        # has to be calculated. Otherwise you can return cached signature
        raise RuntimeError('Undefined base function')

    # -----------------------------------------------------
    # derived functions that do not need to be redefined
    #
    def sig_file(self):
        if self._sigfile is None:
            self._sigfile = os.path.join(env.exec_dir, '.sos', '.runtime',
                                         f'{self.__class__.__name__}_{textMD5(self.target_name())}.file_info')
        return self._sigfile

    def remove_sig(self):
        if self.sig_file() and os.path.isfile(self.sig_file()):
            os.remove(self.sig_file())

    def write_sig(self):
        '''Write .sig file with signature'''
        # path to file
        with open(self.sig_file(), 'w') as sig:
            sig.write(f'{self.target_name()}\t{self.target_signature()}\n')

    def __repr__(self):
        return f'{self.__class__.__name__}("{self.target_name()}")'

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, self.__class__) and self.target_signature() == obj.target_signature()

class sos_variable(BaseTarget):
    '''A target for a SoS variable.'''
    def __init__(self, var):
        super(sos_variable, self).__init__()
        self._var = var

    def target_exists(self, mode='any'):
        return self._var in env.sos_dict

    def target_name(self):
        return self._var

    def target_signature(self, mode='any'):
        return textMD5(self._var)

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return isinstance(other, sos_variable) and self._var == other._var

    def __format__(self, format_spec):
        # handling special !q conversion flag
        if format_spec and format_spec[0] == 'R':
            return self._var.__format__(format_spec[1:])
        else:
            return str(self).__format__(format_spec)


class env_variable(BaseTarget):
    '''A target for an environmental variable.'''
    def __init__(self, var):
        super(env_variable, self).__init__()
        self._var = var

    def target_exists(self, mode='any'):
        return self._var in os.environ

    def target_name(self):
        return self._var

    def target_signature(self, mode='any'):
        return textMD5(repr(os.environ[self._var]))

    def __eq__(self, other):
        return isinstance(other, env_variable) and self._var == other._var

    def __hash__(self):
        return hash(repr(self))

    def __format__(self, format_spec):
        # handling special !q conversion flag
        if format_spec and format_spec[0] == 'R':
            return self._var.__format__(format_spec[1:])
        else:
            return str(self).__format__(format_spec)


class sos_step(BaseTarget):
    '''A target for a step of sos.'''
    def __init__(self, step_name):
        super(sos_step, self).__init__()
        self._step_name = 'default_' + step_name if step_name.isdigit() else step_name

    def target_exists(self, mode='any'):
        # the target exists only if it has been executed?
        # which is indicated by a variable
        return '__completed__' in env.sos_dict and self._step_name in env.sos_dict['__completed__']

    def target_name(self):
        return self._step_name

    def target_signature(self, mode='any'):
        return textMD5(f'sos_step({self._step_name})')

    def write_sig(self):
        pass

    def __eq__(self, other):
        return isinstance(other, sos_step) and self._step_name == other._step_name

    def __hash__(self):
        return hash(repr(self))

    def __format__(self, format_spec):
        # handling special !q conversion flag
        if format_spec and format_spec[0] == 'R':
            return self._step_name.__format__(format_spec[1:])
        else:
            return str(self).__format__(format_spec)


# class bundle(BaseTarget):
#     '''a bundle of other targets'''
#     def __init__(self, *args):
#         super(bundle, self).__init__()
#         self._targets = []
#         for arg in args:
#             if isinstance(arg, str):
#                 self._targets.append(os.path.expanduser(arg))
#             elif isinstance(arg, Iterable):
#                 # in case arg is a Generator, check its type will exhaust it
#                 arg = list(arg)
#                 if not all(isinstance(x, str) for x in arg):
#                     raise RuntimeError('Invalid filename: {}'.format(arg))
#                 self._targets.extend(arg)
#             else:
#                 raise RuntimeError('Unrecognized file: {}'.format(arg))
#
#     def target_exists(self, mode='any'):
#         return all(file_target(x).target_exists(mode) for x in self._targets)
#
#     def target_name(self):
#         return repr(self._targets)
#
#     def target_signature(self, mode='any'):
#         return textMD5(repr(self._targets))
#
class dynamic(BaseTarget):
    '''A dynamic executable that only handles input files when
    it is available. This target is handled directly with its `resolve`
    function called by the executor. '''
    def __init__(self, target):
        self._target = target

    def target_name(self):
        return self._target

    def resolve(self):
        return self._target

    def __format__(self, format_spec):
        # handling special !q conversion flag
        if format_spec and format_spec[0] == 'R':
            return sos_targets(self._target).__format__(format_spec[1:])
        else:
            return str(self).__format__(format_spec)


class remote(BaseTarget):
    '''A remote target is not tracked and not translated during task execution'''
    def __init__(self, *targets):
        super(remote, self).__init__()
        self.__unresolvable_object__ = True
        if len(targets) == 1:
            self._target = targets[0]
        else:
            # multi-item targets
            self._target = targets
        if isinstance(self._target, Sequence) and not isinstance(self._target, str):
            self.__flattenable__ = True

    def target_name(self):
        if isinstance(self._target, str):
            return file_target(self._target).target_name()
        elif isinstance(self._target, BaseTarget):
            return self._target.target_name()
        else:
            return repr(self._target)

    def target_exists(self, mode='any'):
        return True

    def target_signature(self, mode='any'):
        return textMD5(self.target_name())

    def resolve(self):
        return self._target

    def flatten(self):
        return [remote(x) for x in self._target]

    def __format__(self, format_spec):
        # handling special !q conversion flag
        if format_spec and format_spec[0] == 'R':
            return sos_targets(self._target).__format__(format_spec[1:])
        else:
            return str(self).__format__(format_spec)

class executable(BaseTarget):
    '''A target for an executable command.'''

    def __init__(self, cmd, version=None):
        super(executable, self).__init__()
        self._cmd = cmd
        self._md5 = None
        if version is None:
            self._version = ()
        elif isinstance(version, str):
            self._version = (version,)
        else:
            self._version = tuple(version)

    def __eq__(self, other):
        return isinstance(other, executable) and self._cmd == other._cmd and self._version == other._version

    def target_exists(self, mode='any'):
        if mode in ('any', 'target') and shutil.which(shlex.split(self._cmd)[0]):
            if self._version:
                try:
                    output = subprocess.check_output(self._cmd,
                        stderr=subprocess.STDOUT, shell=True, timeout=5).decode()
                except subprocess.TimeoutExpired as e:
                    env.logger.warning(e)
                    return False
                except subprocess.CalledProcessError as e:
                    env.logger.warning(e)
                    return False
                for ver in self._version:
                    if ver in output:
                        return True
                return False
            else:
                return True
        if mode in ('any', 'signature') and os.path.isfile(self.sig_file()):
            return True
        return False

    def target_name(self):
        if self._version:
            return f'{self._cmd} (version={self._version})'
        else:
            return self._cmd

    def target_signature(self, mode='any'):
        if mode != 'target' and self._md5:
            return self._md5
        exe_file = shutil.which(shlex.split(self._cmd)[0])
        if exe_file is None or not os.path.isfile(exe_file):
            self._md5 = None
        else:
            self._md5 = fileMD5(exe_file)
        return self._md5

    def __hash__(self):
        return hash(repr(self))

    def __format__(self, format_spec):
        # handling special !q conversion flag
        if format_spec and format_spec[0] == 'R':
            return self._cmd.__format__(format_spec[1:])
        else:
            return str(self).__format__(format_spec)

class path(type(Path())):
    '''A regular target for files.
    '''
    CONVERTERS = {
        'u': os.path.expanduser,
        'e': lambda x: x.replace(' ', '\\ '),
        'a': lambda x: os.path.abspath(os.path.expanduser(x)),
        'l': lambda x: os.path.realpath(os.path.expanduser(x)),
        'd': os.path.dirname,
        'b': os.path.basename,
        'n': lambda x: os.path.splitext(x)[0],
        'x': lambda x: os.path.splitext(x)[1],
        'q': (lambda x: subprocess.list2cmdline([x])) if sys.platform == 'win32' else quote,
        'p': lambda x: ('/' if len(x) > 1 and x[1] == ':' else '') + x.replace('\\', '/').replace(':', '/'),
        'r': repr,
        's': str,
        # these are handled elsewhere
        ',': lambda x: x,
        '!': lambda x: x,
        'R': lambda x: x,
        }

    def __new__(cls, *args, **kwargs):
        if cls is Path:
            cls = WindowsPath if os.name == 'nt' else PosixPath

        return cls._from_parts(args).expanduser()

    def is_external(self):
        try:
            return os.path.relpath(self.fullname(), env.exec_dir).startswith('..')
        except Exception:
            # under windows the file might be on different volume
            return True

    def fullname(self):
        return os.path.abspath(os.path.expanduser(str(self)))

    def __fspath__(self):
        return self.fullname()

    def __eq__(self, other):
        return os.path.abspath(self.fullname()) == os.path.abspath((other
            if isinstance(other, file_target) else path(other)).fullname())

    def __add__(self, part):
        if isinstance(part, (str, path)):
            return path(str(self) + str(part))
        else:
            raise ValueError(f"Cannot concatenate path to {part} of type {type(part).__name__}: expect a string or path")

    def __format__(self, format_spec):
        # handling special !q conversion flag
        obj = str(self)
        for i,c in enumerate(format_spec):
            if c in self.CONVERTERS:
                obj = self.CONVERTERS[c](obj)
            else:
                # other defined format
                return obj.__format__(format_spec[i:])
        return obj

    def __lt__(self, other):
        return str(self)  < str(other)

    def __hash__(self):
        return hash(repr(self))


class file_target(path, BaseTarget):
    '''A regular target for files.
    '''
    def __init__(self, *args):
        # this is path segments 
        super(file_target, self).__init__(*args)
        if len(args) == 1 and isinstance(args[0], file_target):
            self._md5 = args[0]._md5
            self._attachments = args[0]._attachments
        else:
            self._md5 = None
            self._attachments = []

    def target_exists(self, mode='any'):
        try:
            if mode in ('any', 'target') and self.expanduser().exists():
                return True
            elif mode == 'any' and Path(str(self.expanduser()) + '.zapped').exists():
                return True
            elif mode == 'signature' and Path(self.sig_file()).exists():
                return True
            return False
        except Exception as e:
            env.logger.debug(f"Invalid file_target {self}: {e}")
            return False

    def target_name(self):
        return str(self)

    def __fspath__(self):
        return super(file_target, self).__fspath__()

    def sig_file(self):
        if self._sigfile is not None:
            return self._sigfile
        # If the output path is outside of the current working directory
        fullname = str(self.expanduser().resolve())
        name_md5 = textMD5(fullname)

        if self.is_external():
            self._sigfile = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', name_md5 + '.file_info')
        else:
            self._sigfile = os.path.join(env.exec_dir, '.sos', '.runtime', name_md5 + '.file_info')
        return self._sigfile

    def target_signature(self, mode='any'):
        '''Return file signature'''
        if mode == 'target':
            self._md5 = fileMD5(self.fullname())
        if self._md5 is not None:
            return self._md5
        if os.path.isfile(self.sig_file()) and (not os.path.isfile(self.fullname()) or os.path.getmtime(self.sig_file()) > os.path.getmtime(self.fullname())):
            with open(self.sig_file()) as md5:
                try:
                    line = md5.readline()
                    _, _, _, m = line.rsplit('\t', 3)
                    return m.strip()
                except Exception:
                    pass
        elif os.path.isfile(self.fullname() + '.zapped'):
            with open(self.fullname() + '.zapped') as md5:
                try:
                    line = md5.readline()
                    _, _, _, m = line.rsplit('\t', 3)
                    return m.strip()
                except Exception:
                    pass
        self._md5 = fileMD5(self.fullname())
        return self._md5
    #
    # file_target - specific functions. Not required by other targets
    #
    def add(self, filename):
        '''add related files to the same signature'''
        self._attachments.append(os.path.abspath(os.path.expanduser(filename)))

    def remove(self, mode='both'):
        if mode in ('both', 'target') and os.path.isfile(self.fullname()):
            os.remove(self.fullname())
        if mode in ('both', 'signature') and os.path.isfile(self.sig_file()):
            os.remove(self.sig_file())

    def size(self):
        if self.exists():
            return os.path.getsize(self.fullname())
        elif os.path.isfile(self.sig_file()):
            with open(self.sig_file()) as md5:
                line = md5.readline()
                _, _, s, _ = line.rsplit('\t', 3)
                return int(s.strip())
        elif os.path.isfile(self.fullname() + '.zapped'):
            with open(self.fullname() + '.zapped') as md5:
                line = md5.readline()
                _, _, s, _ = line.rsplit('\t', 3)
                return int(s.strip())
        else:
            raise RuntimeError(f'{self} or its signature does not exist.')

    def mtime(self):
        if self.exists():
            return os.path.getmtime(self.fullname())
        elif os.path.isfile(self.sig_file()):
            with open(self.sig_file()) as md5:
                line = md5.readline()
                _, t, _, _ = line.rsplit('\t', 3)
                return t.strip()
        elif os.path.isfile(self.fullname() + '.zapped'):
            with open(self.fullname() + '.zapped') as md5:
                line = md5.readline()
                _, t, _, _ = line.rsplit('\t', 3)
                return t.strip()
        else:
            raise RuntimeError(f'{self} or its signature does not exist.')

    def write_sig(self):
        '''Write .file_info file with signature'''
        # path to file
        with open(self.sig_file(), 'w') as md5:
            md5.write(
                f'{self.fullname()}\t{os.path.getmtime(self.fullname())}\t{os.path.getsize(self.fullname())}\t{self.target_signature()}\n')
            for f in self._attachments:
                md5.write(f'{f}\t{os.path.getmtime(f)}\t{os.path.getsize(f)}\t{fileMD5(f)}\n')

    def validate(self):
        '''Check if file matches its signature'''
        if not os.path.isfile(self.sig_file()):
            return False
        with open(self.sig_file()) as md5:
            for line in md5:
                f, _, _, m = line.rsplit('\t', 3)
                if not os.path.isfile(f):
                    return False
                if fileMD5(f) != m.strip():
                    env.logger.debug(f'MD5 mismatch {f}')
                    return False
        return True

    def __hash__(self):
        return hash(repr(self))


class paths(Sequence, os.PathLike):
    '''A collection of targets'''
    def __init__(self, *args):
        self._paths = []
        for arg in args:
            self.__append__(arg)

        for t in self._paths:
            if isinstance(t, paths):
                raise RuntimeError(f"Nested paths {t} were introduced by {args}")
            if not isinstance(t, path):
                raise RuntimeError(f"Unrecognized path {t}")

    def paths(self):
        return [str(x) for x in self._paths]

    def __append__(self, arg):
        if isinstance(arg, paths):
            self._paths.extend(arg._paths)
        elif isinstance(arg, str):
            self._paths.append(path(arg))
        elif isinstance(arg, sos_targets):
            if not all(isinstance(x, file_target) for x in arg._targets):
                raise ValueError(f'Cannot convert a sos_targets object {arg} with non-file target to paths')
            self._paths.extend([path(str(x)) for x in arg._targets])
        elif isinstance(arg, file_target):
            self._paths.append(path(str(arg)))
        elif isinstance(arg, Iterable):
            for t in list(arg):
                self.__append__(t)
        elif arg is not None:
            self._paths.append(path(str(arg)))

    def __getstate__(self):
        return self._paths

    def __setstate__(self, paths):
        self._paths = paths

    def __len__(self):
        return len(self._paths)

    def __getitem__(self, i):
        return self._paths[i]

    def __fspath__(self):
        if len(self._paths) == 1:
            return self._paths[0].__fspath__()
        elif len(self._paths) == 0:
            raise ValueError(f"Cannot treat an empty paths as single path")
        else:
            raise ValueError(f'Cannot treat an paths object {self} with more than one paths as a single path')

    def __format__(self, format_spec):
        if ',' in format_spec:
            fmt_spec = format_spec.replace(',', '')
            return ','.join(x.__format__(fmt_spec) for x in self._paths)
        else:
            return ' '.join(x.__format__(format_spec) for x in self._paths)

    def __deepcopy__(self, memo):
        return paths(deepcopy(self._paths))

    def __getattr__(self, key):
        if len(self._paths) == 1:
            return getattr(self._paths[0], key)
        elif len(self._paths) == 0:
            raise AttributeError(f"Cannot get attribute {key} from empty target list")
        else:
            raise AttributeError(f'Cannot get attribute {key} from group of {len(self)} targets {self!r}')

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return self._paths == other._paths if isinstance(other, paths) else other

    def __repr__(self):
        return '[' + ', '.join(repr(x) for x in self._paths) + ']'

    def __str__(self):
        return self.__format__('')

class sos_targets(BaseTarget, Sequence, os.PathLike):
    '''A collection of targets'''
    def __init__(self, *args):
        super(BaseTarget, self).__init__()
        self._targets = []
        for arg in args:
            self.__append__(arg)
        for t in self._targets:
            if isinstance(t, sos_targets):
                raise RuntimeError(f"Nested sos_targets {t} were introduced by {args}")
            if not isinstance(t, BaseTarget):
                raise RuntimeError(f"Unrecognized target {t}")

    def __append__(self, arg):
        if isinstance(arg, Undetermined):
            raise RuntimeError("Undetermined cannot be inserted as a target")
        elif isinstance(arg, paths):
            self._targets.extend([file_target(x) for x in arg._paths])
        elif isinstance(arg, (str, path)):
            self._targets.append(file_target(arg))
        elif isinstance(arg, sos_targets):
            self._targets.extend(arg._targets)
        elif isinstance(arg, BaseTarget):
            self._targets.append(arg)
        elif isinstance(arg, Iterable):
            # in case arg is a Generator, check its type will exhaust it
            for t in list(arg):
                self.__append__(t)
        elif arg is not None:
            raise RuntimeError(f'Unrecognized targets {arg} of type {arg.__class__.__name__}')

    def targets(self):
        return [x.target_name() if isinstance(x, file_target) else x for x in self._targets]

    def extend(self, another):
        if isinstance(another, Undetermined):
            self._targets.append(another)
        else:
            self._targets.extend(sos_targets(another)._targets)

    def has_undetermined(self):
        return any(isinstance(x, Undetermined) for x in self._targets)

    def __getstate__(self):
        return self._targets

    def __setstate__(self, targets):
        self._targets = targets

    def __len__(self):
        return len(self._targets)

    def __getitem__(self, i):
        return self._targets[i]

    def __format__(self, format_spec):
        if ',' in format_spec:
            fmt_spec = format_spec.replace(',', '')
            return ','.join(x.__format__(fmt_spec) for x in self._targets)
        else:
            return ' '.join(x.__format__(format_spec) for x in self._targets)

    def target_signature(self, mode='any'):
        if len(self._targets) == 1:
            return self._targets[0].target_signature()
        else:
            raise ValueError(f'No signature for group of targets {self}')

    def target_exists(self, mode='any'):
        if len(self._targets) == 1:
            return self._targets[0].target_exists(mode)
        else:
            raise ValueError(f'Cannot test existense for group of {len(self)} targets {self!r}')
    
    def __deepcopy__(self, memo):
        return sos_targets(deepcopy(self._targets))

    def __getattr__(self, key):
        if len(self._targets) == 1:
            return getattr(self._targets[0], key)
        elif len(self._targets) == 0:
            raise AttributeError(f"Cannot get attribute {key} from empty target list")
        else:
            raise AttributeError(f'Cannot get attribute {key} from group of {len(self)} targets {self!r}')

    def target_name(self):
        if len(self._targets) == 1:
            return self._targets[0].target_name()
        else:
            raise ValueError(f'Cannot get name() for group of targets {self}')

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return self._targets == other._targets if isinstance(other, sos_targets) else other

    def sig_file(self):
        if len(self._targets) == 1:
            return self._targets[0].sig_file()
        else:
            raise ValueError(f'Cannot get sig_file for group of targets {self}')

    def __add__(self, part):
        if len(self._targets) == 1:
            return self._targets[0].__add__(part)
        elif len(self._targets) == 0:
            raise ValueError(f"Cannot add {part} to empty target list")
        else:
            raise ValueError(f'Cannot add {part} to group of {len(self)} targets {self!r}')

    def __fspath__(self):
        if len(self._targets) == 1:
            return self._targets[0].__fspath__()
        elif len(self._targets) == 0:
            raise ValueError(f"Cannot treat an empty sos_targets as single target")
        else:
            raise ValueError(f'Cannot treat an sos_targets object {self} with more than one targets as a single target')

    def __repr__(self):
        return '[' + ', '.join(repr(x) for x in self._targets) + ']'

    def __str__(self):
        return self.__format__('')

class RuntimeInfo:
    '''Record run time information related to a number of output files. Right now only the
    .exe_info files are used.
    '''
    def __init__(self, step_md5, script, input_files=None, output_files=None, dependent_files = None,
        signature_vars = None, sdict=None):
        '''Runtime information for specified output files

        output_files:
            intended output file

        '''
        if sdict is None:
            sdict = env.sos_dict
        self.step_md5 = step_md5
        self.script = script
        # input can only be a list of files
        if not isinstance(input_files, list):
            if input_files is None:
                self.input_files = []
            else:
                raise RuntimeError('Input files must be a list of filenames for runtime signature.')
        else:
            self.input_files = [file_target(x) if isinstance(x, str) else x for x in input_files]

        if dependent_files is None:
            self.dependent_files = []
        elif isinstance(dependent_files, list):
            self.dependent_files = [file_target(x) if isinstance(x, str) else x for x in dependent_files]
        elif isinstance(dependent_files, Undetermined):
            self.dependent_files = dependent_files
        else:
            raise RuntimeError('Dependent files must be a list of filenames or Undetermined for runtime signature.')

        self.external_output = False
        if output_files is None:
            self.output_files = []
        elif isinstance(output_files, list):
            self.output_files = [file_target(x) if isinstance(x, str) else x for x in output_files]
            self.external_output = self.output_files and isinstance(self.output_files[0], file_target) and self.output_files[0].is_external()
        elif isinstance(output_files, Undetermined):
            self.output_files = output_files
        else:
            raise RuntimeError('Output files must be a list of filenames or Undetermined for runtime signature.')

        self.signature_vars = {} if signature_vars is None else {x: sdict[x] if x in sdict else Undetermined() for x in signature_vars}

        sig_vars = [] if signature_vars is None else sorted([x for x in signature_vars if x in sdict and isPrimitive(sdict[x])])
        self.sig_id = textMD5('{} {} {} {} {}'.format(self.script, self.input_files, output_files, self.dependent_files,
            '\n'.join(f'{x}:{stable_repr(sdict[x])}' for x in sig_vars)))

        if self.external_output:
            # global signature
            self.proc_info = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', f'{self.sig_id}.exe_info')
        else:
            self.proc_info = os.path.join(env.exec_dir, '.sos', '.runtime', f'{self.sig_id}.exe_info')

    def __getstate__(self):
        return {'step_md5': self.step_md5,
                'input_files': self.input_files,
                'output_files': self.output_files,
                'dependent_files': self.dependent_files,
                'signature_vars': self.signature_vars,
                'script': self.script,
                'sig_id': self.sig_id,
                'external': self.external_output}

    def __setstate__(self, sdict):
        self.step_md5 = sdict['step_md5']
        self.input_files = sdict['input_files']
        self.output_files = sdict['output_files']
        self.dependent_files = sdict['dependent_files']
        self.signature_vars = sdict['signature_vars']
        self.script = sdict['script']
        self.sig_id = sdict['sig_id']
        self.external_output = sdict['external']
        #
        # the signature might be on a remote machine and has changed location
        if self.external_output:
            self.proc_info = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', f'{self.sig_id}.exe_info')
        else:
            self.proc_info = os.path.join(env.exec_dir, '.sos', '.runtime', f'{self.sig_id}.exe_info')


    def lock(self):
        # we will need to lock on a file that we do not really write to
        # otherwise the lock will be broken when we write to it.
        self._lock = fasteners.InterProcessLock(self.proc_info + '_')
        if not self._lock.acquire(blocking=False):
            self._lock = None
            raise UnavailableLock((self.input_files, self.output_files, self.proc_info))
        else:
            env.logger.trace(f'Lock acquired for output files {short_repr(self.output_files)}')

    def release(self, quiet=False):
        if self._lock:
            try:
                self._lock.release()
                env.logger.trace(f'Lock released for output files {short_repr(self.output_files)}')
            except Exception as e:
                if not quiet:
                    env.logger.warning(f'Unable to release lock for output files {self.output_files}: {e}')
            finally:
                self._lock = None

    def set(self, files, file_type):
        # add signature file if input and output files are dynamic
        env.logger.trace(f'Set {file_type} of signature to {files}')
        if file_type == 'output':
            self.output_files = [file_target(x) for x in files]
        elif file_type == 'depends':
            self.depends_files = [file_target(x) for x in files]
        else:
            raise RuntimeError(f'Invalid signature file type {file_type}')

    def write(self, rebuild=False):
        '''Write signature file with signature of script, input, output and dependent files.
        Because local input and output files can only be determined after the execution
        of workflow. They are not part of the construction.
        '''
        if isinstance(self.output_files, Undetermined) or isinstance(self.dependent_files, Undetermined):
            env.logger.trace('Write signature failed due to undetermined files')
            return False
        env.logger.trace(f'Write signature {self.proc_info}')
        with open(self.proc_info, 'w') as md5:
            md5.write(f'{textMD5(self.script)}\n')
            md5.write('# input\n')
            for f in self.input_files:
                if f.target_exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write(f'{f}\t{f.target_signature()}\n')
                elif not rebuild and f.target_exists('signature'):
                    md5.write(f'{f}\t{f.target_signature()}\n')
                else:
                    env.logger.warning(f'Failed to create signature: input target {f} does not exist')
                    return False
            md5.write('# output\n')
            for f in self.output_files:
                if f.target_exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write(f'{f}\t{f.target_signature()}\n')
                elif not rebuild and f.target_exists('signature'):
                    md5.write(f'{f}\t{f.target_signature()}\n')
                else:
                    env.logger.warning(f'Failed to create signature: output target {f} does not exist')
                    return False
            md5.write('# dependent\n')
            for f in self.dependent_files:
                if f.target_exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write(f'{f}\t{f.target_signature()}\n')
                elif not rebuild and f.target_exists('signature'):
                    md5.write(f'{f}\t{f.target_signature()}\n')
                else:
                    env.logger.warning(f'Failed to create signature: dependent target {f} does not exist')
                    return False
            # context that will be needed for validation
            md5.write('# init context\n')
            for var in sorted(self.signature_vars.keys()):
                # var can be local and not passed as outside environment
                value = self.signature_vars[var]
                if not isinstance(value, Undetermined):
                    try:
                        var_expr = save_var(var, value)
                        if var_expr:
                            md5.write(var_expr)
                    except Exception:
                        env.logger.debug(f'Variable {var} of value {short_repr(value)} is ignored from step signature')
            # context used to return context
            md5.write('# end context\n')
            for var in sorted(self.signature_vars.keys()):
                # var can be local and not passed as outside environment
                if var in env.sos_dict:
                    value = env.sos_dict[var]
                    try:
                        md5.write(save_var(var, value))
                    except Exception:
                        env.logger.debug(f'Variable {var} of value {short_repr(value)} is ignored from step signature')
            md5.write('# step process\n')
            md5.write(self.script)
            md5.write('# step process\n')
            md5.write(self.script)
        # successfully write signature, write in workflow runtime info
        if '__workflow_sig__' in env.sos_dict and os.path.isfile(env.sos_dict['__workflow_sig__']):
            workflow_sig = env.sos_dict['__workflow_sig__']
            with TimeoutInterProcessLock(workflow_sig + '_'):
                with open(workflow_sig, 'a') as wf:
                    wf.write(
                        f'EXE_SIG\tstep={self.step_md5}\tsession={os.path.basename(self.proc_info).split(".")[0]}\n')
                    for f in self.input_files:
                        if isinstance(f, file_target):
                            wf.write(
                                f'IN_FILE\tfilename={f}\tsession={self.step_md5}\tsize={f.size()}\tmd5={f.target_signature()}\n')
                    for f in self.dependent_files:
                        if isinstance(f, file_target):
                            wf.write(
                                f'IN_FILE\tfilename={f}\tsession={self.step_md5}\tsize={f.size()}\tmd5={f.target_signature()}\n')
                    for f in self.output_files:
                        if isinstance(f, file_target):
                            wf.write(
                                f'OUT_FILE\tfilename={f}\tsession={self.step_md5}\tsize={f.size()}\tmd5={f.target_signature()}\n')
        return True

    def validate(self):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        if not self.proc_info or not os.path.isfile(self.proc_info):
            return f'Missing signature file {self.proc_info}'
        env.logger.trace(f'Validating {self.proc_info}')
        #
        # file not exist?
        if isinstance(self.output_files, Undetermined):
            return "Undetermined output files"
        sig_files = self.input_files + self.output_files + self.dependent_files
        for x in sig_files:
            if not x.target_exists('any'):
                return f'Missing target {x}'
        #
        files_checked = {x.target_name():False for x in sig_files if not isinstance(x, Undetermined)}
        res = {'input': [], 'output': [], 'depends': [], 'vars': {}}
        cur_type = 'input'
        with open(self.proc_info) as md5:
            cmdMD5 = md5.readline().strip()   # command
            if textMD5(self.script) != cmdMD5:
                return "Changed command"
            for line in md5:
                if not line.strip():
                    continue
                if line.startswith('#'):
                    if line == '# input\n':
                        cur_type = 'input'
                    elif line == '# output\n':
                        cur_type = 'output'
                    elif line == '# dependent\n':
                        cur_type = 'depends'
                    elif line == '# init context\n':
                        cur_type = 'init context'
                    elif line == '# end context\n':
                        cur_type = 'end context'
                    elif line == '# step process\n':
                        break
                    else:
                        env.logger.trace(f'Unrecognized line in sig file {line}')
                    continue
                # for validation
                if cur_type == 'init context':
                    key, value = load_var(line)
                    if key not in env.sos_dict:
                        return f'Variable {key} not in running environment'
                    try:
                        try:
                            if env.sos_dict[key] != value:
                                return f'Context variable {key} value mismatch: {short_repr(value)} saved, {short_repr(env.sos_dict[key])} current'
                        except Exception as e:
                            env.logger.debug(f"Variable {key} of type {type(value).__name__} cannot be compared: {e}")
                    except Exception as e:
                        env.logger.warning(f'Failed to restore variable {key} from signature: {e}')
                    continue
                # for return context
                elif cur_type == 'end context':
                    try:
                        key, value = load_var(line)
                        res['vars'][key] = value
                    except Exception as e:
                        env.logger.warning(f'Failed to restore variable from signature: {e}')
                    continue
                try:
                    f, m = line.rsplit('\t', 1)
                    if '(' in f and ')' in f:
                        # this part is hard, because this can be a customized target.
                        target_type = f.split('(')[0]
                        target_class = None
                        if target_type in globals():
                            target_class = eval(target_type)
                        else:
                            # check registry
                            for entrypoint in pkg_resources.iter_entry_points(group='sos_targets'):
                                if entrypoint.name.strip() == target_type:
                                    target_class = entrypoint.load()
                                    break
                        if target_class is None:
                            raise ValueError(f'Failed to identify target class {target_type}')
                        # parameter of class?
                        freal = eval(f, {target_type: target_class})
                    else:
                        freal = file_target(f)
                    if freal.target_exists('target'):
                        fmd5 = freal.target_signature('target')
                    elif freal.target_exists('signature'):
                        fmd5 = freal.target_signature()
                    else:
                        return f'File {f} not exist'
                    res[cur_type].append(freal.target_name() if isinstance(freal, file_target) else freal)
                    if fmd5 != m.strip():
                        return f'File has changed {f}'
                    files_checked[freal.target_name()] = True
                except Exception as e:
                    env.logger.debug(f'Wrong md5 line {line} in {self.proc_info}: {e}')
                    continue
        #
        if not all(files_checked.values()):
            return f'No MD5 signature for {", ".join(x for x,y in files_checked.items() if not y)}'
        env.logger.trace(f'Signature matches and returns {res}')
        # validation success, record signature used
        if '__workflow_sig__' in env.sos_dict and os.path.isfile(env.sos_dict['__workflow_sig__']):
            workflow_sig = env.sos_dict['__workflow_sig__']
            with TimeoutInterProcessLock(workflow_sig + '_'):
                with open(workflow_sig, 'a') as wf:
                    wf.write(self.proc_info + '\n')
        return res

