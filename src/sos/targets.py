#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import pickle
import shlex
import shutil
import subprocess
import sys
from collections.abc import Iterable, Sequence
from copy import deepcopy
from pathlib import Path
from shlex import quote
from typing import Union, Dict, Any
import fasteners
import pkg_resources

from .utils import (Error, env, pickleable, short_repr, stable_repr)

from .signatures import target_signatures, step_signatures, workflow_signatures

try:
    from xxhash import xxh64 as hash_md5
except ImportError:
    from hashlib import md5 as hash_md5


__all__ = ['dynamic', 'executable', 'env_variable', 'sos_variable']


class UnknownTarget(Error):
    def __init__(self, target: 'BaseTarget'):
        Error.__init__(self, 'Target unavailable: %s' % target)
        self.target = target


class RemovedTarget(Error):
    def __init__(self, target: 'BaseTarget'):
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


def objectMD5(obj):
    '''Get md5 of an object'''
    try:
        return textMD5(pickle.dumps(obj))
    except:
        return ''


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
            # otherwise, use the first and last 8M
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


class BaseTarget(object):
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

    def create_placeholder(self):
        pass

    # -----------------------------------------------------
    # derived functions that do not need to be redefined
    #
    def remove_sig(self):
        target_signatures.remove(self)

    def write_sig(self):
        '''Write .sig file with signature'''
        raise RuntimeError('Undefined base function')

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
        if step_name.isdigit():
            self._step_name = 'default_' + step_name.lstrip('0')
        elif '_' in step_name and step_name.rsplit('_', 1)[-1].isdigit():
            n, i = step_name.rsplit('_', 1)
            self._step_name = f'{n}_{int(i)}'
        else:
            self._step_name = step_name

    def target_exists(self, mode='any'):
        # the target exists only if it has been executed?
        # which is indicated by a variable
        return '__completed__' in env.sos_dict and (self._step_name in env.sos_dict['__completed__'] or
                                                    self._step_name in [x.rsplit('_', 1)[0] for x in env.sos_dict['__completed__']])

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


def collapseuser(path):
    home = os.path.expanduser('~')
    if path == home:
        return '~'
    elif path.startswith(home + os.sep):
        return '~' + path[len(home):]
    else:
        return path


class path(type(Path())):
    '''A regular target for files.
    '''
    CONVERTERS = {
        'u': os.path.expanduser,
        'U': collapseuser,
        'e': lambda x: x.replace(' ', '\\ '),
        'a': lambda x: os.path.abspath(os.path.expanduser(x)),
        'l': lambda x: os.path.realpath(os.path.expanduser(x)),
        'd': lambda x: os.path.dirname(x) or '.',
        'b': os.path.basename,
        'n': lambda x: os.path.splitext(x)[0],
        'x': lambda x: os.path.splitext(x)[1],
        'q': (lambda x: subprocess.list2cmdline([x])) if sys.platform == 'win32' else quote,
        'p': lambda x: ('/' if len(x) > 1 and x[1] == ':' else '') + x.replace('\\', '/').replace(':/', '/').replace(':', '/'),
        'r': repr,
        's': str,
        # these are handled elsewhere
        ',': lambda x: x,
        '!': lambda x: x,
        'R': lambda x: x,
    }

    # def __new__(cls, *args, **kwargs):
    #    if cls is Path:
    #        cls = WindowsPath if os.name == 'nt' else PosixPath
    #
    #    return cls._from_parts(args).expanduser()

    def _init(self, template=None):
        super(path, self)._init(template)
        if not (self._drv or self._root) and self._parts and self._parts[0][:1] == '~':
            expanded = self.expanduser()
            self._parts = expanded._parts
            self._drv = expanded._drv
            self._root = expanded._root

    def is_external(self):
        try:
            return os.path.relpath(self.fullname(), env.exec_dir).startswith('..')
        except Exception:
            # under windows the file might be on different volume
            return True

    def fullname(self):
        return os.path.abspath(str(self))

    def __fspath__(self):
        return self.fullname()

    def __eq__(self, other):
        return os.path.abspath(self.fullname()) == os.path.abspath((other
                                                                    if isinstance(other, file_target) else path(other)).fullname())

    def __add__(self, part):
        if isinstance(part, (str, path)):
            return path(str(self) + str(part))
        else:
            raise ValueError(
                f"Cannot concatenate path to {part} of type {type(part).__name__}: expect a string or path")

    def __format__(self, format_spec):
        # handling special !q conversion flag
        obj = str(self)
        for i, c in enumerate(format_spec):
            if c in self.CONVERTERS:
                obj = self.CONVERTERS[c](obj)
            else:
                # other defined format
                return obj.__format__(format_spec[i:])
        return obj

    def __lt__(self, other):
        return str(self) < str(other)

    def __hash__(self):
        return hash(repr(self))

    def zap(self):
        zap_file = self + '.zapped'
        if not self.exists() and zap_file.is_file():
            return
        if not self.exists() or not self.is_file():
            raise FileNotFoundError(str(self))
        with open(zap_file, 'w') as md5:
            md5.write(
                f'{self.resolve()}\t{os.path.getmtime(self)}\t{os.path.getsize(self)}\t{fileMD5(self)}\n')
        self.unlink()


class file_target(path, BaseTarget):
    '''A regular target for files.
    '''

    def __init__(self, *args):
        # this is path segments
        super(file_target, self).__init__(*args)
        if len(args) == 1 and isinstance(args[0], file_target):
            self._md5 = args[0]._md5
        else:
            self._md5 = None

    def _init(self, template=None):
        super(file_target, self)._init(template)
        self._md5 = None

    def create_placeholder(self):
        # create an empty placeholder file
        env.logger.debug(f'Create placeholder target {self}')
        self.touch()
        workflow_signatures.write('placeholder', 'file_target', str(self))

    def target_exists(self, mode='any'):
        try:
            if mode in ('any', 'target') and self.exists():
                return True
            elif mode == 'any' and (self + '.zapped').exists():
                return True
            elif mode == 'signature' and target_signatures.get(self):
                return True
            return False
        except Exception as e:
            env.logger.debug(f"Invalid file_target {self}: {e}")
            return False

    def target_name(self):
        return str(self)

    def __fspath__(self):
        return super(file_target, self).__fspath__()

    def target_signature(self, mode='any'):
        '''Return file signature'''
        if mode == 'target':
            self._md5 = fileMD5(self)
            target_signatures.set(self, os.path.getmtime(
                self), os.path.getsize(self), self._md5)
        if self._md5 is not None:
            return self._md5
        sig = target_signatures.get(self)
        if sig is not None:
            self._md5 = sig.md5
            return self._md5
        self._md5 = fileMD5(self)
        target_signatures.set(self, os.path.getmtime(
            self), os.path.getsize(self), self._md5)
        return self._md5

    def remove(self, mode='both'):
        if mode in ('both', 'target') and self.is_file():
            self.unlink()
        if mode in ('both', 'signature'):
            target_signatures.remove(self)

    def size(self):
        if self.exists():
            return os.path.getsize(self)
        sig = target_signatures.get(self)
        if sig:
            return sig.size
        if (self + '.zapped').is_file():
            with open(self + '.zapped') as md5:
                line = md5.readline()
                _, _, s, _ = line.rsplit('\t', 3)
                return int(s.strip())
        else:
            raise RuntimeError(f'{self} or its signature does not exist.')

    def mtime(self):
        if self.exists():
            return os.path.getmtime(self)
        sig = target_signatures.get(self)
        if sig:
            return sig.mtime
        if (self + '.zapped').is_file():
            with open(self + '.zapped') as md5:
                line = md5.readline()
                _, t, _, _ = line.rsplit('\t', 3)
                return t.strip()
        else:
            raise RuntimeError(f'{self} or its signature does not exist.')

    def write_sig(self):
        '''Write signature to sig store'''
        # path to file
        if not self.exists() and (self + '.zapped').exists():
            with open(self + '.zapped') as md5:
                line = md5.readline()
                _, mtime, size, md5 = line.rsplit('\t', 3)
                target_signatures.set(self, mtime, size, md5.strip())
            return
        if not self._md5:
            self._md5 = fileMD5(self)
        target_signatures.set(self,
                              os.path.getmtime(self), os.path.getsize(self), self._md5)

    def validate(self):
        '''Check if file matches its signature'''
        sig = target_signatures.get(self)
        if not sig or not os.path.isfile(self) or sig.size != os.path.getsize(self):
            return False
        if sig.mtime == os.path.getmtime(self):
            return True
        if not self._md5:
            self._md5 = fileMD5(self)
        return self._md5 == sig.md5

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, file_target) and str(self) == str(obj)


class paths(Sequence, os.PathLike):
    '''A collection of targets'''

    def __init__(self, *args):
        self._paths = []
        for arg in args:
            self.__append__(arg)

        for t in self._paths:
            if isinstance(t, paths):
                raise RuntimeError(
                    f"Nested paths {t} were introduced by {args}")
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
                raise ValueError(
                    f'Cannot convert a sos_targets object {arg} with non-file target to paths')
            self._paths.extend([path(str(x)) for x in arg._targets])
        elif isinstance(arg, file_target):
            self._paths.append(path(str(arg)))
        elif isinstance(arg, Iterable):
            for t in list(arg):
                self.__append__(t)
        elif arg is not None:
            self._paths.append(path(str(arg)))

    def zap(self):
        for p in self._paths:
            p.zap()

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
            raise ValueError(
                f'Cannot treat an paths object {self} with more than one paths as a single path')

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
            raise AttributeError(
                f"Cannot get attribute {key} from empty target list")
        else:
            raise AttributeError(
                f'Cannot get attribute {key} from group of {len(self)} targets {self!r}')

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

    def __init__(self, *args, undetermined: Union[bool, str]=None):
        super(BaseTarget, self).__init__()
        self._targets = []
        if isinstance(undetermined, (bool, str)):
            self._undetermined = undetermined
        else:
            self._undetermined = not bool(args)
        for arg in args:
            self.__append__(arg)
        for t in self._targets:
            if isinstance(t, sos_targets):
                raise RuntimeError(
                    f"Nested sos_targets {t} were introduced by {args}")
            if not isinstance(t, BaseTarget):
                raise RuntimeError(f"Unrecognized target {t}")

    def is_external(self):
        if not self.determined():
            return False
        return all(x.is_external() for x in self._targets if isinstance(x, file_target))

    def determined(self):
        return self._targets or not self._undetermined

    def __append__(self, arg):
        if isinstance(arg, paths):
            self._targets.extend([file_target(x) for x in arg._paths])
        elif isinstance(arg, (str, path)):
            self._targets.append(file_target(arg))
        elif isinstance(arg, sos_targets):
            if arg.determined() and not self.determined():
                self._undetermined = False
            self._targets.extend(arg._targets)
        elif isinstance(arg, BaseTarget):
            self._targets.append(arg)
        elif isinstance(arg, Iterable):
            # in case arg is a Generator, check its type will exhaust it
            for t in list(arg):
                self.__append__(t)
        elif arg is not None:
            raise RuntimeError(
                f'Unrecognized targets {arg} of type {arg.__class__.__name__}')

    def targets(self, file_only=False):
        if file_only:
            return [x.target_name() for x in self._targets if isinstance(x, file_target)]
        else:
            return [x.target_name() if isinstance(x, file_target) else x for x in self._targets]

    def extend(self, another):
        self._targets.extend(sos_targets(another)._targets)

    def zap(self):
        for target in self._targets:
            if isinstance(target, file_target):
                target.zap()
            else:
                env.logger.debug(f'Ignore non-file target {target}')

    def __getstate__(self):
        return self._targets

    def __setstate__(self, targets) -> None:
        self._targets = targets
        self._undetermined = True

    def __len__(self):
        return len(self._targets)

    def __getitem__(self, i):
        return self._targets[i]

    def target_signature(self, mode='any'):
        if len(self._targets) == 1:
            return self._targets[0].target_signature()
        else:
            raise ValueError(f'No signature for group of targets {self}')

    def target_exists(self, mode='any'):
        if len(self._targets) == 1:
            return self._targets[0].target_exists(mode)
        else:
            raise ValueError(
                f'Cannot test existense for group of {len(self)} targets {self!r}')

    def __deepcopy__(self, memo):
        return sos_targets(deepcopy(self._targets), undetermined=self._undetermined)

    def __getattr__(self, key):
        if len(self._targets) == 1:
            return getattr(self._targets[0], key)
        elif len(self._targets) == 0:
            raise AttributeError(
                f"Cannot get attribute {key} from empty target list")
        else:
            raise AttributeError(
                f'Cannot get attribute {key} from group of {len(self)} targets {self!r}')

    def target_name(self):
        if len(self._targets) == 1:
            return self._targets[0].target_name()
        else:
            raise ValueError(f'Cannot get name() for group of targets {self}')

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, other):
        return self._targets == other._targets if isinstance(other, sos_targets) else other

    def __add__(self, part):
        if len(self._targets) == 1:
            return self._targets[0].__add__(part)
        elif len(self._targets) == 0:
            raise ValueError(f"Cannot add {part} to empty target list")
        else:
            raise ValueError(
                f'Cannot add {part} to group of {len(self)} targets {self!r}')

    def __fspath__(self):
        if len(self._targets) == 1:
            return self._targets[0].__fspath__()
        elif len(self._targets) == 0:
            raise ValueError(
                f"Cannot treat an empty sos_targets as single target")
        else:
            raise ValueError(
                f'Cannot treat an sos_targets object {self} with more than one targets as a single target')

    def __repr__(self):
        return ('[' + ', '.join(repr(x) for x in self._targets) + ']') if self.determined() else 'TBD'

    def __str__(self):
        return self.__format__('') if self.determined() else 'TBD'

    def __format__(self, format_spec):
        if not self.determined():
            return 'TBD'
        if ',' in format_spec:
            fmt_spec = format_spec.replace(',', '')
            return ','.join(x.__format__(fmt_spec) for x in self._targets)
        else:
            return ' '.join(x.__format__(format_spec) for x in self._targets)


class InMemorySignature:
    def __init__(self, input_files: sos_targets, output_files: sos_targets,
                 dependent_files: sos_targets, signature_vars: set=set(),
                 sdict: dict={}):
        '''Runtime information for specified output files
        '''
        if not sdict:
            sdict = env.sos_dict
        self.input_files = sos_targets(
            [x for x in input_files._targets if not isinstance(x, sos_step)])
        self.dependent_files = sos_targets(
            [x for x in dependent_files._targets if not isinstance(x, sos_step)])
        self.output_files = sos_targets(
            [x for x in output_files._targets if not isinstance(x, sos_step)])
        self.signature_vars = signature_vars
        # signatures that exist before execution and might change during execution
        self.init_signature = {x: deepcopy(sdict[x]) for x in sorted(
            signature_vars) if x in sdict and not callable(sdict[x]) and pickleable(sdict[x], x)}

    def write(self, rebuild=False):
        if not self.output_files.determined() or not self.dependent_files.determined():
            env.logger.trace(
                'Write signature failed due to undetermined files')
            return False
        input_sig = {}
        for f in self.input_files:
            try:
                if f.target_exists('any'):
                    # this calculates file MD5
                    f.write_sig()
                input_sig[str(f)] = f.target_signature()
            except Exception as e:
                env.logger.debug(
                    f'Failed to create signature: input target {f} does not exist')
                return False
        output_sig = {}
        for f in self.output_files:
            try:
                if f.target_exists('any'):
                    # this calculates file MD5
                    f.write_sig()
                output_sig[str(f)] = f.target_signature()
            except Exception as e:
                env.logger.debug(
                    f'Failed to create signature: output target {f} does not exist')
                return False
        dependent_sig = {}
        for f in self.dependent_files:
            try:
                if f.target_exists('any'):
                    # this calculates file MD5
                    f.write_sig()
                dependent_sig[str(f)] = f.target_signature()
            except Exception as e:
                env.logger.debug(
                    f'Failed to create signature: dependent target {f} does not exist')
                return False
        init_context_sig = {var: objectMD5(self.init_signature[var]) for var in self.init_signature if pickleable(
            self.init_signature[var], var)}
        end_context = {var: env.sos_dict[var] for var in self.signature_vars if var in env.sos_dict and pickleable(
            env.sos_dict[var], var)}

        return {
            'input': input_sig,
            'output': output_sig,
            'depends': dependent_sig,
            'init_context_sig': init_context_sig,
            'end_context': end_context
        }

    def validate(self, signature):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        if not signature:
            return 'Empty signature'
        # file not exist?
        if not self.output_files.determined():
            return "Undetermined output files"
        sig_files = self.input_files._targets + self.output_files._targets + \
            self.dependent_files._targets
        for x in sig_files:
            if not x.target_exists('any'):
                return f'Missing target {x}'
        #
        files_checked = {x.target_name(): False for x in sig_files}
        res = {'input': [], 'output': [], 'depends': [], 'vars': {}}
        cur_type = 'input'
        # old signature
        if 'init_context' in signature:
            for key, value in signature['init_context'].items():
                if key not in env.sos_dict:
                    return f'Variable {key} not in running environment'
                try:
                    if env.sos_dict[key] != value:
                        return f'Context variable {key} value mismatch: {short_repr(value)} saved, {short_repr(env.sos_dict[key])} current'
                except Exception as e:
                    env.logger.debug(
                        f"Variable {key} of type {type(value).__name__} cannot be compared: {e}")
        elif 'init_context_sig' in signature:
            for key, value in signature['init_context_sig'].items():
                if key not in env.sos_dict:
                    return f'Variable {key} not in running environment'
                try:
                    if objectMD5(env.sos_dict[key]) != value:
                        return f'ID of context variable {key} mismatch: {short_repr(env.sos_dict[key])} does not match id {value}'
                except Exception as e:
                    env.logger.debug(f"Variable {key} cannot be compared: {e}")

        res['vars'].update(signature['end_context'])
        #
        for cur_type in ['input', 'output', 'depends']:
            for f, m in signature[cur_type].items():
                try:
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
                            raise ValueError(
                                f'Failed to identify target class {target_type}')
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
                    res[cur_type].append(freal.target_name() if isinstance(
                        freal, file_target) else freal)
                    if fmd5 != m.strip():
                        return f'File has changed {f}'
                    files_checked[freal.target_name()] = True
                except Exception as e:
                    env.logger.debug(
                        f'Wrong md5 in signature: {e}')
        #
        if not all(files_checked.values()):
            return f'No MD5 signature for {", ".join(x for x,y in files_checked.items() if not y)}'
        return res


class RuntimeInfo(InMemorySignature):
    '''Record run time information related to a number of output files. Right now only the
    .exe_info files are used.
    '''

    def __init__(self, step_md5: str, script: str, input_files: sos_targets, output_files: sos_targets,
                 dependent_files: sos_targets, signature_vars: set=set(), sdict: dict={}):
        '''Runtime information for specified output files
        '''
        if not sdict:
            sdict = env.sos_dict
        self.step_md5 = step_md5
        self.script = script
        super(RuntimeInfo, self).__init__(input_files, output_files,
                                          dependent_files, signature_vars)

        # if all output files are external
        self.external_sig = self.output_files.is_external() and self.input_files.is_external()

        self.sig_id = textMD5(
            f'{self.script} {self.input_files} {self.output_files} {self.dependent_files} {stable_repr(self.init_signature)}')

    def __getstate__(self):
        return {'step_md5': self.step_md5,
                'input_files': self.input_files,
                'output_files': self.output_files,
                'dependent_files': self.dependent_files,
                'signature_vars': self.signature_vars,
                'init_signature': self.init_signature,
                'script': self.script,
                'sig_id': self.sig_id,
                'external': self.external_sig}

    def __setstate__(self, sdict: Dict[str, Any]):
        self.step_md5 = sdict['step_md5']
        self.input_files = sdict['input_files']
        self.output_files = sdict['output_files']
        self.dependent_files = sdict['dependent_files']
        self.signature_vars = sdict['signature_vars']
        self.init_signature = sdict['init_signature']
        self.script = sdict['script']
        self.sig_id = sdict['sig_id']
        self.external_sig = sdict['external']

    def lock(self):
        # we will need to lock on a file that we do not really write to
        # otherwise the lock will be broken when we write to it.
        self._lock = fasteners.InterProcessLock(
            os.path.join(env.temp_dir,  self.sig_id + '.lock'))
        if not self._lock.acquire(blocking=False):
            self._lock = None
            raise UnavailableLock(
                (self.input_files, self.output_files, self.sig_id))
        else:
            env.logger.trace(
                f'Lock acquired for output files {short_repr(self.output_files)}')

    def release(self, quiet=False):
        if self._lock:
            try:
                self._lock.release()
                env.logger.trace(
                    f'Lock released for output files {short_repr(self.output_files)}')
            except Exception as e:
                if not quiet:
                    env.logger.warning(
                        f'Unable to release lock for output files {self.output_files}: {e}')
            finally:
                self._lock = None

    def set(self, files: sos_targets, file_type: str):
        # add signature file if input and output files are dynamic
        env.logger.trace(f'Set {file_type} of signature to {files}')
        if file_type == 'output':
            self.output_files = files
        elif file_type == 'depends':
            self.depends_files = files
        else:
            raise RuntimeError(f'Invalid signature file type {file_type}')

    def write(self, rebuild=False):
        '''Write signature file with signature of script, input, output and dependent files.
        Because local input and output files can only be determined after the execution
        of workflow. They are not part of the construction.
        '''
        ret = super(RuntimeInfo, self).write()
        if ret is False:
            return ret

        env.logger.trace(f'Write signature {self.sig_id}')
        step_signatures.set(self.sig_id, ret, self.external_sig)
        # successfully write signature, write in workflow runtime info
        for f in self.input_files:
            if isinstance(f, file_target):
                workflow_signatures.write('input_file', self.step_md5,
                                          f'{{"filename":{str(f)!r},"size":{f.size()},"md5":{f.target_signature()!r}}}')
        for f in self.dependent_files:
            if isinstance(f, file_target):
                workflow_signatures.write('dependent_file', self.step_md5,
                                          f'{{"filename":{str(f)!r},"size":{f.size()},"md5":{f.target_signature()!r}}}')
        for f in self.output_files:
            if isinstance(f, file_target):
                workflow_signatures.write('output_file', self.step_md5,
                                          f'{{"filename":{str(f)!r},"size":{f.size()},"md5":{f.target_signature()!r}}}')
        return True

    def validate(self):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        env.logger.trace(f'Validating {self.sig_id}')
        #
        # file not exist?
        if not self.output_files.determined():
            return "Undetermined output files"
        sig_files = self.input_files._targets + self.output_files._targets + \
            self.dependent_files._targets
        for x in sig_files:
            if not x.target_exists('any'):
                return f'Missing target {x}'
        #
        sig = step_signatures.get(self.sig_id, self.external_sig)
        if not sig:
            return f"No signature found for {self.sig_id}"
        return super(RuntimeInfo, self).validate(sig)
