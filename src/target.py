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
import hashlib
import shlex
import shutil
import fasteners
import pkg_resources
from collections import Iterable

from .utils import env, Error, short_repr, stable_repr, save_var, load_var, isPrimitive
from .sos_eval import Undetermined

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
        Error.__init__(self, ('Failed to obtain a lock for input {} and output {}. It is likely ' +
            'that these files are protected by another SoS process or concurrant task ' +
            'that is generating the same set of files.').format(short_repr(signature[0]),
                short_repr(signature[1])))
        self.input = signature[0]
        self.output = signature[1]
        self.sig_file = signature[2]

#
# Runtime signature
#
def textMD5(text):
    '''Get md5 of a piece of text'''
    m = hashlib.md5()
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
    md5 = hashlib.md5()
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
        sys.exit('Failed to read {}: {}'.format(filename, e))
    return md5.hexdigest()


class BaseTarget:
    '''A base class for all targets (e.g. a file)'''
    def __init__(self):
        self._sigfile = None

    def exists(self, mode='any'):
        # mode should be 'any', 'target', or 'signature'
        raise RuntimeError('Undefined base function')

    def name(self):
        # name of the target, which should be able to differentiate
        # this object with other targets of the same type.
        raise RuntimeError('Undefined base function')

    def signature(self, mode='any'):
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
            self._sigfile = os.path.join(env.exec_dir, '.sos', '.runtime', '{}_{}.file_info'.format(self.__class__.__name__,
                textMD5(self.name())))
        return self._sigfile

    def remove_sig(self):
        if self.sig_file() and os.path.isfile(self.sig_file()):
            os.remove(self.sig_file())

    def write_sig(self):
        '''Write .sig file with signature'''
        # path to file
        with open(self.sig_file(), 'w') as sig:
            sig.write('{}\t{}\n'.format(self.name(), self.signature()))

    def __repr__(self):
        return '{}("{}")'.format(self.__class__.__name__, self.name())

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, self.__class__) and self.signature() == obj.signature()

class sos_variable(BaseTarget):
    '''A target for a SoS variable.'''
    def __init__(self, var):
        super(sos_variable, self).__init__()
        self._var = var

    def exists(self, mode='any'):
        return self._var in env.sos_dict

    def name(self):
        return self._var

    def signature(self, mode='any'):
        return textMD5(self._var)

class env_variable(BaseTarget):
    '''A target for an environmental variable.'''
    def __init__(self, var):
        super(env_variable, self).__init__()
        self._var = var

    def exists(self, mode='any'):
        return self._var in os.environ

    def name(self):
        return self._var

    def signature(self, mode='any'):
        return textMD5(repr(os.environ[self._var]))

class sos_step(BaseTarget):
    '''A target for a step of sos.'''
    def __init__(self, step_name):
        super(sos_step, self).__init__()
        self._step_name = step_name

    def exists(self, mode='any'):
        # the target exists only if it has been executed?
        # which is indicated by a variable
        return '__completed__' in env.sos_dict and self._step_name in env.sos_dict['__completed__']

    def name(self):
        return self._step_name

    def signature(self, mode='any'):
        return textMD5('sos_step({})'.format(self._step_name))

    def write_sig(self):
        pass

def bundle(BaseTarget):
    '''a bundle of other targets'''
    def __init__(self, *args):
        super(sos_step, self).__init__()
        self._targets = []
        for arg in args:
            if isinstance(arg, str):
                self._targets.append(os.path.expanduser(arg))
            elif isinstance(arg, Iterable):
                # in case arg is a Generator, check its type will exhaust it
                arg = list(arg)
                if not all(isinstance(x, str) for x in arg):
                    raise RuntimeError('Invalid filename: {}'.format(arg))
                self._targets.extend(arg)
            else:
                raise RuntimeError('Unrecognized file: {}'.format(arg))

    def exists(self, mode='any'):
        return all(FileTarget(x).exists(mode) for x in self._targets)

    def name(self):
        return repr(self._targets)

    def signature(self, mode='any'):
        return textMD5(repr(self._targets))

class dynamic(BaseTarget):
    '''A dynamic executable that only handles input files when
    it is available. This target is handled directly with its `resolve`
    function called by the executor. '''
    def __init__(self, target):
        self._target = target

    def name(self):
        return self._target

    def resolve(self):
        return self._target

class remote(BaseTarget):
    '''A remote target is not tracked and not translated during task execution'''
    def __init__(self, target):
        super(remote, self).__init__()
        self.__unresolvable_object__ = True
        self._target = target

    def name(self):
        return FileTarget(self._target).name() if isinstance(self._target, str) else self._target.name()

    def exists(self, mode='any'):
        return True

    def signature(self, mode='any'):
        return textMD5(self.name())

    def resolve(self):
        return self._target

class local(BaseTarget):
    '''A local target that counters remote() in -r mode'''
    def __init__(self, target):
        super(local, self).__init__()
        self._target = target

    def name(self):
        return FileTarget(self._target).name() if isinstance(self._target, str) else self._target.name()

    def exists(self, mode='any'):
        return False

    def signature(self, mode='any'):
        return textMD5(self.name())

    def resolve(self):
        return self._target

class executable(BaseTarget):
    '''A target for an executable command.'''

    def __init__(self, cmd, version=[]):
        super(executable, self).__init__()
        self._cmd = cmd
        if isinstance(version, str):
            self._version = (version,)
        else:
            self._version = tuple(version)

    def exists(self, mode='any'):
        if mode in ('any', 'target') and shutil.which(shlex.split(self._cmd)[0]):
            if self._version:
                import subprocess
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

    def name(self):
        if self._version:
            return '{} (version={})'.format(self._cmd, self._version)
        else:
            return self._cmd

    def signature(self, mode='any'):
        if mode != 'target' and hasattr(self, '_md5'):
            return self._md5
        exe_file = shutil.which(shlex.split(self._cmd)[0])
        if exe_file is None or not os.path.isfile(exe_file):
            self._md5 = None
        else:
            self._md5 = fileMD5(exe_file)
        return self._md5

class FileTarget(BaseTarget):
    '''A regular target for files.
    '''
    def __init__(self, filename):
        super(FileTarget, self).__init__()
        self._filename = os.path.expanduser(filename)
        self._md5 = None
        self._attachments = []

    def exists(self, mode='any'):
        if mode in ('any', 'target') and os.path.isfile(self.fullname()):
            return True
        elif mode == 'any' and os.path.isfile(self.fullname() + '.zapped'):
            return True
        elif mode == 'signature' and os.path.isfile(self.sig_file()):
            return True
        return False

    def name(self):
        return self._filename

    # redefine sig_file because of special request to store sig files
    # in different folders
    def __repr__(self):
        return self.name()

    def sig_file(self):
        if self._sigfile is not None:
            return self._sigfile
        # If the output path is outside of the current working directory
        fullname = os.path.abspath(self.name())
        name_md5 = textMD5(fullname)

        if self.is_external():
            self._sigfile = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', name_md5 + '.file_info')
        else:
            self._sigfile = os.path.join(env.exec_dir, '.sos', '.runtime', name_md5 + '.file_info')
        return self._sigfile

    def signature(self, mode='any'):
        '''Return file signature'''
        if mode == 'target':
            self._md5 = fileMD5(self.fullname())
        if self._md5 is not None:
            return self._md5
        if os.path.isfile(self.sig_file()):
            with open(self.sig_file()) as md5:
                try:
                    line = md5.readline()
                    _, _, _, m = line.rsplit('\t', 3)
                    return m.strip()
                except:
                    pass
        elif os.path.isfile(self.fullname() + '.zapped'):
            with open(self.fullname() + '.zapped') as md5:
                try:
                    line = md5.readline()
                    _, _, _, m = line.rsplit('\t', 3)
                    return m.strip()
                except:
                    pass                
        self._md5 = fileMD5(self.fullname())
        return self._md5
    #
    # FileTarget - specific functions. Not required by other targets
    #
    def add(self, filename):
        '''add related files to the same signature'''
        self._attachments.append(os.path.abspath(os.path.expanduser(filename)))

    def remove(self, mode='both'):
        if mode in ('both', 'target') and os.path.isfile(self.fullname()):
            os.remove(self.fullname())
        if mode in ('both', 'signature') and os.path.isfile(self.sig_file()):
            os.remove(self.sig_file())

    def is_external(self):
        try:
            return os.path.relpath(self.fullname(), env.exec_dir).startswith('..')
        except:
            # under windows the file might be on different volume
            return True

    def fullname(self):
        return os.path.abspath(self.name())

    def size(self):
        if os.path.isfile(self._filename):
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
            raise RuntimeError('{} or its signature does not exist.'.format(self._filename))

    def mtime(self):
        if os.path.isfile(self._filename):
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
            raise RuntimeError('{} or its signature does not exist.'.format(self._filename))

    def __eq__(self, other):
        return os.path.abspath(self.fullname()) == os.path.abspath(other.fullname())

    def write_sig(self):
        '''Write .file_info file with signature'''
        # path to file
        with open(self.sig_file(), 'w') as md5:
            md5.write('{}\t{}\t{}\t{}\n'.format(self.fullname(), os.path.getmtime(self.fullname()),
                os.path.getsize(self.fullname()), self.signature()))
            for f in self._attachments:
                md5.write('{}\t{}\t{}\t{}\n'.format(f, os.path.getmtime(f),
                    os.path.getsize(f), fileMD5(f)))

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
                    env.logger.debug('MD5 mismatch {}'.format(f))
                    return False
        return True

class RuntimeInfo:
    '''Record run time information related to a number of output files. Right now only the
    .exe_info files are used.
    '''
    def __init__(self, step_md5, script, input_files=[], output_files=[], dependent_files = [],
        signature_vars = [], sdict=None):
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
            self.input_files = [FileTarget(x) if isinstance(x, str) else x for x in input_files]

        if dependent_files is None:
            self.dependent_files = []
        elif isinstance(dependent_files, list):
            self.dependent_files = [FileTarget(x) if isinstance(x, str) else x for x in dependent_files]
        elif isinstance(dependent_files, Undetermined):
            self.dependent_files = dependent_files
        else:
            raise RuntimeError('Dependent files must be a list of filenames or Undetermined for runtime signature.')

        external_output = False
        if output_files is None:
            self.output_files = []
        elif isinstance(output_files, list):
            self.output_files = [FileTarget(x) if isinstance(x, str) else x for x in output_files]
            external_output = self.output_files and isinstance(self.output_files[0], FileTarget) and self.output_files[0].is_external()
        elif isinstance(output_files, Undetermined):
            self.output_files = output_files
        else:
            raise RuntimeError('Output files must be a list of filenames or Undetermined for runtime signature.')

        self.local_input_files = []
        self.local_output_files = []

        self.signature_vars = {x: sdict[x] if x in sdict else Undetermined() for x in signature_vars}

        sig_vars = sorted([x for x in signature_vars if x in sdict and isPrimitive(sdict[x])])
        self.sig_id = textMD5('{} {} {} {} {}'.format(self.script, self.input_files, output_files, self.dependent_files,
            '\n'.join('{}:{}'.format(x, stable_repr(sdict[x])) for x in sig_vars)))

        if external_output:
            # global signature
            self.proc_info = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', '{}.exe_info'.format(self.sig_id))
        else:
            self.proc_info = os.path.join(env.exec_dir, '.sos', '.runtime', '{}.exe_info'.format(self.sig_id))

    def __getstate__(self):
        return {'step_md5': self.step_md5,
                'input_files': self.input_files,
                'output_files': self.output_files,
                'dependent_files': self.dependent_files,
                'local_input_files': self.local_input_files,
                'local_output_files': self.local_output_files,
                'signature_vars': self.signature_vars,
                'script': self.script}

    def __setstate__(self, dict):
        self.step_md5 = dict['step_md5']
        self.input_files = dict['input_files']
        self.output_files = dict['output_files']
        self.local_input_files = dict['local_input_files']
        self.local_output_files = dict['local_output_files']
        self.dependent_files = dict['dependent_files']
        self.signature_vars = dict['signature_vars']
        self.script = dict['script']
        #
        # the signature might be on a remote machine and has changed location
        self.proc_info = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', '{}.exe_info'.format(
            textMD5('{} {} {} {}'.format(self.script, self.input_files, self.output_files, self.dependent_files))))


    def lock(self):
        # we will need to lock on a file that we do not really write to
        # otherwise the lock will be broken when we write to it.
        self._lock = fasteners.InterProcessLock(self.proc_info + '_')
        if not self._lock.acquire(blocking=False):
            self._lock = None
            raise UnavailableLock((self.input_files, self.output_files, self.proc_info))
        else:
            env.logger.trace('Lock acquired for output files {}'.format(short_repr(self.output_files)))

    def release(self):
        if self._lock:
            self._lock.release()
            env.logger.trace('Lock released for output files {}'.format(short_repr(self.output_files)))

    def set(self, files, file_type):
        # add signature file if input and output files are dynamic
        env.logger.trace('Set {} of signature to {}'.format(file_type, files))
        if file_type == 'output':
            self.output_files = [FileTarget(x) for x in files]
        elif file_type == 'depends':
            self.depends_files = [FileTarget(x) for x in files]
        else:
            raise RuntimeError('Invalid signature file type {}'.format(file_type))

    def write(self, local_input_files, local_output_files, rebuild=False):
        '''Write signature file with signature of script, input, output and dependent files.
        Because local input and output files can only be determined after the execution
        of workflow. They are not part of the construction...
        '''
        if isinstance(self.output_files, Undetermined) or isinstance(self.dependent_files, Undetermined):
            env.logger.trace('Write signature failed due to undetermined files')
            return False
        env.logger.trace('Write signature {}'.format(self.proc_info))
        with open(self.proc_info, 'w') as md5:
            md5.write('{}\n'.format(textMD5(self.script)))
            md5.write('# input\n')
            for f in self.input_files:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                elif not rebuild and f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                else:
                    env.logger.warning('Failed to create signature: input target {} does not exist'.format(f))
                    return False
            md5.write('# output\n')
            for f in self.output_files:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                elif not rebuild and f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                else:
                    env.logger.warning('Failed to create signature: output target {} does not exist'.format(f))
                    return False
            md5.write('# dependent\n')
            for f in self.dependent_files:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                elif not rebuild and f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                else:
                    env.logger.warning('Failed to create signature: dependent target {} does not exist'.format(f))
                    return False
            md5.write('# local input\n')
            for f in [FileTarget(x) if isinstance(x, str) else x for x in local_input_files]:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                elif not rebuild and f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                else:
                    env.logger.warning('Failed to create signature: local input target {} does not exist'.format(f))
                    return False
            md5.write('# local output\n')
            for f in [FileTarget(x) if isinstance(x, str) else x for x in local_output_files]:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                elif not rebuild and f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.signature()))
                else:
                    env.logger.warning('Failed to create signature: local output target {} does not exist'.format(f))
                    return False
            # context that will be needed for validation
            md5.write('# init context\n')
            for var in sorted(self.signature_vars.keys()):
                # var can be local and not passed as outside environment
                value = self.signature_vars[var]
                if not isinstance(value, Undetermined):
                    try:
                        md5.write(save_var(var, value))
                    except Exception as e:
                        env.logger.debug('Variable {} of value {} is ignored from step signature'.format(var, short_repr(value)))
            # context used to return context
            md5.write('# end context\n')
            for var in sorted(self.signature_vars.keys()):
                # var can be local and not passed as outside environment
                if var in env.sos_dict:
                    value = env.sos_dict[var]
                    try:
                        md5.write(save_var(var, value))
                    except Exception as e:
                        env.logger.debug('Variable {} of value {} is ignored from step signature'.format(var, short_repr(value)))
            md5.write('# step process\n')
            md5.write(self.script)
            md5.write('# step process\n')
            md5.write(self.script)
        # successfully write signature, write in workflow runtime info
        if '__workflow_sig__' in env.sos_dict and os.path.isfile(env.sos_dict['__workflow_sig__']):
            workflow_sig = env.sos_dict['__workflow_sig__']
            with fasteners.InterProcessLock(workflow_sig + '_'):
                with open(workflow_sig, 'a') as wf:
                    wf.write('EXE_SIG\tstep={}\tsession={}\n'.format(self.step_md5, os.path.basename(self.proc_info).split('.')[0]))
                    for f in self.input_files:
                        if isinstance(f, FileTarget):
                            wf.write('IN_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.signature()))
                    for f in self.dependent_files:
                        if isinstance(f, FileTarget):
                            wf.write('IN_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.signature()))
                    for f in self.output_files:
                        if isinstance(f, FileTarget):
                            wf.write('OUT_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.signature()))
                    for f in self.local_input_files:
                        if isinstance(f, FileTarget):
                            wf.write('IN_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.signature()))
                    for f in self.local_output_files:
                        if isinstance(f, FileTarget):
                            wf.write('OUT_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.signature()))
        return True

    def validate(self):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        if not self.proc_info or not os.path.isfile(self.proc_info):
            return 'Missing signature file {}'.format(self.proc_info)
        env.logger.trace('Validating {}'.format(self.proc_info))
        #
        # file not exist?
        if isinstance(self.output_files, Undetermined):
            return "Undetermined output files"
        sig_files = self.input_files + self.output_files + self.dependent_files
        for x in sig_files:
            if not x.exists('any'):
                return 'Missing target {}'.format(x)
        #
        files_checked = {x.name():False for x in sig_files if not isinstance(x, Undetermined)}
        res = {'input': [], 'output': [], 'depends': [], 'local_input': [], 'local_output': [], 'vars': {}}
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
                    elif line == '# local input\n':
                        cur_type = 'local_input'
                    elif line == '# local output\n':
                        cur_type = 'local_output'
                    elif line == '# init context\n':
                        cur_type = 'init context'
                    elif line == '# end context\n':
                        cur_type = 'end context'
                    elif line == '# step process\n':
                        break
                    else:
                        env.logger.trace('Unrecognized line in sig file {}'.format(line))
                    continue
                # for validation
                if cur_type == 'init context':
                    key, value = load_var(line)
                    if key not in env.sos_dict:
                        return 'Variable {} not in running environment {}'.format(key)
                    try:
                        try:
                            if env.sos_dict[key] != value:
                                return 'Context variable {} value mismatch: {} saved, {} current'.format(
                                    key, short_repr(value), short_repr(env.sos_dict[key]))
                        except Exception as e:
                            env.logger.warning("Variable {} of type {} cannot be compared: {}".format(key, type(key).__name__, e))
                    except Exception as e:
                        env.logger.warning('Failed to restore variable {} from signature: {}'.format(key, e))
                    continue
                # for return context
                elif cur_type == 'end context':
                    try:
                        key, value = load_var(line)
                        res['vars'][key] = value
                    except Exception as e:
                        env.logger.warning('Failed to restore variable {} from signature: {}'.format(key, e))
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
                            raise ValueError('Failed to identify target class {}'.format(target_type))
                        # parameter of class?
                        freal = eval(f, {target_type: target_class})
                    else:
                        freal = FileTarget(f)
                    if freal.exists('target'):
                        fmd5 = freal.signature('target')
                    elif freal.exists('signature'):
                        fmd5 = freal.signature()
                    else:
                        return 'File {} not exist'.format(f)
                    res[cur_type].append(freal.name() if isinstance(freal, FileTarget) else freal)
                    if fmd5 != m.strip():
                        return 'File has changed {}'.format(f)
                    files_checked[freal.name()] = True
                except Exception as e:
                    env.logger.debug('Wrong md5 line {} in {}: {}'.format(line, self.proc_info, e))
                    continue
        #
        if not all(files_checked.values()):
            return 'No MD5 signature for {}'.format(', '.join(x for x,y in files_checked.items() if not y))
        env.logger.trace('Signature matches and returns {}'.format(res))
        # validation success, record signature used
        if '__workflow_sig__' in env.sos_dict and os.path.isfile(env.sos_dict['__workflow_sig__']):
            workflow_sig = env.sos_dict['__workflow_sig__']
            with fasteners.InterProcessLock(workflow_sig + '_'):
                with open(workflow_sig, 'a') as wf:
                    wf.write(self.proc_info + '\n')
        return res

