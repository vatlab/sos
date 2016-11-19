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
import hashlib
import shutil
import fasteners
from .utils import env, Error, short_repr
from .sos_eval import Undetermined

__all__ = ['dynamic', 'executable', 'env_variable', 'sos_variable']

class UnknownTarget(Error):
    def __init__(self, target):
        Error.__init__(self, 'Unknown target %s' % target)
        self.target = target

class RemovedTarget(Error):
    def __init__(self, target):
        Error.__init__(self, 'Removed target %s' % target)
        self.target = target

class UnavailableLock(Error):
    """Raised when there are errors in prepare mode. Such errors are not raised
    immediately, but will be collected and raised at the end """

    def __init__(self, signature):
        Error.__init__(self, 'Failed to obtain a lock for output %s' % short_repr(signature[0]))
        self.output = signature[0]
        self.sig_file = signature[1]

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
        pass

    def exists(self, mode='any'):
        raise RuntimeError('Undefined base function')

    def sig_file(self):
        raise RuntimeError('Undefined base function')

class FileTarget(BaseTarget):
    '''A regular target for files.
    '''
    def __init__(self, filename):
        self._filename = filename
        self._sig_file = None
        self._md5 = None
        self._attachments = []

    def exists(self, mode='any'):
        if mode in ('any', 'target') and os.path.isfile(self.fullname()):
            return True
        if mode in ('any', 'signature') and os.path.isfile(self.sig_file()):
            return True
        return False

    def add(self, filename):
        '''add related files to the same signature'''
        self._attachments.append(os.path.abspath(os.path.expanduser(filename)))

    def remove(self, mode='both'):
        if mode in ('both', 'target') and os.path.isfile(self.fullname()):
            os.remove(self.fullname())
        if mode in ('both', 'signature') and os.path.isfile(self.sig_file()):
            os.remove(self.sig_file())

    def fullname(self):
        return os.path.expanduser(self._filename)

    def size(self):
        if os.path.isfile(self._filename):
            return os.path.getsize(self.fullname())
        elif os.path.isfile(self.sig_file()):
            with open(self.sig_file()) as md5:
                line = md5.readline()
                _, _, s, _ = line.rsplit('\t', 3)
                return s.strip()
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
        else:
            raise RuntimeError('{} or its signature does not exist.'.format(self._filename))

    def sig_file(self):
        if self._sig_file is not None:
            return self._sig_file
        # If the output path is outside of the current working directory
        fullname = os.path.abspath(self.fullname())
        name_md5 = textMD5(fullname)
        rel_path = os.path.relpath(fullname, env.exec_dir)

        # if this file is not relative to cache, use global signature file
        if rel_path.startswith('../'):
            self._sig_file = os.path.join(os.path.expanduser('~'), '.sos', '.runtime',
                name_md5 + '.file_info')
        else:
            # if this file is relative to cache, use local directory
            self._sig_file = os.path.join('.sos', '.runtime', name_md5 + '.file_info')
        return self._sig_file

    def __eq__(self, other):
        return os.path.abspath(self.fullname()) == os.path.abspath(other.fullname())

    def write_sig(self):
        '''Write .file_info file with signature'''
        # path to file
        with open(self.sig_file(), 'w') as md5:
            self.calc_md5()
            md5.write('{}\t{}\t{}\t{}\n'.format(self.fullname(), os.path.getmtime(self.fullname()),
                os.path.getsize(self.fullname()), self.md5()))
            for f in self._attachments:
                md5.write('{}\t{}\n'.format(f, fileMD5(f)))

    def calc_md5(self):
        if self._md5 is None:
            self._md5 = fileMD5(self.fullname())
        return self._md5

    def md5(self):
        '''Return md5'''
        if self._md5 is not None:
            return self._md5
        if not os.path.isfile(self.sig_file()):
            return self.calc_md5()
        else:
            with open(self.sig_file()) as md5:
                line = md5.readline()
                _, _, _, m = line.rsplit('\t', 3)
                return m.strip()

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

    def __repr__(self):
        return self._filename

class dynamic(BaseTarget):
    '''A dynamic executable that only handles input files when
    it is available.'''
    def __init__(self, target):
        self._target = target

    def exists(self, mode='any'):
        return True

    def resolve(self):
        return self._target

    def __repr__(self):
        return 'dynamic({})'.format(self._target)

class executable(BaseTarget):
    '''A target for an executable command.'''
    _available_commands = set()

    def __init__(self, cmd, version=[], check_command=None):
        self._cmd = cmd
        if isinstance(version, str):
            self._version = [version]
        else:
            self._version = version
        self._check_command = check_command
        if self._version and not self._check_command:
            self._check_command = cmd
        self.sig_file = os.path.join('.sos/.runtime/{}.sig'.format(self.md5()))

    def exists(self, mode='any'):
        if (self._cmd, self._version) in self._available_commands:
            return True
        if mode in ('any', 'target') and shutil.which(self._cmd):
            if self._version:
                import subprocess
                try:
                    output = subprocess.check_output(self._check_command,
                        stderr=subprocess.STDOUT, shell=True, timeout=5).decode()
                except subprocess.TimeoutExpired as e:
                    env.logger.warning(e)
                    return False
                except subprocess.CalledProcessError as e:
                    env.logger.warning(e)
                    return False
                for ver in self._version:
                    if ver in output:
                        self._available_commands.add((self._cmd, self._version))
                        return True
                return False
            else:
                self._available_commands.add((self._cmd, self._version))
                return True
        if mode in ('any', 'signature') and os.path.isfile(self.sig_file):
            return True
        return False

    def fullname(self):
        if self._version:
            return 'command {} (version={})'.format(self._cmd, self._version)
        else:
            return 'command {}'.format(self._cmd)

    def __repr__(self):
        if self._version:
            return 'executable("{}", version={!r})'.format(self._cmd, self._version)
        else:
            return 'executable("{}")'.format(self._cmd)

    def calc_md5(self):
        return textMD5(self._cmd)

    def md5(self):
        return textMD5(self._cmd)

    def write_sig(self):
        '''Write .sig file with signature'''
        # path to file
        with open(self.sig_file, 'w') as md5:
            md5.write('{}\t{}\n'.format(self.fullname(), self.md5()))

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, executable) and self._cmd == obj._cmd and self._version == obj._version

class sos_variable(BaseTarget):
    '''A target for a SoS variable.'''
    def __init__(self, var):
        self._var = var

    def exists(self, mode='any'):
        return self._var in env.sos_dict

    def fullname(self):
        return 'sos_variable {}'.format(self._var)

    def __repr__(self):
        return 'sos_variable("{}")'.format(self._var)

    def calc_md5(self):
        return textMD5(self._var)

    def md5(self):
        return textMD5(self._var)

    def write_sig(self):
        pass

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, sos_variable) and self._var == obj._var


class env_variable(BaseTarget):
    '''A target for an environmental variable.'''
    def __init__(self, var):
        self._var = var

    def exists(self, mode='any'):
        return self._var in os.environ

    def fullname(self):
        return 'env_variable {}'.format(self._var)

    def __repr__(self):
        return 'env_variable("{}")'.format(self._var)

    def calc_md5(self):
        return textMD5(repr(os.environ[self._var]))

    def md5(self):
        return textMD5(repr(os.environ[self._var]))

    def write_sig(self):
        pass

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, sos_variable) and self._var == obj._var

class R_library(BaseTarget):
    '''A target for a R library.'''

    LIB_STATUS_CACHE = {}

    def __init__(self, library, version = None, repos = 'http://cran.us.r-project.org'):
        self._library = library
        self._version = version
        self._repos = repos

    def _install(self, name, version, repos):
        '''Check existence and version match of R library.
        cran and bioc packages are unique yet might overlap with github.
        Therefore if the input name is {repo}/{pkg} the package will be
        installed from github if not available, else from cran or bioc
        '''
        from .pattern import glob_wildcards
        from .sos_eval import interpolate
        import tempfile
        import shlex
        import subprocess

        output_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False).name
        script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.R', delete=False).name
        if len(glob_wildcards('{repo}/{pkg}', [name])['repo']):
            # package is from github
            self._intall('devtools', version, repos)
            install_script = interpolate('''
            options(warn=-1)
            package_repo <- ${name!r}
            package <- basename(package_repo)
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "AVAILABLE"), file="${output_file}")
            } else {
                devtools::install_github(package_repo)
                # if it still does not exist, write the package name to output
                if (require(package, character.only=TRUE, quietly=TRUE)) {
                    write(paste(package, packageVersion(package), "INSTALLED"), file="${output_file}")
                } else {
                    write(paste(package, "NA", "MISSING"), file="${output_file}")
                    quit("no")
                }
            }
            cur_version <- packageVersion(package)
            ''', '${ }', locals())
        else:
            # package is from cran or bioc
            install_script = interpolate('''
            options(warn=-1)
            package <- ${name!r}
            if (require(package, character.only=TRUE, quietly=TRUE)) {
                write(paste(package, packageVersion(package), "AVAILABLE"), file="${output_file}")
            } else {
                install.packages(package, repos="${repos}",
                    quiet=FALSE)
                # if the package still does not exist
                if (!require(package, character.only=TRUE, quietly=TRUE)) {
                    source("http://bioconductor.org/biocLite.R")
                    biocLite(package, suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
                }
                # if it still does not exist, write the package name to output
                if (require(package, character.only=TRUE, quietly=TRUE)) {
                    write(paste(package, packageVersion(package), "INSTALLED"), file="${output_file}")
                } else {
                    write(paste(package, "NA", "MISSING"), file="${output_file}")
                    quit("no")
                }
            }
            cur_version <- packageVersion(package)
            ''', '${ }', locals())
        version_script = ''
        if version is not None:
            version = [version] if isinstance(version, str) else version
            operators = []
            for idx, value in enumerate(version):
                value = str(value)
                if value.endswith('+'):
                    operators.append('>=')
                    version[idx] = value[:-1]
                elif value.endswith('-'):
                    operators.append('<')
                    version[idx] = value[:-1]
                else:
                    operators.append('==')
            # check version and mark version mismatch
            # if current version satisfies any of the
            # requirement the check program quits
            for x, y in zip(version, operators):
                version_script += '''
                if (cur_version {1} {0}) {{
                  quit("no")
                }}
                '''.format(repr(x), y)
            version_script += 'write(paste(package, cur_version, "VERSION_MISMATCH"), file = {})'.\
              format(repr(output_file))
        # temporarily change the run mode to run to execute script
        try:
            with open(script_file, 'w') as sfile:
                sfile.write(install_script + version_script)
            cmd = 'Rscript --default-packages=utils ' + shlex.quote(script_file)
            #
            p = subprocess.Popen(cmd, shell=True)
            ret = p.wait()
            if ret != 0:
                env.logger.warning('Failed to detect or install R library')
                return False
        except Exception as e:
            env.logger.error('Failed to execute script: {}'.format(e))
            return False
        finally:
            os.remove(script_file)

        ret_val = False
        with open(output_file) as tmp:
            for line in tmp:
                lib, version, status = line.split()
                if status.strip() == "MISSING":
                    env.logger.warning('R Library {} is not available and cannot be installed.'.format(lib))
                elif status.strip() == 'AVAILABLE':
                    env.logger.debug('R library {} ({}) is available'.format(lib, version))
                    ret_val = True
                elif status.strip() == 'INSTALLED':
                    env.logger.debug('R library {} ({}) has been installed'.format(lib, version))
                    ret_val = True
                elif status.strip() == 'VERSION_MISMATCH':
                    env.logger.warning('R library {} ({}) does not satisfy version requirement!'.format(lib, version))
                else:
                    raise RuntimeError('This should not happen: {}'.format(line))
        try:
            os.remove(output_file)
        except:
            pass
        return ret_val

    def exists(self, mode='any'):
        if (self._library, self._version) in self.LIB_STATUS_CACHE:
            return self.LIB_STATUS_CACHE[(self._library, self._version)]
        else:
            ret = self._install(self._library, self._version, self._repos)
            self.LIB_STATUS_CACHE[(self._library, self._version)] = ret
            return ret

    def fullname(self):
        return 'R_library {}'.format(self._library)

    def __repr__(self):
        return 'R_library("{}")'.format(self._library)

    def calc_md5(self):
        return textMD5(repr(self._library))

    def md5(self):
        return textMD5(repr(self._library))

    def write_sig(self):
        pass

    def __hash__(self):
        return hash(repr(self))

    def __eq__(self, obj):
        return isinstance(obj, R_library) and self._library == obj._library

class RuntimeInfo:
    '''Record run time information related to a number of output files. Right now only the
    .exe_info files are used.
    '''
    def __init__(self, step_md5, script, input_files=[], output_files=[], dependent_files = [], signature_vars = []):
        '''Runtime information for specified output files

        output_files:
            intended output file

        '''
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

        if isinstance(output_files, list):
            self.output_files = [FileTarget(x) if isinstance(x, str) else x for x in output_files]
        elif isinstance(output_files, Undetermined):
            self.output_files = output_files
        else:
            raise RuntimeError('Output files must be a list of filenames or Undetermined for runtime signature.')

        self.signature_vars = signature_vars

        sig_name = textMD5('{} {} {} {}'.format(self.script, self.input_files, output_files, self.dependent_files))
        info_file = os.path.join('.sos', '.runtime', sig_name)
        if not isinstance(self.output_files, Undetermined) and self.output_files:
            # If the output path is outside of the current working directory
            rel_path = os.path.relpath(os.path.realpath(self.output_files[0].fullname()), env.exec_dir)
            # if this file is not relative to cache, use global signature file
            if rel_path.startswith('../'):
                info_file = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', sig_name.lstrip(os.sep))
        # path to file
        sig_path = os.path.split(info_file)[0]
        self.proc_info = '{}.exe_info'.format(info_file)

        # we will need to lock on a file that we do not really write to
        # otherwise the lock will be broken when we write to it.
        self.lock = fasteners.InterProcessLock(self.proc_info + '_')
        if not self.lock.acquire(blocking=False):
            raise UnavailableLock((self.output_files, self.proc_info))
        else:
            env.logger.trace('Lock acquired for output files {}'.format(short_repr(self.output_files)))

    def __getstate__(self):
        return {'step_md5': self.step_md5,
                'proc_info': self.proc_info,
                'input_files': self.input_files,
                'output_files': self.output_files,
                'dependent_files': self.dependent_files,
                'signature_vars': self.signature_vars,
                'script': self.script}

    def __setstate__(self, dict):
        self.step_md5 = dict['step_md5']
        self.proc_info = dict['proc_info']
        self.input_files = dict['input_files']
        self.output_files = dict['output_files']
        self.dependent_files = dict['dependent_files']
        self.signature_vars = dict['signature_vars']
        self.script = dict['script']

    def release(self):
        self.lock.release()
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

    def write(self):
        '''Write signature file with signature of script, input, output and dependent files.'''
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
                    md5.write('{}\t{}\n'.format(f, f.md5()))
                elif f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.md5()))
                else:
                    return False
            md5.write('# output\n')
            for f in self.output_files:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.md5()))
                elif f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.md5()))
                else:
                    return False
            md5.write('# dependent\n')
            for f in self.dependent_files:
                if f.exists('target'):
                    # this calculates file MD5
                    f.write_sig()
                    md5.write('{}\t{}\n'.format(f, f.md5()))
                elif f.exists('signature'):
                    md5.write('{}\t{}\n'.format(f, f.md5()))
                else:
                    return False
            md5.write('# context\n')
            for var in sorted(self.signature_vars):
                # var can be local and not passed as outside environment
                if var in env.sos_dict:
                    value = env.sos_dict[var]
                    if isinstance(value, (str, bool, int, float, complex, bytes, list, tuple, set, dict)):
                        md5.write('{} = {!r}\n'.format(var, value))
                    else:
                        env.logger.debug('Variable {} of value {} is ignored from step signature'.format(var, value))
            md5.write('# step process\n')
            md5.write(self.script)
        # successfully write signature, write in workflow runtime info
        workflow_sig = env.sos_dict['__workflow_sig__']
        with fasteners.InterProcessLock(workflow_sig + '_'):
            with open(workflow_sig, 'a') as wf:
                wf.write('EXE_SIG\tstep={}\tsession={}\n'.format(self.step_md5, os.path.basename(self.proc_info).split('.')[0]))
                for f in self.input_files:
                    if isinstance(f, FileTarget):
                        wf.write('IN_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.md5()))
                for f in self.dependent_files:
                    if isinstance(f, FileTarget):
                        wf.write('IN_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.md5()))
                for f in self.output_files:
                    if isinstance(f, FileTarget):
                        wf.write('OUT_FILE\tfilename={}\tsession={}\tsize={}\tmd5={}\n'.format(f, self.step_md5, f.size(), f.md5()))
        return True

    def validate(self):
        '''Check if ofiles and ifiles match signatures recorded in md5file'''
        if not self.proc_info or not os.path.isfile(self.proc_info):
            env.logger.trace('Fail because of no signature file {}'.format(self.proc_info))
            return False
        env.logger.trace('Validating {}'.format(self.proc_info))
        #
        # file not exist?
        if isinstance(self.output_files, Undetermined):
            env.logger.trace('Fail because of undetermined output files.')
            return False
        sig_files = self.input_files + self.output_files + self.dependent_files
        for x in sig_files:
            if not x.exists('any'):
                env.logger.trace('Missing target {}'.format(x))
                return False
        #
        if '__hard_target__' in env.sos_dict:
            for x in self.output_files:
                if not x.exists('target'):
                    env.logger.trace('Missing real target {}'.format(x))
                    return False
        #
        files_checked = {x.fullname():False for x in sig_files if not isinstance(x, Undetermined)}
        res = {'input': [], 'output': [], 'depends': [], 'vars': {}}
        cur_type = 'input'
        with open(self.proc_info) as md5:
            cmdMD5 = md5.readline().strip()   # command
            if textMD5(self.script) != cmdMD5:
                env.logger.trace('Fail because of command change')
                return False
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
                    elif line == '# context\n':
                        cur_type = 'context'
                    elif line == '# step process\n':
                        break
                    else:
                        env.logger.trace('Unrecognized line in sig file {}'.format(line))
                    continue
                if cur_type == 'context':
                    key, value = line.split('=', 1)
                    try:
                        res['vars'][key.strip()] = eval(value.strip())
                    except Exception as e:
                        env.logger.warning('Variable {} with value {} cannot be restored from signature'.format(key, value.strip()))
                    continue
                try:
                    f, m = line.rsplit('\t', 1)
                    if '(' in f and ')' in f:
                        freal = eval(f)
                    else:
                        freal = FileTarget(f)
                    if freal.exists('target'):
                        fmd5 = freal.calc_md5()
                    elif freal.exists('signature'):
                        env.logger.info('Validate with signature of non-existing target {}'.format(freal))
                        fmd5 = freal.md5()
                    else:
                        env.logger.trace('File {} not exist'.format(f))
                        return False
                    res[cur_type].append(freal.fullname() if isinstance(freal, FileTarget) else freal)
                    if fmd5 != m.strip():
                        env.logger.trace('MD5 mismatch {}: {} / {}'.format(f, fmd5, m.strip()))
                        return False
                    files_checked[freal.fullname()] = True
                except Exception as e:
                    env.logger.trace('Wrong md5 line {} in {}: {}'.format(line, self.proc_info, e))
                    continue
        #
        if not all(files_checked.values()):
            env.logger.trace('No MD5 signature for {}'.format(', '.join(x for x,y in files_checked.items() if not y)))
            return False
        env.logger.trace('Signature matches and returns {}'.format(res))
        # validation success, record signature used
        workflow_sig = env.sos_dict['__workflow_sig__']
        with fasteners.InterProcessLock(workflow_sig + '_'):
            with open(workflow_sig, 'a') as wf:
                wf.write(self.proc_info + '\n')
        return res

