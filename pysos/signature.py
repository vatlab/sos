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
from .utils import env
from .sos_eval import Undetermined

__all__ = []

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


class FileInfo:
    def __init__(self, filename):
        self.filename = filename

    def first5(self):
        try:
            lines = []
            with open(self.filename, 'r') as ifile:
                for i in range(5):
                    lines.append(ifile.readline())
            return ''.join(lines)
        except:
            return ''

    def gtf(self):
        return self.first5()

    def gff(self):
        return self.first5()

    def describe(self):
        # return description of files
        fn, ext = os.path.splitext(self.filename)
        ext = ext.lower()
        if not ext:
            return ''
        if hasattr(self, ext[1:]):
            return getattr(self, ext[1:])()
        else:
            try:
                return self.first5()
            except:
                return ''

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
        return {'input': self.input_files, 'output': self.output_files, 'depends': self.dependent_files}

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

