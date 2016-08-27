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
            self.sig_file = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', sig_name.lstrip(os.sep) + '.file_info')
        else:
            # if this file is relative to cache, use local directory
            self.sig_file = os.path.join('.sos', '.runtime', rel_path + '.file_info')
        self._md5 = []

    def exist(self):
        return os.path.isfile(self.sig_file)

    def remove(self):
        if os.path.isfile(self.sig_file):
            os.remove(self.sig_file)

    def md5(self):
        if self._md5:
            return self._md5[0]
        with open(self.sig_file) as md5:
            line = md5.readline()
            f, m = line.rsplit('\t', 1)
            return m.strip()

    def add(self, filename):
        '''add related files to the same signature'''
        self.filenames.append(os.path.abspath(os.path.expanduser(filename)))

    def write(self):
        '''Write .file_info file with signature'''
        # path to file
        sig_path = os.path.split(self.sig_file)[0]
        if not os.path.isdir(sig_path):
            try:
                os.makedirs(sig_path)
            except Exception as e:
                raise RuntimeError('Failed to create runtime directory {}: {}'.format(sig_path, e))
        with open(self.sig_file, 'w') as md5:
            for filename in self.filenames:
                self._md5.append(fileMD5(filename))
                md5.write('{}\t{}\n'.format(filename, self._md5[-1]))

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
    def __init__(self, script, input_files=[], output_files=[], dependent_files = [], index=None):
        '''Runtime information for specified output files
        index:
            in case of partial output, output files can be the same form (dynamic) so we need index to differntiate

        output_files:
            intended output file

        '''
        self.script = script
        # input can only be a list of files
        if not isinstance(input_files, list):
            if input_files is None:
                self.input_files = []
            else:
                raise RuntimeError('Input files must be a list of filenames for runtime signature.')
        else:
            self.input_files = input_files

        if dependent_files is None:
            self.dependent_files = []
        elif not isinstance(dependent_files, (list, Undetermined)):
            raise RuntimeError('Dependent files must be a list of filenames or Undetermined for runtime signature.')
        else:
            self.dependent_files = dependent_files

        if not isinstance(output_files, (list, Undetermined)):
            raise RuntimeError('Dependent files must be a list of filenames or Undetermined for runtime signature.')
        else:
            self.output_files = output_files
        
        self.index = index
        self.sig_files = []

        if isinstance(self.output_files, Undetermined) or not self.output_files:
            sig_name = 'Dynamic_' + textMD5('{} {} {} {} {}'.format(self.script, self.input_files, output_files, self.dependent_files, self.index))
        else:
            sig_name = os.path.realpath(os.path.expanduser(self.output_files[0]) + '_' + \
                textMD5('{} {} {} {} {}'.format(self.script, self.input_files, self.output_files, self.dependent_files, self.index)))
        #
        # If the output path is outside of the current working directory
        rel_path = os.path.relpath(sig_name, env.exec_dir)
        # if this file is not relative to cache, use global signature file
        if rel_path.startswith('../'):
            info_file = os.path.join(os.path.expanduser('~'), '.sos', '.runtime', sig_name.lstrip(os.sep))
        else:
            # if this file is relative to cache, use local directory
            info_file = os.path.join('.sos', '.runtime', rel_path)
        # path to file
        sig_path = os.path.split(info_file)[0]
        if not os.path.isdir(sig_path):
            try:
                os.makedirs(sig_path)
            except Exception as e:
                raise RuntimeError('Failed to create runtime directory {}: {}'.format(sig_path, e))
        self.proc_info = '{}.exe_info'.format(info_file)

    def set(self, files, file_type):
        # add signature file if input and output files are dynamic
        env.logger.trace('Set {} of signature to {}'.format(file_type, files))
        if file_type == 'output':
            self.output_files = files
        elif file_type == 'depends':
            self.depends_files = files
        else:
            raise RuntimeError('Invalid signature file type {}'.format(file_type))

    def write(self):
        '''Write signature file with signature of script, input, output and dependent files.'''
        if isinstance(self.output_files, Undetermined) or isinstance(self.dependent_files, Undetermined):
            env.logger.trave('Write signature failed due to undetermined files')
            return False
        env.logger.trace('Write signature {}'.format(self.proc_info))
        with open(self.proc_info, 'w') as md5:
            md5.write('{}\n'.format(textMD5(self.script)))
            md5.write('# input\n')
            for f in self.input_files:
                if os.path.isfile(os.path.expanduser(f)):
                    sig = FileSignature(os.path.expanduser(f))
                    # this calculates file MD5
                    sig.write()
                    md5.write('{}\t{}\n'.format(f, sig.md5()))
                else:
                    sig = FileSignature(os.path.expanduser(f))
                    if sig.exist():
                        md5.write('{}\t{}\n'.format(f, sig.md5()))
                    else:
                        return False
            md5.write('# output\n')
            for f in self.output_files:
                if os.path.isfile(os.path.expanduser(f)):
                    sig = FileSignature(os.path.expanduser(f))
                    sig.write()
                    md5.write('{}\t{}\n'.format(f, sig.md5()))
                else:
                    sig = FileSignature(os.path.expanduser(f))
                    if sig.exist():
                        md5.write('{}\t{}\n'.format(f, sig.md5()))
                    else:
                        return False
            md5.write('# dependent\n')
            for f in self.dependent_files:
                if os.path.isfile(os.path.expanduser(f)):
                    sig = FileSignature(os.path.expanduser(f))
                    sig.write()
                    md5.write('{}\t{}\n'.format(f, sig.md5()))
                else:
                    sig = FileSignature(os.path.expanduser(f))
                    if sig.exist():
                        md5.write('{}\t{}\n'.format(f, sig.md5()))
                    else:
                        return False
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
        self.sig_files = self.input_files + self.output_files + self.dependent_files
        if not all(os.path.isfile(x) or FileSignature(x).exist() for x in self.sig_files):
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
                if not line.strip():
                    continue
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
                    if os.path.isfile(freal):
                        fmd5 = fileMD5(freal)
                    else:
                        sig = FileSignature(freal)
                        if sig.exist():
                            fmd5 = sig.md5()
                        else:
                            env.logger.debug('File {} not exist'.format(f))
                            return False
                    res[cur_type].append(f)
                except Exception as e:
                    env.logger.debug('Wrong md5 line {} in {}: {}'.format(line, self.proc_info, e))
                    continue
                if fmd5 != m.strip():
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

