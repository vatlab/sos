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
from sos.utils import env
from sos.target import BaseTarget, textMD5

class R_library(BaseTarget):
    '''A target for a R library.'''

    LIB_STATUS_CACHE = {}

    def __init__(self, library, version = None, repos = 'http://cran.us.r-project.org'):
        super(R_library, self).__init__()
        self._library = library
        self._version = version
        self._repos = repos

    def _install(self, name, version, repos):
        '''Check existence and version match of R library.
        cran and bioc packages are unique yet might overlap with github.
        Therefore if the input name is {repo}/{pkg} the package will be
        installed from github if not available, else from cran or bioc
        '''
        from sos.pattern import glob_wildcards
        from sos.sos_eval import interpolate
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

    def name(self):
        return self._library

    def md5(self):
        # we are supposed to get signature of the library, but we cannot
        return textMD5(repr(self._library))
