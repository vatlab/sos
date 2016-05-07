#!/usr/bin/env python2.7
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
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
import filecmp
import shutil
from setuptools import setup
from distutils import log
from distutils.command.install import install

try:
    import json
    try:
        from jupyter_client.kernelspec import install_kernel_spec
    except ImportError:
        from IPython.kernel.kernelspec import install_kernel_spec
    from IPython.utils.tempdir import TemporaryDirectory
    from IPython.paths import get_ipython_dir, locate_profile
    ipython = True
except:
    log.info('\nJupyter kernel of SoS is not installed because no Jupyter installation is found.')
    ipython = False

kernel_json = {
    "argv":         ["python", "-m", "pysos.kernel", "-f", "{connection_file}"],
    "display_name": "SoS",
    "language":     "sos",
}

class InstallWithConfigurations(install):
    def run(self):
        # Regular installation
        install.run(self)

        # copy sos.vim to .vim
        vim_dir = os.path.expanduser('~/.vim/syntax')
        vim_file = os.path.join(vim_dir, 'sos.vim')
        if not os.path.isdir(vim_dir):
            os.makedirs(vim_dir)
        shutil.copy('misc/sos.vim', vim_file)
        log.info('\nsos.vim is copied to ~/.vim/syntax. Use')
        log.info('   set syntax=sos')
        log.info('in vim/gvim/mvim to enable syntax highlighting, or add')
        log.info('   autocmd BufNewFile,BufRead *.sos set syntax=sos')
        log.info('to your ~/.vimrc (~/.gvimrc) file')

        if not ipython:
            return
        #
        # copy ipython magic to ~/.ipython/extensions
        ext_dir = os.path.join(get_ipython_dir(), 'extensions')
        ext_file = os.path.join(ext_dir, 'sos_magic.py')
        if not os.path.isdir(ext_dir):
            os.makedirs(ext_dir)
        prof_dir = os.path.join(get_ipython_dir(), 'profile_sos')
        prof_file = os.path.join(prof_dir, 'ipython_config.py')
        if not os.path.isdir(prof_dir):
            os.makedirs(prof_dir)
        #
        shutil.copy('misc/sos_magic.py', ext_file)
        shutil.copy('misc/sos_ipython_profile.py', prof_file)
        #
        log.info('\nsos ipython magic is installed to {}. Use '.format(ext_dir))
        log.info('    ipython --profile sos')
        log.info('to use a profile with pre-loaded sos magic, or use')
        log.info('    %load_ext sos_magic')
        log.info('to load it in an ipython interactive shell, or insert sos_magic to')
        log.info('    c.InteractiveShellApp.extensions')
        log.info('of {}/ipython_config.py to load it automatically to your default profile.'.format(prof_dir))
        #
        # Now write the kernelspec
        with TemporaryDirectory() as td:
            os.chmod(td, 0o755)  # Starts off as 700, not user readable
            shutil.copy('misc/sos_codemirror.js', os.path.join(td, 'kernel.js'))
            with open(os.path.join(td, 'kernel.json'), 'w') as f:
                json.dump(kernel_json, f, sort_keys=True)
            try:
                install_kernel_spec(td, 'sos', user=self.user, replace=True)
                log.info('\nJupyter kernel named "sos" is installed. You can use it by running')
                log.info('    jupyter qtconsole --kernel sos')
                log.info('or start a jupyter notebook with command')
                log.info('    jupyter notebook')
                log.info('and create a new file with kernel SoS.')
            except:
                log.error("\nWARNING: Could not install SoS Kernel as %s user." % self.user)

exec(open('pysos/_version.py').read())

setup(name = "sos",
    version = __version__,
    description = 'Script of Scripts (SoS): a lightweight workflow system for the creation of readable workflows',
    author = 'Bo Peng',
    url = 'https://github.com/BoPeng/SOS',
    author_email = 'bpeng@mdanderson.org',
    maintainer = 'Bo Peng',
    maintainer_email = 'bpeng@mdanderson.org',
    license = 'GPL3',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        ],
    packages = ['pysos'],
    cmdclass={'install': InstallWithConfigurations},
    install_requires=[
          'psutil',
          'pyyaml',
          'docker-py',
          'pycurl',
          'blessings',
          'pygments',
          # for jupyter notebook format conversion
          'nbformat',
      ],
    entry_points='''
[pygments.lexers]
sos = pysos.sos_show:SoS_Lexer
''',
    scripts = [
        'sos',
        'sos-runner',
        ]
)
