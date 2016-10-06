#!/usr/bin/env python
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

import sys, os
import shutil
from setuptools import setup
from distutils import log
from setuptools.command.install import install

# obtain version of SoS
with open('pysos/_version.py') as version:
    exec(version.read())

kernel_json = {
    "argv":         ["python", "-m", "pysos.kernel", "-f", "{connection_file}"],
    "display_name": "SoS",
    "language":     "sos",
}

def patch_spyder():
    '''spyder does not recognize .sos extension. The only change we need to make it work is to
    change the following code from

        ALL_LANGUAGES = {
                     'Python': ('py', 'pyw', 'python', 'ipy'),

    to

        ALL_LANGUAGES = {
                     'Python': ('py', 'pyw', 'python', 'ipy', 'sos'),

    in spyderlib.utils.sourcecode.py, and add 

        (_("SoS files"), ('.sos', )),

    to

        EDIT_FILETYPES = (

    in spyderlib.config.py.
    '''
    try:
        from spyderlib.utils import sourcecode
        src_file = sourcecode.__file__
        with open(src_file, encoding='utf-8') as src:
            content = src.read()
        with open(src_file, 'w', encoding='utf-8') as src:
            src.write(content.replace("'Python': ('py', 'pyw', 'python', 'ipy')",
                "'Python': ('py', 'pyw', 'python', 'ipy', 'sos')")
                .replace(
                '''CELL_LANGUAGES = {'Python': ('#%%', '# %%', '# <codecell>', '# In[')}''',
                '''CELL_LANGUAGES = {'Python': ('#%%', '# %%', '# <codecell>', '# In[', '%cell')}'''))
        #
        from spyderlib import config
        src_file = config.__file__
        with open(src_file, encoding='utf-8') as src:
            content = src.read()
        with open(src_file, 'w', encoding='utf-8') as src:
            src.write(content.replace('''
    (_("Cython/Pyrex files"), ('.pyx', '.pxd', '.pxi')),
    (_("C files"), ('.c', '.h')),''', '''
    (_("Cython/Pyrex files"), ('.pyx', '.pxd', '.pxi')),
    (_("SoS files"), ('.sos', )),
    (_("C files"), ('.c', '.h')),'''))
        #
        log.info('\nAllow spyder to accept .sos as input format.')
    except Exception as e:
        log.info('Failed to patch spyder to accept .sos file: {}'.format(e))

class InstallWithConfigurations(install):
    def run(self):
        # Regular installation
        install.do_egg_install(self)

        # copy sos.vim to .vim
        vim_syntax_dir = os.path.expanduser('~/.vim/syntax')
        vim_syntax_file = os.path.join(vim_syntax_dir, 'sos.vim')
        if not os.path.isdir(vim_syntax_dir):
            os.makedirs(vim_syntax_dir)
        shutil.copy('misc/sos.vim', vim_syntax_file)
        # copy vim-ipython to .vim/ftplugin
        vim_plugin_dir = os.path.expanduser('~/.vim/ftplugin/sos')
        if not os.path.isdir(vim_plugin_dir):
            os.makedirs(vim_plugin_dir)
        shutil.copy('misc/vim-ipython/ipy.vim', os.path.join(vim_plugin_dir, 'ipy.vim'))
        shutil.copy('misc/vim-ipython/vim_ipython.py', os.path.join(vim_plugin_dir, 'vim_ipython.py'))
        #
        # at this point, jupyter and ipython should have been installed.
        import json
        try:
            from jupyter_client.kernelspec import KernelSpecManager as KS
        except ImportError:
            from ipykernel.kernelspec import KernelSpecManager as KS
        from IPython.utils.tempdir import TemporaryDirectory
        from IPython.paths import get_ipython_dir, locate_profile
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
        patch_spyder()
        #
        shutil.copy('misc/sos_magic.py', ext_file)
        shutil.copy('misc/sos_ipython_profile.py', prof_file)
        #
        log.info('\nSoS is installed and configured to use with vim, ipython, spyder, and Jupyter.')
        log.info('Use "set syntax=sos" to enable syntax highlighting.')
        log.info('Use "ipython --profile sos" to start ipython with sos magic.')
        #
        # Now write the kernelspec
        with TemporaryDirectory() as td:
            os.chmod(td, 0o755)  # Starts off as 700, not user readable
            shutil.copy('misc/sos_codemirror.js', os.path.join(td, 'kernel.js'))
            with open(os.path.join(td, 'kernel.json'), 'w') as f:
                json.dump(kernel_json, f, sort_keys=True)
            try:
                KS().install_kernel_spec(td, 'sos', user=self.user, replace=True, prefix=sys.exec_prefix)
                log.info('Use "jupyter notebook" to create or open SoS notebooks.')
            except:
                log.error("\nWARNING: Could not install SoS Kernel as %s user." % self.user)
        log.info('And "sos -h" to start using Script of Scripts.')

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
          # for file lock
          'fasteners',
          'pyyaml',
          'docker-py',
          'blessings',
          'pygments',
          # for jupyter notebook format conversion
          'nbformat',
          'nbconvert>=4.2.0',
          'ipython',
          'notebook',
          # for DAG
          'networkx',
          'pydotplus',
          # for RQ task-execution engine
          'rq-dashboard',
          'rq',
          # for celery task-execution engine
          'celery',
          'flower',
      ],
    entry_points='''
[pygments.lexers]
sos = pysos.converter:SoS_Lexer
''',
    scripts = [
        'sos',
        'sos-runner',
        ]
)
