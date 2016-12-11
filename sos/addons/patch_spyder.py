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

import os
import sys

def patch_spyder2(verbose=False):
    '''Patch spyder to make it work with sos files and sos kernel '''
    try:
        # patch spyderlib/config.py for file extension
        from spyderlib import config
        src_file = config.__file__
        spyderlib_dir = os.path.dirname(src_file)
        if verbose:
            print('Patching {}'.format(src_file))
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
        # patch spyderlib/cli_options.py to add command line option --kernel
        src_file = os.path.join(spyderlib_dir, 'cli_options.py')
        if verbose:
            print('Patching {}'.format(src_file))
        with open(src_file, encoding='utf-8') as src:
            content = src.read()
        with open(src_file, 'w', encoding='utf-8') as src:
            src.write(content.replace('''with Python profiling)")
    options, args = parser.parse_args()''',
            '''with Python profiling)")
    parser.add_option('--kernel', help="Jupyter kernel to start.")
    options, args = parser.parse_args()'''))
        #
        # patch spyderlib.utils.sourcecode.py,
        src_file = os.path.join(spyderlib_dir, 'utils', 'sourcecode.py')
        if verbose:
            print('Patching {}'.format(src_file))
        with open(src_file, encoding='utf-8') as src:
            content = src.read()
        with open(src_file, 'w', encoding='utf-8') as src:
            src.write(content.replace(
                "'Python': ('py', 'pyw', 'python', 'ipy')",
                "'Python': ('py', 'pyw', 'python', 'ipy', 'sos')")
                .replace(
                '''CELL_LANGUAGES = {'Python': ('#%%', '# %%', '# <codecell>', '# In[')}''',
                '''CELL_LANGUAGES = {'Python': ('#%%', '# %%', '# <codecell>', '# In[', '%cell')}''')
            )
        #
        # patch spyderlib/spyder.py
        src_file = os.path.join(spyderlib_dir, 'spyder.py')
        if verbose:
            print('Patching {}'.format(src_file))
        with open(src_file, encoding='utf-8') as src:
            content = src.read()
        with open(src_file, 'w', encoding='utf-8') as src:
            src.write(content.replace(
            '''
    app.exec_()
''',
            r'''
    try:
        if options.kernel == 'sos':
            cfg_file = os.path.expanduser('~/.ipython/profile_default/ipython_config.py')
            has_cfg = os.path.isfile(cfg_file)
            if has_cfg and not os.path.isfile(cfg_file + '.sos_bak'):
                os.rename(cfg_file, cfg_file + '.sos_bak')
            with open(cfg_file, 'w') as cfg:
                cfg.write("c.IPKernelApp.kernel_class =  'sos.jupyter.kernel.SoS_Kernel'\n")
        app.exec_()
    finally:
        if options.kernel == 'sos':
            os.remove(cfg_file)
            if os.path.isfile(cfg_file + '.sos_bak'):
                os.rename(cfg_file + '.sos_bak', cfg_file)
'''))
        #
        print('\nSpyder is successfully patched to accept .sos format and sos kernel.')
        print('Use ')
        print()
        print('    $ spyder --kernel sos')
        print()
        print('to start spyder with sos kernel')
    except Exception as e:
        sys.exit('Failed to patch spyder: {}'.format(e))


def patch_file(src_file, from_text, to_text, verbose):
    with open(src_file, encoding='utf-8') as src:
        content = src.read()
    if to_text in content:
        if verbose:
            print('{} already patched.'.format(src_file))
        return True
    if from_text not in content:
        sys.exit('Failed to patch {}'.format(src_file))
    if verbose:
        print('Patching {}'.format(src_file))
    with open(src_file, 'w', encoding='utf-8') as src:
        src.write(content.replace(from_text, to_text))
    return True

def patch_spyder3(verbose=False):
    '''Patch spyder to make it work with sos files and sos kernel '''
    try:
        # patch spyder/config/utils.py for file extension
        from spyder.config import utils
        src_file = utils.__file__
        spyder_dir = os.path.dirname(os.path.dirname(src_file))
        patch_file(src_file,
        '''
    (_("Cython/Pyrex files"), ('.pyx', '.pxd', '.pxi')),
    (_("C files"), ('.c', '.h')),''',
        '''
    (_("Cython/Pyrex files"), ('.pyx', '.pxd', '.pxi')),
    (_("SoS files"), ('.sos', )),
    (_("C files"), ('.c', '.h')),''',
            verbose=verbose)
        #
        # patch spyder/app/cli_options.py to add command line option --kernel
        patch_file(os.path.join(spyder_dir, 'app', 'cli_options.py'),
            '''help="String to show in the main window title")
    options, args = parser.parse_args()''',
            '''help="String to show in the main window title")
    parser.add_option('--kernel', help="Jupyter kernel to start.")
    options, args = parser.parse_args()''',
            verbose=verbose)
        #
        # patch spyder/utils/sourcecode.py,
        patch_file(os.path.join(spyder_dir, 'utils', 'sourcecode.py'),
                "'Python': ('py', 'pyw', 'python', 'ipy')",
                "'Python': ('py', 'pyw', 'python', 'ipy', 'sos')",
            verbose=verbose)
        patch_file(os.path.join(spyder_dir, 'utils', 'sourcecode.py'),
                '''CELL_LANGUAGES = {'Python': ('#%%', '# %%', '# <codecell>', '# In[')}''',
                '''CELL_LANGUAGES = {'Python': ('#%%', '# %%', '# <codecell>', '# In[', '%cell')}''',
            verbose=verbose)
        #
        # patch spyder/app/mainwindow.py
        patch_file(os.path.join(spyder_dir, 'app', 'mainwindow.py'),
            '''
    app.exec_()
''',
            r'''
    try:
        if options.kernel == 'sos':
            cfg_file = os.path.expanduser('~/.ipython/profile_default/ipython_config.py')
            has_cfg = os.path.isfile(cfg_file)
            if has_cfg and not os.path.isfile(cfg_file + '.sos_bak'):
                os.rename(cfg_file, cfg_file + '.sos_bak')
            with open(cfg_file, 'w') as cfg:
                cfg.write("c.IPKernelApp.kernel_class =  'sos.jupyter.spyder_kernel.SoS_SpyderKernel'\n")
        app.exec_()
    finally:
        if options.kernel == 'sos':
            os.remove(cfg_file)
            if os.path.isfile(cfg_file + '.sos_bak'):
                os.rename(cfg_file + '.sos_bak', cfg_file)
''',
            verbose=verbose)
        #
        print('\nSpyder is successfully patched to accept .sos format and sos kernel.')
        print('Use ')
        print()
        print('    $ spyder --kernel sos')
        print()
        print('to start spyder with sos kernel')
    except Exception as e:
        sys.exit('Failed to patch spyder: {}'.format(e))


def get_patch_spyder_parser():
    import argparse
    parser = argparse.ArgumentParser('patch-spyder',
        description='''Patch spyder to work with SoS kernel''')
    parser.add_argument('-v', '--verbose', action='store_true',
            help='Display details of files to be patched')
    return parser

def patch_spyder(args, unknown_args):
    try:
        from spyderlib import config
        # suppress pyflakes warning
        config
        patch_spyder2(args.verbose)
    except ImportError:
        patch_spyder3(args.verbose)

if __name__ == '__main__':
    parser = get_patch_spyder_parser()
    args, unknown_args = parser.parse_known_args()
    patch_spyder(args, unknown_args)
