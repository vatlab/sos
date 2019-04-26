#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import argparse
import os
import sys
import shutil

_py_ver = sys.version_info
if _py_ver.major == 2 or (_py_ver.major == 3 and
                          (_py_ver.minor, _py_ver.micro) < (6, 0)):
    raise SystemError(
        'sos requires Python 3.6 or higher. Please upgrade your Python {}.{}.{}.'
        .format(_py_ver.major, _py_ver.minor, _py_ver.micro))


def get_install_vim_syntax_parser():
    parser = argparse.ArgumentParser(description='Install vim syntax for sos')
    return parser


def install_vim_syntax(args):

    # copy sos.vim and sos-detect.vim to .vim
    vim_syntax_dir = os.path.join(os.path.expanduser('~'), '.vim', 'syntax')
    vim_syntax_file = os.path.join(vim_syntax_dir, 'sos.vim')
    if not os.path.isdir(vim_syntax_dir):
        os.makedirs(vim_syntax_dir)
    shutil.copy(
        os.path.join(os.path.split(__file__)[0], 'vim', 'sos.vim'),
        vim_syntax_file)
    #
    vim_ftdetect_dir = os.path.join(os.path.expanduser('~'), '.vim', 'ftdetect')
    vim_ftdetect_file = os.path.join(vim_ftdetect_dir, 'sos.vim')
    if not os.path.isdir(vim_ftdetect_dir):
        os.makedirs(vim_ftdetect_dir)
    shutil.copy(
        os.path.join(os.path.split(__file__)[0], 'vim', 'sos-detect.vim'),
        vim_ftdetect_file)
    print(
        'vim syntax file for sos is installed. Use "set syntax=sos" to enable syntax highlighting.'
    )


if __name__ == '__main__':
    parser = get_install_vim_syntax_parser()
    args = parser.parse_args()
    install_vim_syntax(args)
