#!/usr/bin/env python
#
# This file is part of Script of Scripts (sos), a workflow system
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
import argparse
import shutil


def get_install_vim_syntax_parser():
    parser = argparse.ArgumentParser(description='Install vim syntax for sos')
    return parser

def install_vim_syntax(args):
    
    # copy sos.vim and sos-detect.vim to .vim
    vim_syntax_dir = os.path.join(os.path.expanduser('~'), '.vim', 'syntax')
    vim_syntax_file = os.path.join(vim_syntax_dir, 'sos.vim')
    if not os.path.isdir(vim_syntax_dir):
        os.makedirs(vim_syntax_dir)
    shutil.copy(os.path.join(os.path.split(__file__)[0], 'vim', 'sos.vim'), vim_syntax_file)
    #
    vim_ftdetect_dir = os.path.join(os.path.expanduser('~'), '.vim', 'ftdetect')
    vim_ftdetect_file = os.path.join(vim_ftdetect_dir, 'sos.vim')
    if not os.path.isdir(vim_ftdetect_dir):
        os.makedirs(vim_ftdetect_dir)
    shutil.copy(os.path.join(os.path.split(__file__)[0], 'vim', 'sos-detect.vim'), vim_ftdetect_file)
    print('vim syntax file for sos is installed. Use "set syntax=sos" to enable syntax highlighting.')


if __name__ == '__main__':
    parser = get_install_vim_syntax_parser()
    args = parser.parse_args()
    install_vim_syntax(args)

