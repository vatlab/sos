#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import os
import shutil

import pytest

from sos import execute_workflow
from sos.targets import file_target


def test_bash():
    '''Test action bash'''
    script=(r'''
        [0]
        bash:
            echo 'Echo'
    ''')
    with pytest.raises(Exception):
        execute_workflow(script)

def test_sh():
    '''Test action run'''
    script=(r'''
        [0]
        sh:
            echo 'Echo'
    ''')
    with pytest.raises(Exception):
        execute_workflow(script)

def test_csh():
    '''Test action csh'''
    if not shutil.which('csh'):
        return
    execute_workflow(r'''
        [0]
        csh:
            foreach color (red orange yellow green blue)
                echo $color
            end
    ''')

def test_tcsh():
    '''Test action tcsh'''
    if not shutil.which('tcsh'):
        return
    execute_workflow(r'''
        [0]
        tcsh:
            foreach color (red orange yellow green blue)
                echo $color
            end
    ''')

def test_zsh():
    '''Test action zsh'''
    if not shutil.which('zsh'):
        return
    execute_workflow(r'''
        [0]
        zsh:
            echo "Hello World!", $SHELL
    ''')

def test_args():
    '''Test args option of scripts'''
    if os.path.isfile('a.txt'):
        file_target('a.txt').unlink()
    execute_workflow(r'''
        [0]
        sh: args='-n {filename:q}'
            touch a.txt
    ''')
    assert not os.path.exists('a.txt')
