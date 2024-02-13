#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import os
import shutil

import pytest

from sos import execute_workflow


def test_bash():
    '''Test action bash'''
    execute_workflow(r'''
        [0]
        bash:
        echo 'Echo'
        ''')


def test_bash_1():
    with pytest.raises(Exception):
        execute_workflow(r'''
            [0]
            bash:
            echo 'Echo
            ''')


def test_sh():
    '''Test action run'''
    execute_workflow(r'''
        [0]
        sh:
        echo 'Echo'
        ''')


def test_sh_1():
    with pytest.raises(Exception):
        execute_workflow(r'''
            [0]
            sh:
            echo 'Echo
            ''')


@pytest.mark.skipif(not shutil.which('csh'), reason="Needs csh")
def test_csh():
    '''Test action csh'''
    execute_workflow(r'''
        [0]
        csh:
            foreach color (red orange yellow green blue)
                echo $color
            end
    ''')


@pytest.mark.skipif(not shutil.which('tcsh'), reason="Needs tcsh")
def test_tcsh():
    '''Test action tcsh'''
    execute_workflow(r'''
        [0]
        tcsh:
            foreach color (red orange yellow green blue)
                echo $color
            end
    ''')


@pytest.mark.skipif(not shutil.which('zsh'), reason="Needs zsh")
def test_zsh():
    '''Test action zsh'''
    execute_workflow(r'''
        [0]
        zsh:
            echo "Hello World!", $SHELL
    ''')


def test_args(clear_now_and_after):
    '''Test args option of scripts'''
    clear_now_and_after('a.txt')

    execute_workflow(r'''
        [0]
        sh: args='-n {filename:q}'
            touch a.txt
    ''')
    assert not os.path.exists('a.txt')
