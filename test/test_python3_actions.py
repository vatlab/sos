#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos import execute_workflow


def test_python():
    '''Test python command. This might fail if python3 is the
    default interpreter'''
    execute_workflow(r'''
        [0]
        python: expand='${ }'
        a = {'1': 2}
        print(a)
        ''')


def test_python3():
    execute_workflow(r'''
        [0]
        python3: expand='${ }'
        a = {'1', '2'}
        print(a)
        ''')
