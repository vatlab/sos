#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import shutil

import pytest

from sos import execute_workflow


@pytest.mark.skipif(
    not shutil.which('python2.7'),
    reason='Skip test because of no python2.7 installation')
def test_python2():
    execute_workflow(r'''
        [0]
        python2: expand='${ }'
        a = {'1', '2'}
        print a
        ''')
