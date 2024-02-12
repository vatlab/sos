#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil

import pytest

from sos import execute_workflow


@pytest.mark.skipif(not shutil.which('matlab'), reason='Matlab not installed')
def test_matlab(clear_now_and_after):
    '''Test action matlab'''
    clear_now_and_after('/tmp/matlab_example.txt')
    execute_workflow(r'''
        [0]
        matlab:
            f = fopen("/tmp/matlab_example.txt", "w");
            fprintf(f, "A, B, C, D\n");
            fclose(f);
    ''')
    assert os.path.isfile('/tmp/matlab_example.txt')
