#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil

import pytest

from sos import execute_workflow


@pytest.mark.skipif(not shutil.which('julia'), reason="julia not available")
def test_julia(self):
    '''Test action Julia'''
    if os.path.isfile('julia_example.txt'):
        os.remove('julia_example.txt')
    execute_workflow(r'''
        [0]
        julia:
            open("julia_example.txt", "w") do f
                write(f, "A, B, C, D\n")
            end
    ''')
    assert os.path.isfile('julia_example.txt')
