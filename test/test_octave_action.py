#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil

import pytest

from sos import execute_workflow


@pytest.mark.skipif(not shutil.which('octave'), reason='Octave not installed')
def test_octave():
    '''Test action octave'''
    if os.path.isfile('octave_example.txt'):
        os.remove('octave_example.txt')
    execute_workflow(r'''
        [0]
        octave:
            f = fopen("octave_example.txt", "w")
            fprintf(f, "A, B, C, D\n")
            fclose(f)
        ''')
    assert os.path.isfile('octave_example.txt')
