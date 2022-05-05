#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil

import pytest
from sos import execute_workflow


@pytest.mark.skipif(not shutil.which('ruby'), reason='ruby not installed')
def test_ruby(clear_now_and_after):
    '''Test action ruby'''
    clear_now_and_after('sample.txt')
    execute_workflow(r'''
        [10]

        ruby:
        fname = "sample.txt"
        somefile = File.open(fname, "w")
        somefile.puts "Hello file!"
        somefile.close
        ''')
    assert os.path.isfile('sample.txt')
