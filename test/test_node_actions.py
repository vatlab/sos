#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import shutil

import pytest

from sos import execute_workflow


@pytest.mark.skipif(not shutil.which('node'), reason='node not installed')
def test_node():
    '''Test action node'''
    execute_workflow(r'''
        [0]
        node:
        var args = process.argv.slice(2);
        console.log('Hello ' + args.join(' ') + '!');
        ''')
