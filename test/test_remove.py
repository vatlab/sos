#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import subprocess
from pathlib import Path


def assertExists(fdlist):
    for fd in fdlist:
        assert os.path.exists(fd), '{} does not exist'.format(fd)

def assertNonExists(fdlist):
    for fd in fdlist:
        assert not os.path.exists(fd), '{} still exists'.format(fd)

def test_setup(test_workflow):
    assertExists([
        'ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_f1', 'ut_d1/ut_f2',
        'ut_d2/ut_d3/ut_f3'
    ])
    assertExists(
        ['t_f1', 't_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2'])
    # this is the tricky part, directory containing untracked file should remain
    assertExists(['t_d1', 't_d1/ut_f4'])

def test_remove_all_tracked(test_workflow):
    '''test list files'''
    subprocess.call('sos remove . -t -y', shell=True)
    assertExists([
        'ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_f1', 'ut_d1/ut_f2',
        'ut_d2/ut_d3/ut_f3'
    ])
    assertNonExists(['t_d1/t_f2', 't_d2/t_d3/t_f3'])
    # this is the tricky part, directory containing untracked file should remain
    assertExists(['t_d1', 't_f1', 't_d1/ut_f4'])

def test_remove_specific_tracked(test_workflow):
    # note the t_f1, which is under current directory and has to be remove specifically.
    subprocess.call('sos remove t_f1 ut_f1 t_d2 ut_d2 -t -y', shell=True)
    assertExists([
        'ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_f1', 'ut_d1/ut_f2',
        'ut_d2/ut_d3/ut_f3', 't_d1/t_f2', 't_d1', 't_d1/ut_f4'
    ])
    assertNonExists(['t_f1', 't_d2/t_d3/t_f3'])

def test_remove_all_untracked(test_workflow):
    '''test remove all untracked files'''
    subprocess.call('sos remove . -u -y', shell=True)
    assertNonExists(['ut_d1/ut_f2', 't_d1/ut_f4', 'ut_d2/ut_d3/ut_f3'])
    assertExists([
        't_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2', 't_d1', 't_f1'
    ])
    # this is the tricky part, files under the current directory are not removed
    assertExists(['ut_f1'])

def test_remove_specific_untracked(test_workflow):
    # note the t_f1, which is under current directory and has to be remove specifically.
    subprocess.call(
        'sos remove t_f1 ut_f1 ut_d1/ut_f2 t_d1 -u -y', shell=True)
    assertNonExists(['ut_f1', 'ut_d1/ut_f2', 't_d1/ut_f4'])
    assertExists([
        't_d1/t_f2', 't_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2', 't_d1', 't_f1'
    ])
    assertExists(
        ['ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_d2/ut_d3/ut_f3'])

def test_remove_by_age(test_workflow):
    '''test remove by age'''
    subprocess.call('sos remove --age=+1h -y', shell=True)
    # nothing is removed
    assertExists([
        'ut_f1', 'ut_d1/ut_f2', 't_d1/ut_f4', 't_d1/t_f2', 't_d2/t_d3/t_f3',
        't_d2/t_d3', 't_d2', 't_d1', 't_f1', 'ut_d1', 'ut_d2',
        'ut_d2/ut_d3', 'ut_d2/ut_d3/ut_f3'
    ])
    #
    subprocess.call('sos remove -t --age=-1h -y', shell=True)
    assertNonExists(['t_d1/t_f2', 't_d2/t_d3/t_f3'])
    assertExists([
        'ut_f1', 'ut_d1/ut_f2', 't_d1/ut_f4', 't_f1', 't_d2/t_d3', 't_d2',
        't_d1', 'ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_d2/ut_d3/ut_f3'
    ])
    #
    subprocess.call('sos remove -u --age=-1h -y', shell=True)
    assertExists([
        'ut_f1', 't_f1', 't_d2/t_d3', 't_d2', 't_d1', 'ut_d1', 'ut_d2',
        'ut_d2/ut_d3'
    ])

def test_remove_by_size(test_workflow):
    '''test remove by size'''
    subprocess.call('sos remove --size=+10M -y', shell=True)
    # nothing is removed
    assertExists([
        'ut_f1', 'ut_d1/ut_f2', 't_d1/ut_f4', 't_d1/t_f2', 't_d2/t_d3/t_f3',
        't_d2/t_d3', 't_d2', 't_d1', 't_f1', 'ut_d1', 'ut_d2',
        'ut_d2/ut_d3', 'ut_d2/ut_d3/ut_f3'
    ])
    #
    subprocess.call('sos remove -t --size=-1M -y', shell=True)
    assertNonExists(['t_d1/t_f2', 't_d2/t_d3/t_f3'])
    assertExists([
        'ut_f1', 'ut_d1/ut_f2', 't_d1/ut_f4', 't_f1', 't_d2/t_d3', 't_d2',
        't_d1', 'ut_d1', 'ut_d2', 'ut_d2/ut_d3', 'ut_d2/ut_d3/ut_f3'
    ])
    #
    subprocess.call('sos remove -u --size=-1M -y', shell=True)
    assertExists([
        'ut_f1', 't_f1', 't_d2/t_d3', 't_d2', 't_d1', 'ut_d1', 'ut_d2',
        'ut_d2/ut_d3'
    ])

def test_remove_all(test_workflow):
    '''Test remove all specified files'''
    subprocess.call('sos remove ut_d1 t_d1 ut_d2/ut_d3 -y', shell=True)
    assertExists(['t_d2/t_d3/t_f3', 't_d2/t_d3', 't_d2', 't_f1'])

def test_remove_placeholders(test_workflow):
    '''Test remaining placeholder files'''
    # let us create a fake placeholder file
    os.remove('t_f1')
    os.remove('t_d1/t_f2')
    subprocess.call('sos dryrun test.sos', shell=True)
    assert not os.path.isfile('t_f1')
    assert not os.path.isfile('t_d1/t_f2')
    Path('t_f1').touch()
    Path('t_d1/t_f2').touch()
    #
    subprocess.call('sos remove -p', shell=True)
    assert not os.path.isfile('t_f1')
    assert not os.path.isfile('t_d1/t_f2')
