#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os

import pytest
from sos import execute_workflow
# if the test is imported under sos/test, test interacive executor
from sos.workflow_executor import Base_Executor


def test_plain_target(clear_now_and_after):
    '''Test sos run -t filename'''
    clear_now_and_after('t_a.txt')
    execute_workflow(
        '''
        [A]
        output: 't_a.txt'
        _output.touch()
        ''', targets=['t_a.txt'])
    assert os.path.isfile('t_a.txt')


def test_target_in_named(clear_now_and_after):
    '''Test sos run -t filename'''
    clear_now_and_after('t_na.txt')
    execute_workflow(
        '''
        [A]
        output: res='t_na.txt'
        _output.touch()
        ''',
                targets=['t_na.txt'])
    assert os.path.isfile('t_na.txt')


def test_named_output_as_target(clear_now_and_after):
    '''Test sos run -t named_output'''
    clear_now_and_after('t_no.txt')
    execute_workflow(
        '''
        [A]
        output: res='t_no.txt'
        _output.touch()
        ''', targets=['res'])
    assert os.path.isfile('t_no.txt')


def test_provides_target(clear_now_and_after):
    '''Test sos run -t filename with exact match'''
    clear_now_and_after('t_pa.txt')
    execute_workflow(
        '''
        [A: provides="t_pa.txt"]
        _output.touch()
        ''', targets=['t_pa.txt'])
    assert os.path.isfile('t_pa.txt')
    #
    clear_now_and_after('t_pa.txt')
    execute_workflow(
        '''
        [A: provides="t_pa.txt"]
        output: 't_pa.txt'
        _output.touch()
        ''',
        targets=['t_pa.txt'])
    assert os.path.isfile('t_pa.txt')
    #
    clear_now_and_after('t_pa.txt')
    execute_workflow(
        '''
        [A: provides="t_pa.txt"]
        output: pa='t_pa.txt'
        _output.touch()
        ''',
        targets=['t_pa.txt'])
    assert os.path.isfile('t_pa.txt')
    #
    clear_now_and_after('t_pa.txt')
    with pytest.raises(Exception):
        execute_workflow(
            '''
            [A: provides="t_pa.txt"]
            output: pa='t_pa_none.txt'
            _output.touch()
            ''',
            targets=['t_pa.txt'])


def test_provides_pattern(clear_now_and_after):
    '''Test sos run -t filename with pattern matching'''
    clear_now_and_after('t_ma.txt')
    execute_workflow('''
        [A: provides="{filename}.txt"]
        _output.touch()
        ''', targets=['t_ma.txt'])
    assert os.path.isfile('t_ma.txt')
