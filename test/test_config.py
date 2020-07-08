#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import getpass
import subprocess

from sos import execute_workflow
from sos._version import __version__
from sos.utils import env, load_config_files

# if the test is imported under sos/test, test interacive executor


def test_command_line():
    '''Test command line arguments'''
    assert subprocess.call(
        'sos config -h',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    assert subprocess.call(
        'sos config -g --get',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    assert subprocess.call(
        'sos config',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    assert subprocess.call(
        'sos config --get',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    assert subprocess.call(
        'sos config -g --set a 5',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    assert subprocess.call(
        'sos config --get a',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    assert subprocess.call(
        'sos config -g --unset a',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0

def test_config_set():
    '''Test interpolation of config'''
    assert subprocess.call(
        'sos config --set cut 0.5 -c myconfig.yml',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    load_config_files('myconfig.yml')
    assert env.sos_dict['CONFIG']['cut'] == 0.5
    #
    assert subprocess.call(
        'sos config --set cut1 0.5 2 3 -c myconfig.yml',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    load_config_files('myconfig.yml')
    assert env.sos_dict['CONFIG']['cut1'] == [0.5, 2, 3]
    #
    assert subprocess.call(
        'sos config --set cut2 a3 -c myconfig.yml',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    load_config_files('myconfig.yml')
    assert env.sos_dict['CONFIG']['cut2'] == 'a3'
    #
    assert subprocess.call(
        'sos config --set cut3 a b c -c myconfig.yml',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    load_config_files('myconfig.yml')
    assert env.sos_dict['CONFIG']['cut3'] == ['a', 'b', 'c']
    #
    assert subprocess.call(
        '''sos config --set cut4 "{'A': 123}" -c myconfig.yml''',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    load_config_files('myconfig.yml')
    assert env.sos_dict['CONFIG']['cut4'] == {'A': 123}

def test_interpolate():
    '''Test interpolation of config'''
    assert subprocess.call(
        '''sos config --set me '{user_name}@my' -c myconfig.yml''',
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0
    load_config_files('myconfig.yml')
    assert env.sos_dict['CONFIG']['me'] == f'{getpass.getuser().lower()}@my'


def test_global_vars(config_factory):
    '''Test SoS defined variables'''
    execute_workflow("[0]", options={'mode': 'dryrun'})

    assert env.sos_dict['SOS_VERSION'] == __version__
    assert isinstance(env.sos_dict['CONFIG'], dict)

    cfg = config_factory('my_config: 5')

    execute_workflow("[0]", options={'config_file': cfg})
    assert env.sos_dict['CONFIG']['my_config'] == 5
