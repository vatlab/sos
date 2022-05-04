#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import pytest
from sos import execute_workflow
from sos.targets import file_target


def test_py_module(clear_now_and_after):
    '''Test target Py_Module'''
    clear_now_and_after('report.md')
    execute_workflow(r'''
[10]
depends: Py_Module('tabulate', autoinstall=True)
from tabulate import tabulate
table = [["Sun",696000,1989100000],["Earth",6371,5973.6],
    ["Moon",1737,73.5],["Mars",3390,641.85]]
report: output='report.md', expand=True
    {tabulate(table)}
''')
    assert file_target('report.md').target_exists('target')
    with open('report.md') as rep:
        assert rep.read().strip() == '''
-----  ------  -------------
Sun    696000     1.9891e+09
Earth    6371  5973.6
Moon     1737    73.5
Mars     3390   641.85
-----  ------  -------------
'''.strip()


def test_py_module_with_version():
    '''Test target Py_Module'''
    with pytest.raises(Exception):
        execute_workflow(r'''
[10]
depends: Py_Module('tabulate', version='2.0', autoinstall=True)
''')

        #
    with pytest.raises(Exception):
        execute_workflow(r'''
[10]
depends: Py_Module('tabulate>=2.0')
''')
        #

    with pytest.raises(Exception):
        execute_workflow(r'''
[10]
depends: Py_Module('tabulate==20.0')
''')
        #
    execute_workflow(r'''
[10]
depends: Py_Module('tabulate<2.0')
''')


def test_upgrade_py_module():
    '''Test upgrade py module #1246'''
    # first install tabulate == 0.7.5
    execute_workflow(r'''
[10]
depends: Py_Module('tabulate==0.7.5', autoinstall=True)
''')
    # test should pass
    execute_workflow(r'''
[10]
depends: Py_Module('tabulate'), Py_Module('tabulate==0.7.5')
''')
    # test for newer version should fail

    with pytest.raises(Exception):
        execute_workflow(r'''
[10]
depends: Py_Module('tabulate==0.8.3')
''')

    # auto install should work
    execute_workflow(r'''
[10]
depends: Py_Module('tabulate==0.8.3', autoinstall=True)
''')

    # test for old version should fail
    execute_workflow(r'''
        [10]
        depends: Py_Module('tabulate'), Py_Module('tabulate==0.8.3')
        ''')
    # test for old version should fail
    with pytest.raises(Exception):
        execute_workflow(r'''
        [10]
        depends: Py_Module('tabulate==0.7.5')
        ''')
