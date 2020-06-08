#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import os
import shutil
import unittest
import pytest

from sos import execute_workflow
from sos.parser import SoS_Script
from sos.targets import file_target
from sos.utils import env
from sos.workflow_executor import Base_Executor


class TestActions(unittest.TestCase):

    def setUp(self):
        env.reset()
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            file_target(f).unlink()

    def touch(self, files):
        '''create temporary files'''
        if isinstance(files, str):
            files = [files]
        #
        for f in files:
            with open(f, 'w') as tmp:
                tmp.write('test')
        #
        self.temp_files.extend(files)

    def test_bash(self):
        '''Test action bash'''
        script = execute_workflow(r'''
            [0]
            bash:
                echo 'Echo'
        ''')
        with pytest.raises(Exception):
            execute_workflow(script)

    def test_sh(self):
        '''Test action run'''
        script = execute_workflow(r'''
            [0]
            sh:
                echo 'Echo'
        ''')
        with pytest.raises(Exception):
            execute_workflow(script)

    def test_csh(self):
        '''Test action csh'''
        if not shutil.which('csh'):
            return
        script = execute_workflow(r'''
            [0]
            csh:
                foreach color (red orange yellow green blue)
                    echo $color
                end
        ''')

    def test_tcsh(self):
        '''Test action tcsh'''
        if not shutil.which('tcsh'):
            return
        script = execute_workflow(r'''
            [0]
            tcsh:
                foreach color (red orange yellow green blue)
                    echo $color
                end
        ''')

    def test_zsh(self):
        '''Test action zsh'''
        if not shutil.which('zsh'):
            return
        script = execute_workflow(r'''
            [0]
            zsh:
                echo "Hello World!", $SHELL
        ''')

    def test_args(self):
        '''Test args option of scripts'''
        if os.path.isfile('a.txt'):
            file_target('a.txt').unlink()
        script = execute_workflow(r'''
            [0]
            sh: args='-n {filename:q}'
                touch a.txt
        ''')

        assert not os.path.exists('a.txt')


if __name__ == '__main__':
    unittest.main()
