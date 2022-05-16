#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil
import sys

import pytest
from sos import execute_workflow


@pytest.mark.skipif(
    not shutil.which('singularity'),
    reason='Skip test because docker is not installed.')
def test_bash_in_singularity():
    '''Test action bash in singularity environment'''
    execute_workflow(r'''
        [0]
        #run:  container='shub://singularityhub/ubuntu'
        run:  container='docker://ubuntu', engine='singularity'
        echo 'Echo'
        ''')


@pytest.mark.skipif(
    not shutil.which('singularity') or sys.platform == 'win32' or
    'TRAVIS' in os.environ or 'APPVEYOR' in os.environ,
    reason='Skip test because docker is not installed.')
def test_singularity_build_linux_image(self):
    '''Test action singularity build, needs authentication so no travis test'''
    execute_workflow(r'''
        singularity_build: dest='lolcow.simg', sudo=True, notest=True
        Bootstrap: docker
        From: ubuntu:16.04
        %post
        apt-get -y update
        apt-get -y install fortune cowsay lolcat
        %environment
        export LC_ALL=C
        export PATH=/usr/games:$PATH
        %runscript
        fortune | cowsay | lolcat
        ''')


@pytest.mark.skipif(
    not shutil.which('singularity') or sys.platform == 'win32' or
    'TRAVIS' in os.environ or 'APPVEYOR' in os.environ,
    reason='Skip test because docker is not installed.')
def test_singularity_build_from_shub(self):
    # needs authentication and cannot test on travis
    execute_workflow(r'''
        singularity_build(src='shub://GodloveD/lolcow', dest='lolcow_shub.simg', sudo=True, notest=True, force=True)
        ''')


@pytest.mark.skipif(
    not shutil.which('singularity') or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_singularity_build_from_docker(self):
    execute_workflow(r'''
        singularity_build(src='docker://godlovedc/lolcow', dest='lolcow_docker.simg', sudo=True, notest=True, force=True)
        ''')
    # test handling simg with singularity
    execute_workflow(r'''
        run: container='lolcow_docker.simg'
        ls /
        ''')
