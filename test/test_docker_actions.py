#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import signal
import sys
import threading
from contextlib import contextmanager

import pytest

from sos import execute_workflow

try:
    import _thread
except Exception:
    import _dummy_thread as _thread

try:
    from sos.docker.client import SoS_DockerClient
    has_docker = True
except ImportError:
    print('Docker is not available')
    has_docker = False


class TimeoutException(Exception):

    def __init__(self, msg=''):
        self.msg = msg


@contextmanager
def time_limit(seconds, msg=''):
    if sys.platform == 'win32':
        # windows system does not have signal SIGALARM so we will
        # have to use a timer approach, which does not work as well
        # as the signal approach
        def timeout_func():
            #env.logger.error('Timed out for operation {}'.format(msg))
            _thread.interrupt_main()

        timer = threading.Timer(seconds, timeout_func)
        timer.start()
        try:
            yield
        except KeyboardInterrupt:
            # important: KeyboardInterrupt does not interrupt time.sleep()
            # because KeyboardInterrupt is handled by Python interpreter but
            # time.sleep() calls a system function.
            raise TimeoutException("Timed out for operation {}".format(msg))
        finally:
            # if the action ends in specified time, timer is canceled
            timer.cancel()
    else:

        def signal_handler(signum, frame):
            raise TimeoutException("Timed out for option {}".format(msg))

        signal.signal(signal.SIGALRM, signal_handler)
        signal.alarm(seconds)
        try:
            yield
        finally:
            signal.alarm(0)


try:
    with time_limit(2, 'check docker daemon'):
        has_docker = SoS_DockerClient().client is not None
except Exception as e:
    print(
        'Cannot connect to a docker daemon in 2 seconds. Assuming no docker environment.'
    )
    print(e)
    has_docker = False


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_bash_in_docker():
    '''Test action bash in docker environment'''
    execute_workflow(r'''
        [0]
        run:  container='docker://ubuntu'
        echo 'Echo'
        ''')


#    @pytest.mark.skipif(not has_docker or sys.platform != 'win32', reason='Skip test because docker is not installed.')
#    def testBatchScriptInDocker():
#        '''Test action powershell in docker environment'''
#        script = execute_workflow(r'''
#           [0]
#           run:  container='docker://microsoft/windowsservercore'
#           dir c:\
#           ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_sh_in_docker():
    '''Test action sh in docker environment'''
    # test docker
    script = execute_workflow(
        r'''
        [0]
        run: container='docker://ubuntu'
        echo 'Echo
        ''',
        options={'run_mode': 'dryrun'})
    with pytest.raises(Exception):
        execute_workflow(script)


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_docker_build_linux_image():
    '''Test action docker build'''
    execute_workflow(r'''
        [0]
        docker_build:  tag='test/docker_build'
        #
        # Super simple example of a Dockerfile
        #
        FROM ubuntu:latest
        MAINTAINER Andrew Odewahn "odewahn@oreilly.com"

        RUN apt-get update

        WORKDIR /home
        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_docker_build_linux_image_option_label_compress():
    '''Test action docker build'''
    # build with more options
    execute_workflow(r'''
        [0]
        docker_build:  tag='test/docker_build1', label='my_label', compress=True, memory='2G'
        #
        # Super simple example of a Dockerfile
        #
        FROM ubuntu:latest
        MAINTAINER Andrew Odewahn "odewahn@oreilly.com"

        WORKDIR /home
        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32' or 'TRAVIS' in os.environ,
    reason='Skip test because docker is not installed.')
def test_docker_build_linux_image_label_with_space():
    '''Test action docker build'''
    with pytest.raises(Exception):
        execute_workflow(r'''
            [0]
            docker_build:  tag='test/docker_build1', label='my label with space', compress=True, memory='2G'
            #
            # Super simple example of a Dockerfile
            #
            FROM ubuntu:latest
            MAINTAINER Andrew Odewahn "odewahn@oreilly.com"

            WORKDIR /home
            ''')


@pytest.mark.skipif(
    not has_docker or sys.platform != 'win32' or 'APPVEYOR' in os.environ,
    reason='Skip test because docker is not installed.')
def test_docker_build_windows_image():
    '''Test action docker build'''
    execute_workflow(r'''
        [0]
        docker_build:  tag='test/docker_build'
        # Indicates that the windowsservercore image will be used as the base image.
        FROM microsoft/windowsservercore

        # Metadata indicating an image maintainer.
        MAINTAINER someone@microsoft.com

        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_docker_image(fastq_files):
    '''Test docker_image option'''
    execute_workflow(r'''
        import os
        import glob
        [0]
        fastq_files = glob.glob('data/S20_R*.fastq')
        input_volume = os.path.dirname(os.path.abspath(fastq_files[0]))
        output_volume = os.getcwd()

        run: container='docker://compbio/ngseasy-fastqc:1.0-r001',
            volumes=[f"{input_volume}:/input_data", f"{output_volume}:/output_data"]

            ls -l /input_data
            /usr/local/bin/fastqc /input_data/*.fastq --outdir /output_data
        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed, or in travis, which failed for unknown reason'
)
def test_docker_image_from_file():
    '''Test docker_image load from a file.'''
    # image from a saved file
    execute_workflow(r'''
        [0]
        run:   container='docker://blang/busybox-bash'

        [1]
        run:
            docker save blang/busybox-bash > hello.tar
            docker rmi -f blang/busybox-bash

        [2]
        run: container='docker://blang/busybox-bash', docker_file = 'hello.tar'

            echo "a"
        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_docker_script_action():
    '''Test action sh in docker environment'''
    # test docker
    execute_workflow(r'''
        [0]
        script: container='docker://ubuntu', args='{script}'
        echo 'Echo'
        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_port_option():
    '''Test use of option port in action'''
    execute_workflow(r'''
        [0]
        run:  container='ubuntu', port=True
        echo 'Echo'
        ''')
    #
    execute_workflow(r'''
        [0]
        run:  container='ubuntu', port=2345
        echo 'Echo'
        ''')


@pytest.mark.skipif(
    not has_docker or sys.platform == 'win32',
    reason='Skip test because docker is not installed.')
def test_auto_output_mount():
    '''Test use of option port in action'''
    execute_workflow(r'''
        output: '../data/1.txt'
        run: container='ubuntu'
        sh: expand=True
            touch {_output}
        ''')
    assert os.path.isfile('../data/1.txt')
