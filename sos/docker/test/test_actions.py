#!/usr/bin/env python
#
# This file is part of Script of Scripts (SoS), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import unittest
import signal
import sys
import os
import threading
try:
    import _thread
except:
    import _dummy_thread as _thread
from contextlib import contextmanager

from sos.sos_script import SoS_Script
from sos.utils import env
from sos.docker.client import DockerClient
from docker.errors import DockerException
from sos.sos_executor import Base_Executor, ExecuteError
from sos.target import FileTarget

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
        has_docker = DockerClient().client is not None
except (TimeoutException, DockerException) as e:
    print('Cannot connect to a docker daemon in 2 seconds. Assuming no docker environment.')
    has_docker = False

class TestActions(unittest.TestCase):
    def setUp(self):
        self.olddir = os.getcwd()
        try:
            # this only works with nose, but is also
            # only needed by nose
            os.chdir(os.path.dirname(__file__))
        except:
            pass
        env.reset()
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            FileTarget(f).remove('both')
        os.chdir(self.olddir)

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
    
    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testBashInDocker(self):
        '''Test action bash in docker environment'''
        script = SoS_Script(r'''
[0]
bash:  docker_image='ubuntu'
echo 'Echo'
''')
        wf = script.workflow()
        Base_Executor(wf).run()


    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testShInDocker(self):
        '''Test action sh in docker environment'''
        # test docker
        script = SoS_Script(r'''
[0]
sh: docker_image='ubuntu'
echo 'Echo
''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        Base_Executor(wf).dryrun()


    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testPythonInDocker(self):
        '''Test action python in docker environment'''
        script = SoS_Script(r'''
[0]
python:  docker_image='python'
a = {'1': 2}
print(a)
''')
        wf = script.workflow()
        Base_Executor(wf).run()


    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testPythonsInDocker(self):
        '''Test action pythons in docker environment'''
        script = SoS_Script(r'''
[0] 
python3: docker_image='python'
#!/usr/bin/env python3
a = {'1', '2'}
print(a)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testPerlInDocker(self):
        '''Test action perl in docker environment'''
        script = SoS_Script(r'''
[0]
perl: docker_image='ubuntu'
use strict;
use warnings;

print "hi NAME\n";
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testRubyInDocker(self):
        '''Test action ruby in docker environment'''
        script = SoS_Script(r'''
[0]
ruby: docker_image='ruby'
line1 = "Cats are smarter than dogs";
line2 = "Dogs also like meat";

if ( line1 =~ /Cats(.*)/ )
  puts "Line1 contains Cats"
end
if ( line2 =~ /Cats(.*)/ )
  puts "Line2 contains  Dogs"
end
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testNodeInDocker(self):
        '''Test action node in docker environment'''
        script = SoS_Script(r'''
[0]
node: docker_image='node'

var args = process.argv.slice(2);
console.log('Hello ' + args.join(' ') + '!');
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        #
        script = SoS_Script(r'''
[0]
JavaScript: docker_image='node'

var args = process.argv.slice(2);
console.log('Hello ' + args.join(' ') + '!');
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testRInDocker(self):
        '''Test action R in docker environment'''
        script = SoS_Script(r'''
[0]
R: docker_image='r-base'
nums = rnorm(25, mean=100, sd=15)
mean(nums)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testDockerBuild(self):
        '''Test action docker build'''
        script = SoS_Script(r'''
[0]
docker_build:  tag='test/docker_build'
#
# Super simple example of a Dockerfile
#
FROM ubuntu:latest
MAINTAINER Andrew Odewahn "odewahn@oreilly.com"

RUN apt-get update
RUN apt-get install -y python python-pip wget
RUN pip install Flask

WORKDIR /home
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testDockerImage(self):
        '''Test docker_image option'''
        script = SoS_Script(r'''
[0]
fastq_files = glob.glob('data/*.fastq')
input_volume = os.path.dirname(fastq_files[0])
output_volume = os.getcwd()

run: docker_image='compbio/ngseasy-fastqc:1.0-r001', 
    volumes=["${input_volume}:/input_data", "${output_volume}:/output_data"]

    ls -l /input_data
    /usr/local/bin/fastqc /input_data/*.fastq --outdir /output_data
''')
        wf = script.workflow()
        Base_Executor(wf).run()

    @unittest.skipIf(not has_docker, 'Skip test because docker is not installed.')
    def testDockerImageFromFile(self):
        '''Test docker_image load from a file.'''
        # image from a saved file
        script = SoS_Script(r'''
[0]
run:   docker_image='blang/busybox-bash'

[1]
run:
    docker save blang/busybox-bash > hello.tar
    docker rmi -f blang/busybox-bash

[2]
run: docker_image='blang/busybox-bash', docker_file = 'hello.tar'

    echo "a"
''')
        wf = script.workflow()
        Base_Executor(wf).run()

if __name__ == '__main__':
    unittest.main()
