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

import os
import unittest
import time
import shutil
import signal
import sys
import threading
try:
    import _thread
except:
    import _dummy_thread as _thread
from contextlib import contextmanager

from pysos.sos_script import SoS_Script
from pysos.utils import env
from pysos.sos_eval import  Undetermined
from pysos.actions import DockerClient
from docker.errors import DockerException
from pysos.sos_executor import Base_Executor, ExecuteError
from pysos.target import FileTarget

import socket
def internet_on(host='8.8.8.8', port=53, timeout=3):
    '''Test if internet is connected '''
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except Exception as e:
        print(e)
        return False

with_network = internet_on()

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
        env.reset()
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            FileTarget(f).remove('both')

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

    def testSoSAction(self):
        '''Test sos_action decorator'''
        script = SoS_Script(r"""
from pysos.actions import SoS_Action

@SoS_Action(run_mode='run')
def func_run():
    return 1

@SoS_Action(run_mode=['run', 'prepare'])
def func_both():
    return 1

[0: shared=('b', 'c')]
b=func_run()
c=func_both()
""")
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertTrue(isinstance(env.sos_dict['b'], Undetermined))
        self.assertEqual(env.sos_dict['c'], 1)
        #
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], 1)
        self.assertEqual(env.sos_dict['c'], 1)

    def testGetOutput(self):
        '''Test utility function get_output'''
        script = SoS_Script(r"""
[0: shared='ret']
ret = get_output('echo blah')
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['ret'], 'blah\n')
        #
        script = SoS_Script(r"""
[0: shared='ret']
ret = get_output('echo blah', show_command=True)
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['ret'], '$ echo blah\nblah\n')
        #
        script = SoS_Script(r"""
[0: shared='ret']
ret = get_output('echo blah', show_command=True, prompt='% ')
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['ret'], '% echo blah\nblah\n')
        #
        script = SoS_Script(r"""
[0]
get_output('catmouse')
""")
        wf = script.workflow()
        # should fail in prepare mode
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        #
        script = SoS_Script(r"""
[0]
ret = get_output('cat -h')
""")
        wf = script.workflow()
        # this should give a warning and return false
        self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testCheckCommand(self):
        '''Test action check_command'''
        script = SoS_Script(r"""
[0]
check_command('cat')
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        #
        script = SoS_Script(r"""
[0]
check_command('catmouse')
""")
        wf = script.workflow()
        # should fail in prepare mode
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        wf = script.workflow()
        #
        # should also check command with option
        script = SoS_Script(r"""
[0]
check_command('ls -l')
""")
        wf = script.workflow()
        # this should pass
        Base_Executor(wf).run()
        #
        script = SoS_Script(r"""
[0]
fail_if(check_command('cat -h') != 0, 'command return non-zero')
""")
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        # check check_command is the command is stuck
        script = SoS_Script(r"""
[0]
check_command('sleep 4')
""")
        wf = script.workflow()
        # this should pass
        Base_Executor(wf).prepare()
        #
        script = SoS_Script(r"""
[0]
fail_if(check_command('sleep 4') != 0, 'Command time out')
""")
        wf = script.workflow()
        # this should pass
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        # test reading this file
        script = SoS_Script(r"""
[0]
check_command('cat test_actions.py', 'abcde' + 'fgh')
""")
        wf = script.workflow()
        # should raise an error
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #
        script = SoS_Script(r"""
check_command('cat test_actions.py', 'testSearchOutput')
""")
        wf = script.workflow()
        Base_Executor(wf).run()


    def testFailIf(self):
        '''Test action fail if'''
        self.touch('a.txt')
        script = SoS_Script(r"""
[0]
input: 'a.txt'
fail_if(len(input) == 1)
""")
        wf = script.workflow()
        # should fail in prepare mode
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        script = SoS_Script(r"""
[0]
input: 'a.txt', 'b.txt'
fail_if(len(input) == 1)
""")
        wf = script.workflow()

    def testWarnIf(self):
        '''Test action fail if'''
        script = SoS_Script(r"""
[0]
warn_if(input is None, 'Expect to see a warning message')
""")
        wf = script.workflow()
        # should see a warning message.
        Base_Executor(wf).prepare()
        #self.assertRaises(ExecuteError, Base_Executor(wf).run)
        script = SoS_Script(r"""
[0]
input: 'a.txt', 'b.txt'
warn_if(len(input) == 1)
""")
        wf = script.workflow()

    def testStopIf(self):
        '''Test action stop_if'''
        script = SoS_Script(r'''
[0: shared='result']
rep = range(20)
result = []
input: for_each='rep'

stop_if(_rep > 10)
result.append(_rep)
''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertEqual(env.sos_dict['result'], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    def testRun(self):
        '''Test action run'''
        script = SoS_Script(r'''
[0]
run:
echo 'Echo'
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        script = SoS_Script(r'''
[0]
run:
echo 'Echo
''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testBash(self):
        '''Test action bash'''
        script = SoS_Script(r'''
[0]
bash:
echo 'Echo'
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        script = SoS_Script(r'''
[0]
bash:
echo 'Echo
''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
    
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


    def testSh(self):
        '''Test action run'''
        script = SoS_Script(r'''
[0]
sh:
echo 'Echo'
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        script = SoS_Script(r'''
[0]
sh:
echo 'Echo
''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testArgs(self):
        '''Test args option of scripts'''
        FileTarget('a.txt').remove('both')
        script = SoS_Script(r'''
[0]
sh: args='-n'
    touch a.txt
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertFalse(os.path.exists('a.txt'))


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
        # this should give us a warning if RAM is less than 4G
        Base_Executor(wf).prepare()


    def testCsh(self):
        '''Test action csh'''
        if not shutil.which('csh'):
            return
        script = SoS_Script(r'''
[0]
csh:
     foreach color (red orange yellow green blue)
        echo $color
     end
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # no test for docker because standard distributions do not have csh

    def testTcsh(self):
        '''Test action tcsh'''
        if not shutil.which('tcsh'):
            return
        script = SoS_Script(r'''
[0]
tcsh:
     foreach color (red orange yellow green blue)
        echo $color
     end
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # no test for docker because standard distributions do not have tcsh

    def testZsh(self):
        '''Test action zsh'''
        if not shutil.which('zsh'):
            return
        script = SoS_Script(r'''
[0]
zsh:
echo "Hello World!", $SHELL
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # cannot test in docker because no first-tier repository
        # provides zsh.

    def testPython(self):
        '''Test python command. This might fail if python3 is the
        default interpreter'''
        script = SoS_Script(r'''
[0]
python:
a = {'1': 2}
print(a)
''')
        wf = script.workflow()
        Base_Executor(wf).run()

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

    def testPython3(self):
        script = SoS_Script(r'''
[0]
python3:
a = {'1', '2'}
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

    def testPerl(self):
        '''Test action ruby'''
        script = SoS_Script(r'''
[0]
perl:
use strict;
use warnings;

print "hi NAME\n";
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


    def testRuby(self):
        '''Test action ruby'''
        if not shutil.which('ruby'):
            return True
        script = SoS_Script(r'''
[0]
ruby:
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

    def testNode(self):
        '''Test action ruby'''
        if not shutil.which('node'):
            return
        script = SoS_Script(r'''
[0]
node:
var args = process.argv.slice(2);
console.log('Hello ' + args.join(' ') + '!');
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        #
        script = SoS_Script(r'''
[0]
JavaScript:
var args = process.argv.slice(2);
console.log('Hello ' + args.join(' ') + '!');
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

    def testR(self):
        '''Test action JavaScript'''
        if not shutil.which('R'):
            return 
        script = SoS_Script(r'''
[0]
R:
nums = rnorm(25, mean=100, sd=15)
mean(nums)
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

    @unittest.skipIf(not with_network, 'Skip test because of no internet connection')
    def testCheckRLibrary(self):
        '''Test action check_R_library'''
        if not shutil.which('R'):
            return 
        script = SoS_Script(r'''
[0]
check_R_library('edgeR')
''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        script = SoS_Script(r'''
[0]
check_R_library('stephens999/ashr')
''')
        wf = script.workflow()
        Base_Executor(wf).run()

        script = SoS_Script(r'''
[0]
check_R_library('edgeRRRR')
''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).prepare)

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
    volumes=['${input_volume}:/input_data', '${output_volume}:/output_data']

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

    @unittest.skipIf(not with_network, 'Skip test because of no internet connection')
    def testDownload(self):
        '''Test download of resources'''
        if not os.path.isdir('tmp'):
            os.makedirs('tmp')
        #
        for name in ['hapmap_ASW_freq.ann', 'hapmap_ASW_freq-hg18_20100817.DB.gz', 'hapmap_CHB_freq.ann',
                'vt_quickStartGuide.tar.gz']:
            if os.path.isfile(os.path.join('tmp', name)):
                os.remove(os.path.join('tmp', name))
        # test decompress tar.gz file
        script = SoS_Script(r'''
[0]
download(['http://bioinformatics.mdanderson.org/Software/VariantTools/repository/snapshot/vt_quickStartGuide.tar.gz'],
    dest_dir='tmp', decompress=True)
''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertTrue(os.path.isfile('tmp/snapshot.proj'))
        self.assertTrue(os.path.isfile('tmp/snapshot_genotype.DB'))
        #
        # testing the download of single file
        #
        script = SoS_Script(r'''
[0]
download: dest_file='tmp/test.ann'
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/hapmap_ASW_freq.ann
''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertTrue(os.path.isfile('tmp/test.ann'))
        # test option dest_dir
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp'
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/hapmap_ASW_freq.ann
''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertTrue(os.path.isfile('tmp/hapmap_ASW_freq.ann'))
        #
   
        # this will take a while
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/non-existing.gz
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/hapmap_ASW_freq.ann
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/hapmap_ASW_freq-hg18_20100817.DB.gz
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/hapmap_CHB_freq.ann	
''')
        start = time.time()
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).prepare)
        self.assertTrue(os.path.isfile('tmp/hapmap_ASW_freq-hg18_20100817.DB'))
        self.assertGreater(time.time() - start, 3)
        # this will be fast
        start = time.time()
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).prepare)
        self.assertLess(time.time() - start, 3)
        # 
        # test decompress tar.gz file
        script = SoS_Script(r'''
[0]
GATK_URL = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/'
download: dest_dir='tmp', decompress=True
    ${GATK_URL}/1000G_omni2.5.hg19.sites.vcf.idx.gz
    ${GATK_URL}/1000G_omni2.5.hg19.sites.vcf.idx.gz.md5
''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        #
        shutil.rmtree('tmp')

    def testPandoc(self):
        '''Test action pandoc'''
        if not shutil.which('pandoc'):
            return
        script = SoS_Script(r'''
[10]

report:
## Some random figure

Generated by matplotlib


[100]
# generate report
output: 'myreport.html'
pandoc(output=_output[0], to='html')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('myreport.html'))
        #
        FileTarget('myreport.html').remove('both')


    def testOptionRunMode(self):
        '''Testing run mode option for action'''
        FileTarget('a.txt').remove('both')
        #
        script = SoS_Script(r'''
[10]
run: run_mode='prepare'
    touch a.txt

''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertTrue(os.path.isfile('a.txt'))
        os.remove('a.txt')
        #
        script = SoS_Script(r'''
[10]
run: run_mode='prepare'
    touch a.txt

''')
        wf = script.workflow()
        Base_Executor(wf).prepare()
        self.assertTrue(os.path.isfile('a.txt'))
        FileTarget('a.txt').remove('both')

    def testSoSRun(self):
        '''Test action sos_run with keyword parameters'''
        for f in ['0.txt', '1.txt']:
            FileTarget(f).remove('both')
        script = SoS_Script(r'''
[A]
parameter: num=5
sh:
    touch ${num}.txt

[batch]
for k in range(2):
    sos_run('A', num=k)
''')
        env.verbosity=3
        wf = script.workflow('batch')
        Base_Executor(wf).run()
        for f in ['0.txt', '1.txt']:
            self.assertTrue(FileTarget(f).exists())
            FileTarget(f).remove('both')

if __name__ == '__main__':
    unittest.main()
