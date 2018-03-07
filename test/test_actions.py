#!/usr/bin/env python
#
# This file is part of Script of Scripts (SoS), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
import sys
import glob
import unittest
import shutil

from sos.parser import SoS_Script
from sos.utils import env
from sos.eval import  Undetermined
from sos.workflow_executor import Base_Executor, ExecuteError
from sos.targets import file_target

import socket
def internet_on(host='8.8.8.8', port=80, timeout=3):
    '''Test if internet is connected '''
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except Exception as e:
        print(e)
        return False

with_network = internet_on()

def multi_attempts(fn):
    def wrapper(*args, **kwargs):
        for n in range(4):
            try:
                fn(*args, **kwargs)
                break
            except Exception:
                if n > 1:
                    raise
    return wrapper

class TestActions(unittest.TestCase):
    def setUp(self):
        env.reset()
        self.temp_files = []

    def tearDown(self):
        for f in self.temp_files:
            file_target(f).remove('both')

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

    def testAcceptableArgs(self):
        '''test acceptable args of options'''
        script = SoS_Script(r"""
run: unrecog=1
    echo 'a'
""")
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)

    def testGetOutput(self):
        '''Test utility function get_output'''
        script = SoS_Script(r"""
[0: shared='ret']
ret = get_output('echo blah')
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        # use strip because there would be \r\n under windows
        self.assertEqual(env.sos_dict['ret'].strip(), 'blah')
        #
        script = SoS_Script(r"""
[0: shared='ret']
ret = get_output('echo blah', show_command=True)
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        self.assertEqual([x.strip() for x in env.sos_dict['ret'].splitlines()],
            ['$ echo blah', 'blah'])
        #
        script = SoS_Script(r"""
[0: shared='ret']
ret = get_output('echo blah', show_command=True, prompt='% ')
""")
        wf = script.workflow()
        # should be ok
        Base_Executor(wf).run()
        self.assertEqual([x.strip() for x in env.sos_dict['ret'].splitlines()],
            ['% echo blah', 'blah'])
        #
        script = SoS_Script(r"""
[0]
get_output('catmouse')
""")
        wf = script.workflow()
        # should fail in dryrun mode
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

    def testFailIf(self):
        '''Test action fail if'''
        self.touch('a.txt')
        script = SoS_Script(r"""
[0]
input: 'a.txt'
fail_if(len(input) == 1)
""")
        wf = script.workflow()
        # should fail in dryrun mode
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
        Base_Executor(wf).run(mode='dryrun')
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
        Base_Executor(wf).run()
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
        #
        # Under windows, echo 'Echo is perfectly OK
        if sys.platform == 'win32':
            return
        script = SoS_Script(r'''
[0]
run:
echo 'Echo
''')
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)

    def testRunWithShebang(self):
        script = SoS_Script(r'''
[0]
run:
    #!/usr/bin/env python
    print('Echo')
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

    @multi_attempts
    @unittest.skipIf(not with_network or 'TRAVIS' in os.environ, 'Skip test because of no internet connection or in travis test')
    def testDownload(self):
        '''Test download of resources'''
        if not os.path.isdir('tmp'):
            os.makedirs('tmp')
        #
        for name in os.listdir('tmp'):
            if os.path.isfile(os.path.join('tmp', name)):
                os.remove(os.path.join('tmp', name))
        # test decompress tar.gz file
        script = SoS_Script(r'''
[0]
download(['ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.gz'],
    dest_dir='tmp', decompress=True)
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isdir('tmp/pcre-8.41'))
        self.assertTrue(os.path.isfile('tmp/pcre-8.41/pcre_get.c'))
        #
        # testing the download of single file
        #
        script = SoS_Script(r'''
[0]
download: dest_file='tmp/pcre-8.41.zip.sig'
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip.sig
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('tmp/pcre-8.41.zip.sig'))
        # test option dest_dir
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp'
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip.sig
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('tmp/pcre-8.41.zip.sig'))
        #
        # this will take a while
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True

    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/non-existing.gz
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip.sig
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.gz
''')
        #start = time.time()
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        self.assertTrue(os.path.isfile('tmp/pcre-8.41/pcre_get.c'))
        #self.assertGreater(time.time() - start, 3)
        # this will be fast
        #start = time.time()
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        #self.assertLess(time.time() - start, 3)
        #
        # test decompress tar.gz, .zip and .gz files
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.gz
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.bz2
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip.sig
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # run in build mode
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.gz
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.tar.bz2
    ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.41.zip.sig
''')
        wf = script.workflow()
        Base_Executor(wf, config={'sig_mode': 'build'}).run()
        #
        shutil.rmtree('tmp')


    def testPandoc(self):
        '''Test action pandoc'''
        if not shutil.which('pandoc'):
            return
        script = SoS_Script(r'''
[10]

report: output='report.md'
## Some random figure

Generated by matplotlib


[100]
# generate report
output: 'myreport.html'
pandoc(input='report.md', output=_output[0])
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('myreport.html'))
        #
        file_target('myreport.html').remove('both')
        # pandoc with specified input.
        script = SoS_Script(r'''
[10]
report: output='a.md'
## Some random figure

Generated by matplotlib


[100]
# generate report
output: 'myreport.html'
pandoc(input='a.md', output=_output[0])
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('myreport.html'))
        file_target('myreport.html').remove()
        #
        # another case is no output
        script = SoS_Script(r'''
[10]
report: output='a.md'
## Some random figure

Generated by matplotlib


[100]
# generate report
pandoc(input='a.md')
''')
        wf = script.workflow()
        Base_Executor(wf).run()

        # test acceptance of a list of input filenames
        #
        script = SoS_Script(r'''
[10]
report: output='default_10.md'
A_10

[20]
report: output='default_20.md'
A_20

[100]
# generate report
pandoc(input=['default_10.md', 'default_20.md'], output='output.html')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for f in ['default_10.md', 'default_20.md', 'output.html']:
            self.assertTrue(file_target(f).target_exists())
            file_target(f).remove()



    def testReport(self):
        '''Test action report'''
        script = SoS_Script(r'''
[A]
parameter: num=5
report: output='report.txt', expand=True
    touch {num}.txt

''')
        # output to a file
        file_target('report.txt').remove('both')
        wf = script.workflow()
        # run twice
        Base_Executor(wf, args=['--num', '7']).run()
        Base_Executor(wf, args=['--num', '5']).run()
        with open('report.txt') as report:
            self.assertEqual(report.read(), 'touch 5.txt\n\n')
        # test overwrite
        file_target('report.txt').remove('both')
        script = SoS_Script(r'''
[A]
report: output='report.txt', expand=True
    {step_name}

[A_10]
report: output='report.txt', expand=True
    {step_name}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # output to a file
        # run twice
        Base_Executor(wf).run()
        with open('report.txt') as report:
            self.assertEqual(report.read(), 'A_10\n\n')

        #
        # test input from another file
        file_target('report.txt').remove()
        script = SoS_Script(r'''
[A_1]
run: output='a.txt'
    echo something > a.txt
report(input='a.txt', output='out.txt')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for name in ('a.txt', 'out.txt'):
            with open(name) as report:
                self.assertEqual(report.read().strip(), 'something')
            file_target(name).remove()
        #
        script = SoS_Script(r'''
[A_1]
run: output='a.txt'
    echo something > a.txt

[A_2]
run: output='b.txt'
    echo something else > b.txt

[A_3]
report(input=['a.txt', 'b.txt'], output='out.txt')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for name in ('a.txt', 'b.txt', 'out.txt'):
            self.assertTrue(file_target(name).target_exists())
            file_target(name).remove()
        #
        # test report to other types of output: path
        script = SoS_Script(r'''
[A_1]
report: output=path('a.txt')
    something
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # test report to other types of output: paths
        script = SoS_Script(r'''
[A_1]
report: output=paths('a.txt')
    something
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # test report to other types of output: file_target
        script = SoS_Script(r'''
[A_1]
output: 'a.txt'
report: output=_output[0]
    something
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # test report to other types of output: sos_targets
        script = SoS_Script(r'''
[A_1]
output: 'a.txt'
report: output=_output
    something
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        #
        script = SoS_Script(r'''
[A_1]
output: 'a.txt', 'b.txt'
report: output=_output
    something 
''')
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)

    def testOptionWorkdir(self):
        '''Test option workdir of tasks'''
        if not os.path.isdir('temp_wdr'):
            os.mkdir('temp_wdr')
        with open(os.path.join('temp_wdr', 'a.txt'), 'w') as tmp:
            tmp.write('hello')
        script = SoS_Script(r'''
[A_1]
run: workdir='temp_wdr'
  cp -f a.txt a2.txt
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(file_target(os.path.join('temp_wdr', 'a2.txt')).target_exists())
        with open(os.path.join('temp_wdr', 'a.txt')) as tmp:
            self.assertEqual('hello', tmp.read())

    def testActionScript(self):
        '''Test action script'''
        script = SoS_Script(r'''
[A_1]
script: interpreter='python'
  with open('something.txt', 'w') as tmp:
    tmp.write('something')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(file_target('something.txt').target_exists())
        with open('something.txt') as tmp:
            self.assertEqual('something', tmp.read())

    def testRegenerateReport(self):
        '''Testing the regeneration of report once is needed. The problem
        here is the 'input' parameter of report.'''
        script = SoS_Script(r'''
[A_1]
output: 'a1.txt', 'a1.md'
run:
    echo 'a1' >> a1.txt

report: output='a1.md'
    a1

[A_2]
output: 'a2.txt', 'a2.md'
run:
    echo 'a2' >> a2.txt
report: output='a2.md'
    a2

[A_3]
input: 'a1.md', 'a2.md'
output: 'out.md'
report:     input=['a1.md', 'a2.md'], output='out.md'
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        with open('a1.md') as a:
            self.assertEqual(a.read(), 'a1\n\n')
        with open('a2.md') as a:
            self.assertEqual(a.read(), 'a2\n\n')
        with open('out.md') as a:
            self.assertEqual(a.read(), 'a1\n\na2\n\n')
        for name in ('a1.md', 'a2.md', 'out.md'):
            file_target(name).remove('target')
        wf = script.workflow()
        Base_Executor(wf).run()


    def testActiveActionOption(self):
        '''Test the active option of actions'''
        # disallow
        self.assertRaises(Exception, SoS_Script, '''
[1]
rep = range(5)
input: for_each = 'rep'
# ff should change and be usable inside run
ff = f"{_rep}.txt"
run:  expand=True, active=1,2
echo {ff}
touch temp/{ff}
''')
        #
        for active, result in [
            ('0', ['temp/0.txt']),
            ('-1', ['temp/4.txt']),
            ('(1,2)', ['temp/1.txt', 'temp/2.txt']),
            ('[2,3]', ['temp/2.txt', 'temp/3.txt']),
            ('(0,2,4)', ['temp/0.txt', 'temp/2.txt', 'temp/4.txt']),
            ('slice(1,None)', ['temp/1.txt', 'temp/2.txt', 'temp/3.txt', 'temp/4.txt']),
            ('slice(1,-2)', ['temp/1.txt', 'temp/2.txt']),
            ('slice(None,None,2)', ['temp/0.txt', 'temp/2.txt', 'temp/4.txt']),
            ('True', ['temp/0.txt', 'temp/1.txt', 'temp/2.txt', 'temp/3.txt', 'temp/4.txt']),
            ('False', []),
            ]:
            if os.path.isdir('temp'):
                shutil.rmtree('temp')
            os.mkdir('temp')
            # test first iteration
            script = SoS_Script(('''
[1]
rep = range(5)
input: for_each = 'rep'
# ff should change and be usable inside run
ff = f"{_rep}.txt"
run:  expand=True, active=%s
echo {ff}
touch temp/{ff}
''' % active).replace('/', os.sep))
            wf = script.workflow()
            env.config['sig_mode'] = 'force'
            env.config['wait_for_task'] = True
            Base_Executor(wf).run()
            files = list(glob.glob(os.path.join('temp', '*.txt')))
            self.assertEqual(sorted(files), sorted([x.replace('/', os.sep) for x in result]))
            #
            # test last iteration
            shutil.rmtree('temp')



if __name__ == '__main__':
    unittest.main()
