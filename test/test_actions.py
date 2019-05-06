#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import glob
import os
import shutil
import socket
import sys
import time
import unittest

from sos.parser import SoS_Script
from sos.targets import file_target
from sos.utils import env
from sos.workflow_executor import Base_Executor


def internet_on(host='8.8.8.8', port=80, timeout=3):
    '''Test if internet is connected '''
    try:
        socket.setdefaulttimeout(timeout)
        socket.create_connection(("www.google.com", 80))
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
            if file_target(f).exists():
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
        self.assertRaises(Exception, Base_Executor(wf).run)
        #
        #
        script = SoS_Script(r"""
[0]
ret = get_output('cat -h')
""")
        wf = script.workflow()
        # this should give a warning and return false
        self.assertRaises(Exception, Base_Executor(wf).run)

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
        self.assertRaises(Exception, Base_Executor(wf).run)
        script = SoS_Script(r"""
[0]
input: 'a.txt', 'b.txt'
fail_if(len(input) == 2)
""")
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)

    def testDelayedFailIf(self):
        # test fail_if of killing another running substep
        script = SoS_Script(r"""
import time

[10]
time.sleep(8)

[20]
input: None
time.sleep(2)
fail_if(True)
""")
        st = time.time()
        wf = script.workflow()
        self.assertRaises(Exception,
                          Base_Executor(wf, config={
                              'max_procs': 3
                          }).run)
        self.assertTrue(time.time() - st >= 8,
                'Test test should fail only after step 10 is completed')

    def testDelayedFailIfFromNestedWorkflow(self):
        # test fail_if of killing another running substep
        script = SoS_Script(r"""
import time

[default]
sos_run('a')

[a_10]
time.sleep(8)

[a_20]
input: None
time.sleep(2)
fail_if(True)
""")
        st = time.time()
        wf = script.workflow()
        self.assertRaises(Exception,
                          Base_Executor(wf, config={
                              'max_procs': 3
                          }).run)
        self.assertTrue(time.time() - st >= 8,
                'Test test should fail only after step 10 is completed')

    def testWarnIf(self):
        '''Test action fail if'''
        script = SoS_Script(r"""
[0]
warn_if(input is None, 'Expect to see a warning message')
""")
        wf = script.workflow()
        # should see a warning message.
        Base_Executor(wf).run(mode='dryrun')
        #self.assertRaises(Exception, Base_Executor(wf).run)
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
        self.assertEqual(env.sos_dict['result'],
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        # stop_if should not be treated as error so the previously
        # generated output file will be removed
        for rep in range(2):
            file = f'test_stop_if_{rep}.txt'
            if os.path.isfile(file):
                os.remove(file)
        script = SoS_Script(r'''
[10]
rep = range(2)
input: for_each='rep'
output: f'test_stop_if_{_rep}.txt'

_output.touch()
stop_if(_rep == 1, no_output=True)

[20]
assert(step_input.contains('test_stop_if_0.txt'))
assert(not step_input.contains('test_stop_if_1.txt'))

''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('test_stop_if_0.txt'))
        self.assertFalse(os.path.isfile('test_stop_if_1.txt'))
        #
        # stop_if should not be treated as error so the previously
        # generated output file will not be removed
        for rep in range(2):
            file = f'test_stop_if_{rep}.txt'
            if os.path.isfile(file):
                os.remove(file)
        script = SoS_Script(r'''
[10]
rep = range(2)
input: for_each='rep'
output: f'test_stop_if_{_rep}.txt'

_output.touch()
stop_if(_rep == 1)

[20]
assert(step_input.contains('test_stop_if_0.txt'))
assert(step_input.contains('test_stop_if_1.txt'))

''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('test_stop_if_0.txt'))
        self.assertTrue(os.path.isfile('test_stop_if_1.txt'))

    def testSkipIf(self):
        '''Test action stop_if'''
        script = SoS_Script(r'''
[0: shared='result']
rep = range(20)
result = []
input: for_each='rep'

skip_if(_rep > 10)
result.append(_rep)
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['result'],
                         [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

        # stop_if should not be treated as error so the previously
        # generated output file will be removed
        for rep in range(2):
            file = f'test_stop_if_{rep}.txt'
            if os.path.isfile(file):
                os.remove(file)
        script = SoS_Script(r'''
[10]
rep = range(2)
input: for_each='rep'
output: f'test_stop_if_{_rep}.txt'

_output.touch()
skip_if(_rep == 1)

[20]
assert(step_input.contains('test_stop_if_0.txt'))
assert(not step_input.contains('test_stop_if_1.txt'))

''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('test_stop_if_0.txt'))
        self.assertFalse(os.path.isfile('test_stop_if_1.txt'))
        #

    def testDoneIf(self):
        'Test action done_if'
        for rep in range(2):
            file = f'test_done_if_{rep}.txt'
            if os.path.isfile(file):
                os.remove(file)
        script = SoS_Script(r'''
[10]
rep = range(2)
input: for_each='rep'
output: f'test_done_if_{_rep}.txt'

_output.touch()
done_if(_rep == 1)

[20]
assert(step_input.contains('test_done_if_0.txt'))
assert(step_input.contains('test_done_if_1.txt'))

''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('test_done_if_0.txt'))
        self.assertTrue(os.path.isfile('test_done_if_1.txt'))

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
        self.assertRaises(Exception, Base_Executor(wf).run)

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
    @unittest.skipIf(
        True, 'Skip test because of no internet connection or in travis test')
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
download(['http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus-20170912.DB.gz'],
    dest_dir='tmp', decompress=True)
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('tmp/CancerGeneCensus-20170912.DB'))
        #
        # testing the download of single file
        #
        script = SoS_Script(r'''
[0]
download: dest_file='tmp/refgene.ppp'
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/resource/refgene.pkl
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('tmp/refgene.ppp'))
        # test option dest_dir
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp'
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/resource/refgene.pkl
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('tmp/refgene.pkl'))
        #

    @unittest.skipIf(
        not with_network or 'TRAVIS' in os.environ,
        'Skip test because of no internet connection or in travis test')
    def testDownloadMissingFile(self):
        # this will take a while
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True, max_jobs=2
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/resource/non-existing.gz
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus-20170912.DB.gz
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus.ann
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/DGV-hg38_20160831.ann
''')
        #start = time.time()
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)
        self.assertTrue(os.path.isfile('tmp/CancerGeneCensus.ann'))
        #self.assertGreater(time.time() - start, 3)
        # this will be fast
        #start = time.time()
        wf = script.workflow()
        self.assertRaises(Exception, Base_Executor(wf).run)
        #self.assertLess(time.time() - start, 3)
        #

    @unittest.skipIf(
        not with_network or 'TRAVIS' in os.environ,
        'Skip test because of no internet connection or in travis test')
    def testDownloadLargeFile(self):
        # test decompress tar.gz, .zip and .gz files
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus-20170912.DB.gz
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus.ann
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/DGV-hg38_20160831.ann
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # run in build mode
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus-20170912.DB.gz
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/CancerGeneCensus.ann
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/DGV-hg38_20160831.ann
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
        file_target('myreport.html').unlink()
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
        file_target('myreport.html').unlink()
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
            file_target(f).unlink()

    def testReport(self):
        '''Test action report'''
        script = SoS_Script(r'''
[A]
parameter: num=5
report: output='report.txt', expand=True
    touch {num}.txt

''')
        # output to a file
        if file_target('report.txt').exists():
            file_target('report.txt').unlink()
        wf = script.workflow()
        # run twice
        Base_Executor(wf, args=['--num', '7']).run()
        Base_Executor(wf, args=['--num', '5']).run()
        with open('report.txt') as report:
            self.assertEqual(report.read(), 'touch 5.txt\n\n')
        # test overwrite
        if file_target('report.txt').exists():
            file_target('report.txt').unlink()
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
        if file_target('report.txt').exists():
            file_target('report.txt').unlink()
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
            file_target(name).unlink()
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
            file_target(name).unlink()
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
        self.assertTrue(
            file_target(os.path.join('temp_wdr', 'a2.txt')).target_exists())
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
            if file_target(name).exists():
                file_target(name).unlink()
        wf = script.workflow()
        Base_Executor(wf).run()

    def testActiveActionOption(self):
        '''Test the active option of actions'''
        # disallow
        self.assertRaises(
            Exception, SoS_Script, '''
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
            ('slice(1,None)',
             ['temp/1.txt', 'temp/2.txt', 'temp/3.txt', 'temp/4.txt']),
            ('slice(1,-2)', ['temp/1.txt', 'temp/2.txt']),
            ('slice(None,None,2)', ['temp/0.txt', 'temp/2.txt', 'temp/4.txt']),
            ('True', [
                'temp/0.txt', 'temp/1.txt', 'temp/2.txt', 'temp/3.txt',
                'temp/4.txt'
            ]),
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
            Base_Executor(wf).run()
            files = list(glob.glob(os.path.join('temp', '*.txt')))
            self.assertEqual(
                sorted(files), sorted([x.replace('/', os.sep) for x in result]))
            #
            # test last iteration
            shutil.rmtree('temp')


if __name__ == '__main__':
    unittest.main()
