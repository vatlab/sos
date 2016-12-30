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

from sos.sos_script import SoS_Script
from sos.utils import env
from sos.sos_eval import  Undetermined
from sos.sos_executor import Base_Executor, ExecuteError
from sos.target import FileTarget

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
            except:
                if n > 1:
                    raise
    return wrapper

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
from sos.actions import SoS_Action

@SoS_Action(run_mode='run')
def func_run():
    return 1

@SoS_Action(run_mode=['run', 'dryrun'])
def func_both():
    return 1

[0: shared=('b', 'c')]
b=func_run()
c=func_both()
""")
        wf = script.workflow()
        Base_Executor(wf).dryrun()
        self.assertTrue(isinstance(env.sos_dict['b'], Undetermined))
        self.assertTrue(isinstance(env.sos_dict['c'], Undetermined))
        #
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['b'], 1)
        self.assertEqual(env.sos_dict['c'], 1)

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
        Base_Executor(wf).dryrun()
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
sh: args='-n ${filename!q}'
    touch a.txt
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertFalse(os.path.exists('a.txt'))


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


    def testPython3(self):
        script = SoS_Script(r'''
[0]
python3:
a = {'1', '2'}
print(a)
''')
        wf = script.workflow()
        Base_Executor(wf).run()


    def testPython2(self):
        script = SoS_Script(r'''
[0]
python2:
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


    @multi_attempts
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
        Base_Executor(wf).run()
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
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('tmp/test.ann'))
        # test option dest_dir
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp'
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/annoDB/hapmap_ASW_freq.ann
''')
        wf = script.workflow()
        Base_Executor(wf).run()
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
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        self.assertTrue(os.path.isfile('tmp/hapmap_ASW_freq-hg18_20100817.DB'))
        self.assertGreater(time.time() - start, 3)
        # this will be fast
        start = time.time()
        wf = script.workflow()
        self.assertRaises(ExecuteError, Base_Executor(wf).run)
        self.assertLess(time.time() - start, 3)
        # 
        # test decompress tar.gz file
        script = SoS_Script(r'''
[0]
download: dest_dir='tmp', decompress=True
    http://bioinformatics.mdanderson.org/Software/VariantTools/repository/programs/SKAT_0.82.tar.gz
''')
        wf = script.workflow()
        Base_Executor(wf).run()
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
pandoc(output=_output[0])
''')
        wf = script.workflow()
        Base_Executor(wf, config={'report_output': 'report.md'}).run()
        self.assertTrue(os.path.isfile('myreport.html'))
        #
        FileTarget('myreport.html').remove('both')
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
        FileTarget('myreport.html').remove()
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
        Base_Executor(wf, config={'report_output': '${step_name}.md'}).run()
        for f in ['default_10.md', 'default_20.md', 'output.html']:
            self.assertTrue(FileTarget(f).exists())
            FileTarget(f).remove()


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

    def testReport(self):
        '''Test action report'''
        script = SoS_Script(r'''
[A]
parameter: num=5
report:
    touch ${num}.txt

''')
        # output to a file
        FileTarget('report.txt').remove('both')
        wf = script.workflow()
        # run twice
        Base_Executor(wf, args=['--num', '7'], config={'report_output': 'report.txt'}).run()
        Base_Executor(wf, args=['--num', '5'], config={'report_output': 'report.txt'}).run()
        with open('report.txt') as report:
            self.assertEqual(report.read(), 'touch 7.txt\n\ntouch 5.txt\n\n')
        # test overwrite
        FileTarget('report.txt').remove('both')
        script = SoS_Script(r'''
[A]
report: output='report.txt'
    ${step_name}

[A_10]
report: output='report.txt'
    ${step_name}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        # output to a file
        # run twice
        Base_Executor(wf).run()
        with open('report.txt') as report:
            self.assertEqual(report.read(), 'A_10\n\n')

        # test filenames interpolation
        #
        FileTarget('report_1.txt').remove('both')
        FileTarget('report_2.txt').remove('both')
        script = SoS_Script(r'''
[A_1]
report:
    report_${step_name.split('_')[-1]}.txt

[A_2]
report:
    report_${step_name.split('_')[-1]}.txt

''')
        wf = script.workflow()
        Base_Executor(wf,  config={'report_output': 'report_${step_name.split("_")[-1]}.txt'}).run()
        # output to a file
        # run twice
        Base_Executor(wf).run()
        for name in ('report_1.txt', 'report_2.txt'):
            with open(name) as report:
                self.assertEqual(report.read(), name + '\n\n')
            FileTarget(name).remove()
        #
        # test input from another file
        FileTarget('report.txt').remove()
        script = SoS_Script(r'''
[A_1]
sh: output='a.txt'
    echo 'something' > a.txt
report(input='a.txt', output='out.txt')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for name in ('a.txt', 'out.txt'):
            with open(name) as report:
                self.assertEqual(report.read().strip(), 'something')
            FileTarget(name).remove()
        #
        script = SoS_Script(r'''
[A_1]
sh: output='a.txt'
    echo 'something' > a.txt

[A_2]
sh: output='b.txt'
    echo 'something else' > b.txt

[A_3]
report(input=['a.txt', 'b.txt'], output='out.txt')
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for name in ('a.txt', 'b.txt', 'out.txt'):
            self.assertTrue(FileTarget(name).exists())
            FileTarget(name).remove()

    def testRegenerateReport(self):
        '''Testing the regeneration of report once is needed. The problem
        here is the 'input' parameter of report.'''
        script = SoS_Script(r'''
[A_1]
output: 'a1.txt'
run:
    echo 'a1' >> a1.txt

report: output='a1.md'
    a1

[A_2]
output: 'a2.txt'
run:
    echo 'a2' >> a2.txt
report: output='a2.md'
    a2

[A_3]
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
            FileTarget(name).remove('target')
        wf = script.workflow()
        Base_Executor(wf).run()


if __name__ == '__main__':
    unittest.main()
