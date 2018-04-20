#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil
import subprocess
import sys
import unittest

from sos.eval import interpolate
from sos.parser import SoS_Script
from sos.targets import file_target, path, paths, sos_targets
from sos.utils import env
from sos.workflow_executor import Base_Executor


class TestTarget(unittest.TestCase):
    def setUp(self):
        env.reset()
        subprocess.call('sos remove -s', shell=True)
        # self.resetDir('~/.sos')
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

    def testTargetFormat(self):
        '''Test string interpolation of targets'''
        for target, fmt, res in [
                ('a.txt', '', 'a.txt'),
                (sos_targets('a.txt'), '', 'a.txt'),
                (sos_targets(['a.txt']), '', 'a.txt'),
                (sos_targets('/a/b/a.txt'), 'b', 'a.txt'),
                (sos_targets('a b.txt'), 'q', ("'a b.txt'", '"a b.txt"')),
                (sos_targets('a b.txt'), 'x', ".txt"),
        ]:
            if isinstance(res, str):
                self.assertEqual(
                    interpolate('{{target:{}}}'.format(fmt), globals(), locals()), res,
                    "Interpolation of {}:{} should be {}".format(target, fmt, res))
            else:
                self.assertTrue(
                    interpolate('{{target:{}}}'.format(fmt), globals(), locals()) in res,
                    "Interpolation of {}:{} should be one of {}".format(target, fmt, res))

    def testIterTargets(self):
        '''Test iterator interface of targets'''
        t = sos_targets('1', '2', ['3', '4'])
        self.assertEqual(len(t), 4)
        for idx, i in enumerate(t):
            self.assertEqual(str(i), str(idx + 1))

    def resetDir(self, dirname):
        if os.path.isdir(os.path.expanduser(dirname)):
            shutil.rmtree(os.path.expanduser(dirname))
        os.mkdir(os.path.expanduser(dirname))

    def testShared(self):
        '''Test option shared'''
        script = SoS_Script(r"""
parameter: res = 1

[0]
res = 2

[1]
res = 3
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 1)
        #
        script = SoS_Script(r"""
parameter: res = 1

[0: shared='res']
res = 2

[1]
res = 3
""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 2)
        #
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[0: shared='a']
res = 2

[1: shared='res']
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 3)
        self.assertEqual(env.sos_dict['a'], 30)
        # test multiple vars
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=('res', 'a')]
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 3)
        self.assertEqual(env.sos_dict['a'], 5)
        #
        # test expression
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared={'res': 'res + 6', 'c': 'a'}]
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 9)
        self.assertEqual(env.sos_dict['c'], 5)
        # test mixed vars and mapping
        script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=['res', {'c': 'a'}]]
res = 3
a = 5

""")
        wf = script.workflow()
        Base_Executor(wf).run()
        self.assertEqual(env.sos_dict['res'], 3)
        self.assertEqual(env.sos_dict['c'], 5)

#    def testSectionOptionWorkdir(self):
#        '''Test section option workdir'''
#        script = SoS_Script(r"""
#
# [1: workdir='tmp']
# run:
#    touch 'a.txt'
# """)
#        wf = script.workflow()
#        Base_Executor(wf).run()
#        self.assertTrue(os.path.isfile('tmp/a.txt'))
#        shutil.rmtree('tmp')

    def testDependsExecutable(self):
        '''Testing target executable.'''
        script = SoS_Script('''
[0]
depends: executable('python --version', '3.')
run:
    touch a.txt
''')
        wf = script.workflow()
        file_target('a.txt').remove('both')
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('a.txt'))
        file_target('a.txt').remove('both')

    @unittest.skipIf(sys.platform == 'win32', 'Windows executable cannot be created with chmod.')
    def testOutputExecutable(self):
        '''Testing target executable.'''
        # change $PATH so that lls can be found at the current
        # directory.
        os.environ['PATH'] += os.pathsep + '.'
        script = SoS_Script('''
[0]
output: executable('lls')
run:
    touch lls
    chmod +x lls
''')
        wf = script.workflow()
        file_target('lls').remove('both')
        env.config['sig_mode'] = 'force'
        Base_Executor(wf).run()
        # test validation
        env.config['sig_mode'] = 'default'
        Base_Executor(wf).run()
        file_target('lls').remove('both')

    def testDependsEnvVariable(self):
        '''Testing target env_variable.'''
        file_target('a.txt').remove('both')
        if sys.platform == 'win32':
            script = SoS_Script('''
[0]
depends: env_variable('AA')
output:  'a.txt'
run:
    echo %AA% > a.txt
''')
        else:
            script = SoS_Script('''
[0]
depends: env_variable('AA')
output:  'a.txt'
run:
    echo $AA > a.txt
''')
        wf = script.workflow()
        os.environ['AA'] = 'A1'
        Base_Executor(wf).run()
        with open('a.txt') as at:
            self.assertEqual(at.read().strip(), 'A1')
        # test validation
        Base_Executor(wf).run()
        # now if we change var, it should be rerun
        os.environ['AA'] = 'A2'
        Base_Executor(wf).run()
        # allow one second variation
        with open('a.txt') as at:
            self.assertEqual(at.read().strip(), 'A2')
        file_target('a.txt').remove('both')

    @unittest.skipIf(sys.platform == 'win32', 'Windows executable cannot be created with chmod.')
    def testProvidesExecutable(self):
        '''Testing provides executable target.'''
        # change $PATH so that lls can be found at the current
        # directory.
        os.environ['PATH'] += os.pathsep + '.'
        file_target('lls').remove('both')
        script = SoS_Script('''
[lls: provides=executable('lkls')]
run:
    touch lkls
    chmod +x lkls

[c]
depends: executable('lkls')

''')
        wf = script.workflow('c')
        Base_Executor(wf).run()
        file_target('lkls').remove('both')

    def testSharedVarInPairedWith(self):
        self.touch(['1.txt', '2.txt'])
        script = SoS_Script('''
[work_1: shared = {'data': 'step_output'}]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run: expand=True
  touch {_output}

[work_2]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}', paired_with = ['data']
output: expand_pattern('{_name}.out2')
run: expand=True
  touch {_data[0]} {_output}
''')
        wf = script.workflow()
        Base_Executor(wf).run()
        for file in ('1.out', '2.out', '1.out2', '2.out2'):
            file_target(file).remove('both')

    def testSharedVarInForEach(self):
        self.touch(['1.txt', '2.txt'])
        script = SoS_Script('''
[work_1: shared = {'data': 'step_output'}]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run: expand=True
  touch {_output}

[work_2]
input: "1.txt", "2.txt", group_by = 'single', for_each = 'data, data',  pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out2')
run: expand=True
  touch {_data} {_output}

''')
        wf = script.workflow()
        Base_Executor(wf).run()

    def testRemovedDepends(self):
        '''Test a case where a dependent file has signature, but
        gets removed before the next run.'''
        script = SoS_Script('''
[tet: provides='a.txt']
run:
    echo "something" > a.txt

[20]
depends: 'a.txt'
output: 'b.txt'
run:
    cat a.txt > b.txt
''')
        wf = script.workflow()
        # this should be ok.
        Base_Executor(wf).run()
        # now let us remove a.txt (but the signature is still there)
        os.remove('a.txt')
        os.remove('b.txt')
        Base_Executor(wf).run()

    def testSoSStep(self):
        '''Test target sos_step'''
        for file in ['t1.txt', 't2.txt', '5.txt', '10.txt', '20.txt']:
            file_target(file).remove('both')
        script = SoS_Script('''
[t1]
run:
    touch t1.txt

[t2: provides='t2.txt']
depends: sos_step('t1')
run:
    touch t2.txt

[5]
run:
    touch 5.txt

[10]
depends: sos_step('t2')
run:
    touch 10.txt

[20]
depends: sos_step('t1')
task:
run:
    touch 20.txt
''')
        wf = script.workflow()
        env.config['sig_mode'] = 'force'
        # this should be ok.
        Base_Executor(wf).run()
        for file in ['t1.txt', 't2.txt', '5.txt', '10.txt', '20.txt']:
            self.assertTrue(file_target(file).target_exists(), file + ' should exist')
            file_target(file).remove('both')

    def testZap(self):
        '''Test zap'''
        with open('testzap.txt', 'w') as sf:
            sf.write('some text')
        path('testzap.txt').zap()
        self.assertTrue(os.path.isfile('testzap.txt.zapped'))
        self.assertFalse(os.path.isfile('testzap.txt'))
        # re-zap is ok
        file_target('testzap.txt').zap()
        self.assertTrue(os.path.isfile('testzap.txt.zapped'))
        self.assertFalse(os.path.isfile('testzap.txt'))
        # non-existent file
        os.remove('testzap.txt.zapped')
        self.assertRaises(FileNotFoundError, path('testzap.txt').zap)
        #
        with open('testzap.txt', 'w') as sf:
            sf.write('some text')
        with open('testzap1.txt', 'w') as sf:
            sf.write('some text')
        paths('testzap.txt', 'testzap1.txt').zap()
        self.assertTrue(os.path.isfile('testzap.txt.zapped'))
        self.assertFalse(os.path.isfile('testzap.txt'))
        self.assertTrue(os.path.isfile('testzap1.txt.zapped'))
        self.assertFalse(os.path.isfile('testzap1.txt'))
        #
        os.remove('testzap.txt.zapped')
        os.remove('testzap1.txt.zapped')
        with open('testzap.txt', 'w') as sf:
            sf.write('some text')
        with open('testzap1.txt', 'w') as sf:
            sf.write('some text')
        sos_targets(['testzap.txt', 'testzap1.txt']).zap()
        self.assertTrue(os.path.isfile('testzap.txt.zapped'))
        self.assertFalse(os.path.isfile('testzap.txt'))
        self.assertTrue(os.path.isfile('testzap1.txt.zapped'))
        self.assertFalse(os.path.isfile('testzap1.txt'))

    def testZapRun(self):
        '''Test run with zapped input files'''
        with open('zap1.txt', 'w') as sf:
            sf.write('seomething')
        script = SoS_Script('''\
[1]
input: 'zap1.txt'
output: "zap2.txt"
run:
  echo asd>zap2.txt
_input.zap()
''')
        wf = script.workflow()
        env.config['sig_mode'] = 'force'
        Base_Executor(wf).run()
        self.assertTrue(os.path.isfile('zap1.txt.zapped'))
        self.assertFalse(os.path.isfile('zap1.txt'))
        self.assertTrue(os.path.isfile('zap2.txt'))
        # can run again
        env.config['sig_mode'] = 'default'
        Base_Executor(wf).run()
        # now if we remove target
        os.remove('zap2.txt')
        self.assertRaises(Exception, Base_Executor(wf).run)


if __name__ == '__main__':
    unittest.main()
