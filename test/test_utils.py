#!/usr/bin/env python3
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

# these functions are normally not available but can be imported 
# using their names for testing purposes
from pysos.utils import env, logger, WorkflowDict, ProgressBar, text_repr
from pysos.pattern import extract_pattern, expand_pattern
from pysos.sos_eval import interpolate, SoS_eval, InterpolationError
from pysos.actions import downloadURL
from pysos.sos_script import SoS_Script
from pysos.sos_executor import Sequential_Executor

class TestUtils(unittest.TestCase):
    def setUp(self):
        env.reset()

    def testLogger(self):
        '''Test logging level'''
        for verbosity in ['0', '1', '2', '3', '4']:
            env.verbosity = verbosity
            logger.trace('Verbosity {}:trace message with ``empahsized text`` in between'.format(env.verbosity))
            logger.debug('Verbosity {}:debug message with ``empahsized text`` in between'.format(env.verbosity))
            logger.info('Verbosity {}:info message with ``empahsized text`` in between'.format(env.verbosity))
            logger.warning('Verbosity {}:warning message with ``empahsized text`` in between'.format(env.verbosity))
            logger.error('Verbosity {}:error message with ``empahsized text`` in between'.format(env.verbosity))
        # log
        if os.path.isfile('test.log'):
            os.remove('test.log')
        env.logfile = 'test.log'
        for verbosity in ['0', '1', '2', '3', '4']:
            env.verbosity = verbosity
            logger.trace('Verbosity {}:trace message with ``empahsized text`` in between'.format(env.verbosity))
            logger.debug('Verbosity {}:debug message with ``empahsized text`` in between'.format(env.verbosity))
            logger.info('Verbosity {}:info message with ``empahsized text`` in between'.format(env.verbosity))
            logger.warning('Verbosity {}:warning message with ``empahsized text`` in between'.format(env.verbosity))
            logger.error('Verbosity {}:error message with ``empahsized text`` in between'.format(env.verbosity))
        # log file should not have any color codes
        with open('test.log') as logfile:
            line_count = 0
            for line in logfile:
                line_count += 1
                self.assertFalse('\033[' in line)
            # 4 lines for all logging level (logging level of logfile is fixed to DEBUG)
            self.assertEqual(line_count, 20)
        os.remove('test.log')

    def testInterpolation(self):
        '''Test string interpolation'''
        env.sos_dict = WorkflowDict({
            'os': os,
            'a': 100,
            'b': 20,
            'c': ['file1', 'file2', 'file3'],
            'd': {'a': 'file1', 'b':'file2'},
            'e': set([1, 'a']),
            'var1': 1/2.,
            'var2': [1, 2, 3.1],
            'file': 'a/b.txt',
            'files2': ['a/b.txt', 'c.d.txt'],
            'file_ws': ['d i r/f .txt']
            })
        for sigil in ('${ }', '{ }', '[ ]', '%( )', '[[ ]]', '%( )s'):
            l, r = sigil.split(' ')
            for expr, result in [
                ('{0}1{1}', '1'),
                ('{0}a{1}', '100'),
                ('{0}a+b{1}', '120'),
                ('{0}a+b*5{1}', '200'),
                ('{0}a+b*5{1} is 200', '200 is 200'),
                ('{0}a+b*5{1} and {0}a{1}', '200 and 100'),
                ('Pre {0}a+b*5{1} and {0}a{1} after', 'Pre 200 and 100 after'),
                ('Nested {0}a+b*{0}b//2{1}{1}', 'Nested 300'),
                ('Format {0}a:.5f{1}', 'Format 100.00000'),
                ('{0}var2[:2]{1}', '1 2'),
                ('{0}var2[1:]{1}', '2 3.1'),
                # nested
                ('Nested {0}a:.{0}4+1{1}f{1}', 'Nested 100.00000'),
                # deep nested
                ('Triple Nested {0}a:.{0}4+{0}5//5{1}{1}f{1}', 'Triple Nested 100.00000'),
                # nested invalid
                ('Nested invalid {0}"{0}a-{1}"{1}', 'Nested invalid {}a-{}'.format(l, r)),
                ('Nested valid {0}"{0}a{1}-"{1}', 'Nested valid 100-'),
                #
                ('Dict {0}d{1}', ['Dict a b', 'Dict b a']),
                ('set {0}e{1}', ['set 1 a', 'set a 1']),
                ('Fmt {0}var1:.2f{1}', 'Fmt 0.50'),
                ('Fmt {0}var2:.2f{1}', 'Fmt 1.00 2.00 3.10'),
                ('LC {0}[x*2 for x in var2]{1}', 'LC 2 4 6.2'),
                ('LC {0}[x*2 for x in var2]:.2f{1}', 'LC 2.00 4.00 6.20'),
                #
                # [['a':'b', 'c':'d']['a']] works because
                #     ['a':'b', 'c':'d']a
                # is invalid so SoS does not consider ['a'] as nested expression
                #
                ('Newline {0}{{"a": "b", \n"c": "d"}}["a"]{1}', 'Newline b'),
                #
                # string literal within interpolated expression
                ('Literal {0}"{0} {1}"{1}', 'Literal ' + sigil),
                #
                ("{0}' {{}} '.format(a){1}", ' 100 '),
                #
                ("{0}os.path.basename(file){1}", 'b.txt'),
                ('{0}os.path.basename(file_ws[0]){1}', 'f .txt'),
                #
                # ! conversion
                ('{0}file!r{1}', "'a/b.txt'"),
                ('{0}file!s{1}', "a/b.txt"),
                ('''{0}"a'b"!r{1}''', '"a\'b"'),
                ('''{0}'a"b'!r{1}''', "'a\"b'"),
                #
                # !q conversion (added by SoS)
                ('{0}file_ws[0]!q{1}', "'d i r/f .txt'"),
                #
                # !, joined by ,
                ('{0}var2!r,{1}', "1,2,3.1"),
                ('{0}c!r,{1}', "'file1','file2','file3'"),
                ('{0}c!,{1}', "file1,file2,file3"),
                #
                # full name by 'a'
                ('{0}"test_utils.py"!a{1}', os.path.abspath('test_utils.py')),
                ('{0}"~/test_utils.py"!a{1}', os.path.expanduser('~/test_utils.py')),
                ('{0}"~/test_utils.py"!e{1}', os.path.expanduser('~/test_utils.py')),
                ('{0}"test/test_utils.py"!b{1}', "test_utils.py"),
            ]:
                #print('Interpolating "{}" with sigal "{}"'.format(expr.format(l, r).replace('\n', r'\n'), sigil))
                if isinstance(result, str):
                    self.assertEqual(interpolate(expr.format(l, r), sigil=sigil), result)
                else:
                    # for cases when the order of output is not guaranteed
                    self.assertTrue(interpolate(expr.format(l, r), sigil=sigil) in result)
        #
        # locals should be the one passed to the expression
        self.assertTrue('file_ws' in interpolate('${globals().keys()}'))
        # 5:.5.0f does not work.
        self.assertRaises(InterpolationError, interpolate, '${5:.${4/2.}f}')

    def testEval(self):
        '''Test the evaluation of SoS expression'''
        env.sos_dict = WorkflowDict({
            'interpolate': interpolate,
            'a': 100,
            'b': 'file name',
            'c': ['file1', 'file2', 'file 3'],
            'd': {'a': 'file1', 'b':'file2'},
        })
        for expression, result in [
            ('''"This is ${a+100}"''', 'This is 200'),
            ('''"${a+100}" + "${a//100}"''', "2001"),
            ('''"${c[1]}"''', 'file2'),
            ('''"${c[1:]}"''', 'file2 file 3'),
            ('''"${d}"''', ['a b', 'b a']),
            ('''"${d}"*2''', ['a ba b', 'b ab a']),
            ('''"${d['a']}"''', 'file1'),
            ('''"${b!q}"''', "'file name'"),
            ('''"${b!r}"''', "'file name'"),
            ('''"${c!q}"''', "file1 file2 'file 3'"),
            ]:
            if isinstance(result, str):
                self.assertEqual(SoS_eval(expression), result)
            else:
                self.assertTrue(SoS_eval(expression) in result)
        #
        # interpolation will only happen in string
        self.assertRaises(SyntaxError, SoS_eval, '''${a}''')

    def testWorkflowDict(self):
        '''Test workflow dict with attribute access'''
        env.reset()
        d = WorkflowDict()
        d['a'] = 1
        self.assertEqual(d['a'], 1)
        #
        d['a'] += 1
        self.assertEqual(d['a'], 2)

    def testPatternMatch(self):
        '''Test snake match's pattern match facility'''
        res = extract_pattern('{a}-{b}.txt', ['file-1.txt', 'file-ab.txt'])
        self.assertEqual(res['a'], ['file', 'file'])
        self.assertEqual(res['b'], ['1', 'ab'])
        res = extract_pattern('{a}-{b}.txt', ['file--ab--cd.txt'])
        self.assertEqual(res['a'], ['file--ab-'])
        self.assertEqual(res['b'], ['cd'])
        res = extract_pattern('{path}/{to}/{file}.txt', ['/tmp//1.txt'])
        self.assertEqual(res['path'], [None])
        self.assertEqual(res['to'], [None])
        self.assertEqual(res['file'], [None])
        res = extract_pattern('{path}/{to}/{file}.txt', ['/tmp/test/1.txt.txt'])
        self.assertEqual(res['path'], ['/tmp'])
        self.assertEqual(res['to'], ['test'])
        self.assertEqual(res['file'], ['1.txt'])
        # expand_pattern
        env.sos_dict = WorkflowDict({
            'os': os,
            'a': 100,
            'b': 'file name',
            'c': ['file1', 'file2', 'file 3'],
            'd': {'a': 'file1', 'b':'file2'},
        })
        self.assertEqual(expand_pattern('{b}.txt'), ['file name.txt'])
        self.assertEqual(expand_pattern('{c}.txt'), ['file1.txt', 'file2.txt', 'file 3.txt'])
        self.assertEqual(expand_pattern('{a}_{c}.txt'), ['100_file1.txt', '100_file2.txt', '100_file 3.txt'])

    def testProgressBar(self):
        '''Test progress bar'''
        env.verbosity = 1
        prog = ProgressBar('test', 100)
        for i in range(100):
            prog.update(i)
        prog.done()
        prog = ProgressBar('test', 100)
        for i in range(20):
            prog.progress(5)
        prog.done()
        #
        script = SoS_Script('''

[1]
[2]
[3]
[4]
[5]
''')
        wf = script.workflow()
        Sequential_Executor(wf).run()
        # progress bar with nested workflow
        script = SoS_Script('''
import time
time.sleep(0.5)
[sub_1]
[sub_2]
[sub_3]
[sub_4]
[a_1]
[a_2]
[a_3]
sos_run('sub')
[a_4]
[a_5]
''')
        wf = script.workflow('a')
        Sequential_Executor(wf).run()


    def testDownload(self):
        '''Test download file'''
        downloadURL('https://github.com/bpeng2000/SOS/wiki/SoS_March2016.pdf', 'tmp/SoS_March2016.pdf', index=0)
        self.assertTrue(os.path.isfile('tmp/SoS_March2016.pdf'))
        os.remove('tmp/SoS_March2016.pdf')

    def testTextRepr(self):
        '''Test text_repr'''
        for text in ['asdf g', 'a \\ng', r'a\nb']:
            self.assertEqual(text_repr(text), repr(text))
        self.assertEqual(text_repr(r'''a
\nb'''), "'a\\n\\\\nb'")
        self.assertEqual(text_repr(r"""a
\nb'"""), '"a\\n\\\\nb\'"')
        self.assertEqual(text_repr(r"""a
\nb''"""), '"a\\n\\\\nb\'\'"')
        self.assertEqual(text_repr(r"""a
\nb'''"""), '"a\\n\\\\nb\'\'\'"')
        self.assertEqual(text_repr(r"""a
'''\nb'''"""), '"a\\n\'\'\'\\\\nb\'\'\'"')
        self.assertEqual(text_repr(r'''a
b
\nc'''), "r'''a\nb\n\\nc'''")

if __name__ == '__main__':
    unittest.main()
