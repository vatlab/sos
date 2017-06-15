#!/usr/bin/env python3
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
import unittest
import cProfile
import timeit
import textwrap

# these functions are normally not available but can be imported 
# using their names for testing purposes
from sos.utils import env, logger, WorkflowDict, stable_repr
from sos.pattern import extract_pattern, expand_pattern
from sos.sos_eval import interpolate, SoS_eval, InterpolationError, accessed_vars, \
    Undetermined, on_demand_options
from sos.sos_script import SoS_Script
from sos.sos_executor import Base_Executor, analyze_section
from sos.target import executable, remote

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

class TestUtils(unittest.TestCase):
    def setUp(self):
        env.reset()

    def testLogger(self):
        '''Test logging level'''
        for verbosity in [0, 1, 2, 3, 4]:
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
        import logging
        logging.shutdown()
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
            'file_ws': ['d i r/f .txt'],
            'remote': remote
            })
        for sigil in ('${ }', '{ }', '[ ]', '%( )', '[[ ]]', '%( )s', '# #', '` `'):
            l, r = sigil.split(' ')
            for expr, result, nested, exclude in [
                ('{0}1{1}', '1', False, []),
                ('{0}a{1}', '100', False, []),
                ('{0}a+b{1}', '120', False, []),
                ('{0}a+b*5{1}', '200', False, []),
                ('{0}a+b*5{1} is 200', '200 is 200', False, []),
                ('{0}a+b*5{1} and {0}a{1}', '200 and 100', False, []),
                ('Pre {0}a+b*5{1} and {0}a{1} after', 'Pre 200 and 100 after', False, []),
                ('Nested {0}a+b*{0}b//2{1}{1}', 'Nested 300', True, []),
                ('Format {0}a:.5f{1}', 'Format 100.00000', False, []),
                ('{0}var2[:2]{1}', '1 2', False, []),
                ('{0}var2[1:]{1}', '2 3.1', False, []),
                # nested
                ('Nested {0}a:.{0}4+1{1}f{1}', 'Nested 100.00000', True, []),
                # deep nested
                ('Triple Nested {0}a:.{0}4+{0}5//5{1}{1}f{1}', 'Triple Nested 100.00000', True, []),
                # nested invalid
                ('Nested invalid {0}"{0}a-{1}"{1}', 'Nested invalid {}a-{}'.format(l, r), True, []),
                ('Nested valid {0}"{0}a{1}-"{1}', 'Nested valid 100-', True, []),
                #
                ('Dict {0}d{1}', ['Dict a b', 'Dict b a'], False, []),
                ('set {0}e{1}', ['set 1 a', 'set a 1'], False, []),
                ('Fmt {0}var1:.2f{1}', 'Fmt 0.50', False, []),
                ('Fmt {0}var2:.2f{1}', 'Fmt 1.00 2.00 3.10', False, []),
                ('LC {0}[x*2 for x in var2]{1}', 'LC 2 4 6.2', False, []),
                ('LC {0}[x*2 for x in var2]:.2f{1}', 'LC 2.00 4.00 6.20', False, []),
                #
                # [['a':'b', 'c':'d']['a']] works because
                #     ['a':'b', 'c':'d']a
                # is invalid so SoS does not consider ['a'] as nested expression
                #
                ('Newline {0}{{"a": "b", \n"c": "d"}}["a"]{1}', 'Newline b', False, []),
                #
                # string literal within interpolated expression
                ('Literal {0}"{0} {1}"{1}', 'Literal ' + sigil, True, []),
                # this case does not work because {} would become '' with sigil {}
                ("{0}' {{}} '.format(a){1}", ' 100 ', False, ['{ }']),
                #
                ("{0}os.path.basename(file){1}", 'b.txt', False, []),
                ('{0}os.path.basename(file_ws[0]){1}', 'f .txt', False, []),
                #
                # ! conversion
                ('{0}file!r{1}', "'a/b.txt'", False, []),
                ('{0}file!s{1}', "a/b.txt", False, []),
                ('''{0}"a'b"!r{1}''', '"a\'b"', False, []),
                ('''{0}'a"b'!r{1}''', "'a\"b'", False, []),
                #
                # !q conversion (added by SoS)
                ('{0}file_ws[0]!q{1}', '"d i r/f .txt"' if sys.platform == 'win32' else "'d i r/f .txt'", False, []),
                ('{0}file_ws[0]!e{1}', "d\\ i\\ r/f\\ .txt", False, []),
                #
                # !, joined by ,
                ('{0}var2!r,{1}', "1, 2, 3.1", False, []),
                ('{0}c!r,{1}', "'file1', 'file2', 'file3'", False, []),
                ('{0}c!,{1}', "file1, file2, file3", False, []),
                #
                ('{0}10000:,{1}', "10,000", False, []),
                #
                # full name by 'a'
                ('{0}"test_utils.py"!a{1}', os.path.abspath('test_utils.py'), False, []),
                ('{0}"a/b/c/test_utils.py"!b{1}', 'test_utils.py', False, []),
                ('{0}"a/b/c/test_utils.py"!x{1}', '.py', False, []),
                ('{0}"a/b/c/test_utils.py"!d{1}', 'a/b/c', False, []),
                ('{0}"a/b/c/test_utils.py"!dd{1}', 'a/b', False, []),
                ('{0}"a/b/c/test_utils.py"!ddb{1}', 'b', False, []),
                ('{0}"a/b/c/test_utils.py"!n{1}', 'a/b/c/test_utils', False, []),
                ('{0}"a/b/c/test_utils.py"!bn{1}', 'test_utils', False, []),
                ('{0}remote("a/b/c/test_utils.py")!Rbn{1}', 'test_utils', False, []),
                ('{0}"~/test_utils.py"!a{1}', os.path.abspath(os.path.expanduser('~/test_utils.py')), False, []),
                ('{0}"~/test_utils.py"!u{1}', os.path.expanduser('~/test_utils.py'), False, []),
                ('{0}"test/test_utils.py"!b{1}', "test_utils.py", False, []),
                (r'''{0}{1}''', '', False, []),
            ]:
                if l == r and nested:
                    continue
                if sigil in exclude:
                    continue
                if isinstance(result, str):
                    self.assertEqual(interpolate(expr.format(l, r), sigil=sigil), result,
                            'Failed to interpolate {} with sigil {} and expected result {}'.format(expr, sigil, result))
                else:
                    # for cases when the order of output is not guaranteed
                    self.assertTrue(interpolate(expr.format(l, r), sigil=sigil) in result,
                            'Failed to interpolate {} with sigil {} and expected result {}'.format(expr, sigil, result))
        #
        # locals should be the one passed to the expression
        self.assertTrue('file_ws' in interpolate('${globals().keys()}', '${ }'))
        # 5:.5.0f does not work.
        self.assertRaises(InterpolationError, interpolate, '${5:.${4/2.}f}', '${ }')

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
            ('''"${b!q}"''', '"file name"' if sys.platform == "win32" else "'file name'"),
            ('''"${b!r}"''', "'file name'"),
            ('''"${c!q}"''', 'file1 file2 "file 3"' if sys.platform == 'win32' else "file1 file2 'file 3'"),
            ]:
            if isinstance(result, str):
                self.assertEqual(SoS_eval(expression, '${ }'), result,
                        "Test {} with expected result {}".format(expression, result))
            else:
                self.assertTrue(SoS_eval(expression, '${ }') in result,
                        "Test {} with expected result {}".format(expression, result))
        #
        # interpolation will only happen in string
        self.assertRaises(SyntaxError, SoS_eval, '''${a}''', '${ }')

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

    def testAccessedVars(self):
        '''Test accessed vars of a SoS expression or statement.'''
        self.assertEqual(accessed_vars('''a = 1''', '${ }'), {'a'})
        self.assertEqual(accessed_vars('''a = b + 2.0''', '${ }'), {'a', 'b'})
        self.assertEqual(accessed_vars('''a = "C"''', '${ }'), {'a'})
        self.assertEqual(accessed_vars('''a = "C" + "${D}"''', '${ }'), {'a', 'D'})
        self.assertEqual(accessed_vars('''a = 1 + "${D + 20:f}" ''', '${ }'), {'a', 'D'})
        self.assertEqual(accessed_vars('''k, "a.txt", "b.txt", skip=True ''', '${ }'), {'k', 'skip', 'True'})
        # this is a complicated case because the actual variable depends on the
        # result of an expression... However, in the NO-evaluation case, this is
        # the best we can do.
        self.assertEqual(accessed_vars('''c + "${D + '${E}'}" ''', '${ }'), {'c', 'D', 'E'})

    def testProgressBar(self):
        '''Test progress bar'''
        env.verbosity = 1
        #
        script = SoS_Script('''

[1]
[2]
[3]
[4]
[5]
''')
        wf = script.workflow()
        Base_Executor(wf).run()


    def testTextRepr(self):
        # the " as the last character can lead to problems...
        script = SoS_Script('''
run:
    echo "Hi, This is from bash"''')
        wf = script.workflow()
        Base_Executor(wf).run()
        #
        for text in ('"""a"""', '"b"',
            r'"""\na\\nb"""', r"'''a\nb'''",
            """ "a'\\"='" """): 
            script = SoS_Script(r'''
a = 1
python:
   with open('tmp.txt', 'w') as tmp:
      tmp.write({} + '{}')
k = """b"""'''.format(text, '${a}')
)
            wf = script.workflow()
            Base_Executor(wf).run()
            with open('tmp.txt') as tmp:
                self.assertEqual(tmp.read(), eval(text) + '1')
        os.remove('tmp.txt')


    def performanceTestInterpolation(self):
        '''Testing the performance of string interpolation. This test
        will not be executed automatically.'''
        setup_stmt = textwrap.dedent('''
        from sos.utils import env, WorkflowDict
        from sos.sos_eval import interpolate
        env.sos_dict = WorkflowDict({
            'interpolate': interpolate,
            'a': 100,
            'b': 'file name',
            'c': ['file1', 'file2', 'file 3'],
        })
        ''')
        stmt = "interpolate('${a}_${b}_${c}.txt', sigil='${ }')"
        ni_stmt = '''eval("'{}_{}_{}.txt'.format(a, b, c)", env.sos_dict._dict)'''
        exec(setup_stmt)
        cProfile.run(stmt)
        # 
        # original implementation takes ~12s on mac mini
        # No interpolation (python .format()) takes ~2.28 s
        #
        # 9/1/2016: using directly substitution reduced time to interpolate
        #   simple replacement (${var}) to about 6.06s.
        #
        print('Interpolating {} times take {} seconds'.format(
            100000,
            timeit.timeit(stmt, setup=setup_stmt, number=100000)))
        # comparing to non-interpolation
        print('No interpolating {} times take {} seconds'.format(
            100000,
            timeit.timeit(ni_stmt, setup=setup_stmt, number=100000)))

    def testAnalyzeSection(self):
        '''Test analysis of sections (statically)'''
        script = SoS_Script('''
g1 = 'a'
g2 = 1
parameter: p1 = 5
parameter: infiles = 'a.txt'

[A_1: shared='b']
b = p1 + 2
input:  infiles
output: None
c = 5

[A_2]
b = [1, 2, 3]
input: for_each='b'
depends: 'some.txt', executable('ls')
import time
import random

r = random.randint(1, 5)
time.sleep(r)

[A_3]
input: None
print(p1)

''')
        wf = script.workflow('A')
        for section in wf.sections:
            res = analyze_section(section)
            if section.names[0][1] == '1':
                self.assertTrue(isinstance(res['step_input'], Undetermined))
                self.assertEqual(res['step_depends'], [])
                self.assertEqual(res['step_output'], [])
                self.assertEqual(res['environ_vars'], {'p1', 'infiles'})
                self.assertEqual(res['signature_vars'], {'c'})
                self.assertEqual(res['changed_vars'], {'b'})
            elif section.names[0][1] == '2':
                self.assertEqual(res['step_input'], [])
                self.assertEqual(res['step_depends'], ['some.txt', executable('ls')])
                self.assertTrue(isinstance(res['step_output'], Undetermined))
                # for_each will not be used for DAG
                self.assertEqual(res['environ_vars'], {'for_each'})
                self.assertEqual(res['signature_vars'], {'import', 'r', 'time', 'random'})
                self.assertEqual(res['changed_vars'], set())

    def testOnDemandOptions(self):
        '''Test options that are evaluated upon request.'''
        options = on_demand_options(
            {'a': '"est"', 'b': 'c', 'c': 'e + 2'}, '${ }'
        )
        env.sos_dict = WorkflowDict({
            'e': 10,
            })
        self.assertEqual(options['a'], 'est')
        self.assertRaises(KeyError, options.__getitem__, 'd')
        self.assertRaises(ValueError, options.__getitem__, 'b')
        self.assertEqual(options['c'], 12)
        #
        env.sos_dict.set('e', 20)
        self.assertEqual(options['c'], 22)
        env.sos_dict.set('c', 200)
        self.assertEqual(options['b'], 200)

    def testStableRepr(self):
        self.assertEqual(stable_repr({1, 2, '3', '1'}), "{'1', '3', 1, 2}")
        self.assertEqual(stable_repr({1 : 2, 3:4}), "{1:2, 3:4}")
        self.assertEqual(stable_repr([1, 3, 4]), "[1, 3, 4]")

if __name__ == '__main__':
    unittest.main()
