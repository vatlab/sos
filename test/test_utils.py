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

from sostestcase import SoS_TestCase
import shlex

from pysos.utils import *

class TestUtils(SoS_TestCase):
    def testLogger(self):
        'Test logging level'
        for verbosity in ['0', '1', '2', '3']:
            env.verbosity = verbosity
            env.logger.trace('Verbosity {}:trace message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.debug('Verbosity {}:debug message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.info('Verbosity {}:info message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.warning('Verbosity {}:warning message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.error('Verbosity {}:error message with ``empahsized text`` in between'.format(env.verbosity))
        # log
        if os.path.isfile('test.log'):
            os.remove('test.log')
        env.logfile = 'test.log'
        for verbosity in ['0', '1', '2', '3']:
            env.verbosity = verbosity
            env.logger.trace('Verbosity {}:trace message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.debug('Verbosity {}:debug message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.info('Verbosity {}:info message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.warning('Verbosity {}:warning message with ``empahsized text`` in between'.format(env.verbosity))
            env.logger.error('Verbosity {}:error message with ``empahsized text`` in between'.format(env.verbosity))
        # log file should not have any color codes
        with open('test.log') as logfile:
            line_count = 0
            for line in logfile:
                line_count += 1
                self.assertFalse('\033[' in line)
            # 4 lines for all logging level (logging level of logfile is fixed to DEBUG)
            self.assertEqual(line_count, 16)

    def testInterpolation(self):
        lvar = {
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
            }
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
                ('Nested {0}a+b*{0}b/2{1}{1}', 'Nested 300'),
                ('Format {0}a:.5f{1}', 'Format 100.00000'),
                # nested 
                ('Nested {0}a:.{0}4+1{1}f{1}', 'Nested 100.00000'),
                # deep nested
                ('Triple Nested {0}a:.{0}4+{0}5/5{1}{1}f{1}', 'Triple Nested 100.00000'),
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
            ]:
                #print('Interpolating "{}" with sigal "{}"'.format(expr.format(l, r).replace('\n', r'\n'), sigil))
                if isinstance(result, str):
                    self.assertEqual(interpolate(expr.format(l, r), globals(), lvar, sigil=sigil), result)
                else:
                    # for cases when the order of output is not guaranteed
                    self.assertTrue(interpolate(expr.format(l, r), globals(), lvar, sigil=sigil) in result)
        #
        # locals should be the one passed to the expression
        self.assertTrue('file_ws' in interpolate('${locals().keys()}', globals(), lvar))

if __name__ == '__main__':
    unittest.main()
