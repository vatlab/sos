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
from sos.eval import SoS_eval, accessed_vars, Undetermined, on_demand_options
from sos.parser import SoS_Script
from sos.workflow_executor import Base_Executor, analyze_section
from sos.targets import executable, remote, sos_targets

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
        self.assertEqual(accessed_vars('''a = 1'''), {'a'})
        self.assertEqual(accessed_vars('''a = b + 2.0'''), {'a', 'b'})
        self.assertEqual(accessed_vars('''a = "C"'''), {'a'})
        self.assertEqual(accessed_vars('''a = "C" + f"{D}"'''), {'a', 'D'})
        self.assertEqual(accessed_vars('''a = 1 + f"{D + 20:f}" '''), {'a', 'D'})
        self.assertEqual(accessed_vars('''k, "a.txt", "b.txt", skip=True, par=f(something) '''), {'k', 'f', 'something', '__NULLFUNC__'})
        # this is a complicated case because the actual variable depends on the
        # result of an expression... However, in the NO-evaluation case, this is
        # the best we can do.
        self.assertEqual(accessed_vars('''c + f"{D + 1}" '''), {'c', 'D'})

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
        # windows does not have #! mechanism so the python code
        # would incorrectly be executed as bat
        if sys.platform == 'win32':
            return
        for text in ('"""a"""', '"b"',
            r'"""\na\\nb"""', r"'''a\nb'''",
            """ "a'\\"='" """):
            script = SoS_Script(r'''
a = 1
run: expand=True
   #!/usr/bin/env python
   with open('tmp.txt', 'w') as tmp:
      tmp.write({} + '{}')
k = """b"""'''.format(text, '{a}')
)
            wf = script.workflow()
            Base_Executor(wf).run()
            with open('tmp.txt') as tmp:
                self.assertEqual(tmp.read(), eval(text) + '1')
        os.remove('tmp.txt')


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

[A_4]
input: None
task:
python: expand=True
   print(f'{output}')

[A_5]
task:
   print(f'{_output}')
''')
        wf = script.workflow('A')
        for section in wf.sections:
            res = analyze_section(section)
            if section.names[0][1] == '1':
                self.assertTrue(isinstance(res['step_input'], Undetermined))
                self.assertEqual(res['step_depends'], sos_targets())
                self.assertEqual(res['step_output'], sos_targets())
                self.assertEqual(res['environ_vars'], {'b', 'p1', 'infiles'})
                self.assertEqual(res['signature_vars'], {'c'})
                self.assertEqual(res['changed_vars'], {'b'})
            elif section.names[0][1] == '2':
                self.assertEqual(res['step_input'], sos_targets())
                self.assertEqual(res['step_depends'], sos_targets('some.txt', executable('ls')))
                self.assertTrue(isinstance(res['step_output'], Undetermined))
                # for_each will not be used for DAG
                self.assertEqual(res['environ_vars'], {'b', 'for_each', 'executable'})
                self.assertEqual(res['signature_vars'], {'r', 'time', 'random'})
                self.assertEqual(res['changed_vars'], set())
            elif section.names[0][1] == '4':
                self.assertTrue('output' in res['signature_vars'])
            elif section.names[0][1] == '5':
                self.assertTrue('output' not in res['signature_vars'])

    def testOnDemandOptions(self):
        '''Test options that are evaluated upon request.'''
        options = on_demand_options(
            {'a': '"est"', 'b': 'c', 'c': 'e + 2'} )
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
