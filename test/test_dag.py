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
import time
import glob
import unittest
import shutil
from io import StringIO

from pysos import SoS_Script
from pysos.dag import SoS_DAG
from pysos.utils import env
from pysos.sos_eval import Undetermined
from pysos.sos_executor import Sequential_Executor, Interactive_Executor, ExecuteError
from pysos.sos_script import ParsingError


import matplotlib.pyplot as plt


class TestDAG(unittest.TestCase):
    def assertDAG(self, dag, content):
        out = StringIO()
        dag.write_dot(out)
        dot = out.getvalue()
        self.assertEqual(sorted([x.strip() for x in dot.split('\n') if x.strip()]),
            sorted([x.strip() for x in content.split('\n') if x.strip()]))

    def testSimpleDag(self):
        '''Test DAG with simple dependency'''
        # basica case
        # 1 -> 2 -> 3 -> 4
        script = SoS_Script('''
[A_1]

[A_2]

[A_3]

[A_4]

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''strict digraph "" {
A_2;
A_4;
A_1;
A_3;
A_2 -> A_3;
A_1 -> A_2;
A_3 -> A_4;
}
''')
        # basica case
        # 1 -> 2 -> 3 -> 4
        script = SoS_Script('''
[A_1]

[A_2]

[A_3]
input: 'c.txt'

[A_4]

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''strict digraph "" {
A_2;
A_4;
A_1;
A_3;
A_1 -> A_2;
A_3 -> A_4;
}
''')


        #
        # 1 -> 2 -> 3 -> 4
        #
        script = SoS_Script('''
[A_1]
input: 'a.txt'
output: 'b.txt'

[A_2]
input: 'b.txt'
output: 'c.txt'

[A_3]
input: 'c.txt'
output: 'd.txt'

[A_4]
input: 'd.txt'
output: 'e.txt'

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''strict digraph "" {
A_2;
A_4;
A_1;
A_3;
A_2 -> A_3;
A_1 -> A_2;
A_3 -> A_4;
}
''')
        #
        # 1 -> 2
        # 3 -> 4 (3 does not have any input)
        #
        script = SoS_Script('''
[B_1]
input: 'a.txt'
output: 'b.txt'

[B_2]
input: 'b.txt'
output: 'c.txt'

[B_3]
input: []
output: 'd.txt'

[B_4]
input: 'd.txt'
output: 'e.txt'

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''strict digraph "" {
B_2;
B_4;
B_1;
B_3;
B_1 -> B_2;
B_3 -> B_4;
}
''')
        #
        # 1 -> 2
        # 3 -> 4 (3 depends on something else)
        #
        script = SoS_Script('''
[B_1]
input: 'a.txt'
output: 'b.txt'

[B_2]
input: 'b.txt'
output: 'c.txt'

[B_3]
input: 'a1.txt'
output: 'd.txt'

[B_4]
input: 'd.txt'
output: 'e.txt'

        ''')

        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''strict digraph "" {
B_1;
B_4;
B_2;
B_3;
B_1 -> B_2;
B_3 -> B_4;
}
''')
        #
        # (1) -> 2
        # (1) -> 3 -> 4
        #
        # 2 and 3 depends on the output of 1
        script = SoS_Script('''
[C_1]
input: 'a.txt'
output: 'b.txt'

[C_2]
input: 'b.txt'
output: 'c.txt'

[C_3]
input:  'b.txt'
output: 'd.txt'

[C_4]
depends: 'd.txt'
output: 'e.txt'

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''
strict digraph "" {
C_1;
C_4;
C_2;
C_3;
C_1 -> C_2;
C_1 -> C_3;
C_3 -> C_4;
}
''')

    def testUndetermined(self):
        '''Test DAG with undetermined input.'''
        #
        # input of step 3 is undertermined so
        # it depends on all its previous steps.
        script = SoS_Script('''
[C_1]
input: 'a.txt'
output: 'b.txt'

[C_2]
input: 'b.txt'
output: 'c.txt'

[C_3]
input:  dynamic('*.txt')
output: 'd.txt'

[C_4]
depends: 'd.txt'
output: 'e.txt'

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        dag.show_nodes()
        dag.write_dot('a.dot')
        self.assertDAG(dag,
'''
strict digraph "" {
C_1;
C_4;
C_2;
C_3;
C_1 -> C_2;
C_2 -> C_3;
C_3 -> C_4;
}
''')
        #
        # output of step
        #
        script = SoS_Script('''
[C_1]
input: 'a.txt'
output: 'b.txt'

[C_2]
input: 'b.txt'
output: 'c.txt'

[C_3]
input:  dynamic('*.txt')

[C_4]
depends: 'd.txt'
output: 'e.txt'

        ''')
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''
strict digraph "" {
C_1;
C_4;
C_2;
C_3;
C_1 -> C_2;
C_2 -> C_3;
C_3 -> C_4;
}
''')

    def testAuxiliarySteps(self):
        script = SoS_Script('''
[K: provides='{name}.txt']
output: '${name}.txt'

sh:
    touch '${name}.txt'

[C_2]
input: 'b.txt'
output: 'c.txt'

sh:
    touch c.txt

[C_3]
input: 'a.txt'

        ''')
        # a.txt exists and b.txt does not exist
        with open('a.txt', 'w') as atfile:
            atfile.write('garbage')
        if os.path.isfile('b.txt'):
            os.remove('b.txt')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        #dag.write_dot('a.dot')
        #dag.show_nodes()
        self.assertDAG(dag,
'''
strict digraph "" {
K;
C_2;
C_3;
K -> C_2;
}
''')

if __name__ == '__main__':
    unittest.main()
