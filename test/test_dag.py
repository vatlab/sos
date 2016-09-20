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
from pysos.signature import FileTarget


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
        for filename in ('a.txt', 'a1.txt'):
            with open(filename, 'w') as tmp:
                tmp.write('hey')
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
input: 'a.txt'

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
        for filename in ('a.txt', 'a1.txt'):
            os.remove(filename)

    def testUndetermined(self):
        '''Test DAG with undetermined input.'''
        #
        for filename in ('a.txt', 'd.txt'):
            with open(filename, 'w') as tmp:
                tmp.write('hey')
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
        for filename in ('a.txt', 'd.txt'):
            os.remove(filename)

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

    def testCycle(self):
        '''Test cycle detection of DAG'''
        #
        #  A.txt --> B.txt
        #
        #  B.txt --> C.txt
        #
        #  C.txt --> A.txt
        #
        script = SoS_Script('''
[A_1]
input: 'A.txt'
output: 'B.txt'

[A_2]
output: 'C.txt'

[A_3]
output: 'A.txt'
        ''')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow()
        self.assertRaises(RuntimeError, Sequential_Executor(wf).prepare)

    def testLongChain(self):
        '''Test long make file style dependencies.'''
        #
        for f in ['A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
            FileTarget(f).remove('both')
        #
        #  A1 <- B1 <- B2 <- B3 
        #   |
        #   |
        #  \/
        #  A2 <- B2 <- C1 <- C2 <- C4
        #                    C3
        #
        script = SoS_Script('''
[A_1]
input: 'B1.txt'
output: 'A1.txt'
sh:
    touch A1.txt

[A_2]
depends:  'B2.txt'
sh:
    touch A2.txt

[B1: provides='B1.txt']
depends: 'B2.txt'
sh:
    touch B1.txt

[B2: provides='B2.txt']
depends: 'B3.txt', 'C1.txt'
sh:
    touch B2.txt

[B3: provides='B3.txt']
sh:
    touch B3.txt

[C1: provides='C1.txt']
depends: 'C2.txt', 'C3.txt'
sh:
    touch C1.txt

[C2: provides='C2.txt']
depends: 'C4.txt'
sh:
    touch C2.txt

[C3: provides='C3.txt']
depends: 'C4.txt'
sh:
    touch C3.txt

[C4: provides='C4.txt']
sh:
    touch C4.txt

        ''')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''
strict digraph "" {
A_2;
C3;
C4;
B2;
B3;
A_1;
C1;
B1;
C2;
C3 -> C1;
C4 -> C3;
C4 -> C2;
B2 -> A_2;
B2 -> B1;
B3 -> B2;
A_1 -> A_2;
C1 -> B2;
B1 -> A_1;
C2 -> C1;
}
''')
        Sequential_Executor(wf).run(dag)
        for f in ['A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
            t = FileTarget(f)
            self.assertTrue(t.exists())
            t.remove('both')
    
    def testTarget(self):
        '''Test executing only part of a workflow.'''
        #
        for f in ['A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
            FileTarget(f).remove('both')
        #
        #  A1 <- B1 <- B2 <- B3 
        #   |
        #   |
        #  \/
        #  A2 <- B2 <- C1 <- C2 <- C4
        #                    C3
        #
        script = SoS_Script('''
[A_1]
input: 'B1.txt'
output: 'A1.txt'
sh:
    touch A1.txt

[A_2]
depends:  'B2.txt'
sh:
    touch A2.txt

[B1: provides='B1.txt']
depends: 'B2.txt'
sh:
    touch B1.txt

[B2: provides='B2.txt']
depends: 'B3.txt', 'C1.txt'
sh:
    touch B2.txt

[B3: provides='B3.txt']
sh:
    touch B3.txt

[C1: provides='C1.txt']
depends: 'C2.txt', 'C3.txt'
sh:
    touch C1.txt

[C2: provides='C2.txt']
depends: 'C4.txt'
sh:
    touch C2.txt

[C3: provides='C3.txt']
depends: 'C4.txt'
sh:
    touch C3.txt

[C4: provides='C4.txt']
sh:
    touch C4.txt

        ''')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow()
        #
        # test 1, we only need to generate target 'B1.txt'
        dag = Sequential_Executor(wf).prepare(targets=['B1.txt'])
        # note that A2 is no longer mentioned
        self.assertDAG(dag,
'''
strict digraph "" {
B3;
B1;
B1;
C2;
C4;
B2;
C1;
C3;
B3 -> B2;
C2 -> C1;
C4 -> C3;
C4 -> C2;
B2 -> B1;
B2 -> B1;
C1 -> B2;
C3 -> C1;
}
''')
        Sequential_Executor(wf).run(dag)
        for f in ['A1.txt', 'A2.txt']:
            self.assertFalse(FileTarget(f).exists())
        for f in ['C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
            t = FileTarget(f)
            self.assertTrue(t.exists())
            t.remove('both')
        #
        # test 2, we would like to generate two files
        dag = Sequential_Executor(wf).prepare(targets=['B2.txt', 'C2.txt'])
        # note that A2 is no longer mentioned
        self.assertDAG(dag,
'''
strict digraph "" {
C1;
C2;
B3;
B2;
B2;
C4;
C3;
C1 -> B2;
C1 -> B2;
C3 -> C1;
B3 -> B2;
B3 -> B2;
C4 -> C3;
C4 -> C2;
C2 -> C1;
}
''')
        Sequential_Executor(wf).run(dag)
        for f in ['A1.txt', 'B1.txt', 'A2.txt']:
            self.assertFalse(FileTarget(f).exists())
        for f in ['C2.txt', 'B2.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
            t = FileTarget(f)
            self.assertTrue(t.exists())
            t.remove('both')
        # 
        # test 3, generate two separate trees
        #
        dag = Sequential_Executor(wf).prepare(targets=['B3.txt', 'C2.txt'])
        # note that A2 is no longer mentioned
        self.assertDAG(dag,
'''
strict digraph "" {
B3;
C2;
C4;
C4 -> C2;
}
''')
        Sequential_Executor(wf).run(dag)
        for f in ['A1.txt', 'B1.txt', 'A2.txt', 'B2.txt', 'C1.txt', 'C3.txt']:
            self.assertFalse(FileTarget(f).exists())
        for f in ['C2.txt', 'B3.txt', 'C4.txt']:
            t = FileTarget(f)
            self.assertTrue(t.exists())
            t.remove('both')

    def testPatternReuse(self):
        '''Test repeated use of steps that use pattern and produce different files.'''
        #
        for f in ['A1.txt', 'A2.txt', 'B1.txt', 'B1.txt.p', 'B2.txt', 'B2.txt.p']:
            FileTarget(f).remove('both')
        #
        #  A1 <- P <- B1
        #  A1 <- P <- B2
        #  A2  
        #
        script = SoS_Script('''
[A_1]
input: 'B1.txt.p', 'B2.txt.p'
output: 'A1.txt'
sh:
    touch A1.txt

[A_2]
sh:
    touch A2.txt

[B1: provides='B1.txt']
sh:
    touch B1.txt

[B2: provides='B2.txt']
sh:
    touch B2.txt

[P: provides='{filename}.p']
input: filename
sh:
    touch ${output}
''')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''
strict digraph "" {
B2;
A_2;
P;
B1;
P;
A_1;
B2 -> P;
P -> A_1;
B1 -> P;
P -> A_1;
A_1 -> A_2;
}
''')
        Sequential_Executor(wf).run(dag)
        for f in ['A1.txt', 'A2.txt', 'B1.txt', 'B1.txt.p', 'B2.txt', 'B2.txt.p']:
            t = FileTarget(f)
            self.assertTrue(t.exists())
            t.remove('both')


    def testParallelExecution(self):
        '''Test basic parallel execution'''
        '''
        A1 <- None
        A2 <- B2
        '''
        for f in ['A1.txt', 'B2.txt', 'A2.txt']:
            FileTarget(f).remove('both')
        script = SoS_Script('''
[A_1]
output: 'A1.txt'
sh:
    sleep 3
    touch A1.txt

[A_2]
input:  'B2.txt'
output: 'A2.txt'
sh:
    sleep 3
    touch A2.txt

[B: provides='B2.txt']
output: 'B2.txt'
sh:
    touch B2.txt


''')
        # the workflow should call step K for step C_2, but not C_3
        wf = script.workflow()
        dag = Sequential_Executor(wf).prepare()
        self.assertDAG(dag,
'''
strict digraph "" {
B;
A_2;
A_1;
B -> A_2;
}
''')
        env.max_jobs = 4
        st = time.time()
        Sequential_Executor(wf).run(dag)
        self.assertLess(time.time() - st, 4)
        for f in ['A1.txt', 'B2.txt', 'A2.txt']:
            FileTarget(f).remove('both')


if __name__ == '__main__':
    unittest.main()
