#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import textwrap
from io import StringIO

import pytest

from sos import execute_workflow
from sos.parser import SoS_Script
from sos.targets import file_target
from sos.utils import env
# if the test is imported under sos/test, test interacive executor
from sos.workflow_executor import Base_Executor


def assertDAG(dag, content):
    '''Compare the content of dag in a file, to one of more DAG (strings).'''
    if isinstance(dag, str):
        with open(dag) as d:
            # only get the first DAG
            dot = 'strict' + d.read().split('strict')[1]
    else:
        out = StringIO()
        dag.save(out)
        dot = out.getvalue()

    def sorted_dot(dot):
        return sorted([
            x.strip()
            for x in dot.split('\n')
            if x.strip() and not 'digraph' in x
        ])

    if isinstance(content, str):
        assert sorted_dot(dot) == sorted_dot(content)
    else:
        assert sorted_dot(dot) in [sorted_dot(x) for x in content]


def get_initial_dag(test):
    test = textwrap.dedent(test)
    script = SoS_Script(test)
    wf = script.workflow()
    dag = Base_Executor(wf).initialize_dag()
    return dag


def test_simple_dag(temp_factory):
    '''Test DAG with simple dependency'''
    temp_factory('a.txt', 'a1.txt')
    # basica case
    # 1 -> 2 -> 3 -> 4

    assertDAG(
        get_initial_dag('''
            [A_1]

            [A_2]

            [A_3]

            [A_4]'''),
        textwrap.dedent('''strict digraph "" {
            A_2;
            A_4;
            A_1;
            A_3;
            A_2 -> A_3;
            A_1 -> A_2;
            A_3 -> A_4;
            }
            '''))
    # basica case
    # 1 -> 2 -> 3 -> 4
    assertDAG(
        get_initial_dag('''
            [A_1]

            [A_2]

            [A_3]
            input: 'a.txt'

            [A_4]'''),
        textwrap.dedent('''strict digraph "" {
            A_2;
            A_4;
            A_1;
            A_3;
            A_1 -> A_2;
            A_3 -> A_4;
            }
            '''))
    #
    # 1 -> 2 -> 3 -> 4
    #
    assertDAG(
        get_initial_dag('''
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
            '''),
        textwrap.dedent('''strict digraph "" {
            A_2;
            A_4;
            A_1;
            A_3;
            A_2 -> A_3;
            A_1 -> A_2;
            A_3 -> A_4;
            }
            '''))
    #
    # 1 -> 2
    # 3 -> 4 (3 does not have any input)
    #
    assertDAG(
        get_initial_dag('''
            [B_1]
            input: 'a.txt'
            output: 'b.txt'

            [B_2]
            input: 'b.txt'
            output: 'c.txt'

            [B_3]
            input: None
            output: 'd.txt'

            [B_4]
            input: 'd.txt'
            output: 'e.txt'
            '''),
        textwrap.dedent('''strict digraph "" {
            B_2;
            B_4;
            B_1;
            B_3;
            B_1 -> B_2;
            B_3 -> B_4;
            }
            '''))
    #
    # 1 -> 2
    # 3 -> 4 (3 depends on something else)
    #
    assertDAG(
        get_initial_dag('''
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
            '''),
        textwrap.dedent('''strict digraph "" {
            B_1;
            B_4;
            B_2;
            B_3;
            B_1 -> B_2;
            B_3 -> B_4;
            }
            '''))
    #
    # (1) -> 2
    # (1) -> 3 -> 4
    #
    # 2 and 3 depends on the output of 1
    assertDAG(
        get_initial_dag('''
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
            '''),
        textwrap.dedent('''
            strict digraph "" {
            C_1;
            C_4;
            C_2;
            C_3;
            C_1 -> C_2;
            C_1 -> C_3;
            C_3 -> C_4;
            }
            '''))


def test_undetermined(temp_factory):
    '''Test DAG with undetermined input.'''
    #
    temp_factory('a.txt', 'd.txt')
    # input of step 3 is undertermined so
    # it depends on all its previous steps.
    assertDAG(
        get_initial_dag('''
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
        # dag.show_nodes()
        # dag.save('a.dot')
        ,
        textwrap.dedent('''
            strict digraph "" {
            C_1;
            C_4;
            C_2;
            C_3;
            C_1 -> C_2;
            C_2 -> C_3;
            C_3 -> C_4;
            }
            '''))
    #
    # output of step
    #
    assertDAG(
        get_initial_dag('''
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
            '''),
        textwrap.dedent('''
            strict digraph "" {
            C_1;
            C_4;
            C_2;
            C_3;
            C_1 -> C_2;
            C_2 -> C_3;
            C_3 -> C_4;
            }
            '''))


def test_auxiliary_steps(temp_factory, clear_now_and_after):
    graph = textwrap.dedent(('''
        [K: provides='{name}.txt']
        output: f"{name}.txt"

        run: expand=True
        touch '{name}.txt'

        [C_2]
        input: 'b.txt'
        output: 'c.txt'

        run:
        touch c.txt

        [C_3]
        input: 'a.txt'

        '''))
    # a.txt exists and b.txt does not exist
    temp_factory('a.txt')
    clear_now_and_after('b.txt')
    # the workflow should call step K for step C_2, but not C_3
    dag = get_initial_dag(graph)
    #
    # Ticket 363:
    #
    # we have two possibilities here, one is to ignore a.txt,
    # and one is to regenerate a.txt because it is not generated
    # by sos (without signature)        #
    #
    dag.show_nodes()
    assertDAG(
        dag,
        textwrap.dedent('''
        strict digraph "" {
        "K (b.txt)";
        C_3;
        C_2;
        "K (b.txt)" -> C_2;
        }
        '''))


def test_cycle():
    '''Test cycle detection of DAG'''
    #
    #  A.txt --> B.txt
    #
    #  B.txt --> C.txt
    #
    #  C.txt --> A.txt
    #
    script = SoS_Script(
        textwrap.dedent('''
    [A_1]
    input: 'A.txt'
    output: 'B.txt'

    [A_2]
    output: 'C.txt'

    [A_3]
    output: 'A.txt'
    '''))
    # the workflow should call step K for step C_2, but not C_3
    wf = script.workflow()
    pytest.raises(RuntimeError, Base_Executor(wf).initialize_dag)


def test_long_chain(clear_now_and_after):
    '''Test long make file style dependencies.'''
    #
    clear_now_and_after('A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt',
                        'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt')

    #
    #  A1 <- B1 <- B2 <- B3
    #   |
    #   |
    #  \/
    #  A2 <- B2 <- C1 <- C2 <- C4
    #                    C3
    #
    # the workflow should call step K for step C_2, but not C_3
    #env.verbosity = 4
    assertDAG(
        get_initial_dag('''
        [A_1]
        input: 'B1.txt'
        output: 'A1.txt'
        run:
            touch A1.txt

        [A_2]
        depends:  'B2.txt'
        output: 'A2.txt'
        run:
            touch A2.txt

        [B1: provides='B1.txt']
        depends: 'B2.txt'
        run:
            touch B1.txt

        [B2: provides='B2.txt']
        depends: 'B3.txt', 'C1.txt'
        run:
            touch B2.txt

        [B3: provides='B3.txt']
        run:
            touch B3.txt

        [C1: provides='C1.txt']
        depends: 'C2.txt', 'C3.txt'
        run:
            touch C1.txt

        [C2: provides='C2.txt']
        depends: 'C4.txt'
        run:
            touch C2.txt

        [C3: provides='C3.txt']
        depends: 'C4.txt'
        run:
            touch C3.txt

        [C4: provides='C4.txt']
        run:
            touch C4.txt

        '''),
        textwrap.dedent('''
        strict digraph "" {
        "C4 (C4.txt)";
        "B1 (B1.txt)";
        "C1 (C1.txt)";
        "C2 (C2.txt)";
        "C3 (C3.txt)";
        A_1;
        "B2 (B2.txt)";
        "B3 (B3.txt)";
        A_2;
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        "C4 (C4.txt)" -> "C3 (C3.txt)";
        "B1 (B1.txt)" -> A_1;
        "C1 (C1.txt)" -> "B2 (B2.txt)";
        "C2 (C2.txt)" -> "C1 (C1.txt)";
        "C3 (C3.txt)" -> "C1 (C1.txt)";
        A_1 -> A_2;
        "B2 (B2.txt)" -> "B1 (B1.txt)";
        "B2 (B2.txt)" -> A_2;
        "B3 (B3.txt)" -> "B2 (B2.txt)";
        }
        '''))


def test_target(clear_now_and_after):
    '''Test executing only part of a workflow.'''
    #
    clear_now_and_after('A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt',
                        'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt')
    #
    #  A1 <- B1 <- B2 <- B3
    #   |
    #   |
    #  \/
    #  A2 <- B2 <- C1 <- C2 <- C4
    #                    C3
    #
    script = SoS_Script(
        textwrap.dedent('''
    [A1]
    input: 'B1.txt'
    output: 'A1.txt'
    run:
        touch A1.txt

    [A2]
    depends:  'B2.txt'
    run:
        touch A2.txt

    [B1: provides='B1.txt']
    depends: 'B2.txt'
    run:
        touch B1.txt

    [B2: provides='B2.txt']
    depends: 'B3.txt', 'C1.txt'
    run:
        touch B2.txt

    [B3: provides='B3.txt']
    run:
        touch B3.txt

    [C1: provides='C1.txt']
    depends: 'C2.txt', 'C3.txt'
    run:
        touch C1.txt

    [C2: provides='C2.txt']
    depends: 'C4.txt'
    run:
        touch C2.txt

    [C3: provides='C3.txt']
    depends: 'C4.txt'
    run:
        touch C3.txt

    [C4: provides='C4.txt']
    run:
        touch C4.txt

    '''))
    # the workflow should call step K for step C_2, but not C_3
    wf = script.workflow(use_default=False)
    #
    # test 1, we only need to generate target 'B1.txt'
    dag = Base_Executor(wf).initialize_dag(targets=['B1.txt'])
    # note that A2 is no longer mentioned
    assertDAG(
        dag,
        textwrap.dedent('''
        strict digraph "" {
        "B3 (B3.txt)";
        "C4 (C4.txt)";
        "C2 (C2.txt)";
        "C1 (C1.txt)";
        "B1 (B1.txt)";
        "B2 (B2.txt)";
        "C3 (C3.txt)";
        "B3 (B3.txt)" -> "B2 (B2.txt)";
        "C4 (C4.txt)" -> "C3 (C3.txt)";
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        "C2 (C2.txt)" -> "C1 (C1.txt)";
        "C1 (C1.txt)" -> "B2 (B2.txt)";
        "B2 (B2.txt)" -> "B1 (B1.txt)";
        "C3 (C3.txt)" -> "C1 (C1.txt)";
        }
        '''))
    #
    # test 2, we would like to generate two files
    dag = Base_Executor(wf).initialize_dag(targets=['B2.txt', 'C2.txt'])
    # note that A2 is no longer mentioned
    assertDAG(
        dag,
        textwrap.dedent('''
        strict digraph "" {
        "C4 (C4.txt)";
        "B2 (B2.txt)";
        "C3 (C3.txt)";
        "B3 (B3.txt)";
        "C2 (C2.txt)";
        "C1 (C1.txt)";
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        "C4 (C4.txt)" -> "C3 (C3.txt)";
        "C3 (C3.txt)" -> "C1 (C1.txt)";
        "B3 (B3.txt)" -> "B2 (B2.txt)";
        "C2 (C2.txt)" -> "C1 (C1.txt)";
        "C1 (C1.txt)" -> "B2 (B2.txt)";
        }
        '''))
    #
    # test 3, generate two separate trees
    #
    dag = Base_Executor(wf).initialize_dag(targets=['B3.txt', 'C2.txt'])
    # note that A2 is no longer mentioned
    assertDAG(
        dag,
        textwrap.dedent('''
        strict digraph "" {
        "B3 (B3.txt)";
        "C2 (C2.txt)";
        "C4 (C4.txt)";
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        }
        '''))


def test_pattern_reuse(clear_now_and_after):
    '''Test repeated use of steps that use pattern and produce different files.'''
    #
    clear_now_and_after('A1.txt', 'A2.txt', 'B1.txt', 'B1.txt.p', 'B2.txt',
                        'B2.txt.p')
    #
    #  A1 <- P <- B1
    #  A1 <- P <- B2
    #  A2
    #
    # the workflow should call step K for step C_2, but not C_3
    assertDAG(
        get_initial_dag('''
        [A_1]
        input: 'B1.txt.p', 'B2.txt.p'
        output: 'A1.txt'
        run:
            touch A1.txt

        [A_2]
        output: 'A2.txt'
        run:
            touch A2.txt

        [B1: provides='B1.txt']
        run:
            touch B1.txt

        [B2: provides='B2.txt']
        run:
            touch B2.txt

        [P: provides='{filename}.p']
        input: filename
        run: expand=True
            touch {_output}
        '''),
        textwrap.dedent('''
        strict digraph "" {
        "P (B2.txt.p)";
        "B1 (B1.txt)";
        "B2 (B2.txt)";
        A_2;
        A_1;
        "P (B1.txt.p)";
        "P (B2.txt.p)" -> A_1;
        "B1 (B1.txt)" -> "P (B1.txt.p)";
        "B2 (B2.txt)" -> "P (B2.txt.p)";
        A_1 -> A_2;
        "P (B1.txt.p)" -> A_1;
        }
        '''))


def test_parallel_execution(clear_now_and_after):
    '''Test basic parallel execution
    A1 <- None
    A2 <- B2
    '''
    clear_now_and_after('A1.txt', 'B2.txt', 'A2.txt')
    # the workflow should call step K for step C_2, but not C_3
    assertDAG(
        get_initial_dag('''
        [A_1]
        output: 'A1.txt'
        run:
            sleep 0
            touch A1.txt

        [A_2]
        input:  'B2.txt'
        output: 'A2.txt'
        run:
            sleep 0
            touch A2.txt

        [B: provides='B2.txt']
        output: 'B2.txt'
        run:
            touch B2.txt

        '''),
        textwrap.dedent('''
        strict digraph "" {
        A_1;
        A_2;
        "B (B2.txt)";
        "B (B2.txt)" -> A_2;
        }
        '''))
    env.max_jobs = 4
    #env.verbosity = 4
    # the process is slower after switching to spawn mode


def test_shared_dependency(clear_now_and_after):
    #
    # shared variable should introduce additional dependency
    #
    clear_now_and_after('A1.txt')
    #
    # A1 introduces a shared variable ss, A3 depends on ss but not A2
    #
    script = SoS_Script(
        textwrap.dedent('''
    [A_1: shared='ss']
    ss = 'A1'

    [A_2]
    input: None

    run:
        sleep 0

    [A_3]
    input: None
    import time
    time.sleep(0)
    with open(f"{ss}.txt", 'w') as tmp:
        tmp.write('test')

    '''))
    wf = script.workflow('A')
    dag = Base_Executor(wf).initialize_dag()
    assertDAG(
        dag,
        textwrap.dedent('''
        strict digraph "" {
        A_3;
        A_1;
        A_2;
        A_1 -> A_3;
        }
        '''))
    env.max_jobs = 3


def test_literal_connection(clear_now_and_after):
    '''Testing the connection of steps with by variables.'''
    clear_now_and_after('A1.txt')
    #
    # A1 introduces a shared variable ss, A3 depends on ss but not A2
    #
    script = SoS_Script(
        textwrap.dedent('''
    [A_1: shared='p']
    run:
        touch 'A1.txt'

    p = 'A1.txt'

    [A_2]
    input: None

    run:
        sleep 0

    [A_3]
    input: p
    depends: sos_variable('p')

    run:
        sleep 0

    [A_4]
    input: p
    depends: sos_variable('p')
    run:
        sleep 0

    [A_5]
    input: dynamic(p)
    depends: sos_variable('p')
    '''))
    wf = script.workflow('A')
    dag = Base_Executor(wf).initialize_dag()
    assertDAG(
        dag,
        textwrap.dedent('''
        strict digraph "" {
        A_1;
        A_4;
        A_2;
        A_3;
        A_5;
        A_1 -> A_4;
        A_1 -> A_3;
        A_1 -> A_5;
        }
        '''))
    env.max_jobs = 3


def test_variable_target():
    '''Test dependency caused by variable usage.'''
    script = SoS_Script(
        textwrap.dedent(r'''
    [A: shared='b']
    b = 1

    [C: shared={'c':'k'}]
    k = 2

    [all: shared='p']
    depends: sos_variable('c'), sos_variable('b')

    p = c + b

    '''))
    wf = script.workflow('all')
    Base_Executor(wf).run()
    assert env.sos_dict['p'] == 3


def test_reverse_shared_variable(clear_now_and_after):
    '''Test shared variables defined in auxiliary steps'''
    clear_now_and_after('a.txt')
    script = SoS_Script(
        textwrap.dedent(r'''
    [A: shared='b', provides='a.txt']
    b = 1
    run:
    touch a.txt

    [B_1]
    depends: 'a.txt'

    [B_2]
    print(b)

    '''))
    wf = script.workflow('B')
    Base_Executor(wf).run()
    assert env.sos_dict['b'] == 1


def test_chained_depends(temp_factory):
    '''Test chain dependent'''
    temp_factory('a.bam', 'a.bam.bai', 'a.vcf')
    script = SoS_Script(
        textwrap.dedent(r'''
    # this step provides variable `var`
    [index: provides='{filename}.bam.bai']
    input: f"{filename}.bam"
    run: expand=True
    echo "Generating {_output}"
    touch {_output}

    [call: provides='{filename}.vcf']
    input:   f"{filename}.bam"
    depends: f"{_input}.bai"
    run: expand=True
    echo "Calling variants from {_input} with {_depends} to {_output}"
    touch {_output}
    '''))
    # if file_target('a.bam.bai').exists():
    #     file_target('a.bam.bai').unlink()
    # if file_target('a.vcf').exists():
    #     file_target('a.vcf').unlink()
    # self.touch('a.bam')
    Base_Executor(script.workflow()).run(targets=['a.vcf'])
    # for f in ('a.vcf', 'a.bam', 'a.bam.bai'):
    #     if file_target(f).exists():
    #         file_target(f).unlink()


@pytest.mark.skipif(True, reason='This test is failing')
def test_output_of_dag(clear_now_and_after):
    '''Test output of dag'''
    #
    #for f in ['A1.txt', 'A2.txt', 'C2.txt', 'B2.txt', 'B1.txt', 'B3.txt', 'C1.txt', 'C3.txt', 'C4.txt']:
    #    if file_target(f).exists():
    #        file_target(f).unlink()
    #
    #  A1 <- B1 <- B2 <- B3
    #   |
    #   |
    #  \/
    #  A2 <- B2 <- C1 <- C2 <- C4
    #                    C3
    #
    script = SoS_Script(
        textwrap.dedent('''
    [A_1]
    input: 'B1.txt'
    output: 'A1.txt'
    run:
        touch A1.txt

    [A_2]
    depends:  'B2.txt'
    run:
        touch A2.txt

    [B1: provides='B1.txt']
    depends: 'B2.txt'
    run:
        touch B1.txt

    [B2: provides='B2.txt']
    depends: 'B3.txt', 'C1.txt'
    run:
        touch B2.txt

    [B3: provides='B3.txt']
    run:
        touch B3.txt

    [C1: provides='C1.txt']
    depends: 'C2.txt', 'C3.txt'
    run:
        touch C1.txt

    [C2: provides='C2.txt']
    depends: 'C4.txt'
    run:
        touch C2.txt

    [C3: provides='C3.txt']
    depends: 'C4.txt'
    run:
        touch C3.txt

    [C4: provides='C4.txt']
    run:
        touch C4.txt

    '''))
    # the workflow should call step K for step C_2, but not C_3
    wf = script.workflow(use_default=False)
    #
    # test 1, we only need to generate target 'B1.txt'
    Base_Executor(
        wf, config={
            'output_dag': 'test_outofdag1.dot',
            'trace_existing': True
        }).initialize_dag(targets=['B1.txt'])
    # note that A2 is no longer mentioned
    assertDAG(
        'test_outofdag1.dot',
        textwrap.dedent('''
        strict digraph "" {
        "B3 (B3.txt)";
        "C4 (C4.txt)";
        "C2 (C2.txt)";
        "C1 (C1.txt)";
        "B1 (B1.txt)";
        "B2 (B2.txt)";
        "C3 (C3.txt)";
        "B3 (B3.txt)" -> "B2 (B2.txt)";
        "C4 (C4.txt)" -> "C3 (C3.txt)";
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        "C2 (C2.txt)" -> "C1 (C1.txt)";
        "C1 (C1.txt)" -> "B2 (B2.txt)";
        "B2 (B2.txt)" -> "B1 (B1.txt)";
        "C3 (C3.txt)" -> "C1 (C1.txt)";
        }
        '''))
    # test 2, we would like to generate two files
    Base_Executor(
        wf, config={
            'output_dag': 'test_outofdag2.dot',
            'trace_existing': True
        }).initialize_dag(targets=['B2.txt', 'C2.txt'])
    # note that A2 is no longer mentioned
    assertDAG(
        'test_outofdag2.dot',
        textwrap.dedent('''
        strict digraph "" {
        "C4 (C4.txt)";
        "B2 (B2.txt)";
        "C3 (C3.txt)";
        "B3 (B3.txt)";
        "C2 (C2.txt)";
        "C1 (C1.txt)";
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        "C4 (C4.txt)" -> "C3 (C3.txt)";
        "C3 (C3.txt)" -> "C1 (C1.txt)";
        "B3 (B3.txt)" -> "B2 (B2.txt)";
        "C2 (C2.txt)" -> "C1 (C1.txt)";
        "C1 (C1.txt)" -> "B2 (B2.txt)";
        }
        '''))
    # test 3, generate two separate trees
    #
    Base_Executor(
        wf, config={
            'output_dag': 'test_outofdag3.dot',
            'trace_existing': True
        }).initialize_dag(targets=['B3.txt', 'C2.txt'])
    # note that A2 is no longer mentioned
    assertDAG(
        'test_outofdag3.dot',
        textwrap.dedent('''
        strict digraph "" {
        "B3 (B3.txt)";
        "C2 (C2.txt)";
        "C4 (C4.txt)";
        "C4 (C4.txt)" -> "C2 (C2.txt)";
        }
        '''))
    clear_now_and_after('C2.txt', 'B3.txt', 'C4.txt', 'test.dot', 'test_2.dot')


@pytest.mark.skipif(True, reason='This test is failing')
def test_step_with_multiple_output(clear_now_and_after):
    '''Test addition of steps with multiple outputs. It should be added only once'''
    script = SoS_Script(
        textwrap.dedent('''
    [test_1: provides=['{}.txt'.format(i) for i in range(10)]]
    output: ['{}.txt'.format(i) for i in range(10)]
    run:
    touch {output}

    [test_2: provides=['{}.txt'.format(i) for i in range(10, 20)]]
    depends: ['{}.txt'.format(i) for i in range(10)]
    output: ['{}.txt'.format(i) for i in range(10, 20)]
    run:
    touch {output}

    [default]
    depends: ['{}.txt'.format(i) for i in range(10, 20)]
    '''))
    wf = script.workflow()
    Base_Executor(wf, config={'output_dag': 'test.dot'}).initialize_dag()
    with open('test.dot') as dot:
        lc = len(dot.readlines())
    assert lc == 6
    clear_now_and_after('test.dot')


def test_auxiliary_sos_step():
    '''Testing the use of sos_step with auxiliary step. #736'''
    execute_workflow('''
        [default]
        depends: '1.txt'

        [A_1]
        print("Hi")


        [C_1: provides = "1.txt"]
        depends: sos_step("A_1")
        run:
        touch 1.txt
        ''')


def test_forward_style_depend(clear_now_and_after, temp_factory):
    '''Test the execution of forward-style workflow with undtermined dependency'''
    clear_now_and_after('a.txt.bak')
    temp_factory('a.txt')
    execute_workflow('''
        [10]
        input: 'a.txt'
        output: f"{_input}.bak"
        run: expand=True
            cp {_input} {_output}

        [20]
        depends: "a.txt.bak"
        run: expand=True
            ls {_depends}
        ''')
    assert file_target('a.txt.bak').target_exists()


def test_sos_step_miniworkflow(clear_now_and_after):
    '''Test the addition of mini forward workflows introduced by sos_step'''
    script = SoS_Script(
        textwrap.dedent('''
    [a_1]
    print(step_name)

    [a_2]
    print(step_name)
    [a_20]
    print(step_name)

    [b_1]
    print(step_name)

    [b_2]
    print(step_name)

    [b_20]
    depends: sos_step('c')
    print(step_name)

    [c_1]
    print(step_name)

    [c_2]
    print(step_name)

    [c_20]
    print(step_name)

    [default]
    depends: sos_step('a'), sos_step('b')
    '''))
    wf = script.workflow()
    Base_Executor(wf, config={'output_dag': 'test.dot'}).run()
    # note that A2 is no longer mentioned
    assertDAG(
        'test.dot',
        textwrap.dedent('''
        strict digraph "" {
        default;
        a_1;
        a_2;
        a_20;
        b_1;
        b_2;
        b_20;
        c_1;
        c_2;
        c_20;
        a_1 -> a_2;
        a_2 -> a_20;
        a_20 -> default;
        b_1 -> b_2;
        b_2 -> b_20;
        b_20 -> default;
        c_1 -> c_2;
        c_2 -> c_20;
        c_20 -> b_20;
        }
        '''))
    clear_now_and_after('test.dot')


def test_compound_workflow(clear_now_and_after):
    '''Test the DAG of compound workflow'''
    script = SoS_Script(
        textwrap.dedent('''
    [A_1]
    [A_2]
    [B]
    '''))
    wf = script.workflow('A+B')
    dag = Base_Executor(wf).initialize_dag()
    assertDAG(
        dag,
        textwrap.dedent('''strict digraph "" {
        A_1;
        A_2;
        B;
        A_1 -> A_2;
        A_2 -> B;
        }'''))
    # with empty depends
    script = SoS_Script(
        textwrap.dedent('''
    [A_1]
    [A_2]
    [B]
    depends:
    '''))
    wf = script.workflow('A+B')
    dag = Base_Executor(wf).initialize_dag()
    assertDAG(
        dag,
        textwrap.dedent('''strict digraph "" {
        A_1;
        A_2;
        B;
        A_1 -> A_2;
        A_2 -> B;
        }'''))
    clear_now_and_after('a.txt')
    script = SoS_Script(
        textwrap.dedent('''
    [A_1]
    [A_2]
    [C]
    output: 'a.txt'
    _output.touch()

    [B]
    depends: 'a.txt'
    '''))
    # with more depends
    wf = script.workflow('A+B')
    dag = Base_Executor(wf).initialize_dag()
    assertDAG(
        dag,
        textwrap.dedent('''strict digraph "" {
        A_1;
        A_2;
        B;
        "C (a.txt)";
        A_1 -> A_2;
        A_2 -> B;
        "C (a.txt)" -> B;
        }'''))


def test_provides_sos_variable():
    '''Test provides non-filename targets #1341'''
    execute_workflow('''
        [count: provides=sos_variable('numNotebooks')]
        numNotebooks = 1

        [default]
        depends: sos_variable('numNotebooks')
        print(f"There are {numNotebooks} notebooks in this directory")
        ''')


def test_multi_named_output(clear_now_and_after):
    '''Test DAG built from multiple named_output #1166'''
    clear_now_and_after('a.txt', 'b.txt', 'test_named_output.dot')

    execute_workflow(
        '''
        [A]
        output: A='a.txt', B='b.txt'
        _output.touch()

        [default]
        input: named_output('A'), named_output('B')
        ''',
        options={'output_dag': 'test_named_output.dot'})
    # note that A2 is no longer mentioned
    assertDAG('test_named_output.dot', [
        '''
        strict digraph "" {
        default;
        "A (B)";
        "A (B)" -> default;
        }
        ''', '''
        strict digraph "" {
        default;
        "A (A)";
        "A (A)" -> default;
        }
        '''
    ])
