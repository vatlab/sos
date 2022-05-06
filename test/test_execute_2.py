import os
import subprocess

import pytest
from sos import execute_workflow
from sos._version import __version__
from sos.parser import SoS_Script
from sos.utils import env
# if the test is imported under sos/test, test interacive executor
from sos.workflow_executor import Base_Executor


def test_for_each_nested_list(temp_factory):
    """Test for_each option of input"""
    # test dictionary format of for_each
    temp_factory("p1.txt", "p2.txt", "p3.txt", "a.txt", "b.txt")
    # test same-level for loop and parameter with nested list
    execute_workflow(
        r"""
        [0: shared=['processed']]
        files = ['a.txt', 'b.txt']

        input: files, for_each={'par':[(1, 2), (1, 3), (2, 3)], 'res': ['p1.txt', 'p2.txt', 'p3.txt']}
        output: res

        processed = (par, res)
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["step_processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]


def test_removal_of_output_from_failed_step(clear_now_and_after):
    """Test the removal of output files if a step fails #1055"""
    clear_now_and_after("failed.csv", "result.csv")
    script = """
[sub: provides='{file}.csv']
sh: expand=True
  touch {_output}
  eco "something wrong"

[step]
depends: 'failed.csv'
path('result.csv').touch()
"""
    with pytest.raises(Exception):
        execute_workflow(script, workflow='step')
    # rerun should still raise
    with pytest.raises(Exception):
        execute_workflow(script, workflow='step')

    assert not os.path.isfile("failed.csv")
    assert not os.path.isfile("result.csv")


def test_error_handling_of_substeps(clear_now_and_after):
    clear_now_and_after(
        [f"test_{i}.txt" for i in range(10)],
        [f"test_{i}.bak" for i in range(10)],
        "test_30.txt",
    )

    script = r"""
import time

[10]
input: for_each=dict(i=range(10)), concurrent=False
output: f'test_{i}.txt'

_output.touch()

fail_if(_index == 5, 'fail at 5')

[20]
output: _input.with_suffix('.bak')
_output.touch()

[30]
input: None
output: 'test_30.txt'
time.sleep(2)
_output.touch()
"""

    #
    # default mode
    #
    with pytest.raises(Exception):
        execute_workflow(script)

    for i in range(10):
        if i == 5:
            assert not os.path.isfile(f"test_{i}.txt")
            assert not os.path.isfile(f"test_{i}.bak")
        else:
            assert os.path.isfile(f"test_{i}.txt")
            # the following step is not executed
            assert not os.path.isfile(f"test_{i}.bak")
    # but the another branch continues
    assert os.path.isfile("test_30.txt")
    #
    # ignore mode
    #
    clear_now_and_after(
        [f"test_{i}.txt" for i in range(10)],
        [f"test_{i}.bak" for i in range(10)],
        "test_30.txt",
    )

    execute_workflow(script, options={"error_mode": "ignore"})
    for i in range(6, 10):
        if i == 5:
            assert not os.path.isfile(f"test_{i}.txt")
            assert not os.path.isfile(f"test_{i}.bak")
        else:
            assert os.path.isfile(f"test_{i}.txt")
            assert os.path.isfile(f"test_{i}.bak")
    assert os.path.isfile("test_30.txt")
    #
    # abort mode
    #
    clear_now_and_after(
        [f"test_{i}.txt" for i in range(10)],
        [f"test_{i}.bak" for i in range(10)],
        "test_30.txt",
    )

    with pytest.raises(Exception):
        execute_workflow(script, options={"error_mode": "abort"})
    for i in range(10):
        if i < 5:
            assert os.path.isfile(f"test_{i}.txt")
            assert not os.path.isfile(f"test_{i}.bak")
        else:
            assert not os.path.isfile(f"test_{i}.txt")
            assert not os.path.isfile(f"test_{i}.bak")
    assert not os.path.isfile("test_30.txt")


def test_parallel_nestedworkflow(clear_now_and_after):
    # 1375
    clear_now_and_after([f'{i+1}.txt' for i in range(5)])
    clear_now_and_after([f'{i+1}.out' for i in range(5)])
    execute_workflow(r"""
        [global]
        parameter: num = [x+1 for x in range(5)]

        [1]
        input: for_each = 'num'
        output: f'{_num}.txt'
        bash: expand = True
            touch {_output}

        [2]
        parameter: a = 1
        output: f'{_input:n}.out'
        sos_run('a' if a >=1 else 'b')

        [a,b]
        output: f'{_input:n}.out'
        bash: expand = True
            touch {_output}
        """)


def test_for_each_as_target_property_nested_list(temp_factory):
    """Test for_each option of input"""
    # test same-level for loop and parameter with nested list
    temp_factory("a.txt", "b.txt", "p1.txt", "p2.txt", "p3.txt")
    execute_workflow(
        r"""
        [0: shared=['processed']]
        files = ['a.txt', 'b.txt']
        processed = None

        input: files, for_each={'par':[(1, 2), (1, 3), (2, 3)], 'res': ['p1.txt', 'p2.txt', 'p3.txt']}
        output: _input.res

        processed = (_input.par, _input.res)
        """,
        options={"run_mode": "dryrun"},
    )

    assert env.sos_dict["step_processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]


def test_rerun_with_zap(clear_now_and_after):
    clear_now_and_after([f"zapped_example_{i}.txt.zapped" for i in range(3)])
    clear_now_and_after([f"zapped_example_{i}.bak" for i in range(3)])

    script = '''
        [step_10]
        input: for_each={'i': range(3)}
        output: f'zapped_example_{i}.txt'
        sh: expand=True
            echo "hello" > {_output}

        [step_20]
        input: group_by=1
        output: _input.with_suffix('.bak')
        sh: expand=True
            cp {_input} {_output}

        _input.zap()
    '''
    execute_workflow(script)
    execute_workflow(script)


def test_for_each_as_target_property_same_level_loop(temp_factory):
    """Test for_each option of input"""
    # test same-level for loop and parameter with nested list
    temp_factory("a.txt", "b.txt", "p1.txt", "p2.txt", "p3.txt")
    execute_workflow(
        r"""
        [0: shared=['processed']]
        files = ['a.txt', 'b.txt']
        par = [(1, 2), (1, 3), (2, 3)]
        res = ['p1.txt', 'p2.txt', 'p3.txt']

        input: files, for_each='par,res'
        output: res, group_by=1

        print([x._dict for x in step_input._groups])
        processed = (_input._par, _input._res)
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["step_processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]
