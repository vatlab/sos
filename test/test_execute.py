#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import glob
import os
import subprocess
import sys
import textwrap
import time

import pytest

from sos import execute_workflow
from sos._version import __version__
from sos.parser import SoS_Script
from sos.targets import file_target, sos_targets
from sos.utils import env
# if the test is imported under sos/test, test interacive executor
from sos.workflow_executor import Base_Executor


def multi_attempts(fn):

    def wrapper(*args, **kwargs):
        for n in range(4):
            try:
                fn(*args, **kwargs)
                break
            except Exception:
                if n > 1:
                    raise

    return wrapper


def test_command_line():
    """Test command line arguments"""
    with open("test_cl.sos", "w") as cl:
        cl.write("""\
#!/usr/bin/env sos-runner
#fileformat=SOS1.0

[L]
a =1
""")
    result = subprocess.check_output(
        "sos --version", stderr=subprocess.STDOUT, shell=True).decode()
    assert result.startswith("sos {}".format(__version__))
    assert (subprocess.call(
        "sos", stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        shell=True) == 0)
    assert (subprocess.call(
        "sos -h",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True) == 0)
    assert (subprocess.call(
        "sos run -h",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    #
    assert (subprocess.call(
        "sos run test_cl -w -W",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 1)
    assert (subprocess.call(
        "sos-runner -h",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    assert (subprocess.call(
        "sos dryrun -h",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    assert (subprocess.call(
        "sos dryrun test_cl",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    assert (subprocess.call(
        "sos dryrun test_cl.sos",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    assert (subprocess.call(
        "sos dryrun test_cl L",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    assert (subprocess.call(
        "sos-runner test_cl L",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    # script help
    assert (subprocess.call(
        "sos-runner test_cl -h",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)
    assert (subprocess.call(
        "sos convert -h",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0)


def test_func_def(temp_factory):
    """Test defintion of function that can be used by other steps"""
    temp_factory("aa.txt", "ab.txt")
    execute_workflow(
        r"""
        def myfunc(a):
            sum(range(5))
            return ['a' + x for x in a]

        [0: shared={'test':'step_input'}]
        input: myfunc(['a.txt', 'b.txt'])
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["test"] == ["aa.txt", "ab.txt"]


def test_input(temp_factory):
    """Test input specification"""
    temp_factory("test_input.txt", "test_input1.txt")
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        input: '*.txt'
        output: [x + '.res' for x in _input]
        """,
        options={"run_mode": "dryrun"},
    )
    assert file_target("test_input.txt.res").resolve() in env.sos_dict["res"]
    assert file_target("test_input1.txt.res").resolve() in env.sos_dict["res"]


def test_for_each(temp_factory):
    """Test for_each option of input"""
    temp_factory("a.txt", "b.txt", "a.pdf")
    execute_workflow(r"""
        [0: shared=['counter', 'all_names', 'all_loop']]
        files = ['a.txt', 'b.txt']
        names = ['a', 'b', 'c']
        c = ['1', '2']
        counter = 0
        all_names = ''
        all_loop = ''

        input: 'a.pdf', files, group_by='single', paired_with='names', for_each='c'

        all_names += str(_names[0]) + " "
        all_loop += str(_c) + " "

        counter = counter + 1
        """)
    assert env.sos_dict["counter"] == 6
    assert env.sos_dict["all_names"] == "a b c a b c "
    assert env.sos_dict["all_loop"] == "1 1 1 2 2 2 "


def test_for_each_same_level(temp_factory):
    """Test for_each option of input"""
    temp_factory("a.txt", "b.txt", "a.pdf")
    # test same-level for loop and parameter with nested list
    execute_workflow(
        r"""
        [0: shared=['processed']]
        files = ['a.txt', 'b.txt']
        par = [(1, 2), (1, 3), (2, 3)]
        res = ['p1.txt', 'p2.txt', 'p3.txt']
        processed = []

        input: files, for_each='par,res'
        output: res, group_by=1

        processed.append((_par, _res))
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]


def test_for_each_dataframe(temp_factory):
    """Test for_each option of input"""
    temp_factory("a.txt", "b.txt", "a.pdf")
    # test for each for pandas dataframe
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

        input: for_each='data'
        output: f"{_data['A']}_{_data['B']}_{_data['C']}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["1_2_Hello.txt", "2_4_World.txt"]


def test_for_each_dict(temp_factory):
    """Test for_each option of input"""
    # test dictionary format of for_each
    temp_factory("a.txt", "b.txt", "a.pdf")
    execute_workflow(r"""
        [0: shared=['counter', 'all_names', 'all_loop']]
        files = ['a.txt', 'b.txt']
        names = ['a', 'b', 'c']
        counter = 0
        all_names = ''
        all_loop = ''

        input: 'a.pdf', files, group_by='single', paired_with='names', for_each={'c':  ['1', '2']}

        all_names += str(_names[0]) + " "
        all_loop += c + " "

        counter = counter + 1
        """)
    assert env.sos_dict["counter"] == 6
    assert env.sos_dict["all_names"] == "a b c a b c "
    assert env.sos_dict["all_loop"] == "1 1 1 2 2 2 "


def test_for_each_multikey_dict(temp_factory):
    """Test for_each option of input"""
    # test multi-key dictionary format of for_each
    temp_factory("a.txt")
    execute_workflow(r"""
        import itertools
        [0: shared=['counter', 'all_names', 'all_loop']]
        parameter: n = [300, 100]
        parameter: p = [50, 200, 100]
        parameter: outfile = ['1', '2', '3', '4', '5', '6']
        counter = 0
        all_names = ''
        all_loop = ''
        input: 'a.txt', group_by='single', for_each={'_n,_p': [(_n,_p) for _n,_p in itertools.product(n,p) if _n > _p]}

        all_names += outfile[_index] + " "
        all_loop += '{} {} '.format(_n, _p)
        counter = counter + 1
        """)
    assert env.sos_dict["counter"] == 4
    assert env.sos_dict["all_names"] == "1 2 3 4 "
    assert env.sos_dict["all_loop"] == "300 50 300 200 300 100 100 50 "


def test_for_each_nested_list(temp_factory):
    """Test for_each option of input"""
    # test dictionary format of for_each
    temp_factory("p1.txt", "p2.txt", "p3.txt", "a.txt", "b.txt")
    # test same-level for loop and parameter with nested list
    execute_workflow(
        r"""
        [0: shared=['processed']]
        files = ['a.txt', 'b.txt']
        processed = []

        input: files, for_each={'par':[(1, 2), (1, 3), (2, 3)], 'res': ['p1.txt', 'p2.txt', 'p3.txt']}
        output: res

        processed.append((par, res))
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]


def test_for_each_pandas(temp_factory):
    """Test for_each option of input"""
    # test for each for pandas dataframe
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        input: for_each={'data': pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])}
        output: f"{data['A']}_{data['B']}_{data['C']}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["1_2_Hello.txt", "2_4_World.txt"]


def test_for_each_series():
    """Test for_each option of input"""
    # support for pands Series and Index types
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])
        input: for_each={'A': data['A']}
        output: f"a_{A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["a_1.txt", "a_2.txt"]


def test_for_each_dataframe_row():
    """Test for_each option of input"""
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])
        data.set_index('C', inplace=True)
        input: for_each={'A': data.index}
        output: f"{A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["Hello.txt", "World.txt"]


def test_for_each_each_series():
    # test for each of Series
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(0, 1, 'Ha'), (1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

        data.set_index('A', inplace=True)
        data = data.tail(2)
        input: for_each={'A': data['B']}
        output: f"{A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["2.txt", "4.txt"]


def test_for_each_iterable():
    # test iterable
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(0, 1, 'Ha'), (1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

        data.set_index('A', inplace=True)
        data = data.tail(2)
        input: for_each={'A,B': zip(data['B'],data['C'])}
        output: f"{A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["2.txt", "4.txt"]


def test_for_each_as_target_property(temp_factory):
    """Test for_each option of input"""
    temp_factory("a.txt", "b.txt", "a.pdf")
    execute_workflow(r"""
        [0: shared=['counter', 'all_names', 'all_loop']]
        files = ['a.txt', 'b.txt']
        names = ['a', 'b', 'c']
        c = ['1', '2']
        counter = 0
        all_names = ''
        all_loop = ''

        input: 'a.pdf', files, group_by='single', paired_with='names', for_each='c'

        all_names += str(_input._names) + " "
        all_loop += str(_input._c) + " "

        counter = counter + 1
        """)
    assert env.sos_dict["counter"] == 6
    assert env.sos_dict["all_names"] == "a b c a b c "
    assert env.sos_dict["all_loop"] == "1 1 1 2 2 2 "


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
        processed = []

        input: files, for_each='par,res'
        output: res, group_by=1

        print([x._dict for x in step_input._groups])
        processed.append((_input._par, _input._res))
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]


def test_for_each_as_target_property_pandas():
    # test for each for pandas dataframe
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

        input: for_each='data'
        output: f"{_input._data['A']}_{_input._data['B']}_{_input._data['C']}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["1_2_Hello.txt", "2_4_World.txt"]


def test_for_each_as_target_property_pandas_direct_value(temp_factory):
    """Test for_each option of input"""
    # test for each for pandas dataframe
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        input: for_each={'data': pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])}
        output: f"{_input.data['A']}_{_input.data['B']}_{_input.data['C']}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["1_2_Hello.txt", "2_4_World.txt"]
    #


def test_for_each_as_target_property_dict(temp_factory):
    """Test for_each option of input"""
    # test dictionary format of for_each
    temp_factory("a.txt", "b.txt", "a.pdf")
    execute_workflow(r"""
        [0: shared=['counter', 'all_names', 'all_loop']]
        files = ['a.txt', 'b.txt']
        names = ['a', 'b', 'c']
        counter = 0
        all_names = ''
        all_loop = ''

        input: 'a.pdf', files, group_by='single', paired_with='names', for_each={'c':  ['1', '2']}

        all_names += str(_input._names) + " "
        all_loop += _input.c + " "

        counter = counter + 1
        """)
    assert env.sos_dict["counter"] == 6
    assert env.sos_dict["all_names"] == "a b c a b c "
    assert env.sos_dict["all_loop"] == "1 1 1 2 2 2 "


def test_for_each_as_target_property_multidict(temp_factory):
    """Test for_each option of input"""
    # test multi-key dictionary format of for_each
    temp_factory("a.txt")
    execute_workflow(r"""
        import itertools
        [0: shared=['counter', 'all_names', 'all_loop']]
        parameter: n = [300, 100]
        parameter: p = [50, 200, 100]
        parameter: outfile = ['1', '2', '3', '4', '5', '6']
        counter = 0
        all_names = ''
        all_loop = ''
        input: 'a.txt', group_by='single', for_each={'_n,_p': [(_n,_p) for _n,_p in itertools.product(n,p) if _n > _p]}

        all_names += outfile[_index] + " "
        all_loop += '{} {} '.format(_input._n, _input._p)
        counter = counter + 1
        """)
    assert env.sos_dict["counter"] == 4
    assert env.sos_dict["all_names"] == "1 2 3 4 "
    assert env.sos_dict["all_loop"] == "300 50 300 200 300 100 100 50 "


def test_for_each_as_target_property_nested_list(temp_factory):
    """Test for_each option of input"""
    # test same-level for loop and parameter with nested list
    temp_factory("a.txt", "b.txt", "p1.txt", "p2.txt", "p3.txt")
    execute_workflow(
        r"""
        [0: shared=['processed']]
        files = ['a.txt', 'b.txt']
        processed = []

        input: files, for_each={'par':[(1, 2), (1, 3), (2, 3)], 'res': ['p1.txt', 'p2.txt', 'p3.txt']}
        output: _input.res

        processed.append((_input.par, _input.res))
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["processed"] == [
        ((1, 2), "p1.txt"),
        ((1, 3), "p2.txt"),
        ((2, 3), "p3.txt"),
    ]


def test_for_each_as_target_property_index_type(temp_factory):
    """Test for_each option of input"""
    # support for pands Series and Index types
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])
        input: for_each={'A': data['A']}
        output: f"a_{_input.A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["a_1.txt", "a_2.txt"]


def test_for_each_as_target_property_dataframe_keys(temp_factory):
    """Test for_each option of input"""
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])
        data.set_index('C', inplace=True)
        input: for_each={'A': data.index}
        output: f"{_input.A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["Hello.txt", "World.txt"]


def test_for_each_as_target_property_series(temp_factory):
    """Test for_each option of input"""
    # test for each of Series
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(0, 1, 'Ha'), (1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

        data.set_index('A', inplace=True)
        data = data.tail(2)
        input: for_each={'A': data['B']}
        output: f"{_input.A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["2.txt", "4.txt"]


def test_for_each_as_target_property_iterables(temp_factory):
    """Test for_each option of input"""
    # test iterable
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        import pandas as pd
        data = pd.DataFrame([(0, 1, 'Ha'), (1, 2, 'Hello'), (2, 4, 'World')], columns=['A', 'B', 'C'])

        data.set_index('A', inplace=True)
        data = data.tail(2)
        input: for_each={'A,B': zip(data['B'],data['C'])}
        output: f"{_input.A}.txt"
        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["res"] == ["2.txt", "4.txt"]


def test_for_each_with_full_context(temp_factory):
    # test iterable
    execute_workflow(r"""
        input: for_each = [
            {'A': 1, 'B': 2},
            {'A': 5, 'B': 7}
            ]

        print(f'A={A}, B={B}')
        """)


def test_group_by_with_no_input():
    """Test group_by with no input file"""
    execute_workflow(r"""
        [0]
        input: group_by=2
        """)


def test_paired_with(temp_factory, clear_now_and_after):
    """Test option paired_with """
    temp_factory("a.txt", "b.txt")
    clear_now_and_after("a.txt1", "b.txt2")
    #
    # string input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]

        input: files, paired_with='vars', group_by=1
        output: f"{_input}{_vars[0]}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
    #
    # list input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]
        vars2 = ['a', 'b']

        input: files, paired_with=('vars', 'vars2'), group_by=1
        output: f"{_input}{_vars[0]}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
    #
    # dict input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        input: files, paired_with={'var': [1,2], 'var2': ['a', 'b']}, group_by=1
        output: f"{_input}{var[0]}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
    #
    # test paired_with paths
    clear_now_and_after("1.txt", "2.txt")
    execute_workflow(r"""
        [1]
        input: 'a.txt', 'b.txt', paired_with = dict(bgen=paths(['1.txt', '2.txt'])), group_by=1
        output: _input.bgen
        _output.touch()
        """)
    for ofile in ["1.txt", "2.txt"]:
        assert file_target(ofile).target_exists("target")


def test_paired_with_as_target_property(temp_factory):
    """Test option paired_with with values accessed by individual target """
    temp_factory("a.txt", "b.txt")
    for ofile in ["a.txt1", "b.txt2"]:
        if file_target(ofile).exists():
            file_target(ofile).unlink()
    #
    # string input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]

        input: files, paired_with='vars', group_by=1
        output: f"{_input}{_input._vars}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # list input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]
        vars2 = ['a', 'b']

        input: files, paired_with=('vars', 'vars2'), group_by=1
        output: f"{_input}{_input._vars}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # dict input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        input: files, paired_with={'var': [1,2], 'var2': ['a', 'b']}, group_by=1
        output: f"{_input}{_input.var}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()


def test_group_with(temp_factory):
    """Test option group_with """
    temp_factory("a.txt", "b.txt")
    for ofile in ["a.txt1", "b.txt2"]:
        if file_target(ofile).exists():
            file_target(ofile).unlink()
    #
    # string input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]

        input: files, group_with='vars', group_by=1
        output: f"{_input}{_vars}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # list input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1]
        vars2 = ['a']

        input: files, group_with=('vars', 'vars2'), group_by=2
        output: f"{_input[0]}{_vars}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # dict input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        input: files, group_with={'var': [1], 'var2': ['a']}, group_by=2
        output: f"{_input[0]}{var}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()


def test_output_group_with(temp_factory):
    """Test option group_with in output statement"""
    temp_factory("a.txt", "b.txt")
    for ofile in ["a.txt1", "b.txt2"]:
        if file_target(ofile).exists():
            file_target(ofile).unlink()
    #
    # string input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]

        input: files, group_by=1
        output: f"{_input}.bak", group_with=dict(_vars=vars[_index])
        run: expand=True
            touch {_output}

        [1]
        assert(_vars == _index + 1)
        """)
    for ofile in ["a.txt.bak", "b.txt.bak"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # list input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1]
        vars2 = ['a']

        input: files, group_by=2
        output: f"{_input[0]}1", group_with=('vars', 'vars2')
        run: expand=True
            touch {_output}

        [1]
        assert(_vars == 1)
        assert(_input._vars2 == 'a')
        """)
    for ofile in ["a.txt1"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # dict input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        input: files, group_by=2
        output: f"{_input[0]}.bak",  group_with={'var': [1], 'var2': ['a']}
        run: expand=True
            touch {_output}

        [1]
        assert(var == 1)
        assert(var2 == 'a')
        """)
    for ofile in ["a.txt.bak"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()


def test_group_with_as_target_property(temp_factory):
    """Test option group_with """
    temp_factory("a.txt", "b.txt")
    for ofile in ["a.txt1", "b.txt2"]:
        if file_target(ofile).exists():
            file_target(ofile).unlink()
    #
    # string input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1, 2]

        input: files, group_with='vars', group_by=1
        output: f"{_input}{_input._vars}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1", "b.txt2"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # list input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        vars = [1]
        vars2 = ['a']

        input: files, group_with=('vars', 'vars2'), group_by=2
        output: f"{_input[0]}{_input._vars}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()
    #
    # dict input
    execute_workflow(r"""
        [0]
        files = ['a.txt', 'b.txt']
        input: files, group_with={'var': [1], 'var2': ['a']}, group_by=2
        output: f"{_input[0]}{_input.var}"
        run: expand=True
            touch {_output}
        """)
    for ofile in ["a.txt1"]:
        assert file_target(ofile).target_exists("target")
        file_target(ofile).unlink()


def test_input_pattern(temp_factory):
    """Test option pattern of step input """
    # env.verbosity = 4
    temp_factory("a-20.txt", "b-10.txt")
    execute_workflow(
        r"""
        [0: shared=['base', 'name', 'par', '_output']]

        files = ['a-20.txt', 'b-10.txt']
        input: files, pattern=['{name}-{par}.txt', '{base}.txt']
        output: ['{}-{}-{}.txt'.format(x,y,z) for x,y,z in zip(_base, _name, _par)]

        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["base"] == ["a-20", "b-10"]
    assert env.sos_dict["name"] == ["a", "b"]
    assert env.sos_dict["par"] == ["20", "10"]
    assert env.sos_dict["_output"] == ["a-20-a-20.txt", "b-10-b-10.txt"]


def test_input_pattern_as_target_property(temp_factory):
    """Test option pattern of step input """
    # env.verbosity = 4
    temp_factory("a-20.txt", "b-10.txt")
    execute_workflow(
        r"""
        [0: shared=['_output']]

        files = ['a-20.txt', 'b-10.txt']
        input: files, pattern=['{name}-{par}.txt', '{base}.txt']
        output: [f'{x._base}-{x._name}-{x._par}.txt' for x in _input]

        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["_output"] == ["a-20-a-20.txt", "b-10-b-10.txt"]


def test_output_pattern(temp_factory):
    """Test option pattern of step output"""
    # env.verbosity = 4
    temp_factory("a-20.txt", "b-10.txt")
    execute_workflow(
        r"""
        [0: shared=['base', 'name', 'par', '_output']]

        files = ['a-20.txt', 'b-10.txt']
        input: files, pattern=['{name}-{par}.txt', '{base}.txt']
        output: expand_pattern('{base}-{name}-{par}.txt'), expand_pattern('{par}.txt')

        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["base"] == ["a-20", "b-10"]
    assert env.sos_dict["name"] == ["a", "b"]
    assert env.sos_dict["par"] == ["20", "10"]
    assert env.sos_dict["_output"] == [
        "a-20-a-20.txt",
        "b-10-b-10.txt",
        "20.txt",
        "10.txt",
    ]


def test_output_from_input(temp_factory):
    """Test deriving output files from input files"""
    temp_factory("a.txt", "b.txt")
    execute_workflow(
        r"""
        [0: shared={'counter':'counter', 'step':'step_output'}]
        files = ['a.txt', 'b.txt']
        counter = 0

        input: files, group_by='single'
        output: _input[0] + '.bak'
        _output.touch()
        counter += 1
        """,
        options={"run_mode": "run"},
        config={"sig_mode": "force"},
    )
    assert env.sos_dict["counter"] == 2
    assert env.sos_dict["step"] == ["a.txt.bak", "b.txt.bak"]


def test_depends_from_input(temp_factory, clear_now_and_after):
    """Test deriving dependent files from input files"""
    temp_factory("a.txt", "b.txt")
    clear_now_and_after("a.txt.bak", "b.txt.bak")
    execute_workflow(
        r"""
        [bak: provides='{file}.bak']
        _output.touch()

        [0: shared={'counter': 'counter', 'step':'step_depends'}]
        counter = 0

        input: 'a.txt', 'b.txt', group_by='single'
        depends: _input[0] + '.bak'
        counter += 1
        """,
        options={"run_mode": "run"},
        config={"sig_mode": "force"},
    )
    assert env.sos_dict["counter"] == 2
    assert env.sos_dict["step"] == ["a.txt.bak", "b.txt.bak"]
    assert env.sos_dict["step"].groups == [["a.txt.bak"], ["b.txt.bak"]]


def test_local_namespace(temp_factory):
    """Test if steps are well separated."""
    # interctive mode behave differently
    temp_factory("a.txt")
    script = textwrap.dedent(r"""
    [1]
    a = 1

    [2]
    # this should fail because a is defined in another step
    print(a)

    """)
    with pytest.raises(Exception):
        execute_workflow(script)
    # however, alias should be sent back
    execute_workflow(
        r"""
        [1: shared={'shared': 'step_output'}]
        input: 'a.txt'
        output: 'b.txt'

        [2: shared={'tt':'step_output'}]
        print(shared)

        output: [x + '.res' for x in shared]

        """,
        options={"run_mode": "dryrun"},
    )
    assert env.sos_dict["shared"] == ["b.txt"]
    assert env.sos_dict["tt"] == ["b.txt.res"]
    #
    # this include other variables set in the step
    execute_workflow(
        r"""
        [1: shared={'shared':'c', 'd':'d'}]
        input: 'a.txt'
        output: 'b.txt'

        c = 'c.txt'
        d = 1

        [2: shared={'d': 'e'}]
        # this should fail because a is defined in another step
        print(shared)

        output: shared

        e = d + 1

        """,
        options={"run_mode": "dryrun"},
    )
    # I would like to disallow accessing variables defined
    # in other cases.
    assert env.sos_dict["shared"] == "c.txt"
    assert env.sos_dict["d"] == 2


def test_dynamic_output(temp_factory):
    """Testing dynamic output"""
    #
    temp_factory(dir="temp")
    #
    execute_workflow("""
        [10: shared={'test':'step_output'}]
        ofiles = []
        output: dynamic(ofiles)

        for i in range(4):
            ff = 'temp/something{}.html'.format(i)
            ofiles.append(ff)
            with open(ff, 'w') as h:
                h.write('a')
        """)
    assert env.sos_dict["test"] == [
        "temp/something{}.html".format(x) for x in range(4)
    ]


def test_dynamic_input(temp_factory):
    """Testing dynamic input"""
    #
    temp_factory(dir="temp")
    #
    script = SoS_Script("""
import os
[1]

from pathlib import Path
for i in range(5):
    Path(os.path.join('temp', f'test_{i}.txt')).touch()

[10: shared={'test':'step_output'}]
input: dynamic(os.path.join('temp', '*.txt')), group_by='single'
output: dynamic(os.path.join('temp', '*.txt.bak'))

run: expand=True
    touch {_input}.bak
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["test"], (
        sos_targets([
            os.path.join("temp", "test_{}.txt.bak".format(x)) for x in range(5)
        ]) ==
        f"Expecting {[os.path.join('temp', 'test_{}.txt.bak'.format(x)) for x in range(5)]} observed {env.sos_dict['test']}"
    )
    # this time we use th existing signature
    Base_Executor(wf).run()
    assert env.sos_dict["test"], (
        sos_targets([
            os.path.join("temp", "test_{}.txt.bak".format(x)) for x in range(5)
        ]) ==
        f"Expecting {[os.path.join('temp', 'test_{}.txt.bak'.format(x)) for x in range(5)]} observed {env.sos_dict['test']}"
    )


def test_use_of_runmode(temp_factory):
    """Test the use of run_mode variable in SoS script"""
    #
    temp_factory(dir="temp")
    env.config["sig_mode"] = "ignore"
    execute_workflow("""
        [1: shared={'res': '_output'}]
        import random
        for i in range(3):
            with open(f"temp/test_{random.randint(1, 100000)}.txt", 'w') as res:
                res.write(str(i))

        """)
    # we should have 9 files
    files = glob.glob(os.path.join("temp", "*.txt"))
    assert len(files) == 3


def test_action_before_input():
    """Testing the execution of actions before input directive
    (variables such as _index should be made available)."""
    execute_workflow(
        """
        [0]
        run('echo "A"')
        input:
        """,
        options={"run_mode": "dryrun"},
    )


def test_duplicate_io_files(temp_factory):
    """Test interpretation of duplicate input/output/depends"""
    temp_factory(dir="temp")
    # Test duplicate input
    os.system("touch temp/1.txt")
    execute_workflow("""
        [1]
        input: ['temp/1.txt' for x in range(5)]
        run: expand=True
        touch temp/{len(_input)}.input
                """)
    assert os.path.isfile("temp/5.input")
    # Test duplicate output
    script = textwrap.dedent("""
    [1]
    output: ['temp/2.txt' for x in range(5)]
    run: expand=True
    touch temp/2.txt
    touch temp/{len(_output)}.output
    """)
    with pytest.raises(Exception):
        execute_workflow(script)
    # Test duplicate depends
    script = textwrap.dedent("""
    [1]
    input: 'temp/1.txt'
    depends: ['temp/2.txt' for x in range(5)]
    output: 'temp/3.txt'
    run: expand=True
    touch temp/3.txt
    touch temp/{len(_depends)}.depends
    """)
    with pytest.raises(Exception):
        execute_workflow(script)


def test_output_in_loop(temp_factory):
    """Test behavior of {_output} when used in loop"""
    temp_factory(dir="temp")
    env.config["sig_mode"] = "ignore"
    execute_workflow("""
        [default]
        s = [x for x in range(5)]
        output_files = ['temp/{}.txt'.format(x) for x in range(5)]
        input: for_each = ['s'], concurrent=False
        output: output_files[_index]
        run: active = 0
        rm -f temp/out.log
        run: expand=True
        echo {step_output} >> temp/out.log
        touch {step_output}
                """)
    # output should have 1, 2, 3, 4, 5, respectively
    with open("temp/out.log") as out:
        assert len(out.read().split()) == 5
    #
    temp_factory(dir="temp")
    env.config["sig_mode"] = "ignore"
    execute_workflow("""
    [default]
    s = [x for x in range(5)]
    output_files = ['temp/{}.txt'.format(x) for x in range(5)]
    input: for_each = ['s'], concurrent=False
    output: output_files[_index]
    run: active = 0
    rm -f temp/out.log
    run: expand=True
    echo {step_output} >> temp/out.log
    touch {step_output}
    """)
    with open("temp/out.log") as out:
        assert len(out.read().split()) == 5


@multi_attempts
def test_execution_lock():
    """Test execution lock of two processes"""
    with open("lock.sos", "w") as lock:
        lock.write(r"""
import time
[A_1]
output: 'a.txt'
with open('a.txt', 'w') as txt:
txt.write('A1\n')

# A1 and A2 are independent
[A_2]
input: None
output: 'b.txt'
with open('b.txt', 'w') as txt:
txt.write('A2\n')
    """)
    ret1 = subprocess.Popen("sos run lock -j1", shell=True)
    ret2 = subprocess.Popen("sos run lock -j1", shell=True)
    ret1.wait()
    ret2.wait()
    # two processes execute A_1 and A_2 separately, usually
    # takes less than 5 seconds
    file_target("lock.sos").unlink()


def test_removed_intermediate_files(clear_now_and_after):
    """Test behavior of workflow with removed internediate files"""
    clear_now_and_after("a.txt", "aa.txt")
    script = SoS_Script("""
[10]
output: 'a.txt'
run:
echo "a" > a.txt

[20]
output: 'aa.txt'
run: expand=True
cat {_input} > {_output}
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert file_target("aa.txt").target_exists()
    # rerun should be faster
    Base_Executor(wf).run()
    # if we remove the middle result, it should not matter
    os.remove("a.txt")
    Base_Executor(wf).run()
    #
    # if we remove the final result, it will be rebuilt
    os.remove("aa.txt")
    Base_Executor(wf).run()
    #
    # now we request the generation of target
    file_target("a.txt").unlink()
    file_target("aa.txt").unlink()
    Base_Executor(wf).run()
    #
    file_target("a.txt").unlink()
    file_target("aa.txt").unlink()


def test_stopped_output():
    """test output with stopped step"""
    for file in ["{}.txt".format(a) for a in range(10)]:
        if file_target(file).exists():
            file_target(file).unlink()
    execute_workflow("""
        [test_1]
        input: for_each={'a': range(10)}
        output: f"{a}.txt"

        stop_if(a % 2 == 0, no_output=True)
        run: expand=True
            touch {_output}

        [test_2]
        assert(len(step_input) == 5)
        """)
    for idx in range(10):
        if idx % 2 == 0:
            assert not file_target("{}.txt".format(idx)).target_exists()
        else:
            assert file_target("{}.txt".format(idx)).target_exists()
            file_target(f"{idx}.txt").unlink()


def test_allow_error(clear_now_and_after):
    """Test option allow error"""
    clear_now_and_after("a.txt")
    execute_workflow("""
        [test]
        run:  allow_error=True
            something_wrong

        run:
            touch a.txt
        """)
    assert file_target("a.txt").target_exists()


def test_concurrent_worker(clear_now_and_after):
    """Test the starting of multiple workers #493 """
    with open("test_script.sos", "w") as script:
        script.write("""
[10]
input: for_each={'i': range(1)}

[20]
input: for_each={'i': range(2)}
""")
    subprocess.call("sos run test_script", shell=True)
    clear_now_and_after("test_script.sos")


def test_depends_caused_dependency():
    # test for #674
    for tfile in ("1.txt", "2.txt", "3.txt"):
        if file_target(tfile).exists():
            file_target(tfile).unlink()
    execute_workflow("""
        [1: shared = {'dfile':'_output'}]
        output: '1.txt'
        run:
            echo 1 > 1.txt

        [2: shared = {'ifile':'_output'}]
        output: '2.txt'
        run: expand=True
            echo {_input} > 2.txt

        [3]
        depends: sos_variable('ifile'), sos_variable('dfile'), ifile
        input: dfile
        output: '3.txt'
        run: expand=True
            cat {_input} > {_output}
        """)
    for tfile in ("1.txt", "2.txt", "3.txt"):
        assert file_target(tfile).target_exists()
        if file_target(tfile).exists():
            file_target(tfile).unlink()


def test_concurrent_input_option(temp_factory):
    """Test input option"""
    temp_factory("1.txt", "2.txt")
    execute_workflow("""
        [1]
        n =[str(x) for x in range(2)]
        input: [f'{x+1}.txt' for x in range(2)], paired_with = 'n', concurrent = True
        run: expand = True
        echo {_n} {_input}
        """)


def test_non_existent_dependent_target():
    """Test non existent dependent targets"""
    with pytest.raises(Exception):
        execute_workflow(r"""
            [1]

            [2]
            depends: sos_step('wrong')
            """)
    with pytest.raises(Exception):
        execute_workflow(r"""
            [1]

            [2]
            depends: 'non-existent.txt'
            """)


def test_output_report(clear_now_and_after):
    """Test generation of report"""
    clear_now_and_after("report.html")
    script = SoS_Script(r"""
[1: shared = {'dfile':'_output'}]
output: '1.txt'
run:
echo 1 > 1.txt

[2: shared = {'ifile':'_output'}]
output: '2.txt'
run: expand=True
echo {_input} > 2.txt

[3]
depends: ifile, sos_variable('ifile'), sos_variable('dfile')
input: dfile
output: '3.txt'
run: expand=True
cat {_input} > {_output}
""")
    env.config["output_report"] = "report.html"
    wf = script.workflow()
    Base_Executor(wf).run()
    assert os.path.isfile("report.html")


@pytest.mark.skipif(
    sys.platform == "win32", reason="Graphviz not available under windows")
def test_output_report_with_dag(clear_now_and_after):
    # test dag
    clear_now_and_after("report.html")
    execute_workflow(
        r"""
        [1: shared = {'dfile':'_output'}]
        output: '1.txt'
        run:
            echo 1 > 1.txt

        [2: shared = {'ifile':'_output'}]
        output: '2.txt'
        run: expand=True
            echo {_input} > 2.txt

        [3]
        depends: ifile, sos_variable('ifile'), sos_variable('dfile')
        input: dfile
        output: '4.txt'
        run: expand=True
            cat {_input} > {_output}
        """,
        options={
            "output_report": "report.html",
            "output_dag": "report.dag"
        },
    )
    with open("report.html") as rep:
        content = rep.read()
    assert "Execution DAG" in content


def test_sos_step_with_output():
    """Test checking output of sos_step #981"""
    execute_workflow("""
        [step]
        output: 'a'
        sh:
        touch a

        [default]
        depends: sos_step('step')
        """)


def test_multi_sos_step():
    """Test matching 'a_1', 'a_2' etc with sos_step('a')"""
    for file in ("a_1", "a_2"):
        if file_target(file).exists():
            file_target(file).unlink()
    res = execute_workflow("""
    [a_b_1]
    output: "a_1"
    sh:
    echo whatever > a_1

    [a_b_2]
    output: "a_2"
    sh: expand=True
    cp {_input} {_output}

    [default]
    depends: sos_step('a_b')
    """)
    assert res["__completed__"]["__step_completed__"] == 3
    assert os.path.isfile("a_1")
    assert os.path.isfile("a_2")
    with open("a_1") as a1, open("a_2") as a2:
        assert a1.read() == a2.read()


def test_depends_auxi_and_forward():
    """Test depends on auxiliary, which then depends on a forward-workflow #983"""
    for f in ("a_1", "a_2"):
        if file_target(f).exists():
            file_target(f).unlink()
    res = execute_workflow("""

    [hg_1]
    output: 'a_1'
    sh:
    echo "something" > a_1

    [hg_2]

    [star: provides = "a_2"]
    depends: sos_step('hg')
    sh:
    cp  a_1 a_2

    [default]
    depends: "a_2"
            """)
    assert res["__completed__"]["__step_completed__"] == 4
    assert os.path.isfile("a_1")
    assert os.path.isfile("a_2")
    with open("a_1") as a1, open("a_2") as a2:
        assert a1.read() == a2.read()


def test_depends_auxi_and_single_step_forward():
    """Test depends on auxiliary, which then depends on a single-step forward-workflow"""
    for f in ("a_1", "a_2"):
        if file_target(f).exists():
            file_target(f).unlink()
    res = execute_workflow("""

    [hg_1]
    output: 'a_1'
    sh:
    echo "something" > a_1

    [star: provides = "a_2"]
    depends: sos_step('hg')
    sh:
    cp  a_1 a_2

    [default]
    depends: "a_2"
    """)
    assert res["__completed__"]["__step_completed__"] == 3
    assert os.path.isfile("a_1")
    assert os.path.isfile("a_2")
    with open("a_1") as a1, open("a_2") as a2:
        assert a1.read() == a2.read()


def test_dryrun_placeholder():
    """Test the creation and removal of placeholder files in dryrun mode"""
    if file_target("1.txt").exists():
        file_target("1.txt").unlink()
    execute_workflow(
        """
        a = '1.txt'

        [out: provides=a]
        output: a
        run: expand = True
        touch {a}

        [1]
        depends: a
        """,
        options={"run_mode": "dryrun"},
    )
    # but the file would be removed afterwards
    assert not os.path.isfile("1.txt")


def test_dryrun_in_sos_run(temp_factory):
    """Test dryrun mode with sos_run #1007"""
    temp_factory("1.txt")
    script = SoS_Script(
        textwrap.dedent("""
        [remove]
        run:
        rm 1.txt

        [default]
        sos_run('remove')
        """))
    wf = script.workflow()
    Base_Executor(wf).run(mode="dryrun")
    assert os.path.isfile("1.txt")
    Base_Executor(wf).run(mode="run")
    assert not os.path.isfile("1.txt")


def test_concurrent_with_dynamic_output(clear_now_and_after):
    """Test concurrent steps with dynamic output"""
    douts = glob.glob("*.dout")
    for dout in douts:
        clear_now_and_after(dout)
    execute_workflow("""
        input: for_each={'i': range(3)}, concurrent=True
        output: dynamic('*.dout')
        import random
        path(f'{random.randint(0, 1000000)}.dout').touch()
        """)
    douts = glob.glob("*.dout")
    assert len(douts) == 3


def test_group_by_with_emtpy_input():
    """ Test option group by with empty input #1044"""
    execute_workflow("""
        [1]
        input: group_by=1
        print(_input)
        """)


def test_removal_of_output_from_failed_step(clear_now_and_after):
    """Test the removal of output files if a step fails #1055"""
    clear_now_and_after("failed.csv", "result.csv")
    script = SoS_Script("""
[sub: provides='{file}.csv']
sh: expand=True
touch {_output}
eco "something wrong"

[step]
depends: 'failed.csv'
path('result.csv').touch()
""")
    wf = script.workflow()
    with pytest.raises(Exception):
        Base_Executor(wf).run()
    # rerun should still raise
    with pytest.raises(Exception):
        Base_Executor(wf).run()

    assert not os.path.isfile("failed.csv")
    assert not os.path.isfile("result.csv")


def test_depends_to_concurrent_substep():
    """Testing forward style example"""
    # sos_variable('data') is passed to step [2]
    # but it is not passed to concurrent substep because
    # the variable is not used in the substep. This test
    # should fail at least under windows
    execute_workflow("""
        [1: shared={'data': 'step_output'}]
        output: 'a.txt'
        _output.touch()

        [2]
        depends: sos_variable('data')
        input: for_each={'i': range(2)}, group_by=1, concurrent=True
        print(1)
        """)


def test_pass_of_target_source():
    """Test passing of source information from step_output"""
    execute_workflow("""
        [1]
        output: 'a.txt'
        _output.touch()

        [2]
        assert step_input.labels == ['1']
        """)
    #
    execute_workflow("""
        [1]
        input: for_each={'i': range(2)}
        output: 'a.txt', 'b.txt', group_by=1
        _output.touch()

        [2]
        assert step_input.labels == ['1', '1']
        """)
    #
    file_target("c.txt").touch()
    execute_workflow("""
        [1]
        input: for_each={'i': range(2)}
        output: 'a.txt', 'b.txt', group_by=1
        _output.touch()

        [2]
        input: 'c.txt'
        assert step_input.labels == ['2']
        """)


def test_rerun_with_zap():
    execute_workflow("""
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
        """)
    execute_workflow("""
        [step_10]
        input: for_each={'i': range(3)}
        output: f'zapped_example_{i}.txt'
        sh: expand=True
            echo "hello" > {_output}

        [step_20]
        input: group_by=1
        output: _input.with_suffix('.bak')
        print(_output)
        sh: expand=True
            cp {_input} {_output}

        _input.zap()
        """)
    for i in range(3):
        os.remove(f"zapped_example_{i}.txt.zapped")


def test_return_output_in_step_output():
    """Testing the return of _output as groups of step_output"""
    execute_workflow("""\
        [1]
        input: for_each=dict(i=range(5))
        output: f'a_{i}.txt'
        _output.touch()
        assert(_input.i == i)

        [2]
        assert(len(step_input.groups) == 5)
        assert(len(step_input) == 5)
        assert(step_input.groups[0] == 'a_0.txt')
        assert(step_input.groups[4] == 'a_4.txt')
        #assert(_input.i == _index)
        """)
    #
    # test accumulation of named output
    execute_workflow("""\
        [1]
        input: for_each=dict(i=range(5))
        output: a=f'a_{i}.txt', b=f'b_{i}.txt'
        _output.touch()

        [2]
        assert(len(step_input.groups) == 5)
        assert(len(step_input) == 10)
        assert(step_input.groups[0] == ['a_0.txt', 'b_0.txt'])
        assert(step_input.groups[0].labels == ['a', 'b'])
        assert(step_input.groups[4] == ['a_4.txt', 'b_4.txt'])
        assert(step_input.groups[4].labels == ['a', 'b'])
        """)


def test_output_from():
    """Testing output_from input function"""
    script = SoS_Script("""\
[A]
input: for_each=dict(i=range(5))
output: f'a_{i}.txt'
_output.touch()

[A1]
input: for_each=dict(i=range(4))
output: aa=f'a_{i}.txt'
_output.touch()

[B]
input: output_from('A')
assert(len(step_input.groups) == 5)
assert(len(step_input) == 5)
assert(step_input.labels == ['A']*5)
assert(step_input.groups[0] == 'a_0.txt')
assert(step_input.groups[4] == 'a_4.txt')

[C]
input: K=output_from('A')
assert(len(step_input.groups) == 5)
assert(step_input.labels == ['K']*5)

[D]
input: K=output_from('A', group_by='all')
assert(len(step_input) == 5)
assert(len(step_input.groups) == 1)
assert(step_input.labels == ['K']*5)

[E]
input: output_from('A1', group_by='all')
assert(len(step_input) == 4)
assert(len(step_input.groups) == 1)
assert(step_input.labels == ['aa']*4)

[F]
input: K=output_from('A1', group_by='all')['aa']
assert(len(step_input) == 4)
assert(len(step_input.groups) == 1)
assert(step_input.labels == ['K']*4)

[G_0]
input: for_each=dict(i=range(4))
output: f'g_{i}.txt'
_output.touch()

[G_100]
input: K=output_from(-1, group_by=2)
assert(len(step_input) == 4)
assert(len(step_input.groups) == 2)
assert(step_input.labels == ['K']*4)

[H_0]
input: for_each=dict(i=range(4))
output: f'g_{i}.txt'
_output.touch()

[H_100]
input: K=output_from([-1, 'A1'], group_by=2)
assert(len(step_input) == 8)
assert(len(step_input.groups) == 4)
assert(step_input.labels == ['K']*8)

""")
    for wf in ("B", "C", "D", "E", "F", "G", "H"):
        wf = script.workflow(wf)
        Base_Executor(wf).run()


def test_named_output1336():
    "Test issue 1336"
    execute_workflow("""
        import time

        [10]
        output: bak='A.bak'
        time.sleep(2)
        path('A.bak').touch()


        [20]
        input: named_output('bak')
        output: 'B.txt'
        _output.touch()
        """)


def test_set_variables_to_output():
    """Test assigning variables to _output"""
    execute_workflow("""\
        [10]
        output: 'a.txt'
        _output[0].set(tvar=1)
        _output.set(gvar=2)
        _output.touch()

        [20]
        assert(gvar == 2)
        assert(_input.gvar == 2)
        assert(_input.tvar == 1)
        assert(_input[0].tvar == 1)
        """)
    # if there are substeps
    execute_workflow("""\
        [10]
        input: for_each=dict(i=range(4))
        output: f'a_{i}.txt'
        _output[0].set(tvar=i)
        _output.set(gvar=i)
        _output.touch()

        [20]
        assert(gvar == _index)
        assert(_input.gvar == _index)
        assert(_input.tvar == _index)
        assert(_input[0].tvar == _index)
        """)


def test_step_id_vars():
    """Test variables in a step"""
    execute_workflow("""
        [nested]
        print(f'Workflow {workflow_id}: step name={step_name}')
        print(f'Workflow {workflow_id}: step id={step_id}')
        print(f'Workflow {workflow_id}: workflow id={workflow_id}')
        print(f'Workflow {workflow_id}: master id={master_id}')
        assert step_name == 'nested'
        assert workflow_id != master_id

        [default]
        print(f'Workflow {workflow_id}: step name={step_name}')
        print(f'Workflow {workflow_id}: step id={step_id}')
        print(f'Workflow {workflow_id}: workflow id={workflow_id}')
        print(f'Workflow {workflow_id}: master id={master_id}')
        assert step_name == 'default'
        assert workflow_id == master_id
        sos_run('nested')
        """)


def test_reexecution_of_dynamic_depends(clear_now_and_after):
    """Testing the rerun of steps to verify dependency"""
    clear_now_and_after("a.bam", "a.bam.bai")
    script = """
        [BAI: provides='{filename}.bam.bai']
        _output.touch()

        [BAM]
        output: 'a.bam'
        _output.touch()

        [default]
        input: 'a.bam'
        depends: _input.with_suffix('.bam.bai')
        """
    execute_workflow(script)
    # if we run again, because depends, the step will be re-checked
    os.remove("a.bam")
    res = execute_workflow(script)
    assert res["__completed__"]["__step_completed__"] == 2
    assert res["__completed__"]["__step_skipped__"] == 0
    #
    os.remove("a.bam")
    res = execute_workflow(script, options={"trace_existing": True})
    assert res["__completed__"]["__step_completed__"] == 2
    assert res["__completed__"]["__step_skipped__"] == 1


def test_traced_function(clear_now_and_after):
    clear_now_and_after("a.bam", "a.bam.bai")
    execute_workflow("""
        [BAI: provides='{filename}.bam.bai']
        _output.touch()

        [BAM]
        output: 'a.bam'
        _output.touch()

        [default]
        input: 'a.bam'
        depends: traced(_input.with_suffix('.bam.bai'))
    """)
    # if we run again, because depends, the step will be re-checked
    os.remove("a.bam")
    res = execute_workflow("""
        [BAI: provides='{filename}.bam.bai']
        _output.touch()

        [BAM]
        output: 'a.bam'
        _output.touch()

        [default]
        input: 'a.bam'
        depends: traced(_input.with_suffix('.bam.bai'))
        """)
    assert res["__completed__"]["__step_completed__"] == 2
    assert res["__completed__"]["__step_skipped__"] == 1


def test_error_handling_of_step():
    # test fail_if of killing another running substep
    def cleanup():
        for step in (10, 11):
            if os.path.isfile(f"{step}.txt"):
                os.remove(f"{step}.txt")

    script = textwrap.dedent(r"""
    import time

    [10]
    output: '10.txt'
    time.sleep(8)
    _output.touch()

    [11]
    output: '11.txt'
    _output.touch()

    [20]
    input: None
    time.sleep(2)
    fail_if(True)
    """)
    #
    # default mode
    #
    cleanup()
    st = time.time()
    with pytest.raises(Exception):
        execute_workflow(script)
    assert (time.time() - st >=
            8), "Test test should fail only after step 10 is completed"
    assert os.path.isfile("10.txt")
    assert os.path.isfile("11.txt")
    #
    # ignore mode
    #
    cleanup()
    #
    st = time.time()
    execute_workflow(script, options={"error_mode": "ignore"})
    assert (time.time() - st >=
            8), "Test test should fail only after step 10 is completed"
    assert os.path.isfile("10.txt")
    assert os.path.isfile("11.txt")
    #
    # abort mode
    #
    cleanup()
    #
    with pytest.raises(Exception):
        execute_workflow(script, options={"error_mode": "abort"})
    assert not os.path.isfile("10.txt")
    assert not os.path.isfile("11.txt")


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


def test_error_handling_of_concurrent_substeps():

    def cleanup():
        for i in range(200):
            if os.path.isfile(f"test_{i}.txt"):
                os.remove(f"test_{i}.txt")

    script = textwrap.dedent(r"""
    import time

    [10]
    input: for_each=dict(i=range(200))
    output: f'test_{i}.txt'

    _output.touch()

    fail_if(_index == 5, 'fail at 5')
    fail_if(_index == 10, 'fail at 10')
    """)
    #
    # default mode
    #
    cleanup()
    with pytest.raises(Exception):
        execute_workflow(script)
    for i in (5, 10):
        assert not os.path.isfile(f"test_{i}.txt")
    for i in range(190, 200):
        assert os.path.isfile(f"test_{i}.txt")


def test_error_handling_of_tasks(clear_now_and_after):
    clear_now_and_after(
        [f"test_{i}.txt" for i in (0, 1, 2, 30, 31)],
        [f"test_{i}.bak" for i in (0, 1, 2, 30, 31)],
    )

    script = r"""
import time

[10]
input: for_each=dict(i=range(3)), concurrent=False
output: f'test_{i}.txt'

task: queue='localhost'
fail_if(i == 1, 'fail at 1')
_output.touch()

[20]
output: _input.with_suffix('.bak')
_output.touch()

[30]
input: None
output: 'test_30.txt'
time.sleep(18)
_output.touch()

[31]
output: 'test_31.txt'
_output.touch()
    """
    #
    # default mode
    #

    with pytest.raises(Exception):
        execute_workflow(script)

    for i in (0, 2, 30, 31):
        assert os.path.isfile(f"test_{i}.txt")
        assert not os.path.isfile(f"test_{i}.bak")
    assert not os.path.isfile("test_1.txt")
    #
    #  ignore mode
    #
    clear_now_and_after(
        [f"test_{i}.txt" for i in (0, 1, 2, 30, 31)],
        [f"test_{i}.bak" for i in (0, 1, 2, 30, 31)],
    )

    execute_workflow(script, options={"error_mode": "ignore"})
    for i in (0, 2, 30, 31):
        assert os.path.isfile(f"test_{i}.txt")
    for i in (0, 2):
        assert not os.path.isfile(f"test_{i}.bak")
    assert not os.path.isfile("test_1.txt")
    #
    # abort mode
    #
    clear_now_and_after(
        [f"test_{i}.txt" for i in (0, 1, 2, 30, 31)],
        [f"test_{i}.bak" for i in (0, 1, 2, 30, 31)],
    )

    with pytest.raises(Exception):
        execute_workflow(script, options={"error_mode": "abort"})
    for i in (0, 2):
        assert os.path.isfile(f"test_{i}.txt")
    for i in (30, 31):
        assert not os.path.isfile(f"test_{i}.txt")
    for i in (0, 2):
        assert not os.path.isfile(f"test_{i}.bak")
    assert not os.path.isfile("test_1.txt")


def test_stmt_before_input(clear_now_and_after):
    """Bug #1270, if there is any statement before input, the step will be undetermined"""

    clear_now_and_after("test_1270.txt", "test_1270.out")
    execute_workflow(r"""
        [10]
        with open('test_1270.txt', 'w') as t1:
            t1.write('something')

        input: 'test_1270.txt'
        output: 'test_1270.out'
        with open(_input, 'r') as ifile, open(_output, 'w') as ofile:
            ofile.write(ifile.read())

        """)

    assert os.path.isfile("test_1270.txt")
    assert os.path.isfile("test_1270.out")


def test_param_with_step_no_statement():
    # 1375
    execute_workflow(
        r"""
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
        task: trunk_workers = 1, trunk_size = 1, walltime = '3m', mem = '1G', cores = 1
        bash: expand = True
        echo {a} > {_output}
            """,
        options={"default_queue": "localhost"},
    )


def test_parallel_nestedworkflow():
    # 1375
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


def test_concurrent_running_tasks(script_factory):
    """Test two sos instances running the same task"""
    script = script_factory("""
        [1]
        task:
        import time
        time.sleep(5)
        """)

    ret1 = subprocess.Popen(
        ["sos", "run", script, "-s", "force", "-q", "localhost"])
    ret2 = subprocess.Popen(
        ["sos", "run", script, "-s", "force", "-q", "localhost"])
    ret1.wait()
    ret2.wait()
    assert ret1.returncode == 0
    assert ret2.returncode == 0


def test_interpolation():
    """Test string interpolation during execution"""
    execute_workflow(r"""
        [0: shared='res']
        res = ''
        b = 200
        res += f"{b}"
    """)
    assert env.sos_dict["res"] == "200"


def test_interpolation_1():
    """Test string interpolation during execution"""
    execute_workflow(r"""
        [0: shared='res']
        res = ''
        for b in range(5):
            res += f"{b}"
    """)
    assert env.sos_dict["res"] == "01234"


def test_interpolation_2(temp_factory):
    """Test string interpolation during execution"""
    temp_factory("a_1.txt", "b_2.txt", "c_2.txt")

    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
        output: [f'{x}_{y}_processed.txt' for x,y in zip(name, model)]
        """,
        options={"run_mode": "dryrun"},
    )

    assert env.sos_dict["res"] == [
        "a_1_processed.txt",
        "b_2_processed.txt",
        "c_2_processed.txt",
    ]


def test_interpolation_3(temp_factory):
    """Test string interpolation during execution"""
    temp_factory("a_1.txt", "b_2.txt", "c_2.txt")

    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
        output: [f"{x}_{y}_process.txt" for x,y in zip(name, model)]
        """,
        options={"run_mode": "dryrun"},
    )

    assert env.sos_dict["res"] == [
        "a_1_process.txt",
        "b_2_process.txt",
        "c_2_process.txt",
    ]


def test_interpolation_4(temp_factory):
    """Test string interpolation during execution"""
    temp_factory("a_1.txt", "b_2.txt", "c_2.txt")
    execute_workflow(
        r"""
        [0: shared={'res':'step_output'}]
        def add_a(x):
            return ['a'+_x for _x in x]

        input: 'a_1.txt', 'b_2.txt', 'c_2.txt', pattern='{name}_{model}.txt'
        output: add_a([f"{x}_{y}_process.txt" for x,y in zip(name, model)])

        """,
        options={"run_mode": "dryrun"},
    )

    assert env.sos_dict["res"] == [
        "aa_1_process.txt",
        "ab_2_process.txt",
        "ac_2_process.txt",
    ]


def test_named_output():
    """Testing named_output input function"""
    for wf in ("B", "C", "D"):
        execute_workflow(
            """
            [A]
            input: for_each=dict(i=range(4))
            output: aa=f'a_{i}.txt', bb=f'b_{i}.txt'
            _output.touch()

            [B]
            input: named_output('aa')
            assert(len(step_input.groups) == 4)
            assert(len(step_input) == 4)
            assert(step_input.labels == ['aa']*4)
            assert(step_input.groups[0] == 'a_0.txt')
            assert(step_input.groups[3] == 'a_3.txt')

            [C]
            input: K=named_output('bb')
            assert(len(step_input.groups) == 4)
            assert(len(step_input) == 4)
            assert(step_input.labels == ['K']*4)
            assert(step_input.groups[0] == 'b_0.txt')
            assert(step_input.groups[3] == 'b_3.txt')

            [D]
            input: K=named_output('bb', group_by=2)
            assert(len(step_input.groups) == 2)
            assert(len(step_input) == 4)
            assert(step_input.labels == ['K']*4)
            assert(step_input.groups[1] == ['b_2.txt', 'b_3.txt'])
            """,
            workflow=wf,
        )


def test_auto_provide():
    """Testing steps to provide plain output"""
    execute_workflow("""
        [global]

        a = 'a.txt'

        [b]
        output: a
        _output.touch()

        [default]

        input: a
        output: f"{file_target(a):n}.out"

        _output.touch()
        """)


def test_error_handling_of_missing_input(clear_now_and_after):
    # finishing another step if currently missing
    clear_now_and_after("11.txt", "22.txt")

    st = time.time()

    script = """
        import time

        [10]
        time.sleep(8)

        [11]
        output: '11.txt'
        _output.touch()

        [20]
        input: None
        time.sleep(2)

        [21]
        input: 'no_existent.txt'

        [22]
        output: '22.txt'
        _output.touch()
        """
    with pytest.raises(Exception):
        execute_workflow(script)

    assert (time.time() - st >=
            8), "Test test should fail only after step 10 is completed"
    assert os.path.isfile("11.txt")

    clear_now_and_after("11.txt", "22.txt")

    execute_workflow(script, options={"error_mode": "ignore"})

    assert (time.time() - st >=
            8), "Test test should fail only after step 10 is completed"
    assert not os.path.isfile("22.txt")


def test_assignment_after_input(temp_factory):
    """Testing assignment after input should be usable inside step process."""
    #
    temp_factory(dir="temp")

    execute_workflow(
        """
        [1]
        rep = range(5)
        input:  for_each='rep'
        output: f"temp/{_rep}.txt"

        # ff should change and be usable inside run
        ff = f"{_rep}.txt"
        run: expand=True
        echo {ff}
        touch temp/{ff}
        """,
        options={"sig_mode": "ignore"},
    )


def test_1379(clear_now_and_after):
    clear_now_and_after("a.txt", "b.txt")

    execute_workflow("""
        [multi: provides=['a.txt', 'b.txt']]

        out = 'a.txt', 'b.txt'
        output: out

        import time
        time.sleep(2)
        _output.touch()


        [step_1]
        input: 'a.txt'

        [step_2]
        input: 'b.txt'
        """)


def test_remove_empty_groups_default():
    """Test remove of empty groups"""
    # case 1, default output
    execute_workflow("""
        [10]
        input: for_each=dict(i=range(4))
        output: f'a_{i}.txt'
        _output.touch()
        skip_if(i==2)

        [20]
        assert len(step_input.groups) == 4
        """)


def test_remove_empty_groups_false():
    """Test remove of empty groups"""
    # case 2, use default remove_empty_groups=False
    execute_workflow("""
        [A]
        input: for_each=dict(i=range(4))
        output: f'a_{i}.txt'
        _output.touch()
        skip_if(i==2)

        [default]
        input: output_from('A')
        assert len(step_input.groups) == 4
        """)


def test_remove_empty_groups_true():
    """Test remove of empty groups"""
    # case 3, use remove_empty_groups=True
    execute_workflow("""
        [A]
        input: for_each=dict(i=range(4))
        output: f'a_{i}.txt'
        _output.touch()
        skip_if(i==2)

        [default]
        input: output_from('A', remove_empty_groups=True)
        assert len(step_input.groups) == 3
        """)


def test_remove_empty_groups_named():
    """Test remove of empty groups"""
    # case 4, use named_output
    execute_workflow("""
        [A]
        input: for_each=dict(i=range(4))
        output: A=f'a_{i}.txt'
        _output.touch()
        skip_if(i==2)

        [default]
        input: named_output('A')
        assert len(step_input.groups) == 4
        """)


def test_remove_empty_groups_empty_named():
    """Test remove of empty groups"""
    # case 5, use named_output
    execute_workflow("""
        [A]
        input: for_each=dict(i=range(4))
        output: A=f'a_{i}.txt'
        _output.touch()
        skip_if(i==2)

        [default]
        input: named_output('A', remove_empty_groups=True)
        assert len(step_input.groups) == 3
        """)


def test_multi_depends(clear_now_and_after, temp_factory):
    """Test a step with multiple depdendend steps"""

    clear_now_and_after(
        "dbsnp.vcf",
        "hg19.fa",
        "f1.fastq",
        "f2.fastq",
        "f1.bam",
        "f2.bam",
        "f1.bam.idx",
        "f2.bam.idx",
    )
    temp_factory("f1.fastq", "f2.fastq")

    execute_workflow(
        """
        import time

        [refseq: provides='hg19.fa']
        time.sleep(1)
        _output.touch()

        [dbsnp: provides='dbsnp.vcf']
        _output.touch()

        [align_10]
        depends: 'hg19.fa'
        input: 'f1.fastq', 'f2.fastq', group_by=1, concurrent=True
        output: _input.with_suffix('.bam')
        _output.touch()

        [align_20]
        input: group_by=1, concurrent=True
        output: _input.with_suffix('.bam.idx')
        _output.touch()

        [call_10]
        depends: 'dbsnp.vcf', 'hg19.fa'

        [call_20]
        """,
        workflow="align+call",
    )

    for file in (
            "dbsnp.vcf",
            "hg19.fa",
            "f1.bam",
            "f2.bam",
            "f1.bam.idx",
            "f2.bam.idx",
    ):
        assert os.path.isfile(file)


def test_execute_ipynb():
    """Test extracting and executing workflow from .ipynb files"""
    script = SoS_Script(filename="sample_workflow.ipynb")
    wf = script.workflow()
    Base_Executor(wf).run()


@pytest.mark.skipif(
    True,
    reason="Skip test because travis fails on this test for unknown reason, also due to a bug in psutil under windows",
)
def test_kill_worker(script_factory):
    """Test if the workflow can error out after a worker is killed"""
    import time

    import psutil

    script_file = script_factory("""
        import time

        [1]
        time.sleep(4)

        [2]
        time.sleep(4)
        """)
    ret = subprocess.Popen(["sos", "run", script_file])
    proc = psutil.Process(ret.pid)
    while True:
        children = proc.children(recursive=True)
        if children:
            children[0].terminate()
            break
        time.sleep(0.1)
    ret.wait()

    ret = subprocess.Popen(["sos", "run", script_file])
    proc = psutil.Process(ret.pid)
    while True:
        children = proc.children(recursive=True)
        if children:
            children[0].kill()
            break
        time.sleep(0.1)
    ret.wait()


def test_kill_substep_worker(script_factory):
    """Test if the workflow can error out after a worker is killed"""
    import time

    import psutil

    script_file = script_factory("""
        import time

        [1]
        input: for_each=dict(i=range(4))
        time.sleep(2)
        """)
    ret = subprocess.Popen(["sos", "run", script_file, "-j3"])
    proc = psutil.Process(ret.pid)
    while True:
        children = proc.children(recursive=True)
        print(children)
        if children:
            children[-1].terminate()
            break
        time.sleep(0.1)
    ret.wait()

    ret = subprocess.Popen(["sos", "run", script_file, "-j3"])
    proc = psutil.Process(ret.pid)
    while True:
        children = proc.children(recursive=True)
        if children:
            children[-1].kill()
            break
        time.sleep(0.1)
    ret.wait()


@pytest.mark.skipif(
    True, reason="This test needs to be improved to make it consistent")
def test_kill_task(script_factory):
    """Test if the workflow can error out after a worker is killed"""
    subprocess.call(["sos", "purge", "--all"])
    import time

    import psutil

    script_file = script_factory("""
        [1]
        task:
        import time
        time.sleep(10)
        """)
    ret = subprocess.Popen(
        ["sos", "run", script_file, "-s", "force", "-q", "localhost"])
    proc = psutil.Process(ret.pid)

    while True:
        children = proc.children(recursive=True)
        execute = [x for x in children if "execute" in x.cmdline()]
        if len(execute) >= 1:
            # a bug: if the process is killed too quickly (the signal
            # function is not called), this will fail.
            time.sleep(1)
            execute[0].terminate()
            break
        time.sleep(0.1)
    ret.wait()
    assert ret.returncode != 0

    ret = subprocess.Popen(["sos", "run", script_file, "-q", "localhost"])
    proc = psutil.Process(ret.pid)
    while True:
        children = proc.children(recursive=True)
        execute = [x for x in children if "execute" in x.cmdline()]
        if len(execute) >= 1:
            time.sleep(1)
            execute[0].kill()
            break
        time.sleep(0.1)
    ret.wait()
    assert ret.returncode != 0


@pytest.mark.skipif(
    True, reason="This test needs to be improved to make it consistent")
def test_restart_orphaned_tasks(script_factory):
    """Test restarting orphaned tasks which displays as running at first."""
    import time

    import psutil

    subprocess.call(["sos", "purge", "--all"])

    script_file = script_factory("""
        [1]
        task:
        import time
        time.sleep(12)
        """)
    ret = subprocess.Popen(
        ["sos", "run", script_file, "-s", "force", "-q", "localhost"])
    proc = psutil.Process(ret.pid)

    while True:
        children = proc.children(recursive=True)
        execute = [x for x in children if "execute" in x.cmdline()]
        if len(execute) >= 1:
            # a bug: if the process is killed too quickly (the signal
            # function is not called), this will fail.
            time.sleep(1)
            execute[0].kill()
            break
        time.sleep(0.1)
    proc.kill()
    #
    ret = subprocess.Popen(["sos", "run", script_file, "-q", "localhost"])
    ret.wait()
    assert ret.returncode == 0


def test_chdir_in_step(clear_now_and_after):
    """Test for #1415"""
    clear_now_and_after("test_chdir")
    execute_workflow("""
import os
os.makedirs('test_chdir', exist_ok=True)
os.chdir('test_chdir')

input: for_each=dict(i=range(2))
output: f'a_{i}.txt'
_output.touch()

    """)
    assert os.path.isfile("test_chdir/a_0.txt")
    assert os.path.isfile("test_chdir/a_0.txt")
