#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import os
import shutil
import subprocess

import pytest
from sos import execute_workflow
from sos.parser import SoS_Script
from sos.targets import file_target
from sos.utils import env


def test_progress_bar():
    # progress bar with nested workflow
    execute_workflow(
        '''
import time
time.sleep(0)
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
''',
        workflow='a')


def test_nested_workflow(temp_factory, clear_now_and_after):
    '''Test the creation and execution of combined workfow'''
    clear_now_and_after('a.done')
    temp_factory('a.txt', 'b.txt', 'b.begin')
    execute_workflow(
        '''
        if 'executed' not in locals():
            executed = []
        if 'inputs' not in locals():
            inputs = []

        [a_1: shared=['executed', 'inputs']]
        executed.append(step_name)
        inputs.append(_input)

        [a_2: shared=['executed', 'inputs']]
        executed.append(step_name)
        inputs.append(_input)

        [a_3: shared=['executed', 'inputs']]
        executed.append(step_name)
        inputs.append(_input)

        [a_4: shared=['executed', 'inputs']]
        executed.append(step_name)
        output: 'a.done'
        inputs.append(_input)
        run: expand=True
        touch {_output}

        [b_1: shared=['executed', 'inputs']]
        executed.append(step_name)
        input: 'b.begin'
        inputs.append(_input)

        [b_2: shared=['executed', 'inputs']]
        executed.append(step_name)
        inputs.append(_input)

        [b_3: shared=['executed', 'inputs']]
        executed.append(step_name)
        inputs.append(_input)

        [b_4: shared=['executed', 'inputs']]
        executed.append(step_name)
        output: 'b.txt'
        inputs.append(_input)

        [c: shared=['executed', 'inputs']]
        executed.append(step_name)
        input: 'a.txt'
        output: 'b.txt'
        inputs.append(_input)
        sos_run('a+b', shared=['executed', 'inputs'])
        ''',
        workflow='c',
        options={'sig_mode': 'ignore'})

    # order of execution is not guaranteed
    assert sorted(env.sos_dict['executed']) == sorted(
        ['c', 'a_1', 'a_2', 'a_3', 'a_4', 'b_1', 'b_2', 'b_3', 'b_4'])
    env.sos_dict.pop('executed', None)


def test_looped_nested_workflow(temp_factory, clear_now_and_after):
    # step will be looped
    clear_now_and_after('a.txt.a1', 'a.txt.a1.a2', 'b.txt.a1', 'b.txt.a1.a2')
    temp_factory('a.txt', 'b.txt')

    execute_workflow(
        '''
        if 'executed' not in locals():
            executed = []
        if 'inputs' not in locals():
            inputs = []

        [a_1:shared=['executed', 'inputs']]
        executed.append(step_name)
        output: _input[0] + '.a1'
        inputs.append(_input)
        run: expand=True
        touch {_output}

        [a_2:shared=['executed', 'inputs']]
        executed.append(step_name)
        output: _input[0] + '.a2'
        inputs.append(_input)
        run: expand=True
        touch {_output}

        [c:shared=['executed', 'inputs']]
        executed.append(step_name)
        input: 'a.txt', 'b.txt', group_by='single'
        inputs.append(_input)
        sos_run('a', shared=['executed', 'inputs'])
        ''',
        workflow='c')
    assert env.sos_dict['executed'] == ['c', 'a_1', 'a_2', 'a_1', 'a_2']


def test_single_looped_nested_workflow(temp_factory, clear_now_and_after):
    temp_factory('a.txt', 'b.txt')
    env.sos_dict.pop('executed', None)

    execute_workflow(
        '''
        if 'executed' not in locals():
            executed = []
        [a_1:shared='executed']
        executed.append(step_name)
        [a_2:shared='executed']
        executed.append(step_name)
        [c_0:shared='executed']
        executed.append(step_name)
        [c_1:shared='executed']
        depends: sos_variable('executed')
        executed.append(step_name)
        input: 'a.txt', 'b.txt', group_by='single'
        sos_run('a:2', shared='executed')
        ''',
        workflow='c')
    assert env.sos_dict['executed'] == ['c_0', 'c_1', 'a_2', 'a_2']
    env.sos_dict.pop('executed', None)
    # allow specifying a single step
    # step will be looped
    execute_workflow(
        '''
        if 'executed' not in locals():
            executed = []

        [a_1:shared='executed']
        executed.append(step_name)

        [a_2:shared='executed']
        executed.append(step_name)

        [c_0:shared='executed']
        executed.append(step_name)

        [c_1:shared='executed']
        depends: sos_variable('executed')
        executed.append(step_name)
        input: 'a.txt', 'b.txt', group_by='single'
        sos_run('a:2', shared='executed')
        ''',
        workflow='c')

    assert env.sos_dict['executed'] == ['c_0', 'c_1', 'a_2', 'a_2']


def test_recursive_nested_workflow(temp_factory):
    # recursive subworkflow not allowed
    temp_factory('a.txt', 'b.txt')
    with pytest.raises(Exception):
        execute_workflow(
            '''
            if 'executed' not in locals():
                executed = []
            [a_1:shared='executed']
            executed.append(step_name)
            [a_2:shared='executed']
            executed.append(step_name)
            [c_0:shared='executed']
            executed.append(step_name)
            [c_1:shared='executed']
            executed.append(step_name)
            input: 'a.txt', 'b.txt', group_by='single'
            sos_run('a_2+c', shared='executed')
            ''',
            workflow='c')
    #
    env.sos_dict.pop('executed', None)
    # nested subworkflow is allowed
    execute_workflow(
        '''
        if 'executed' not in locals():
            executed = []
        [a_1:shared='executed']
        executed.append(step_name)
        [a_2:shared='executed']
        executed.append(step_name)
        [a_3:shared='executed']
        executed.append(step_name)
        [b_1:shared='executed']
        executed.append(step_name)
        [b_2:shared='executed']
        executed.append(step_name)
        sos_run('a:1-2', shared='executed')
        [c_0:shared='executed']
        executed.append(step_name)
        [c_1:shared='executed']
        depends: sos_variable('executed')
        executed.append(step_name)
        input: 'a.txt'
        sos_run('a+b', shared='executed')
        ''',
        workflow='c')

    assert env.sos_dict['executed'] == [
        'c_0', 'c_1', 'a_1', 'a_2', 'a_3', 'b_1', 'b_2', 'a_1', 'a_2'
    ]

    env.sos_dict.pop('executed', None)


def test_subworkflow_with_options(temp_factory, clear_now_and_after):
    # nested subworkflow with step option and others
    temp_factory('a.txt', 'b.txt')
    clear_now_and_after('a.done')
    script = '''
        if 'executed' not in locals():
            executed = []
        [a_1:shared='executed']
        executed.append(step_name)

        [a_2:shared='executed']
        executed.append(step_name)

        [a_3:shared='executed']
        executed.append(step_name)

        [b:shared='executed']
        executed.append(step_name)
        input: 'a.txt', 'b.txt', group_by='single'
        sos_run('a:3+a:1', shared='executed')

        [d:shared='executed']
        executed.append(step_name)
        input: 'a.txt', 'b.txt', group_by='single'
        sos_run('a:2', shared='executed')

        [e2_2:shared='executed']
        executed.append(step_name)
        input: 'a.txt', 'b.txt', group_by='single'
        '''
    execute_workflow(script, workflow='b')
    assert env.sos_dict['executed'] == ['b', 'a_3', 'a_1', 'a_3', 'a_1']

    env.sos_dict.pop('executed', None)

    execute_workflow(script, workflow='d')
    assert env.sos_dict['executed'] == ['d', 'a_2', 'a_2']
    env.sos_dict.pop('executed', None)

    execute_workflow(script, workflow='e2')
    assert env.sos_dict['executed'] == ['e2_2']


def test_dynamic_nested_workflow():
    '''Test nested workflow controlled by command line option'''
    script = '''
        if 'executed' not in locals():
            executed = []

        parameter: wf='a'

        [a_1:shared='executed']
        executed.append(step_name)
        [a_2:shared='executed']
        executed.append(step_name)
        [a_3:shared='executed']
        executed.append(step_name)
        [b_1:shared='executed']
        executed.append(step_name)
        [b_2:shared='executed']
        executed.append(step_name)
        [b_3:shared='executed']
        executed.append(step_name)

        [default:shared='executed']
        executed.append(step_name)
        sos_run(wf, shared='executed')
        '''
    execute_workflow(script, args=['--wf', 'b'])
    assert env.sos_dict['executed'] == ['default', 'b_1', 'b_2', 'b_3']
    #
    env.sos_dict.pop('executed', None)
    execute_workflow(script, args=['--wf', 'a'])
    assert env.sos_dict['executed'] == ['default', 'a_1', 'a_2', 'a_3']


def test_sos_run(clear_now_and_after):
    '''Test action sos_run with keyword parameters'''
    clear_now_and_after('0.txt', '1.txt')
    execute_workflow(
        r'''
        [A]
        parameter: num=5

        run: expand=True
            touch {num}.txt

        [batch]
        for k in range(2):
            sos_run('A', num=k)
        ''',
        workflow='batch')

    for f in ['0.txt', '1.txt']:
        assert file_target(f).target_exists()
    #
    clear_now_and_after('0.txt', '1.txt', '5.txt')

    execute_workflow(
        r'''
        [A]
        parameter: num=5
        run: expand=True
            touch {num}.txt

        [batch]
        for num in range(2):
            sos_run('A')
        ''',
        workflow='batch')
    for f in ['0.txt', '1.txt']:
        assert not file_target(f).target_exists()

    assert file_target('5.txt').target_exists()

    clear_now_and_after('10.txt', '11.txt')
    #
    # test parameter shared to send and return vars
    #
    execute_workflow(
        r'''
        [A: shared='k']
        k += 10

        [batch]
        for k in range(2):
            sos_run('A', shared='k')
            run(f"touch {k}.txt")
        ''',
        workflow='batch')
    for f in ['10.txt', '11.txt']:
        assert file_target(f).target_exists()


def test_da_gof_dynamic_nested_workflow(clear_now_and_after):
    #
    # Because we are not sure which workflows would be executed
    # until run time, the DAG should not contain nested workflow
    # until runtime.
    #
    clear_now_and_after('B0.txt', 'B0.txt.p', 'B1.txt', 'B1.txt.p', 'B2.txt',
                        'B2.txt.p')
    #
    #  A1 <- P <- B
    #  A1 <- P <- B
    #  A2
    #
    #  ALL calls A and B with parameter
    #
    execute_workflow(
        '''
        [A_1]
        parameter: num = 2
        input: f"B{num}.txt.p"

        [B: provides='B{num}.txt']
        run: expand=True
        touch 'B{num[0]}.txt'

        [P: provides='{filename}.p']
        input: filename
        run: expand=True
        touch {_output}

        [ALL]

        for i in range(3):
            sos_run('A', num=i)
    ''',
        workflow='ALL')
    for f in ['B0.txt', 'B0.txt.p', 'B1.txt', 'B1.txt.p', 'B2.txt', 'B2.txt.p']:
        assert file_target(f).target_exists()


def test_outcome_oriented_nested_workflow(clear_now_and_after):
    '''test nested workflow triggered by targets'''
    clear_now_and_after('test_15.txt')

    execute_workflow('''
[A: provides='test_{idx}.txt']
_output.touch()

[default]
sos_run(targets='test_15.txt')
    ''')
    assert os.path.isfile('test_15.txt')


def test_passing_vars_to_nested_workflow():
    '''Test if variables can be passed to nested workflows'''
    execute_workflow(r"""
        import time
        import random

        [nested]
        parameter: nested=True
        parameter: seed=1
        print(f'I am nested {nested} with seed {seed}')

        [0]
        reps = range(5)
        input: for_each='reps'
        import random
        nested = _reps
        seed = random.randint(1, 1000)
        print(f'Passing {seed} to {nested}')
        sos_run('nested', nested=nested, seed=seed)
    """)


def test_user_defined_func():
    '''Test the use of user-defined functions in SoS script'''
    execute_workflow(
        r"""
        def myfunc():
            return 'a'

        [1: shared={'test':'_output'}]
        output: myfunc()

        myfunc()

        """,
        options={'run_mode': 'dryrun'})
    assert env.sos_dict['test'] == ['a']
    # User defined function should also work under nested workflows
    # This is difficult because the 'local namespace' is usually
    # not seen inside function definition. The solution now is to
    # use a single workspace.
    execute_workflow(
        r"""
        def myfunc():
            # test if builtin functions (sum and range) can be used here.
            return 'a' + str(sum(range(10)))

        [1: shared={'test':'_output'}]
        output: [myfunc() for i in range(10)][0]

        myfunc()

        """,
        options={'run_mode': 'dryrun'})
    assert env.sos_dict['test'] == ['a45']


def test_config_file_of_nested_workflow(config_factory):
    '''Test passing of configurationg to nested workflow'''
    cfg = config_factory("""{'1':'hi'}""")

    execute_workflow(
        '''
[test_1]
parameter: key = None
print(CONFIG[key])

[default_1]
sos_run('test:1', key = '1')
''',
        options={'config_file': cfg})


def test_error_from_subworkflow():
    '''Test if error from subworkflow is passed to master (#396)'''
    with pytest.raises(Exception):
        execute_workflow('''
            [test_1]
            R:
            set.seed(xxx)

            [default]
            sos_run('test')
            ''')


def test_fun_def(temp_factory):
    '''Test defintion of function that can be used by other steps'''
    temp_factory('aa.txt', 'ab.txt')
    # in nested workflow?
    execute_workflow(
        r"""
        def myfunc(a):
            return ['a' + x for x in a]

        [mse: shared={'test':'_output'}]
        input: myfunc(['a.txt', 'b.txt'])

        [1]
        sos_run('mse')
        """,
        options={'run)mode': 'dryrun'})
    #
    # Names defined in subworkflow is not returned to the master dict
    assert 'test' not in env.sos_dict


def test_search_path(clear_now_and_after):
    '''Test if any action should exit in five seconds in dryrun mode'''
    clear_now_and_after('crazy_path', 'test.yml')

    sos_config_file = os.path.join(
        os.path.expanduser('~'), '.sos', 'config.yml')
    shutil.copy(sos_config_file, 'test.yml')
    #
    subprocess.call(
        'sos config --set sos_path {0}/crazy_path {0}/crazy_path/more_crazy/'
        .format(os.getcwd()),
        shell=True)
    #
    if not os.path.isdir('crazy_path'):
        os.mkdir('crazy_path')
        os.mkdir(os.path.join('crazy_path', 'more_crazy'))
    with open(os.path.join('crazy_path', 'crazy_master.sos'), 'w') as crazy:
        crazy.write('''
[0]
sos_run('cc', source='crazy_slave.sos')

''')
    with open(os.path.join('crazy_path', 'more_crazy', 'crazy_slave.sos'),
              'w') as crazy:
        crazy.write('''
[cc_0]
print('hay, I am crazy')
''')

    script = SoS_Script(filename='crazy_master.sos')
    script.workflow()
    #
    os.remove(sos_config_file)
    shutil.copy('test.yml', sos_config_file)


def test_nested_workdir(clear_now_and_after, purge_tasks):
    '''Test nested runtime option for work directory'''
    clear_now_and_after('tmp', 'tmp1')
    execute_workflow(
        '''
        [step]
        task: workdir='tmp1'
        run:
            touch 'a.txt'

        [default]
        sos_run('step', workdir='tmp')
        ''',
        options={'default_queue': 'localhost'})

    assert os.path.isfile('tmp1/a.txt')


def test_failure_of_nested_workflow(clear_now_and_after):
    '''Test failure of nested workflow #838'''
    clear_now_and_after('a.txt')
    with pytest.raises(Exception):
        execute_workflow('''
            [something]
            input: 'a.txt'

            [default]
            sos_run('something')
            ''')


def test_nested_from_another_file(clear_now_and_after):
    '''Test nested runtime option for work directory'''
    clear_now_and_after('a.txt', 'another.sos')
    with open('another.sos', 'w') as another:
        another.write('''
[whatever]
run:
touch 'a.txt'

''')
    execute_workflow('''
        [default]
        sos_run('whatever', source='another.sos')
        ''')
    assert os.path.isfile(
        'a.txt'
    ), 'a.txt should have been created by nested workflow from another file'


def test_concurrent_sub_workflow():
    '''Test concurrent subworkflow sos_run '''
    import time
    st = time.time()

    execute_workflow(
        '''
        [A]
        parameter: idx=0
        import time
        time.sleep(5)

        [default]
        input: for_each=dict(i=range(6))
        sos_run('A', idx=i)
        ''',
        options={'worker_procs': ['8']})
    assert time.time() - st < 30


def test_sos_multi_workflow():
    '''Test multiple workflows in sos_run '''
    import time
    st = time.time()

    execute_workflow(
        '''
        [B]
        parameter: idx=2
        import time
        time.sleep(idx)

        [A]
        import time
        time.sleep(idx)

        [default]
        input: for_each=dict(i=range(4))
        sos_run(['A', 'B'], idx=i)
    ''',
        options={'worker_procs': ['8']})
    assert time.time() - st < 20


def test_pass_of_args():
    '''Test passing of arguments through sos_run #1164'''
    with pytest.raises(Exception):
        execute_workflow('''
            parameter: b=2
            [A]
            parameter: c=1
            print(c)

            [default]
            sos_run('A', a=2)
        ''')
    #
    execute_workflow('''
        parameter: b=2
        [A]
        parameter: c=1
        print(c)

        [default]
        sos_run('A', b=2)
        ''')
    #
    execute_workflow('''
        parameter: b=2
        [A]
        parameter: c=1
        print(c)

        [default]
        sos_run('A', c=2)
    ''')


def test_nested_dynamic_depends(clear_now_and_after):
    '''Test the execution of nested workflow with dynamic depends'''
    clear_now_and_after('B30.txt', 'B30.txt.p', 'B0.txt.p', 'B0.txt')
    execute_workflow('''
        [A_1]
        parameter: num = 20
        input: f"B{num}.txt"
        depends: _input + '.p'

        [B: provides='B{num}.txt']
        run: expand=True
            touch {_output}

        [P: provides='{filename}.p']
        input: filename
        run: expand=True
            touch {_output}

        [default]
        sos_run('A', num=0)
    ''')


def test_neste_with_both_workflow_and_targets(clear_now_and_after):
    '''Test nested workflow with both workflow and targets'''
    clear_now_and_after('A_out.txt', 'B_out.txt')
    execute_workflow(r'''
        [A_20]
        output: 'A_out.txt'
        _output.touch()

        [B]
        output: 'B_out.txt'
        _output.touch()

        [default]
        sos_run('A', targets='B_out.txt')
        ''')
    assert os.path.isfile('A_out.txt')
    assert os.path.isfile('B_out.txt')
