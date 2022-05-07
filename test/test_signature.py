#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil
import subprocess
import sys
import time

import pytest
from sos import execute_workflow
from sos.parser import SoS_Script
from sos.targets import file_target, sos_targets
from sos.utils import env
from sos.workflow_executor import Base_Executor


def _testSignature(script, steps):
    '''Test recognizing the format of SoS script'''
    #
    # only the first step
    res = execute_workflow(
        script,
        workflow=':0',
        options={
            'sig_mode': 'force',
            'default_queue': 'localhost',
        })
    assert os.path.isfile('temp/a.txt')
    assert os.path.isfile('temp/b.txt')
    ## FIXME
    # assert res['__completed__']['__step_completed__'], steps
    # FIXME
    with open('temp/a.txt') as ta:
        assert ta.read().strip() == 'a.txt'
    with open('temp/b.txt') as tb:
        assert tb.read().strip() == 'b.txt'

    res = execute_workflow(
        script,
        workflow=':0',
        options={
            'sig_mode': 'assert',
            'default_queue': 'localhost',
        })

    assert res['__completed__']['__step_completed__'] == 0
    # all of them
    res = execute_workflow(
        script, options={
            'sig_mode': 'default',
            'default_queue': 'localhost',
        })
    # now, rerun in build mode
    res = execute_workflow(
        script, options={
            'sig_mode': 'build',
            'default_queue': 'localhost',
        })

    assert res['__completed__']['__step_completed__'] == 0
    #
    assert os.path.isfile('temp/c.txt')
    assert os.path.isfile('temp/d.txt')
    with open('temp/c.txt') as tc:
        assert tc.read().strip() == 'a.txt'
    with open('temp/d.txt') as td:
        assert td.read().strip() == 'b.txt'
    assert env.sos_dict['oa'] == sos_targets('temp/c.txt', 'temp/d.txt')
    #
    # now in assert mode, the signature should be there
    res = execute_workflow(
        script, options={
            'sig_mode': 'assert',
            'default_queue': 'localhost',
        })
    assert res['__completed__']['__step_completed__'] == 0

    #
    res = execute_workflow(
        script, options={
            'sig_mode': 'default',
            'default_queue': 'localhost',
        })
    assert res['__completed__']['__step_completed__'] == 0

    #
    # change script a little bit
    res = execute_workflow(
        ' ' * (len(script) - len(script.lstrip())) + '# comment\n' + script,
        options={
            'sig_mode': 'assert',
            'default_queue': 'localhost',
        })

    assert res['__completed__']['__step_completed__'] == 0

    # add some other variable?
    #script = SoS_Script('comment = 1\n' + text)
    #wf = script.workflow()
    #env.config['sig_mode'] = 'assert'
    #self.assertRaises(Exception, Base_Executor(wf).run)


def test_signature(clear_now_and_after, clear_signatures):
    clear_now_and_after('temp/a.txt', 'temp/b.txt')
    _testSignature(
        r"""
        [*_0]
        output: 'temp/a.txt', 'temp/b.txt'
        task:
        run('''echo "a.txt" > temp/a.txt ''')
        run('''echo "b.txt" > temp/b.txt ''')

        [1: shared={'oa':'step_output'}]
        dest = ['temp/c.txt', 'temp/d.txt']
        input: group_by='single', paired_with='dest'
        output: _dest

        run(f" cp {_input} {_dest[0]} ")
        """, 2)


def test_signature1(clear_now_and_after, clear_signatures):
    clear_now_and_after('temp/a.txt', 'temp/b.txt')
    _testSignature(
        r"""
        [*_0]
        output: 'temp/a.txt', 'temp/b.txt'

        task:
        run('''echo "a.txt" > temp/a.txt ''')
        run('''echo "b.txt" > temp/b.txt ''')

        [1: shared={'oa':'step_output'}]
        dest = ['temp/c.txt', 'temp/d.txt']
        input: group_by='single', paired_with='dest'
        output: _dest

        run(f" cp {_input} {_dest[0]} ")
        """, 2)
    # script format


def test_signature2(clear_now_and_after, clear_signatures):
    clear_now_and_after('temp/a.txt', 'temp/b.txt')
    _testSignature(
        r"""
        [*_0]
        output: 'temp/a.txt', 'temp/b.txt'

        run:
        echo "a.txt" > temp/a.txt

        run:

        echo "b.txt" > temp/b.txt

        [1: shared={'oa':'step_output'}]
        dest = ['temp/c.txt', 'temp/d.txt']
        input: group_by='single', paired_with='dest'
        output: _dest

        task:
        run: expand=True
        echo cp {_input} {_dest[0]}
        cp {_input} {_dest[0]}
        """, 2)


def test_signature_with_shared_variable(clear_now_and_after, clear_signatures):
    '''Test restoration of signature from variables.'''
    clear_now_and_after('a.txt')
    # shared
    script = r"""
        [0: shared='a']
        output: 'a.txt'
        run:
            touch a.txt

        a = 5

        [1]
        print(a)
        """
    # alias should also be recovered.
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 2
    # rerun
    env.sos_dict.pop('a')
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1


def test_signature_without_output(clear_now_and_after, clear_signatures):
    # signature without output file
    clear_now_and_after('temp/a.txt', 'temp/b.txt', 'temp')

    _testSignature(
        r"""
        [*_0]
        output: []

        run:
            [ -d temp ] || mkdir temp
            echo "a.txt" > temp/a.txt
            echo "b.txt" > temp/b.txt

        [1: shared={'oa':'step_output'}]
        dest = ['temp/c.txt', 'temp/d.txt']
        input: 'temp/a.txt', 'temp/b.txt', group_by='single', paired_with='dest'
        output: _dest

        run: expand=True
            cp {_input} {_dest[0]}
        """, 2)


def test_reexecution(clear_signatures):
    '''Test -f option of sos run'''
    script = '''
        [0]
        output: 'a.txt'
        run(f"touch {_output}")
        '''
    try:
        # remove existing output if exists
        file_target('a.txt').unlink()
    except Exception:
        pass
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    # now, rerun should be much faster
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0
    # rerun takes less than 1 second
    #
    # force rerun mode
    res = execute_workflow(script, options={'sig_mode': 'ignore'})
    assert res['__completed__']['__step_completed__'] == 1
    # regularly take more than 5 seconds to execute
    try:
        # remove existing output if exists
        os.remove('a.txt')
    except Exception:
        pass


@pytest.mark.skipif(
    sys.platform == 'win32',
    reason='Windows executable cannot execute bash loop.')
def test_signature_after_removal_of_files(clear_signatures,
                                          clear_now_and_after):
    '''test action to_named_path'''
    clear_now_and_after('largefile.txt')
    script = r'''
        [10]

        # generate a file
        output: 'largefile.txt'

        run: expand='${ }'
        for x in {1..1000}
        do
            echo $x >> ${_output}
        done

        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    # rerun, because this is the final target, it has to be
    # re-generated
    os.remove('largefile.txt')
    res = execute_workflow(script)
    assert os.path.isfile('largefile.txt')
    assert res['__completed__']['__step_completed__'] == 1
    #
    # we discard the signature, the step would still be
    # skipped because file signature will be calculated
    # during verification
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0
    #
    # now if we touch the file, it needs to be regenerated
    with open('largefile.txt', 'a') as lf:
        lf.write('something')
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1


@pytest.mark.skipif(
    sys.platform == 'win32',
    reason='Windows executable cannot be created with chmod.')
def test_removal_of_intermediate_files(clear_signatures, clear_now_and_after):
    # if we zap the file, it
    clear_now_and_after('midfile.txt', 'finalfile.txt', 'midfile.txt.zapped')

    script = r'''
        [10]

        # generate a file
        output: 'midfile.txt'

        run: expand='${ }'
        rm -f ${_output}
        for x in {1..1000}
        do
            echo $x >> ${_output}
        done

        [20]
        output: 'finalfile.txt'
        run: expand=True
        cp {_input} {_output}
        echo "MORE" >> {_output}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 2
    #
    # remove middle file, rerun
    os.remove('midfile.txt')
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    assert os.path.isfile('midfile.txt')
    #
    # now if we touch the mid file, it needs to be regenerated
    with open('midfile.txt', 'a') as lf:
        lf.write('something')
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    #
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0
    #
    # if we zap the mid file, it does not need to be rerun
    subprocess.call('sos remove midfile.txt --zap -y', shell=True)
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0


def test_signature_with_parameter(clear_now_and_after):
    '''Test signature'''
    clear_now_and_after('myfile.txt')
    #
    script = r'''
        parameter: gvar = 10

        [10]
        # generate a file
        output: 'myfile.txt'
        # additional comment
        run: expand=True
        echo {gvar} > {_output:q}

        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    with open('myfile.txt') as tmp:
        assert tmp.read().strip() == '10'
    #
    # now if we change parameter, the step should be rerun
    res = execute_workflow(script, args=['--gvar', '20'])
    assert res['__completed__']['__step_completed__'] == 1
    with open('myfile.txt') as tmp:
        assert tmp.read().strip() == '20'
    #
    # do it again, signature should be effective
    res = execute_workflow(script, args=['--gvar', '20'])
    assert res['__completed__']['__step_completed__'] == 0
    with open('myfile.txt') as tmp:
        assert tmp.read().strip() == '20'

    #
    script = r'''
        parameter: gvar = 10

        [10]
        # generate a file
        output: 'myfile.txt'
        # additional comment
        run: expand=True
        echo {gvar} > {_output:q}
        '''
    res = execute_workflow(script)
    with open('myfile.txt') as tmp:
        assert tmp.read().strip() == '10'
    assert res['__completed__']['__step_completed__'] == 1
    #
    # now if we change parameter, the step should be rerun
    res = execute_workflow(script, args=['--gvar', '20'])
    with open('myfile.txt') as tmp:
        assert tmp.read().strip() == '20'
    assert res['__completed__']['__step_completed__'] == 1
    #
    # do it again, signature should be effective
    res = execute_workflow(script, args=['--gvar', '20'])
    assert res['__completed__']['__step_completed__'] == 0
    with open('myfile.txt') as tmp:
        assert tmp.read().strip() == '20'


def test_loop_wise_signature(clear_now_and_after):
    '''Test partial signature'''
    clear_now_and_after([f'myfile_{i}.txt' for i in range(10, 12)])
    #
    script = r'''
        parameter: gvar = 10

        [10]
        tt = [gvar]
        input: for_each='tt'
        output: f"myfile_{_tt}.txt"
        run: expand=True
        echo "DO {_tt}"
        echo {_tt} > {_output:q}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    ts = os.path.getmtime('myfile_10.txt')
    #
    # now we modify the script
    script = r'''
        parameter: gvar = 10

        [10]
        tt = [gvar, gvar + 1]
        input: for_each='tt'
        output: f"myfile_{_tt}.txt"
        run: expand=True
        echo "DO {_tt}"
        echo {_tt} > {_output:q}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0.5
    # this file is not regenerated
    assert ts == os.path.getmtime('myfile_10.txt')
    ts1 = os.path.getmtime('myfile_11.txt')
    #
    # run it again, neither needs to be rerun
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0
    assert ts == os.path.getmtime('myfile_10.txt')
    assert ts1 == os.path.getmtime('myfile_11.txt')
    #
    # change again, the second one is already there.
    script = r'''
        parameter: gvar = 10

        [10]
        tt = [gvar + 1]
        input: for_each='tt'
        output: f"myfile_{_tt}.txt"
        run: expand=True
        echo "DO {_tt}"
        echo {_tt} > {_output:q}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0
    assert ts1 == os.path.getmtime('myfile_11.txt')
    #
    for t in range(10, 12):
        with open('myfile_{}.txt'.format(t)) as tmp:
            assert tmp.read().strip() == str(t)


def test_output_from_signature(temp_factory, clear_now_and_after):
    'Test restoration of output from signature' ''
    temp_factory('1.txt', '2.txt')
    clear_now_and_after('1.out', '2.out', '1.2.out', '2.3.out')
    script = '''
        parameter: K = [2,3]

        [work_1]
        input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
        output: expand_pattern('{_name}.out')
        run: expand=True
        touch {_output}

        [work_2]

        input: group_by = 'single', pattern = '{name}.{ext}', paired_with = ['K']
        output: expand_pattern('{_name}.{_K}.out')
        run: expand=True
        touch {_output}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 2
    # for the second run, output should be correctly constructed
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0


def test_signature_with_vars(temp_factory, clear_now_and_after):
    '''Test revaluation with variable change'''
    temp_factory('a1.out', 'a2.out')
    clear_now_and_after('b1.out', 'b2.out')
    script = '''
        parameter: DB = {'input': ['a1.out'], 'output': ['b1.out']}
        parameter: input_file = DB['input']
        parameter: output_file =  DB['output']

        [2]
        input: input_file, group_by = 1
        output: output_file[_index]
        run: expand=True
        touch {_output}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1
    ts = os.path.getmtime('b1.out')
    #
    script = '''
        parameter: DB = {'input': ['a1.out', 'a2.out'], 'output': ['b1.out', 'b2.out']}
        parameter: input_file = DB['input']
        parameter: output_file =  DB['output']

        [2]
        input: input_file, group_by = 1
        output: output_file[_index]
        run: expand=True
        touch {_output}
        '''
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 0.5
    assert ts == os.path.getmtime('b1.out')


def test_action_signature(clear_signatures, clear_now_and_after):
    '''Test action signature'''
    clear_now_and_after('test_action.txt', 'lc.txt')
    with open('test_action.txt', 'w') as ta:
        ta.write('#something\n')
    script = r'''
        [1]
        input: 'test_action.txt'
        run: output='lc.txt', expand=True, tracked='test_action.txt'
        sleep 5
        wc -l {_input[0]} > lc.txt
        '''
    execute_workflow(script)
    # the second time, should skip
    st = time.time()
    execute_workflow(script)
    # should skip
    assert time.time() - st < 5
    #
    script = r'''
        [1]
        input: 'test_action.txt'
        print('step is changed')
        run: output='lc.txt', expand=True, tracked='test_action.txt'
        sleep 5
        wc -l {_input[0]} > lc.txt
        '''
    st = time.time()
    execute_workflow(script)
    assert time.time() - st < 5

    # force
    execute_workflow(script, options={'sig_mode': 'build'})


def test_signature_with_without_task(clear_now_and_after):
    '''Test the inclusion of task would not trigger rerun'''
    script = r'''
        [1]
        output: 'aa'
        sh:
            echo aa > aa
        '''
    clear_now_and_after('aa')
    res = execute_workflow(script)
    assert res['__completed__']['__step_completed__'] == 1

    script = r'''
        [1]
        output: 'aa'
        task:
        sh:
            echo aa > aa
        '''
    res = execute_workflow(script, config={'default_queue': 'localhost'})
    assert res['__completed__']['__step_completed__'] == 0


def test_signature_with_dynamic_output(clear_signatures, clear_now_and_after):
    '''Test return of output from dynamic output'''
    clear_now_and_after([f'rep_{i}' for i in range(5)])

    script = r'''
        [1: shared={'step1': 'step_output'}]
        input: for_each={'i': range(5)}, concurrent=True
        output: dynamic(f'rep_{i}/*.res')

        import random, os
        os.makedirs(f'rep_{i}', exist_ok=True)
        path(f'rep_{i}/{random.randint(0, 10000)}.res').touch()
        '''
    res = execute_workflow(script)
    files = env.sos_dict['step1']
    assert len(files) == 5
    assert res['__completed__']['__substep_completed__'] == 5
    # rerun
    res = execute_workflow(script)
    files_again = env.sos_dict['step1']
    assert files == files_again
    assert res['__completed__']['__substep_completed__'] == 0


def test_ignore_signature(clear_signatures, clear_now_and_after):
    '''Test ignore signature mode #1028 '''
    clear_now_and_after([f'out_{i}.txt' for i in range(3)])

    execute_workflow(
        r'''
        input: for_each={'i': range(3)}, concurrent=True
        output: f'out_{i}.txt'
        sh: expand=True
        touch {_output}
        ''',
        options={'sig_mode': 'ignore'})


def test_rebuid_signature(clear_signatures, clear_now_and_after):
    '''Test rebuilding signature'''
    clear_now_and_after('a.txt')
    script = r'''
        [A_1]
        output: 'a.txt'
        _output.touch()
        '''
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_completed__'] == 1
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 1
    res = execute_workflow(script, options={'sig_mode': 'build'})
    assert res['__completed__']['__substep_skipped__'] == 1
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 1
    # if a.txt is changed, rebuild will rerun
    with open('a.txt', 'a') as atxt:
        atxt.write('aaa')
    res = execute_workflow(script, options={'sig_mode': 'build'})
    assert res['__completed__']['__substep_skipped__'] == 1
    # rerun?
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 1
    #
    res = execute_workflow(script, options={'sig_mode': 'force'})
    assert res['__completed__']['__substep_completed__'] == 1
    #
    res = execute_workflow(script, options={'sig_mode': 'ignore'})
    assert res['__completed__']['__substep_completed__'] == 1


def test_rebuid_signature_with_substeps(clear_signatures, clear_now_and_after):
    '''Test rebuilding signature'''
    clear_now_and_after([f'a_{i}.txt' for i in range(4)])
    script = r'''
        [A_1]
        input: for_each=dict(i=range(4))
        output: f'a_{i}.txt'
        _output.touch()
        '''
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_completed__'] == 4
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 4
    res = execute_workflow(script, options={'sig_mode': 'build'})
    assert res['__completed__']['__substep_skipped__'] == 4
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 4
    # if a.txt is changed, rebuild will rerun
    for i in range(4):
        with open(f'a_{i}.txt', 'a') as atxt:
            atxt.write('aaa')

    res = execute_workflow(script, options={'sig_mode': 'build'})
    assert res['__completed__']['__substep_skipped__'] == 4
    # rerun?
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 4
    #

    res = execute_workflow(script, options={'sig_mode': 'force'})
    assert res['__completed__']['__substep_completed__'] == 4
    #
    res = execute_workflow(script, options={'sig_mode': 'ignore'})
    assert res['__completed__']['__substep_completed__'] == 4


def test_rebuid_signature_with_tasks(clear_signatures, clear_now_and_after):
    '''Test rebuilding signature'''
    clear_now_and_after([f'a_{i}.txt' for i in range(4)])
    script = r'''
[A_1]
input: for_each=dict(i=range(2))
output: f'a_{i}.txt'
task:
_output.touch()
'''
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_completed__'] == 2

    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 2

    res = execute_workflow(script, options={'sig_mode': 'build'})
    assert res['__completed__']['__substep_skipped__'] == 2

    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 2
    # if a.txt is changed, rebuild will rerun
    for i in range(2):
        with open(f'a_{i}.txt', 'a') as atxt:
            atxt.write('aaa')

    res = execute_workflow(script, options={'sig_mode': 'build'})
    assert res['__completed__']['__substep_skipped__'] == 2
    # rerun?
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 2
    #
    res = execute_workflow(script, options={'sig_mode': 'force'})
    assert res['__completed__']['__substep_completed__'] == 2
    #
    res = execute_workflow(script, options={'sig_mode': 'ignore'})
    assert res['__completed__']['__substep_completed__'] == 2


def test_signature_with_dependency_tracing_and_vars(clear_signatures,
                                                    clear_now_and_after):
    '''Test signature with parameters with option -T #1200'''
    clear_now_and_after('out.txt')
    script = r'''
        [analyze]
        parameter: par=5
        output: 'out.txt'
        print(par)
        _output.touch()

        [default]
        input: 'out.txt'
        '''
    res = execute_workflow(
        script, options={'trace_existing': True}, args=['--par', '10'])
    assert res['__completed__']['__substep_completed__'] == 2
    # change parameter, the step should not be skipped
    res = res = execute_workflow(
        script, options={'trace_existing': True}, args=['--par', '20'])
    assert res['__completed__']['__substep_completed__'] == 2


def test_skip_mode(clear_signatures, temp_factory, clear_now_and_after):
    '''Test skipping mode of signature'''
    clear_now_and_after([f'a_{i}.bak' for i in range(4)])
    temp_factory([f'a_{i}.txt' for i in range(4)])
    #
    script = r'''
        [A_1]
        input: [f'a_{i}.txt' for i in range(4)], group_by=1
        output: _input.with_suffix('.bak')

        with open(_input) as ifile, open(_output, 'w') as ofile:
            ofile.write(ifile.read())
        '''
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_completed__'] == 4
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_skipped__'] == 4
    #
    res = execute_workflow(script, options={'sig_mode': 'skip'})
    assert res['__completed__']['__substep_skipped__'] == 4
    #
    # if result file is changed, skip will still skip
    for i in range(4):
        with open(f'a_{i}.bak', 'a') as bak:
            bak.write('extra')
    #
    res = execute_workflow(script, options={'sig_mode': 'skip'})
    assert res['__completed__']['__substep_skipped__'] == 4
    # now if we change to default mode, will not skip
    res = execute_workflow(script, options={'sig_mode': 'default'})
    assert res['__completed__']['__substep_completed__'] == 4


def test_signature_with_output_vars(clear_signatures, clear_now_and_after):
    '''Test persistence of _output with vars #1355'''
    clear_now_and_after('test_sig_with_vars.txt')
    script = '''
        [10]
        output: f'test_sig_with_vars.txt'
        _output.touch()
        _output.set(seed=1)

        [20]
        print(_input.seed)
            '''
    execute_workflow(script)
    # second time sshould be fine (using signature
    execute_workflow(script)
    #
    # worker...
    for i in range(2):
        if os.path.isfile(f'test_sig_with_vars_{i}.txt'):
            os.remove(f'test_sig_with_vars_{i}.txt')
    script = '''
        [10]
        input: for_each=dict(i=range(2))
        output: f'test_sig_with_vars_{i}.txt'
        _output.touch()
        _output.set(seed=1)

        [20]
        print(_input.seed)
            '''
    execute_workflow(script)
    # second time sshould be fine (using signature
    execute_workflow(script)
    # task
    # worker...
    for i in range(2):
        if os.path.isfile(f'test_sig_with_vars_{i}.txt'):
            os.remove(f'test_sig_with_vars_{i}.txt')
    script = '''
        [10]
        input: for_each=dict(i=range(2))
        output: f'test_sig_with_vars_{i}.txt'
        task: queue='localhost'
        _output.touch()
        _output.set(seed=1)

        [20]
        print(_input.seed)
            '''
    execute_workflow(script)
    # second time sshould be fine (using signature
    execute_workflow(script)


def test_signature_with_local_parameter(clear_signatures):
    '''test for #1372'''
    script = '''
        [default]
        input:
        output: f'constant-file-name.txt'
        parameter: n_mx = int
        bash: expand=True

        echo {n_mx} > constant-file-name.txt
        '''
    execute_workflow(script, args=['--n-mx', '60'])

    assert os.path.isfile('constant-file-name.txt')
    assert open('constant-file-name.txt').read().strip() == '60'

    execute_workflow(script, args=['--n-mx', '80'])
    assert open('constant-file-name.txt').read().strip() == '80'