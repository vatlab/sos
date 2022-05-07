#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import sys
import random

import pytest
from sos.eval import accessed_vars, on_demand_options
from sos.parser import SoS_Script
from sos.pattern import expand_pattern, extract_pattern
from sos.targets import executable, file_target, sos_step, sos_targets
# these functions are normally not available but can be imported
# using their names for testing purposes
from sos.utils import (WorkflowDict, as_fstring, env, get_logger,
                       split_fstring, stable_repr, fileMD5)
from sos.workflow_executor import Base_Executor, analyze_section


def test_logger():
    '''Test logging level'''
    for verbosity in [0, 1, 2, 3, 4]:
        env.verbosity = verbosity
        get_logger().debug(
            'Verbosity {}:debug message with ``empahsized text`` in between'
            .format(env.verbosity))
        get_logger().info(
            'Verbosity {}:info message with ``empahsized text`` in between'
            .format(env.verbosity))
        get_logger().warning(
            'Verbosity {}:warning message with ``empahsized text`` in between'
            .format(env.verbosity))
        get_logger().error(
            'Verbosity {}:error message with ``empahsized text`` in between'
            .format(env.verbosity))


def test_workflow_dict():
    '''Test workflow dict with attribute access'''
    env.reset()
    d = WorkflowDict()
    d['a'] = 1
    assert d['a'] == 1
    #
    d['a'] += 1
    assert d['a'] == 2


def test_pattern_match():
    '''Test snake match's pattern match facility'''
    res = extract_pattern('{a}-{b}.txt', ['file-1.txt', 'file-ab.txt'])
    assert res['a'] == ['file', 'file']
    assert res['b'] == ['1', 'ab']
    res = extract_pattern('{a}-{b}.txt', ['file--ab--cd.txt'])
    assert res['a'] == ['file--ab-']
    assert res['b'] == ['cd']
    res = extract_pattern('{path}/{to}/{file}.txt', ['/tmp//1.txt'])
    assert res['path'] == [None]
    assert res['to'] == [None]
    assert res['file'] == [None]
    res = extract_pattern('{path}/{to}/{file}.txt', ['/tmp/test/1.txt.txt'])
    assert res['path'] == ['/tmp']
    assert res['to'] == ['test']
    assert res['file'] == ['1.txt']
    # expand_pattern
    env.sos_dict = WorkflowDict({
        'os': os,
        'a': 100,
        'b': 'file name',
        'c': ['file1', 'file2', 'file 3'],
        'd': {
            'a': 'file1',
            'b': 'file2'
        },
    })
    assert expand_pattern('{b}.txt') == ['file name.txt']
    assert expand_pattern('{c}.txt') == ['file1.txt', 'file2.txt', 'file 3.txt']
    assert expand_pattern('{a}_{c}.txt') == [
        '100_file1.txt', '100_file2.txt', '100_file 3.txt'
    ]


def test_accessed_vars():
    '''Test accessed vars of a SoS expression or statement.'''
    assert accessed_vars('''a = 1''') == set()
    assert accessed_vars('''a = b + 2.0''') == {'b'}
    assert accessed_vars('''a = "C"''') == set()
    assert accessed_vars('''a = "C" + f"{D}"''') == {'D'}
    assert accessed_vars('''a = 1 + f"{D + 20:f}" ''') == {'D'}
    assert accessed_vars(
        '''k, "a.txt", "b.txt", par=f(something) ''',
        mode='eva') == {'k', 'f', 'something'}
    # this is a complicated case because the actual variable depends on the
    # result of an expression... However, in the NO-evaluation case, this is
    # the best we can do.
    assert accessed_vars('''c + f"{D + 1}" ''') == {'c', 'D'}
    assert accessed_vars('''a, b=2, c=d ''', mode='eva') == {'a', 'd'}


def test_progress_bar():
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


def test_text_repr():
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
    for text in ('"""a"""', '"b"', r'"""\na\\nb"""', r"'''a\nb'''",
                 """ "a'\\"='" """):
        script = SoS_Script(r'''
a = 1
run: expand=True
#!/usr/bin/env python
with open('tmp.txt', 'w') as tmp:
    tmp.write({} + '{}')
k = """b"""'''.format(text, '{a}'))
        wf = script.workflow()
        Base_Executor(wf).run()
        with open('tmp.txt') as tmp:
            assert tmp.read() == eval(text) + '1'
    os.remove('tmp.txt')


def test_analyze_section():
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
    Base_Executor(wf)
    for section in wf.sections:
        res = analyze_section(section)
        if section.names[0][1] == '1':
            assert res['step_input'].undetermined()
            assert res['step_depends'] == sos_targets()
            assert res['step_output'] == sos_targets()
            assert res['environ_vars'] == {'b', 'p1', 'infiles'}
            assert res['signature_vars'] == {'c'}
            assert res['changed_vars'] == {'b'}
        elif section.names[0][1] == '2':
            assert res['step_input'] == sos_targets()
            assert res['step_depends'] == sos_targets('some.txt',
                                                      executable('ls'))
            assert res['step_output'].unspecified()
            # for_each will not be used for DAG
            assert res['environ_vars'] == {'b', 'for_each', 'executable'}
            assert res['signature_vars'] == {'r', 'time', 'random'}
            assert res['changed_vars'] == set()
        elif section.names[0][1] == '4':
            assert 'output' in res['signature_vars']
        elif section.names[0][1] == '5':
            assert 'output' not in res['signature_vars']


def test_on_demand_options(reset_env):
    '''Test options that are evaluated upon request.'''
    options = on_demand_options({'a': '"est"', 'b': 'c', 'c': 'e + 2'})
    env.sos_dict = WorkflowDict({
        'e': 10,
    })
    assert options['a'] == 'est'
    with pytest.raises(KeyError):
        options.__getitem__('d')
    with pytest.raises(ValueError):
        options.__getitem__('b')
    assert options['c'] == 12
    #
    env.sos_dict.set('e', 20)
    assert options['c'] == 22
    env.sos_dict.set('c', 200)
    assert options['b'] == 200


def test_stable_repr():
    assert stable_repr({1, 2, '3', '1'}) == "{'1', '3', 1, 2}"
    assert stable_repr({1: 2, 3: 4}) == "{1:2, 3:4}"
    assert stable_repr([1, 3, 4]) == "[1, 3, 4]"


def test_split_fstring():
    '''Test function to split f-string in pieces '''
    for string, pieces in [
        ('hello world', ['hello world']),
        ('hello world {', None),
        ('{ hello world', None),
        ('hello world {}', None),
        ('hello world {a b}', None),
        ('{{hello world', ['{{hello world']),
        ('hello {{}} world', ['hello {{}} world']),
        ('hello {{ world }}', ['hello {{ world }}']),
        ('hello }} world', ['hello }} world']),
        ('hello {1} world', ['hello ', '1', ' world']),
        ('hello {a+b } }} world', ['hello ', 'a+b ', ' }} world']),
        ('hello {a+b:r} }} world', ['hello ', 'a+b:r', ' }} world']),
        ('hello {{{a+b!r} }} world', ['hello {{', 'a+b!r', ' }} world']),
        ('hello {a+b + {1,2}.pop() } }} world',
         ['hello ', 'a+b + {1,2}.pop() ', ' }} world']),
    ]:
        if pieces is None:
            with pytest.raises(SyntaxError):
                split_fstring(string)
        else:
            assert split_fstring(string) == pieces


def test_as_fstring():
    '''Test as fstring '''
    for string, fstring in [
        ('hello world', 'fr"""hello world"""'),
        ('{{hello world', 'fr"""{{hello world"""'),
        ('hello {{}} world', 'fr"""hello {{}} world"""'),
        ('hello {{ world }}', 'fr"""hello {{ world }}"""'),
        ('hello }} \n world', 'fr"""hello }} \n world"""'),
        (r'hello }} \n world', 'fr"""hello }} \\n world"""'),
        ('hello {1} world', 'fr"""hello {1} world"""'),
        ('hello {a+b } }} world', 'fr"""hello {a+b } }} world"""'),
        ('hello {a+b:r} }} world', 'fr"""hello {a+b:r} }} world"""'),
        ('''hello """ \'\'\' {a+b } }} world''',
         'f\'hello """ \\\'\\\'\\\' {a+b } }} world\''),
        ('hello {{{a+b!r} }} world', 'fr"""hello {{{a+b!r} }} world"""'),
        ('hello {a+b + {1,2}.pop() } }} world',
         'fr"""hello {a+b + {1,2}.pop() } }} world"""'),
        ('''hello {'a'+b !r} }} world''',
         'fr"""hello {\'a\'+b !r} }} world"""'),
        ('''hello """ \'\'\' {'a'+"b" + {"c", "D"}.pop() } }} world''',
         '\'hello """ \\\'\\\'\\\' {0} }} world\'.format(\'a\'+"b" + {"c", "D"}.pop() )'
        ),
    ]:
        assert as_fstring(string) == fstring


def test_analyze_output_from():
    '''Test extracting of from=value option from input'''
    script = SoS_Script('''
[A_1]
input:  output_from('B')

[A_2]
input: something_unknown, sos_groups(output_from(['C1', 'C2']), by=2), group_by=1
''')
    wf = script.workflow('A')
    Base_Executor(wf)
    for section in wf.sections:
        res = analyze_section(section, analysis_type='forward')
        if section.names[0][1] == 1:
            assert res['step_depends'] == sos_targets(sos_step('B'))
        if section.names[0][1] == 2:
            assert res['step_depends'] == sos_targets(
                sos_step('C1'), sos_step('C2'))


def test_file_sig(clear_now_and_after):
    '''test save and validate of file signature'''
    clear_now_and_after('test_sig.txt', 'test_sig.txt.zapped')

    with open('test_sig.txt', 'w') as ts:
        ts.write('ba')
    a = file_target('test_sig.txt')
    a.write_sig()
    assert a.validate()
    #
    a.zap()
    assert a.validate()

    with open('test_sig.txt', 'w') as ts:
        ts.write('bac')

    assert not a.validate()

@pytest.mark.parametrize('fsize', [12354, 33554432, 34605213])
def test_file_md5(fsize, temp_factory):
    '''test save and validate of file signature'''
    fname = 'test_md5.txt'
    temp_factory(fname, size = fsize)

    partial_md5, full_md5 = fileMD5(fname, sig_type='both')
    assert partial_md5 == fileMD5(fname, sig_type='partial')
    assert full_md5 == fileMD5(fname, sig_type='full')
