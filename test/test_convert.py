#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import subprocess
import textwrap

from sos.converter import extract_workflow


def test_script_to_html(temp_factory):
    '''Test sos show script --html'''
    script1 = textwrap.dedent('''
    [0]
    seq = range(3)
    input: for_each='seq'
    output: "test${_seq}.txt"
    print(output)
    ''')
    script2 = textwrap.dedent('''
    [0]
    seq = range(3)
    input: for_each='seq'
    output: "test${_seq}.txt"
    run:			concurrent=True
        echo 'this is test script'

    [10]
    report('this is action report')
    ''')
    temp_factory('temp1.sos', content=script1)
    temp_factory('temp2.sos', content=script2)

    scripts = ['temp1.sos', 'temp2.sos']
    for script_file in scripts:
        assert subprocess.call(
            f'sos convert {script_file} {script_file}.html', shell=True) == 0
        assert subprocess.call(
            f'sos convert {script_file} {script_file}.html --linenos',
            shell=True) == 0
        #
        assert subprocess.call(['sos', 'convert', script_file, '--to',
                                'html']) == 0


def test_extract_workflow(sample_workflow):
    '''Test extract workflow from ipynb file'''
    content = extract_workflow('sample_workflow.ipynb')
    print(content)
    assert content == textwrap.dedent('''\
    #!/usr/bin/env sos-runner
    #fileformat=SOS1.0

    # this comment will be included but not shown in help message
    # because it is for the global
    [global]
    a = 1
    # this comment will become the comment for parameter b
    parameter: b=2
    parameter: c=3
    # this comment will become the comment for parameter d
    parameter: d='d'

    # this is a section comment, will be displayed
    [default]
    print(f'Hello {a}')

    ''')
