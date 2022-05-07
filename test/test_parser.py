#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import subprocess

import pytest
from sos import execute_workflow
from sos.converter import extract_workflow
from sos.parser import ParsingError, SoS_Script
from sos.targets import file_target, path, paths, sos_targets
from sos.utils import ArgumentError, env
from sos.workflow_executor import Base_Executor

section1_sos = """
#!/usr/bin/env sos-runner
#fileformat=SOS1.1

#
# this is a test sos script
#
var1='value1'
var2 = 'value2'
var3 = [var1,
  var2]

# section
# description for workflow section
#
[parameters]
# par1 string
par1 = 'var1'

# par2 list
par2 = ['a', 'b', 'c']

# par3 multiline
par3 = ['a',
	'b']

[*_0]
var0 = '0'

[section_10]
#
#step 10
var1 = 'a'

[section_2 : shared='var3']
var2 = 'a'
input: var1
output: var2

var3 = 'a'

[section_3, *_4 : shared='var4']
output:
	var2,
	var3

print()
var4= 'value4'

[chapter_5]
var4='5'
 """


def test_file_format():
    """Test recognizing the format of SoS script"""
    # file format must be 'fileformat=SOSx.x'
    with pytest.raises(ParsingError):
        SoS_Script("#fileformat=SS2")
    with pytest.raises(ParsingError):
        SoS_Script("#fileformat=SOS1.0beta")
    #
    # parse a larger script with gormat 1.1
    SoS_Script(section1_sos)
    # not the default value of 1.0
    # assert script.format_version, '1.1')


def test_workflows():
    """Test workflows defined in SoS script"""
    script = SoS_Script("""[0]""")
    assert sorted(script.workflows) == [""]
    script = SoS_Script("""[0]\n[1]""")
    assert sorted(script.workflows) == [""]
    script = SoS_Script("""[0]\n[*_1]""")
    assert sorted(script.workflows) == [""]
    script = SoS_Script("""[0]\n[*_1]\n[auxiliary:provides='{a}.txt']""")
    assert sorted(script.workflows) == ["", "auxiliary"]
    script = SoS_Script("""[0]\n[*_1]\n[human_1]""")
    assert sorted(script.workflows) == ["", "human"]
    script = SoS_Script("""[0]\n[*_1]\n[human_1]\n[mouse_2]""")
    assert sorted(script.workflows) == ["", "human", "mouse"]
    script = SoS_Script("""[0]\n[*_1]\n[human_1]\n[mouse_2]\n[s*_2]""")
    assert sorted(script.workflows) == ["", "human", "mouse"]
    # unnamed
    script = SoS_Script("""[0]\n[*_1]\n[human_1]\n[mouse]\n[s*_2]""")
    assert sorted(script.workflows) == ["", "human", "mouse"]
    #
    # workflow name with -
    script = SoS_Script("""[proc-1]\n[test-case_2]""")
    assert sorted(script.workflows) == ["proc-1", "test-case"]
    script.workflow("proc-1")
    script.workflow("proc-1 + test-case:2")


def test_sections():
    """Test section definitions"""
    # bad names
    for badname in ["56_1", "_a", "a_", "1x", "*", "?"]:
        with pytest.raises(ParsingError):
            SoS_Script("[{}]".format(badname))
    # bad options
    for badoption in ["ss"]:
        with pytest.raises(ParsingError):
            SoS_Script("[0:{}]".format(badoption))
    # allowed names
    for name in ["a5", "a_5", "*_0", "a*1_100"]:
        SoS_Script("[{}]".format(name))
    # allowed names with alias
    for name in ["a5 (p1)", "a_5 (something fun)", "*_0 (no way)", "a*1_100"]:
        SoS_Script("[{}]".format(name))
    # duplicate sections
    with pytest.raises(ParsingError):
        SoS_Script("""[1]\n[1]""")
    with pytest.raises(ParsingError):
        SoS_Script("""[1]\n[3]\n[2,1]""")
    with pytest.raises(ParsingError):
        SoS_Script("""[a_1]\n[a_3]\n[*_1]""")
    #
    # no duplicated section header
    SoS_Script("""[a_1]\n[a_3]\n[b*_1]""")
    #
    # global section
    with pytest.raises(ParsingError):
        SoS_Script("""[global, step_10]""")


def test_parameters():
    """Test ending parameters with new line #1311"""
    with pytest.raises(ParsingError):
        SoS_Script("""input: ['a.txt'\n\n'b.txt']\n""")


def test_global_variables():
    """Test definition of variables"""
    # allow definition
    SoS_Script("""a = '1' """)
    SoS_Script("""a = ['a', 'b'] """)
    # but this one has incorrect syntax
    with pytest.raises(ParsingError):
        SoS_Script("""a = 'b  """)

    SoS_Script('''a = """
this is a multi line
string """
''')
    # multi-line list literal, even with newline in between
    SoS_Script("""a = [
'a',

'b'
]
""")
    #
    SoS_Script(section1_sos)
    # not the default value of 1.0
    #
    script = SoS_Script("""
[global]
a = 1

[b]
print(a)
""")
    wf = script.workflow()
    Base_Executor(wf).run()


def test_parameters_section():
    """Test parameters section"""
    # directive not allowed in parameters
    script = SoS_Script(section1_sos)
    wf = script.workflow("chapter:0")

    #
    script = SoS_Script("""
parameter: a = [1, 2]
[0]
""")
    wf = script.workflow()
    assert list(wf.parameters().keys()) == ["a"]
    Base_Executor(wf).run()
    assert env.sos_dict["a"] == [1, 2]
    env.sos_dict.pop("a")
    wf = script.workflow()
    Base_Executor(wf, args=["--a", "3"]).run()
    assert env.sos_dict["a"] == [3]
    env.sos_dict.pop("a")
    wf = script.workflow()
    Base_Executor(wf, args=["--a", "3", "5"]).run()
    assert env.sos_dict["a"] == [3, 5]
    env.sos_dict.pop("a")
    #
    script = SoS_Script("""
# comment
# comment
parameter: a = ['a.txt', 'b.txt']
[0]
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["a"] == ["a.txt", "b.txt"]
    env.sos_dict.pop("a")
    wf = script.workflow()
    Base_Executor(wf, args=["--a", "3"]).run()
    assert env.sos_dict["a"] == ["3"]
    env.sos_dict.pop("a")
    wf = script.workflow()
    Base_Executor(wf, args=["--a", "3", "5"]).run()
    assert env.sos_dict["a"] == ["3", "5"]
    env.sos_dict.pop("a")
    #
    # test parameter using global definition
    script = SoS_Script("""
a="100"

# comment
parameter: b=str(int(a)+1)
[0]
""")
    wf = script.workflow()
    assert list(wf.parameters().keys()) == ["b"]
    Base_Executor(wf).run()
    assert env.sos_dict["b"] == "101"
    env.sos_dict.pop("b")
    #
    env.sos_dict.clear()
    script = SoS_Script("""
a=100

parameter: b=a+1
[0]
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["b"] == 101
    env.sos_dict.pop("b")
    #
    script = SoS_Script("""
a=100

parameter: b=a+1.
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "1000"]).run()
    #
    assert env.sos_dict["b"] == 1000
    env.sos_dict.pop("b")
    #
    # argument has hve a value
    with pytest.raises(ParsingError):
        SoS_Script("""
[0]
parameter: b=

""")
    # if it is a type, must provide value
    script = SoS_Script("""
# comment
parameter: b = int
[0]
""")
    wf = script.workflow()
    with pytest.raises(ArgumentError):
        Base_Executor(wf)
    #
    script = SoS_Script("""
parameter: b = list
[0]
""")
    assert list(wf.parameters().keys()) == ["b"]
    wf = script.workflow()
    with pytest.raises(ArgumentError):
        Base_Executor(wf)
    # also require the type
    script = SoS_Script("""
parameter: b = int
[0]
""")
    wf = script.workflow()
    with pytest.raises(ArgumentError):
        Base_Executor(wf, args=["--b", "file"])
    #
    script = SoS_Script("""
parameter: b = int
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == 5
    env.sos_dict.pop("b")
    # string
    script = SoS_Script("""
parameter: b = str
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == "5"
    env.sos_dict.pop("b")
    # list is ok
    script = SoS_Script("""
parameter: b = list
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == ["5"]
    # bool required
    script = SoS_Script("""
# comment
parameter: b = bool
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b"]).run(mode="dryrun")
    assert env.sos_dict["b"] == True
    env.sos_dict.pop("b")
    Base_Executor(wf, args=["--no-b"]).run(mode="dryrun")
    assert env.sos_dict["b"] == False
    env.sos_dict.pop("b")
    # bool with default True
    script = SoS_Script("""
parameter: b = True
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=[]).run(mode="dryrun")
    assert env.sos_dict["b"] == True
    env.sos_dict.pop("b")
    Base_Executor(wf, args=["--b"]).run(mode="dryrun")
    assert env.sos_dict["b"] == True
    env.sos_dict.pop("b")
    Base_Executor(wf, args=["--no-b"]).run(mode="dryrun")
    assert env.sos_dict["b"] == False
    env.sos_dict.pop("b")
    # bool with default False
    script = SoS_Script("""
parameter: b = False
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=[]).run(mode="dryrun")
    assert env.sos_dict["b"] == False
    env.sos_dict.pop("b")
    Base_Executor(wf, args=["--b"]).run(mode="dryrun")
    assert env.sos_dict["b"] == True
    env.sos_dict.pop("b")
    Base_Executor(wf, args=["--no-b"]).run(mode="dryrun")
    assert env.sos_dict["b"] == False
    env.sos_dict.pop("b")
    #
    # parameters cannot coincide with a readonly global variable
    # are masked by previous definition
    script = SoS_Script("""
a = 4
parameter: a = 5
[0]
""")
    wf = script.workflow()
    with pytest.raises(Exception):
        Base_Executor(wf, args=["--a", 7])
    # assert env.sos_dict['a'], 4)
    #
    # test parameters with dash
    script = SoS_Script("""
parameter: a_b = 5
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a_b", "10"]).run(mode="dryrun")
    assert env.sos_dict["a_b"] == 10
    env.sos_dict.pop("a_b")
    Base_Executor(wf, args=["--a-b", "10"]).run(mode="dryrun")
    assert env.sos_dict["a_b"] == 10
    env.sos_dict.pop("a_b")
    #
    #
    script = SoS_Script("""
parameter: a_b = int
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a_b", "10"]).run(mode="dryrun")
    assert env.sos_dict["a_b"] == 10
    env.sos_dict.pop("a_b")
    Base_Executor(wf, args=["--a-b", "10"]).run(mode="dryrun")
    assert env.sos_dict["a_b"] == 10
    env.sos_dict.pop("a_b")
    #
    # test support for type path, paths, file_target and sos_targets
    script = SoS_Script("""
parameter: path_var = path
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], path)
    #
    script = SoS_Script("""
parameter: path_var = file_target
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], file_target)
    #
    script = SoS_Script("""
parameter: path_var = paths
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt", "b.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], paths)
    #
    #
    script = SoS_Script("""
parameter: path_var = sos_targets
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], sos_targets)

    #
    script = SoS_Script("""
parameter: path_var = path('a.txt')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], path)
    #
    script = SoS_Script("""
parameter: path_var = file_target('a.txt')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], file_target)
    #
    script = SoS_Script("""
parameter: path_var = paths('a.txt')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt", "b.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], paths)
    #
    #
    script = SoS_Script("""
parameter: path_var = sos_targets('a.txt')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--path-var", "a.txt"]).run(mode="run")
    assert isinstance(env.sos_dict["path_var"], sos_targets)

    #
    # Test allow the use of sos keywords as parameters #1041
    script = SoS_Script("""\
[1]
parameter: input = 5
output = 10
python: expand=True
print({input})
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    # multiple parameters
    script = SoS_Script("""
parameter: a_b = int
[0]
parameter: c_b = list
""")
    wf = script.workflow()
    assert sorted(list(wf.parameters().keys())) == ["a_b", "c_b"]


def test_param_in_task():
    """Test specification of parameters in tasks"""
    with pytest.raises(Exception):
        SoS_Script(
            """\
[1]
output: 'a.txt'
task:
parameter: value = 'a'
sh: expand=True
echo {value} >  a.txt
""",)
    # multiple parameters
    script = SoS_Script("""
parameter: a_b = int
[0]
parameter: c_b = list
""")
    wf = script.workflow()
    assert sorted(list(wf.parameters().keys())) == ["a_b", "c_b"]


def test_type_trait_parameter():
    # type trait
    script = SoS_Script("""
parameter: b
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == "5"
    #
    script = SoS_Script("""
parameter: b :str
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == "5"
    #
    script = SoS_Script("""
parameter: b : list
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == ["5"]

    #
    script = SoS_Script("""
parameter: b : list = 5
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == 5
    #
    script = SoS_Script("""
parameter: b : int
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "5"]).run(mode="dryrun")
    assert env.sos_dict["b"] == 5


def test_input_target():
    # test input of targets
    script = SoS_Script("""
parameter: b : file_target
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "file_target"
    #
    script = SoS_Script("""
parameter: b = file_target('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "file_target"
    #
    script = SoS_Script("""
parameter: b : sos_targets
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa", "bbb"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "sos_targets"
    #
    script = SoS_Script("""
parameter: b = sos_targets('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "sos_targets"
    #
    script = SoS_Script("""
parameter: a_b : file_target
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "file_target"
    #
    script = SoS_Script("""
parameter: a_b = file_target('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "file_target"
    #
    script = SoS_Script("""
parameter: a_b : sos_targets
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa", "bbb"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "sos_targets"
    #
    script = SoS_Script("""
parameter: a_b = sos_targets('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "sos_targets"
    #
    #
    #
    #
    script = SoS_Script("""
parameter: b : path
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "path"
    #
    script = SoS_Script("""
parameter: b = path('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "path"
    #
    script = SoS_Script("""
parameter: b : paths
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa", "bbb"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "paths"
    #
    script = SoS_Script("""
parameter: b = paths('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["b"].__class__.__name__ == "paths"
    #
    script = SoS_Script("""
parameter: a_b : path
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "path"
    #
    script = SoS_Script("""
parameter: a_b = path('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "path"
    #
    script = SoS_Script("""
parameter: a_b : paths
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa", "bbb"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "paths"
    #
    script = SoS_Script("""
parameter: a_b = paths('file')
[0]
""")
    wf = script.workflow()
    Base_Executor(wf, args=["--a-b", "aaa"]).run(mode="dryrun")
    assert env.sos_dict["a_b"].__class__.__name__ == "paths"
    #


def test_section_directives():
    """Test directives of sections"""
    # cannot be in the global section
    # multi-line OK
    SoS_Script("""
[0]
input: 'filename',
'filename1'

""")
    # An abusive case with multi-line OK, from first column ok
    SoS_Script("""
[0]
input: 'filename',
'filename1',
filename4,
opt1=value
output:
blah

depends:
'something else'
""")
    # option with expression ok
    SoS_Script("""
[0]
input: 'filename',  'filename2', opt=value==1

""")
    # unrecognized directive, allowed now
    SoS_Script("""
[0]
something: 'filename',  filename2, opt=value==1
""")
    # need commma
    with pytest.raises(ParsingError):
        SoS_Script("""
[0]
input: 'filename'  filename2
""",)
    # can be after action
    SoS_Script("""
[0]
func()
input: 'filename',  'filename2', opt=value==1
""")
    # assignments between directives are allowed
    SoS_Script("""
[0]
input: 'filename',  'filename2', opt=value==1
a = 'some text'
output: 'filename',  'filename2', opt=value==1
""")
    # can have action between directives
    SoS_Script("""
[0]
input: 'filename',  'filename2', opt=value==1
abc
output: 'filename',  'filename2', opt=value==1
""")


def test_script_format():
    """Test putting scripts directly under action"""
    script = SoS_Script("""
[0]
input: 'filename',  'filename2', opt=value==1
R:

open('something')
save.put()

""")
    script = SoS_Script("""
[0]
input: 'filename',  'filename2', opt=value==1
R: concurrent = True

open('something')
save.put()

""")
    script = SoS_Script("""
[0]
input: 'filename',  'filename2', opt=value==1
R: concurrent = True,
workdir = 'someelse else'

open('something')
save.put()

""")
    # test dedent
    script = SoS_Script("""
[0]
python3:
from sos.runtime import logger
logger.warning('I am from a dented text')
if 1:
    logger.warning('Another warning')
""")
    script.workflow()
    # with section head in the script,
    # this will not work even if the embedded
    # python script is perfectly valid.
    with pytest.raises(ParsingError):
        SoS_Script(
            '''
[0]
input: 'filename',  'filename2', opt=value==1
python3:

with open('something') as e:
e.write("""
[section]
""")

''',)
    # scripts with section head-like lines
    script = SoS_Script("""
[0]
R:
some::function(param)
""")
    script.workflow()
    #
    # script with first-line indent
    #
    script = SoS_Script("""
[0]
sh:
echo "a"

sh('echo "b"')
""")
    script.workflow()
    #
    # script with triple quote and format string
    # #1211
    script = SoS_Script('''
[1]
python3: expand = "${ }"

ld = ${'100'}
a =  """doc"""

''')
    wf = script.workflow()
    Base_Executor(wf).run()
    script = SoS_Script('''
[default]
report: expand = "${ }"
ld_file = ${_input['']:r}
""" \'\'\'
{}, ${_output:r},

''')
    wf = script.workflow()
    Base_Executor(wf).run()


def test_input(temp_factory):
    """Test input directive"""
    temp_factory("a.txt", "b.txt", "a.pdf", "a0", "a1")
    script = SoS_Script("""
[0]
files = ['a.txt', 'b.txt']

input: 'a.pdf', files

""")
    wf = script.workflow()
    Base_Executor(wf).run(mode="dryrun")
    #
    # test input types
    script = SoS_Script("""
[0:shared={'i':'_input', 'o':'_output'}]
files = (f"a{i}" for i in range(2))
input: {'a.txt', 'b.txt'}, files
output: (f"a{x}" for x in _input)

""")
    wf = script.workflow()
    Base_Executor(wf).run(mode="dryrun")
    assert sorted(env.sos_dict["i"]) == sos_targets(
        ["a.txt", "a0", "a1", "b.txt"])
    assert sorted(env.sos_dict["o"]) == sos_targets(
        ["aa.txt", "aa0", "aa1", "ab.txt"])


def test_group_by(temp_factory, clear_now_and_after):
    """Test group_by parameter of step input"""
    clear_now_and_after('a.txt', 'b.txt', 'c.txt' 'xx.txt')
    clear_now_and_after([f"b{x+1}.txt" for x in range(6)])
    # group_by = 'all'
    temp_factory([f"a{x}.txt" for x in range(15)])
    #
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: [f'a{x}.txt' for x in range(1, 5)], group_by='all'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt", "a3.txt", "a4.txt")
    ]
    assert env.sos_dict["executed"][0].labels == ["0"] * 4
    # group_by = 'single'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='single'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt"),
        sos_targets("a2.txt"),
        sos_targets("a3.txt"),
        sos_targets("a4.txt"),
    ]
    # group_by = 'pairs'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: [f'a{x}.txt' for x in range(1, 5)], group_by='pairs'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a3.txt"),
        sos_targets("a2.txt", "a4.txt")
    ]
    # group_by = 'pairs2'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: [f'a{x}.txt' for x in range(1, 9)], group_by='pairs2'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt", "a5.txt", "a6.txt"),
        sos_targets("a3.txt", "a4.txt", "a7.txt", "a8.txt"),
    ]
    # group_by = 'pairs3'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: [f'a{x}.txt' for x in range(1, 13)], group_by='pairs3'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt", "a3.txt", "a7.txt", "a8.txt", "a9.txt"),
        sos_targets("a4.txt", "a5.txt", "a6.txt", "a10.txt", "a11.txt",
                    "a12.txt"),
    ]

    # group_by = 'pairwise'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='pairwise'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt"),
        sos_targets("a2.txt", "a3.txt"),
        sos_targets("a3.txt", "a4.txt"),
    ]
    # group_by = 'pairwiseN'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 7)], group_by='pairwise2'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt", "a3.txt", "a4.txt"),
        sos_targets("a3.txt", "a4.txt", "a5.txt", "a6.txt"),
    ], f'obtained {env.sos_dict["executed"]}'

    # group_by = 'combinations'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='combinations'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"], [
        sos_targets("a1.txt", "a2.txt"),
        sos_targets("a1.txt", "a3.txt"),
        sos_targets("a1.txt", "a4.txt"),
        sos_targets("a2.txt", "a3.txt"),
        sos_targets("a2.txt", "a4.txt"),
        sos_targets("a3.txt", "a4.txt"),
    ]

    # group_by = 'combinations3'
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 5)], group_by='combinations3'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets(["a1.txt", "a2.txt", "a3.txt"]),
        sos_targets(["a1.txt", "a2.txt", "a4.txt"]),
        sos_targets(["a1.txt", "a3.txt", "a4.txt"]),
        sos_targets(["a2.txt", "a3.txt", "a4.txt"]),
    ], f'obtained {env.sos_dict["executed"]}'
    # group_by chunks specified as integers
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=3

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt", "a3.txt"),
        sos_targets("a4.txt", "a5.txt", "a6.txt"),
        sos_targets("a7.txt", "a8.txt", "a9.txt"),
    ]
    # group_by chunks specified as integer strings
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='3'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "a2.txt", "a3.txt"),
        sos_targets("a4.txt", "a5.txt", "a6.txt"),
        sos_targets("a7.txt", "a8.txt", "a9.txt"),
    ]
    # number of files should be divisible by group_by
    temp_factory(["a{}.txt".format(x) for x in range(1, 10)])
    execute_workflow(
        """
        [0]

        executed = []
        input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by=4

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    # incorrect value causes an exception
    with pytest.raises(Exception):
        execute_workflow(
            """
            [0]

            executed = []
            input: ['a{}.txt'.format(x) for x in range(1, 10)], group_by='something'

            executed.append(_input)

            """,
            options={'run_mode': "dryrun"})

    #
    # group by label
    file_target("c.txt").touch()
    execute_workflow(
        """
        [A]
        output: 'a.txt'
        _output.touch()

        [B]
        input: for_each={'i': range(2)}
        output: 'b.txt', 'b1.txt', group_by=1
        _output.touch()

        [0: shared='executed']
        executed = []

        input: 'c.txt', output_from(['A', 'B']), group_by='label'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("c.txt"),
        sos_targets("a.txt"),
        sos_targets("b.txt", "b1.txt"),
    ]
    #
    # group_by='pairlabel'
    file_target("c.txt").touch()
    execute_workflow(
        """
        [A]
        output: 'a1.txt', 'a2.txt'
        _output.touch()

        [B]
        output: 'b1.txt', 'b2.txt'
        _output.touch()

        [0: shared='executed']
        executed = []

        input: output_from(['A', 'B']), group_by='pairlabel'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a1.txt", "b1.txt"),
        sos_targets("a2.txt", "b2.txt")
    ]

    # group_by='pairlabel3'
    temp_factory(["c{}.txt".format(x) for x in range(1, 7)])

    execute_workflow(
        """
        [A]
        output: [f'a{x}.txt' for x in range(1, 7)]
        _output.touch()

        [B]
        output: [f'b{x}.txt' for x in range(1, 7)]
        _output.touch()

        [0: shared='executed']
        executed = []

        input: [f'c{x}.txt' for x in range(1, 7)], output_from(['A', 'B']), group_by='pairlabel3'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets(
            "c1.txt",
            "c2.txt",
            "c3.txt",
            "a1.txt",
            "a2.txt",
            "a3.txt",
            "b1.txt",
            "b2.txt",
            "b3.txt",
        ),
        sos_targets(
            "c4.txt",
            "c5.txt",
            "c6.txt",
            "a4.txt",
            "a5.txt",
            "a6.txt",
            "b4.txt",
            "b5.txt",
            "b6.txt",
        ),
    ]
    # group_by='pairlabel3'
    temp_factory(["c{}.txt".format(x) for x in range(1, 7)])

    execute_workflow(
        """
        [A]
        output: 'a1.txt'
        _output.touch()

        [B]
        output: [f'b{x}.txt' for x in range(1, 7)]
        _output.touch()

        [0: shared='executed']
        executed = []

        input: [f'c{x}.txt' for x in range(1, 3)], output_from(['A', 'B']), group_by='pairlabel3'

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("c1.txt", "a1.txt", "b1.txt", "b2.txt", "b3.txt"),
        sos_targets("c2.txt", "a1.txt", "b4.txt", "b5.txt", "b6.txt"),
    ]
    #
    # group by function
    file_target("c.txt").touch()
    execute_workflow(
        """
        [0: shared='executed']
        executed = []

        def grp(x):
            return  [x[:3], x[3:]]

        input: ['a{}.txt'.format(x) for x in range(5)], group_by=grp

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [["a0.txt", "a1.txt", "a2.txt"],
                                        ["a3.txt", "a4.txt"]]

    #
    # group by lambda function
    file_target("c.txt").touch()
    execute_workflow(
        """
        [0: shared='executed']
        executed = []

        input: ['a{}.txt'.format(x) for x in range(6)], group_by=lambda x: zip(x[:3], x[3:])

        executed.append(_input)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [["a0.txt", "a3.txt"],
                                        ["a1.txt", "a4.txt"],
                                        ["a2.txt", "a5.txt"]]


def test_output_group_by(temp_factory):
    """Test group_by parameter of step output"""
    # group_by = 'all'
    temp_factory(["a{}.txt".format(x) for x in range(4)])
    #
    execute_workflow(
        """
        [0: shared='executed']

        executed = []
        input: ['a{}.txt'.format(x) for x in range(4)], group_by=2
        output: ['a{}.txt.bak'.format(x) for x in range(4)], group_by=2

        executed.append(_output)

        """,
        options={'run_mode': "dryrun"})
    assert env.sos_dict["executed"] == [
        sos_targets("a0.txt.bak", "a1.txt.bak"),
        sos_targets("a2.txt.bak", "a3.txt.bak"),
    ]


def test_steps_with_step_name(clear_now_and_after):
    """Test from steps"""
    clear_now_and_after('a.txt')
    execute_workflow("""
        [step_10]

        output: 'a.txt'
        _output.touch()

        [step_20]
        input: output_from(step_name.split('_')[0] + '_10')
        print(_input)
        """)
    #
    execute_workflow("""
        [step_10]

        output: 'a.txt'
        _output.touch()

        [step_20]
        input: output_from(10)
        print(_input)
        """)


def test_section_actions():
    """Test actions of sections"""
    SoS_Script("""
[0]
func('''
multiline
string''', with_option=1
)
""")
    with pytest.raises(ParsingError):
        SoS_Script("""
[0]
func(
""",)


def test_longer_code():
    """Test definition of classes (with intermediate newlines) in step."""
    execute_workflow("""
        # first block
        [0: shared='b']
        class A:
            def __init__(self):
                pass

            # the newline above should be fine because SoS treat this as
            # regular lines
            def __call__(self):
                return 0

        b = A()()

        """)
    assert env.sos_dict["b"] == 0


def test_combined_workflow():
    """Test the creation and execution of combined workfow"""
    script = """
        a0 = 0
        if 'executed' not in locals():
            executed = []
        parameter: a = a0 + 1
        [a_1: shared='executed']
        executed.append(step_name)
        [a_2: shared='executed']
        executed.append(step_name)
        [a_3: shared='executed']
        executed.append(step_name)
        [a_4: shared='executed']
        executed.append(step_name)
        output: 'out_a_4'
        [b_1: shared=['executed', 'input_b1']]
        executed.append(step_name)
        input_b1 = _input
        [b_2: shared='executed']
        executed.append(step_name)
        [b_3: shared='executed']
        executed.append(step_name)
        [b_4: shared='executed']
        executed.append(step_name)
        [c: shared='executed']
        executed.append(step_name)
        [d: shared='executed']
        executed.append(step_name)
        """
    execute_workflow(script, workflow="a+b", options={'run_mode': 'dryrun'})
    assert env.sos_dict["executed"] == [
        "a_1", "a_2", "a_3", "a_4", "b_1", "b_2", "b_3", "b_4"
    ]

    assert env.sos_dict["a"] == 1
    assert env.sos_dict["input_b1"] == ["out_a_4"]
    #
    env.sos_dict.pop("executed", None)
    execute_workflow(
        script, workflow="a: 1-2 + a:4 + b:3-", options={'run_mode': 'dryrun'})
    assert env.sos_dict["executed"] == ["a_1", "a_2", "a_4", "b_3", "b_4"]
    #
    env.sos_dict.pop("executed", None)
    execute_workflow(script, workflow="a+c+d", options={'run_mode': 'dryrun'})
    assert env.sos_dict["executed"] == ["a_1", "a_2", "a_3", "a_4", "c", "d"]


def test_yaml_config(clear_now_and_after):
    """Test config file in yaml format"""
    clear_now_and_after('myconfig.yml', 'config.sos')
    with open("myconfig.yml", "w") as config:
        config.write("""
# Lines beginning with # are skipped when the JSON is parsed, so we can
# put comments into our JSON configuration files
{
    StoreOwner : "John Doe",

    # List of items that we sell
    Fruits: [ "apple", "banana", "pear" ],
    Price: 1.05
}
""")
    with open("config.sos", "w") as sos:
        sos.write("""
[0]
print(CONFIG['StoreOwner'])
print(CONFIG.get('StoreOwner', 'something'))
print(CONFIG.get('StoreOwnerSpouse', 'someone else'))
#print(CONFIG.StoreOwner)
""")
    # run the command
    assert subprocess.call(
        "sos run config.sos -c myconfig.yml",
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        shell=True,
    ) == 0
    # now test the value
    script = SoS_Script(filename="config.sos")
    wf = script.workflow()
    Base_Executor(wf, config={"config_file": "myconfig.yml"}).run()
    assert env.sos_dict["CONFIG"]["Price"] == 1.05
    assert env.sos_dict["CONFIG"]["StoreOwner"] == "John Doe"
    assert env.sos_dict["CONFIG"]["Fruits"] == ["apple", "banana", "pear"]


def test_var_output():
    """Test early appearance of variable output"""
    execute_workflow(
        """
        [0]
        seq = range(3)
        input: for_each='seq'
        output: f"test{_seq}.txt"
        print(_output)
        """,
        options={'run_mode': "dryrun"})


def test_cell():
    """Test ignoring %cell"""
    SoS_Script("""
%cell 1
[step ]
a = 1
""")


def test_overwrite_keyword(clear_now_and_after):
    """Test overwrite sos keyword with user defined one."""
    clear_now_and_after("a.txt")
    #
    execute_workflow("""
def run(script):
    pass

[1]
run:
    touch a.txt
""")
    assert not os.path.isfile("a.txt")
    #
    execute_workflow("""
parameter: run = 5

[1]
run:
    touch a.txt
""")


def test_help_message(sample_workflow):
    """Test help message from ipynb notebook"""
    msg = subprocess.check_output(
        f"sos run {sample_workflow} -h", shell=True).decode()
    assert "this comment will be included but not shown in help" not in msg
    assert msg.count(
        "this comment will become the comment for parameter b") == 1
    assert msg.count(
        "this comment will become the comment for parameter d") == 1
    assert msg.count("--c 3 (as int)") == 1
    assert msg.count("this is a section comment, will be displayed") == 1
    assert "this is a test workflow" not in msg
    assert "this is a cell with another kernel" not in msg
    assert "this comment will not be included in exported workflow" not in msg


def test_help_on_multi_workflow(script_factory):
    """Test help message of sos file (#985)"""
    test_sos = script_factory("""\
        [workflow_a_10,workflow_b]
        [workflow_a_20]
        [default]
        """)
    msg = subprocess.check_output(
        f"sos run {test_sos} -h", shell=True).decode()
    assert "workflow_a" in msg
    assert "workflow_b" in msg
    # when introducing sections
    assert "workflow_a_10, workflow_b" in msg
    assert "default" in msg


def test_parameter_abbreviation(clear_now_and_after):
    """Test potential problem caused by parameter abbreviation #1053"""
    clear_now_and_after("0914.txt")
    execute_workflow(
        """
        [global]
        parameter: name = '0914'

        [1]
        parameter: n = 4
        output: f'{name}.txt'
        print(_output)
        _output.touch()
        """,
        args=["--n", "5"])
    assert os.path.isfile("0914.txt")


def test_named_input(temp_factory, clear_now_and_after):
    """Test named input"""
    clear_now_and_after('c.txt')
    temp_factory('a.txt', content='a.txt' + '\n')
    temp_factory('b.txt', content='b.txt' + '\n')
    execute_workflow("""
        [1]
        input: {'A': 'a.txt', 'B': 'b.txt'}, group_by='pairlabel'
        output: 'c.txt'
        assert _input['A'] == [file_target('a.txt')]
        assert _input['B'] == [file_target('b.txt')]
        with open(_output, 'w') as out:
            out.write(open(_input['A']).read())
            out.write(open(_input['B']).read())
        """)
    assert open("c.txt").read() == "a.txt\nb.txt\n"


def test_named_output_in_depends(clear_now_and_after):
    """Test named_output in depends statement"""
    clear_now_and_after('a.txt')
    execute_workflow("""
        [A]
        output: A='a.txt'
        _output.touch()

        [default]
        depends: named_output('A')
        """)


def test_output_from_in_depends(clear_now_and_after):
    """Test output_from in depends statement"""
    clear_now_and_after('a.txt')
    execute_workflow("""
        [A]
        output: A='a.txt'
        _output.touch()

        [default]
        depends: output_from('A')
        """)


def test_named_output_in_output(clear_now_and_after):
    """Test named_output in output statement"""
    clear_now_and_after('a.txt')
    with pytest.raises(Exception):
        execute_workflow("""
            [A]
            output: Aa='a.txt'
            _output.touch()

            [default]
            output: named_output('Aa')
            """)


def test_output_from_in_output():
    """Test output_from in output statement"""
    with pytest.raises(Exception):
        execute_workflow("""
            [Aa]
            output: A='a.txt'
            _output.touch()

            [default]
            output: output_from('Aa')
            """)


def test_sos_step_in_input():
    """Test sos_step in input statement"""
    with pytest.raises(Exception):
        execute_workflow("""
            [A]
            output: A='a.txt'
            _output.touch()

            [default]
            input: sos_step('A')
            """)


def test_sos_variable_in_input():
    """Test sos_variable in input statement"""
    with pytest.raises(Exception):
        execute_workflow("""
            [A]
            a=1

            [default]
            input: sos_variable('a')
            """)


def test_sos_variable_in_output():
    """Test sos_variable in output statement"""

    with pytest.raises(Exception):
        execute_workflow("""
            [A]
            a = 1

            [default]
            output: sos_variable('a')
            """)


def test_sos_variable_with_keywordargument():
    """Test output_from in output statement"""

    with pytest.raises(Exception):
        execute_workflow("""
            [A]
            a = 1

            [default]
            depends: sos_variable(var='a')
            """)


def test_wide_card_step_name():
    """test resolving step name with *"""
    execute_workflow("""
        [A_1]

        [*_2]
        assert step_name == 'A_2', f'step_name is {step_name}, A_2 expected'
        """)


def test_outfrom_prev_step(clear_now_and_after):
    """Test output_from(-1) from output_from """
    clear_now_and_after('A_1.txt')
    execute_workflow("""
        [A_1]
        output: 'A_1.txt'
        _output.touch()

        [A_2]
        input: output_from(-1)

        [default]
        depends: sos_step('A_2')

        """)


def test_step_from_numeric_step():
    """Test sos_step(index) #1209"""
    execute_workflow("""
[1]

[2]
depends: sos_step('1')
""")

    #
    execute_workflow("""
[1]

[2]
depends: sos_step(1)

""")


def test_depends_on_step_with_unspecified_input(clear_now_and_after):
    clear_now_and_after("A_1.txt", "A_2.txt", "A_3.txt", "A_4.txt")
    #
    execute_workflow("""
        [A_1]
        output: f'{step_name}.txt'
        _output.touch()

        [A_2]
        output: f'{step_name}.txt'
        _output.touch()

        [A_3]
        output: f'{step_name}.txt'
        _output.touch()

        [A_4]
        output: f'{step_name}.txt'
        _output.touch()

        [default]
        depends: sos_step('A_2')
        """)
    assert os.path.isfile("A_1.txt")
    assert os.path.isfile("A_2.txt")
    assert not os.path.isfile("A_3.txt")
    assert not os.path.isfile("A_4.txt")
    #
    for file in ("A_1.txt", "A_2.txt", "A_3.txt", "A_4.txt"):
        if os.path.isfile(file):
            os.remove(file)
    #
    execute_workflow("""
        [A_1]
        output: f'{step_name}.txt'
        _output.touch()

        [A_2]
        output: f'{step_name}.txt'
        _output.touch()

        [A_3]
        output: f'{step_name}.txt'
        _output.touch()

        [A_4]
        output: f'{step_name}.txt'
        _output.touch()

        [default]
        input: output_from('A_2')
        """)
    assert os.path.isfile("A_1.txt")
    assert os.path.isfile("A_2.txt")
    assert not os.path.isfile("A_3.txt")
    assert not os.path.isfile("A_4.txt")
    #
    for file in ("A_1.txt", "A_2.txt", "A_3.txt", "A_4.txt"):
        if os.path.isfile(file):
            os.remove(file)
    #
    execute_workflow("""
        [A_1]
        output: f'{step_name}.txt'
        _output.touch()

        [A_2]
        input: None
        output: f'{step_name}.txt'
        _output.touch()

        [A_3]
        output: f'{step_name}.txt'
        _output.touch()

        [A_4]
        output: f'{step_name}.txt'
        _output.touch()

        [default]
        input: output_from('A_2')
        """)
    assert not os.path.isfile("A_1.txt")
    assert os.path.isfile("A_2.txt")
    assert not os.path.isfile("A_3.txt")
    assert not os.path.isfile("A_4.txt")


def test_output_from_workflow(clear_now_and_after):
    """Test output from workflow"""
    #
    #
    clear_now_and_after("A_1.txt", "A_2.txt")
    execute_workflow(r"""
        [A_1]
        output: f'{step_name}.txt'

        with open(_output, 'w') as out:
            out.write(step_name)
            out.write('\n')

        [A_2]
        output: f'{step_name}.txt'

        with open(_output, 'w') as out:
            out.write(step_name)
            out.write('\n')

        [default]
        input: output_from('A')

        with open(_input) as content:
            assert content.read() == 'A_2\n'

        """)
    # this should execute A_2 and A_1
    assert os.path.isfile("A_1.txt")
    assert os.path.isfile("A_2.txt")


def test_execute_global_section():
    """Global section should be executed only once #1219"""
    execute_workflow(r"""
        [global]
        parameter: A=[1, 2]
        parameter: B=[]
        B.extend(A)

        [default]
        assert B == [1, 2]
        input: for_each=dict(i=range(2))
        assert B == [1, 2]
        """)


def test_para_from_nested_workflow():
    """Test passing arguments to nested workflow #1229"""
    execute_workflow(
        r"""
        [A]
        parameter: param = str
        print(param)

        [B]
        print(step_name)

        [default]
        sos_run('B+A')
        """,
        args=["--param", "2"])


def test_indented_action():
    """Test the use of indented action in script format"""
    #
    # not indented
    execute_workflow(r"""
        [1]
        if True:
            sh:
        echo something

        [2]
        for i in range(2):
            sh:
        echo something
        """)
    # not indented double action
    with pytest.raises(Exception):
        execute_workflow(r"""
            [1]
            if True:
                sh:
            echo something

                python:
            print(1)
        """)
    # not indented Python structure
    with pytest.raises(Exception):
        execute_workflow(r"""
            [1]
            if True:
                sh:
            echo something
            else:
                python:
            print(1)
            """)
    #  indented
    execute_workflow(r"""
[1]
if True:
    sh:
        echo something

[2]
for i in range(2):
    sh:
        echo something
""")

    # indented double action
    execute_workflow(r"""
[1]
if True:
    sh:
        echo something
    python:
        print(1)
""")
    # indented Python structure
    execute_workflow(r"""
[1]
if True:
    sh:
        echo something
else:
    python:
        print(1)
""")
    # indented, nested structure
    execute_workflow(r"""
[1]
for i in range(2):
    if True:
        sh:
            echo something
    else:
        python:
            print(1)
""")
    #
    # wrong nested action
    execute_workflow(r"""
report:
    name: 'a.txt'
""")


def test_task_param_var():
    """Test global parameter passed to task parameters #1281"""
    execute_workflow(
        r"""
        [global]
        parameter: job_size = 60

        [1]
        task: trunk_size = job_size
        bash:
                echo 1
        """,
        options={"default_queue": "localhost"},
    )


def test_task_param_var_to_substep():
    """Test global parameter passed to task parameters in substep #1281"""
    execute_workflow(
        r"""
        [global]
        parameter: job_size = 60

        [1]
        input: for_each=dict(i=range(3))

        task: trunk_size = job_size
        bash: expand=True
                echo {i}
        """,
        options={"default_queue": "localhost"},
    )


def test_empty_parameter():
    # parameter: without content #1283
    execute_workflow(r"""
        parameter:

        print(1)
        """)


def test_sequential_substeps():
    """Test sequential execution of substeps"""
    execute_workflow(r"""
        [10: shared='sum']
        sum = 0
        input: for_each=dict(i=range(4)), concurrent=False
        sum += i
        print(f'sum is {sum} at index {_index}')
        """)
    assert env.sos_dict["sum"] == 6


def test_limited_concurrency():
    """Set concurrent=INT"""
    execute_workflow(r"""
        [10]
        input: for_each=dict(i=range(6)), concurrent=2
        print(i)
        """)


def test_concurrent_substep_with_step_import():
    """ Test concurrent substep with step leval import statement #1354"""
    execute_workflow("""
        import time
        input: for_each=dict(i=range(2))
        time.sleep(0)
        """)


def test_type_hint():
    """We should be able to differentiate type hint and sos action"""
    SoS_Script("""a : list = 5""")
    SoS_Script("""a : list""")
    # action
    SoS_Script("""a : input='filename' """)
    # action
    SoS_Script("""a : expand='${ }' """)


def test_sos_step_in_output():
    """Test sos_step in output statement"""
    with pytest.raises(Exception):
        execute_workflow("""
            [A]
            output: A='a.txt'
            _output.touch()

            [default]
            output: sos_step('A')
        """)


def test_comments(sample_workflow):
    """Test the use of comments in sos script"""
    # extract workflow from ipynb
    wf = extract_workflow(sample_workflow)
    assert "this is a test workflow" not in wf
    assert wf.count("this comment will be included but not shown in help") == 1
    assert wf.count("this comment will become the comment for parameter b") == 1
    assert wf.count("this comment will become the comment for parameter d") == 1
    assert "this is a cell with another kernel" not in wf
    assert "this comment will not be included in exported workflow" not in wf
