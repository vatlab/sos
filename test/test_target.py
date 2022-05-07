#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import copy
import os
import sys

import pytest
from sos import execute_workflow
from sos.eval import interpolate
from sos.parser import SoS_Script
from sos.targets import file_target, path, paths, sos_step, sos_targets
from sos.utils import env
from sos.workflow_executor import Base_Executor


def test_target_label():
    """Test labels of sos_targets"""
    a = sos_targets("a")
    assert a.labels == [""]
    b = sos_targets(["a", "b"])
    assert b.labels == ["", ""]
    c = sos_targets(["a1", "b1"], _source="here")
    assert c.labels == ["here", "here"]
    c.extend(b)
    assert c.labels == ["here", "here", "", ""]
    #
    assert c.select("").labels == ["", ""]
    assert c.select("here").labels == ["here", "here"]
    assert c["here"], [file_target("a1") == file_target("b1")]
    assert isinstance(c["here"], sos_targets)
    #
    # function item
    assert isinstance(c.select(1), sos_targets)
    assert c.select(1).labels == ["here"]
    assert c.select(1) == ["b1"]
    #
    # test slice of groups
    res = sos_targets(a=["a.txt", "b.txt"], b=["c.txt", "d.txt"], group_by=1)
    assert len(res.groups) == 4
    assert res.labels, ["a", "a", "b", "b"]
    res_a = res["a"]
    assert len(res_a) == 2
    assert res_a.labels == ["a", "a"]
    assert len(res_a.groups) == 4
    assert len(res_a.groups[0]) == 1
    assert len(res_a.groups[1]) == 1
    assert len(res_a.groups[2]) == 0
    assert len(res_a.groups[3]) == 0


def test_remove_targets():
    """Test sos_target.remove_targets()"""
    a = sos_targets(sos_step("1"), "a.txt")._group(by=1)
    a.remove_targets(type=sos_step)
    assert len(a) == 1
    assert len(a.groups) == 2
    assert len(a._groups[0]._indexes) == 0
    assert len(a._groups[0]._labels) == 0
    assert a._groups[1]._indexes == [0]
    assert len(a._groups[1]._labels) == 1


def test_sos_targets_signature(clear_now_and_after):
    """Test save and validate signatures of sos_targets"""
    clear_now_and_after('a.txt', 'b.txt')
    with open("a.txt", "w") as a:
        a.write("text1")
    with open("b.txt", "w") as b:
        b.write("text2")
    t = sos_targets("a.txt", "b.txt")
    sig = t.target_signature()
    assert t.validate(sig)
    # variables does not affect signature
    t[0].set("a", 2)
    assert sig == t.target_signature()
    #
    t.set("cc", "another string")
    assert t.validate(sig)


def test_target_set_get():
    """Test set and get attributes from targets"""
    a = file_target("a")
    a.set(b=1)
    assert a.b == 1
    assert a.get("b") == 1
    with pytest.raises(Exception):
        a.set("touch", 1)


def test_target_group_by():
    """Test new option group_by to sos_targets"""
    res = sos_targets(
        "e.txt",
        "f.ext",
        a=["a.txt", "b.txt"],
        b=["c.txt", "d.txt"],
        group_by=1)
    assert len(res.groups) == 6
    assert res.labels == ["", "", "a", "a", "b", "b"]
    #
    res = sos_targets(res, group_by=2)
    assert len(res.groups) == 3


def test_target_paired_with():
    """Test paired_with targets with vars"""
    res = sos_targets(
        "e.txt",
        "f.ext",
        a=["a.txt", "b.txt"],
        b=["c.txt", "d.txt"],
        group_by=1).paired_with("_name", ["e", "f", "a", "b", "c", "d"])
    for i, n in enumerate(["e", "f", "a", "b", "c", "d"]):
        assert res[i]._name == n
    #
    res = copy.deepcopy(res)
    for i, n in enumerate(["e", "f", "a", "b", "c", "d"]):
        assert res[i]._name == n
    #
    # test assert for length difference
    with pytest.raises(Exception):
        sos_targets("e.txt",
                    "f.ext").paired_with("name", ["e", "f", "a", "b", "c", "d"])


def test_target_group_with():
    """Test group_with targets with vars"""
    res = sos_targets(
        "e.txt",
        "f.ext",
        a=["a.txt", "b.txt"],
        b=["c.txt", "d.txt"],
        group_by=2).group_with("name", ["a1", "a2", "a3"])
    for i, n in enumerate(["a1", "a2", "a3"]):
        assert res.groups[i].name == n
    #
    res = copy.deepcopy(res)
    for i, n in enumerate(["a1", "a2", "a3"]):
        assert res.groups[i].name == n
    #
    # test assert for length difference
    with pytest.raises(Exception):

        sos_targets(
            "e.txt", "f.ext", group_by=1).group_with(
                "name",
                ["e", "f", "g"],
            )


def test_group_with_with_no_output():
    execute_workflow(r"""
[10]
input: for_each=dict(i=range(3))
output: group_with=dict(var=['A', 'B', 'C'][_index])
print(i)

[20]
print(f'Input is {_input} {var}')
""")


def test_merging_of_sos_targets():
    """Test merging of multiple sos targets"""
    # merge 0 to 0
    res = sos_targets("a.txt", "b.txt",
                      sos_targets("c.txt", "d.txt", group_by=1))
    assert len(res) == 4
    assert len(res.groups) == 2
    assert res.groups[0] == ["a.txt", "b.txt", "c.txt"]
    assert res.groups[1] == ["a.txt", "b.txt", "d.txt"]
    # merge N to N
    N1 = sos_targets("c.txt", "d.txt", group_by=1)
    N2 = sos_targets("a1.txt", "a2.txt", "a3.txt", "a4.txt", group_by=2)
    res = sos_targets(N1, N2)
    assert len(res) == 6
    assert len(res.groups) == 2
    assert res.groups[0] == ["c.txt", "a1.txt", "a2.txt"]
    assert res.groups[1] == ["d.txt", "a3.txt", "a4.txt"]
    # test N to M
    N1 = sos_targets("c.txt", "d.txt", group_by=1)
    N2 = sos_targets("a1.txt", "a2.txt", "a3.txt", "a4.txt", group_by=1)

    with pytest.raises(Exception):
        sos_targets(N1, N2)
    # merge 1 to N
    N1 = sos_targets("c.txt", "d.txt", group_by="all")
    N2 = sos_targets("a1.txt", "a2.txt", "a3.txt", "a4.txt", group_by=2)
    res = sos_targets(N1, N2)
    assert len(res) == 6
    assert len(res.groups) == 2
    assert res.groups[0] == ["c.txt", "d.txt", "a1.txt", "a2.txt"]
    assert res.groups[1] == ["c.txt", "d.txt", "a3.txt", "a4.txt"]
    # merge N to 1
    res = sos_targets(N2, N1)
    assert len(res) == 6
    assert len(res.groups) == 2
    assert res.groups[0] == ["a1.txt", "a2.txt", "c.txt", "d.txt"]
    assert res.groups[1] == ["a3.txt", "a4.txt", "c.txt", "d.txt"]


def test_target_format():
    """Test string interpolation of targets"""
    for target, fmt, res in [
        ("a.txt", "", "a.txt"),
        (sos_targets("a.txt"), "", "a.txt"),
        (sos_targets(["a.txt"]), "", "a.txt"),
        (sos_targets([r"c:\path\a.txt"]), "p", "/c/path/a.txt"),
        (sos_targets([r"c:path\a.txt"]), "p", "/c/path/a.txt"),
        (sos_targets("/a/b/a.txt"), "b", "a.txt"),
        (sos_targets("a b.txt"), "q", ("'a b.txt'", '"a b.txt"')),
        (sos_targets("a b.txt"), "x", ".txt"),
    ]:
        if isinstance(res, str):
            assert interpolate(
                "{{target:{}}}".format(fmt), globals(),
                locals()) == res, "Interpolation of {}:{} should be {}".format(
                    target, fmt, res)

        else:
            assert interpolate("{{target:{}}}".format(fmt), globals(), locals(
            )) in res, "Interpolation of {}:{} should be one of {}".format(
                target, fmt, res)


def test_iter_targets():
    """Test iterator interface of targets"""
    t = sos_targets("1", "2", ["3", "4"])
    assert len(t) == 4
    for idx, i in enumerate(t):
        assert str(i) == str(idx + 1)


def test_expand_wildcard():
    """test wildcard expansion of sos_targets"""
    a = sos_targets("*.py")
    assert len(a) > 1
    #
    a = paths("*.py")
    assert len(a) > 1


def test_shared():
    """Test option shared"""
    script = SoS_Script(r"""
parameter: res = 1

[0]
res = 2

[1]
res = 3
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["res"] == 1
    #
    env.sos_dict.pop("res", None)
    script = SoS_Script(r"""
parameter: res = 1

[0: shared='res']
res = 2

[1]
res = 3
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["res"] == 2
    #
    env.sos_dict.pop("res", None)
    script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[0: shared='a']
res = 2

[1: shared='res']
res = 3
a = 5

""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["res"] == 3
    assert env.sos_dict["a"] == 30
    # test multiple vars
    env.sos_dict.pop("res", None)
    script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=('res', 'a')]
res = 3
a = 5

""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["res"] == 3
    assert env.sos_dict["a"] == 5
    #
    # test expression
    env.sos_dict.pop("res", None)
    script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared={'res': 'res + 6', 'c': 'a'}]
res = 3
a = 5

""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["res"] == 9
    assert env.sos_dict["c"] == 5
    # test mixed vars and mapping
    env.sos_dict.pop("res", None)
    script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=['res', {'c': 'a'}]]
res = 3
a = 5

""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["res"] == 3
    assert env.sos_dict["c"] == 5
    # test the step_ version of variables
    script = SoS_Script(r"""
parameter: res = 1
parameter: a = 30

[1: shared=['res', {'c': 'sum(step_a)'}]]
input: for_each={'i': range(10)}
a = _index**2

""")
    wf = script.workflow()
    Base_Executor(wf).run()
    assert env.sos_dict["c"] == sum(x**2 for x in range(10))


#    def testSectionOptionWorkdir():
#        '''Test section option workdir'''
#        script = SoS_Script(r"""
#
# [1: workdir='tmp']
# run:
#    touch 'a.txt'
# """)
#        wf = script.workflow()
#        Base_Executor(wf).run()
#        assert os.path.isfile('tmp/a.txt'))
#        shutil.rmtree('tmp')


def test_depends_executable():
    """Testing target executable."""
    script = SoS_Script("""
[0]
depends: executable('python --version', '3.')
run:
touch a.txt
""")
    wf = script.workflow()
    if file_target("a.txt").exists():
        file_target("a.txt").unlink()
    Base_Executor(wf).run()
    assert os.path.isfile("a.txt")
    file_target("a.txt").unlink()


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Windows executable cannot be created with chmod.")
def test_output_executable(clear_now_and_after):
    """Testing target executable."""
    # change $PATH so that lls can be found at the current
    # directory.
    clear_now_and_after('lls')
    os.environ["PATH"] += os.pathsep + "."
    script = SoS_Script("""
[0]
output: executable('lls')
run:
touch lls
chmod +x lls
""")
    wf = script.workflow()
    if file_target("lls").exists():
        file_target("lls").unlink()
    env.config["sig_mode"] = "force"
    Base_Executor(wf).run()
    # test validation
    env.config["sig_mode"] = "default"
    Base_Executor(wf).run()
    if file_target("ls").exists():
        file_target("lls").unlink()


def test_depends_env_variable():
    """Testing target env_variable."""
    if file_target("a.txt").exists():
        file_target("a.txt").unlink()
    if sys.platform == "win32":
        script = SoS_Script("""
[0]
depends: env_variable('AA')
output:  'a.txt'
run:
echo %AA% > a.txt
""")
    else:
        script = SoS_Script("""
[0]
depends: env_variable('AA')
output:  'a.txt'
run:
echo $AA > a.txt
""")
    wf = script.workflow()
    os.environ["AA"] = "A1"
    Base_Executor(wf).run()
    with open("a.txt") as at:
        assert at.read().strip() == "A1"
    # test validation
    Base_Executor(wf).run()
    # now if we change var, it should be rerun
    os.environ["AA"] = "A2"
    Base_Executor(wf).run()
    # allow one second variation
    with open("a.txt") as at:
        assert at.read().strip() == "A2"
    file_target("a.txt").unlink()


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Windows executable cannot be created with chmod.")
def test_provides_executable():
    """Testing provides executable target."""
    # change $PATH so that lls can be found at the current
    # directory.
    os.environ["PATH"] += os.pathsep + "."
    if file_target("lls").exists():
        file_target("lls").unlink()
    script = SoS_Script("""
[lls: provides=executable('lkls')]
run:
touch lkls
chmod +x lkls

[c]
depends: executable('lkls')

""")
    wf = script.workflow("c")
    Base_Executor(wf).run()
    file_target("lkls").unlink()


def test_shared_var_in_paired_with(temp_factory):
    temp_factory(["1.txt", "2.txt"])
    for file in ("1.out", "2.out", "1.out2", "2.out2"):
        if file_target(file).exists():
            file_target(file).unlink()
    script = SoS_Script("""
[work_1: shared = {'data': 'step_output'}]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run: expand=True
touch {_output}

[work_2]
depends: sos_variable('data')
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}', paired_with = dict(data1=data)
output: expand_pattern('{_name}.out2')
run: expand=True
touch {data1[0]} {_output}
""")
    wf = script.workflow()
    Base_Executor(wf).run()
    for file in ("1.out", "2.out", "1.out2", "2.out2"):
        if file_target(file).exists():
            file_target(file).unlink()


def test_shared_var_in_for_each(temp_factory, clear_now_and_after):
    temp_factory("1.txt", "2.txt")
    clear_now_and_after("1.out", "2.out", "1.out2", "2.out2", '2.out_2.out2', '1.out_1.out2', '1.out_2.out2', '2.out_1.out2')
    script = SoS_Script("""
[work_1: shared = {'data': 'step_output'}]
input: "1.txt", "2.txt", group_by = 'single', pattern = '{name}.{ext}'
output: expand_pattern('{_name}.out')
run: expand=True
touch {_output}

[work_2]
depends: sos_variable('data')
input: "1.txt", "2.txt", group_by = 'single', for_each = dict(data=data),  pattern = '{name}.{ext}'
output: expand_pattern('{data}_{_name}.out2')
run: expand=True
touch {_output}

""")
    wf = script.workflow()
    Base_Executor(wf).run()


def test_removed_depends(clear_now_and_after):
    """Test a case where a dependent file has signature, but
    gets removed before the next run."""
    clear_now_and_after('a.txt', 'b.txt')
    script = """
        [tet: provides='a.txt']
        run:
        echo "something" > a.txt

        [20]
        depends: 'a.txt'
        output: 'b.txt'
        run:
        cat a.txt > b.txt
        """
    # this should be ok.
    execute_workflow(script)
    # now let us remove a.txt (but the signature is still there)
    clear_now_and_after('a.txt', 'b.txt')
    execute_workflow(script)


def test_sos_step(clear_now_and_after):
    """Test target sos_step"""
    clear_now_and_after("t1.txt", "t2.txt", "5.txt", "10.txt", "20.txt")

    script = """
        [t1]
        run:
        touch t1.txt

        [t2: provides='t2.txt']
        depends: sos_step('t1')
        run:
        touch t2.txt

        [5]
        run:
        touch 5.txt

        [10]
        depends: sos_step('t2')
        run:
        touch 10.txt

        [20]
        depends: sos_step('t1')
        run:
        touch 20.txt
        """
    # this should be ok.
    execute_workflow(script, options={'sig_mode': 'force'})
    for file in ["t1.txt", "t2.txt", "5.txt", "10.txt", "20.txt"]:
        assert file_target(file).target_exists(), file + " should exist"


def test_zap(clear_now_and_after):
    """Test zap"""
    clear_now_and_after('testzap1.txt.zapped', 'testzap.txt.zapped')
    with open("testzap.txt", "w") as sf:
        sf.write("some text")
    path("testzap.txt").zap()
    assert os.path.isfile("testzap.txt.zapped")
    assert not os.path.isfile("testzap.txt")
    # re-zap is ok
    file_target("testzap.txt").zap()
    assert os.path.isfile("testzap.txt.zapped")
    assert not os.path.isfile("testzap.txt")
    # non-existent file
    os.remove("testzap.txt.zapped")
    with pytest.raises(FileNotFoundError):
        path("testzap.txt").zap()
    #
    with open("testzap.txt", "w") as sf:
        sf.write("some text")
    with open("testzap1.txt", "w") as sf:
        sf.write("some text")
    paths("testzap.txt", "testzap1.txt").zap()
    assert os.path.isfile("testzap.txt.zapped")
    assert not os.path.isfile("testzap.txt")
    assert os.path.isfile("testzap1.txt.zapped")
    assert not os.path.isfile("testzap1.txt")
    #
    os.remove("testzap.txt.zapped")
    os.remove("testzap1.txt.zapped")
    with open("testzap.txt", "w") as sf:
        sf.write("some text")
    with open("testzap1.txt", "w") as sf:
        sf.write("some text")
    sos_targets(["testzap.txt", "testzap1.txt"]).zap()
    assert os.path.isfile("testzap.txt.zapped")
    assert not os.path.isfile("testzap.txt")
    assert os.path.isfile("testzap1.txt.zapped")
    assert not os.path.isfile("testzap1.txt")


def test_zap_run(clear_now_and_after):
    """Test run with zapped input files"""
    clear_now_and_after('zap1.txt.zapped')
    with open("zap1.txt", "w") as sf:
        sf.write("seomething")
    script = """\
        [1]
        input: 'zap1.txt'
        output: "zap2.txt"
        run:
            echo asd > zap2.txt
        _input.zap()
        """
    execute_workflow(script, options={"sig_mode": "force"})

    assert os.path.isfile("zap1.txt.zapped")
    assert not os.path.isfile("zap1.txt")
    assert os.path.isfile("zap2.txt")
    # can run again
    env.config["sig_mode"] = "default"

    execute_workflow(script)
    # now if we remove target
    os.remove("zap2.txt")
    with pytest.raises(Exception):
        execute_workflow(script)


def test_system_resource():
    """Test targtet system_resource"""
    execute_workflow("""\
        [1: shared='a']
        depends: system_resource(mem='1M',disk='1M')
        a = 1
        """)
    assert env.sos_dict["a"] == 1

    #
    with pytest.raises(Exception):
        execute_workflow("""\
            [1: shared='a']
            depends: system_resource(mem='1T')
            a = 1
            """)

    with pytest.raises(Exception):
        execute_workflow("""\
            [1: shared='a']
            depends: system_resource(disk='10P')
            a = 1
        """)


def test_temp_file():
    """Test sos_tempfile target"""
    execute_workflow("""
        filename = sos_tempfile(name='a')

        assert filename == sos_tempfile(name='a')
        assert filename != sos_tempfile(name='b')
        assert filename != sos_tempfile()
        assert sos_tempfile() != sos_tempfile()

        assert sos_tempfile(prefix='some') != sos_tempfile(prefix='some')
        assert sos_tempfile(suffix='.sh') != sos_tempfile(suffix='.sh')
        assert sos_tempfile(dir='.') != sos_tempfile(dir='.')

        assert sos_tempfile('a.txt') == sos_tempfile('a.txt')
        assert sos_tempfile(path('a.txt')) == sos_tempfile(path('a.txt'))
        assert sos_tempfile(file_target('a.txt')) == sos_tempfile(file_target('a.txt'))

    """)


def test_named_path():
    """Test the use of option name of path"""
    execute_workflow(
        """
        import os
        # windows might not have HOME
        if 'HOME' in os.environ:
            assert path('#home') == os.environ['HOME']
            assert 'home' in path.names()
            assert 'home' in path.names('docker')
        """,
        options={
            "config_file": os.path.join(os.path.expanduser("~"), "docker.yml")
        },
    )


@pytest.mark.skipif(
    sys.platform == 'win32', reason='Graphviz not available under windows')
def test_to_named_path_path():
    execute_workflow(
        """
        [10: shared="a"]
        a = path('/root/xxx/whatever').to_named_path(host='docker')
        """,
        options={
            "config_file": os.path.join(os.path.expanduser("~"), "docker.yml"),
        })
    assert env.sos_dict['a'] == '#home/xxx/whatever'
