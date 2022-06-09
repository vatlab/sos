#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil
import subprocess
import sys
import time
from contextlib import contextmanager

import pytest
from sos import execute_workflow
from sos.parser import SoS_Script
from sos.tasks import TaskFile, TaskParams
from sos.utils import env, textMD5
from sos.workflow_executor import Base_Executor

has_docker = sys.platform != "win32"
try:
    if sys.platform != "win32":
        subprocess.check_output(
            "docker ps | grep test_sos", shell=True).decode()
except subprocess.CalledProcessError:
    subprocess.call("sh build_test_docker.sh", shell=True)
    try:
        subprocess.check_output(
            "docker ps | grep test_sos", shell=True).decode()
    except subprocess.CalledProcessError:
        print("Failed to set up a docker machine with sos")
        has_docker = False


@contextmanager
def cd_new(path):
    old_dir = os.getcwd()
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_dir)


def get_tasks():
    from sos.signatures import WorkflowSignatures
    env.exec_dir = os.path.join(
        os.path.expanduser("~"), ".sos", textMD5(os.getcwd()))
    db = WorkflowSignatures()
    conn = db.conn
    #conn = sqlite3.connect(os.path.join(env.exec_dir, "workflow_signatures.db"))
    cur = conn.cursor()
    cur.execute('SELECT DISTINCT id FROM workflows WHERE entry_type = "task"')
    return [x[0] for x in cur.fetchall() if os.path.isfile(
        os.path.join(os.path.expanduser('~'), '.sos', 'tasks', x[0] + '.task'
    ))]


def test_task_file():
    """Test task file handling"""
    for ext in (".pulse", ".out", ".err", ".task", ".sh"):
        filename = os.path.join(
            os.path.expanduser("~"), ".sos", "tasks", "ffffffffffffffff" + ext)
        if os.path.isfile(filename):
            os.remove(filename)
    params = TaskParams(
        name="ffffffffffffffff",
        global_def=None,
        task="b=a",
        sos_dict={"a": 1},
        tags=["b", "a"],
    )
    a = TaskFile("ffffffffffffffff")
    a.save(params)
    assert a.tags == "a b"
    for ext in (".pulse", ".out", ".err", ".sh"):
        with open(
                os.path.join(
                    os.path.expanduser("~"), ".sos", "tasks",
                    "ffffffffffffffff" + ext),
                "w",
        ) as fh:
            fh.write(ext)
    assert not a.has_stdout()
    assert not a.has_stderr()
    a.add_outputs()
    #
    assert a.params.sos_dict["a"] == 1
    assert a.status == "new"
    a.status = "completed"
    assert time.time() - a.last_updated < 2
    assert a.status == "completed"
    #
    # get and reset info
    info = a.info
    a.status = "running"
    assert a.status == "running"
    a.info = info
    assert a.status == "completed"
    assert a.has_stdout()
    #
    a.add_result({"ret_code": 5})
    #
    a.tags = ["ee", "d"]
    assert a.tags == "d ee"
    # a.add_tags(['kk'])
    # assert a.tags.split(), ['d', 'ee', 'kk'])
    #
    assert a.params.sos_dict["a"] == 1
    assert a.params.task == "b=a"
    #

    assert a.stdout == ".out"
    assert a.stderr == ".err"
    assert a.pulse == ".pulse"
    assert a.shell == ".sh"
    assert a.has_stdout()
    assert a.has_stderr()
    assert a.has_pulse()
    assert a.has_shell()
    #
    #
    a.reset()
    assert a.status == "new"
    assert a.stdout == ""
    assert a.stderr == ""
    assert a.signature == {}
    a.add_outputs()
    a.add_result({"ret_code": 5})
    assert a.result["ret_code"] == 5


def test_workdir():
    """Test workdir option for runtime environment"""
    import tempfile

    tdir = tempfile.mkdtemp()
    with open(os.path.join(tdir, "aaa.pp"), "w") as aaa:
        aaa.write("something")
    script = r"""
import os
[0]
task: workdir={0!r}

with open(os.path.join({1!r}, 'result.txt'), 'w') as res:
    for file in os.listdir({1!r}):
        res.write(file + '\n')
""".format(os.path.split(tdir)[0],
           os.path.split(tdir)[1])

    execute_workflow(
        script, options={
            "sig_mode": "force",
            "default_queue": "localhost",
        })
    with open(os.path.join(tdir, "result.txt")) as res:
        content = [x.strip() for x in res.readlines()]
        assert "aaa.pp" in content


def test_concurrency():
    """Test concurrency option for runtime environment"""
    execute_workflow(
        r"""
[0]

repeat = range(4)
input: for_each='repeat'

task:

import time
print('I am {}, waited {} seconds'.format(_index, _repeat + 1))
time.sleep(_repeat + 1)
print('I am {}, done'.format(_index))
""",
        options={
            "default_queue": "localhost",
            "sig_mode": "force"
        },
    )


def test_prepend_path():
    """Test prepend path"""
    import stat

    if not os.path.isdir("temp"):
        os.mkdir("temp")
    if sys.platform == "win32":
        with open(r"temp\temp_cmd.bat", "w") as tc:
            tc.write('echo "a"')
    else:
        with open("temp/temp_cmd", "w") as tc:
            tc.write('echo "a"')
        os.chmod("temp/temp_cmd", stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR)
    #
    # use option env
    execute_workflow(
        script=r"""
import os
[1]
task: env={'PATH': 'temp' + os.pathsep + os.environ['PATH']}
run:
temp_cmd
""",
        options={
            "default_queue": "localhost",
            "sig_mode": "force"
        },
    )
    #
    #
    execute_workflow(
        script=r"""
[1]
task: prepend_path='temp'
run:
temp_cmd
""",
        options={"default_queue": "localhost"},
    )


def test_no_task():
    execute_workflow(
        r"""
        [10]
        task:
        run:
        sleep 0
        """,
        options={
            "default_queue": None,
            'sig_mode': 'force',
        })
    #
    env.config["sig_mode"] = "force"
    execute_workflow(
        r"""
        [10]
        task: queue=None
        run:
        sleep 0
        """,
        options={
            "default_queue": "localhost",
            'sig_mode': 'force',
        })


def test_passing_var_to_task(clear_now_and_after):
    """Test passing used variable to tasks"""
    clear_now_and_after([f"myfile_{i}.txt" for i in range(10, 13)])
    #
    env.config["sig_mode"] = "force"

    execute_workflow(
        r"""
        parameter: gvar = 10

        [10]
        # generate a file
        tt = range(gvar, gvar + 3)
        input: for_each='tt'
        output: f"myfile_{_tt}.txt"
        # additional comment

        # _tt should be used in task
        task:
        run: expand=True
        echo {_tt}_{_index} > {_output:q}

        """,
        config={
            "default_queue": "localhost",
            "sig_mode": "force",
            "max_jobs": 4
        })
    for t in range(10, 13):
        with open(f"myfile_{t}.txt") as tmp:
            assert tmp.read().strip() == str(t) + "_" + str(t - 10)


def test_max_jobs():
    """Test default max number of jobs"""
    execute_workflow(
        r"""
        [10]
        input: for_each=[{'a': range(2)}, {'b': range(3)}]

        task:
        run: expand=True
        echo "a = {a}, b = {b}"
        """,
        options={
            "max_running_jobs": 2,
            "default_queue": "localhost",
        })


def test_kill_and_purge():
    """Test no wait"""
    subprocess.call(["sos", "purge"])
    with open("test_purge.sos", "w") as script:
        script.write(r"""
[10]
task:
run:
echo Try to kill
sleep 20
""")
    subprocess.Popen("sos run test_purge.sos -s force -q localhost", shell=True)
    time.sleep(5)
    subprocess.call(["sos", "kill", "--all"])
    for i in range(20):
        output = subprocess.check_output(["sos", "status", "-v", "1"]).decode()
        if "killed" in output or "aborted" in output or "completed" in output:
            break
        assert i <= 10, "Task should be killed within 10 seconds, got {}".format(
            output)
        time.sleep(1)
    # test purge by status
    subprocess.call(["sos", "purge", "--status", "aborted"])
    assert "killed" not in subprocess.check_output(["sos", "status", "-v",
                                                    "3"]).decode()
    # purge by all is not tested because it is dangerous


def test_concurrent_task(clear_now_and_after):
    """Test submitting tasks from concurrent substeps"""
    clear_now_and_after([f"con_{x}.txt" for x in range(5)])
    execute_workflow(
        """
        [10]
        input: for_each={'i': range(5)}
        output: f'con_{i}.txt'

        task:
        run: expand=True
        echo {i} > {_output}
        """,
        options={
            "sig_mode": "force",
            "default_queue": "localhost"
        })


def test_shared_option(clear_now_and_after):
    """Test shared option of task"""
    clear_now_and_after("a.txt", "a100.txt")
    execute_workflow(
        """
        [10: shared = 'a']
        output: 'a.txt'
        task: shared={'a': 'int(open("a.txt").read())'}
        run:
        echo 100 > a.txt

        [20]
        run: expand=True
        touch a{a}.txt
        """,
        options={
            "sig_mode": "force",
            "default_queue": "localhost"
        })
    assert os.path.isfile("a100.txt")
    # sequence of var or mapping
    clear_now_and_after("a.txt", "a100.txt")

    execute_workflow(
        """
        [10: shared = ['a', 'b']]
        output: 'a.txt'
        task: shared=[{'a': 'int(open("a.txt").read())'}, 'b']
        b = 20
        run:
        echo 100 > a.txt

        [20]
        run: expand=True
        touch a{a}_{b}.txt
        """,
        options={
            "sig_mode": "force",
            "default_queue": "localhost"
        })
    assert os.path.isfile("a100_20.txt")

    execute_workflow(
        """
        [10 (simulate): shared=['rng', 'step_rng']]
        input: for_each={'i': range(5)}
        task: shared='rng'
        print(f"{i}")
        import random
        rng = random.randint(1, 1000)
        """,
        options={"default_queue": "localhost"})

    var = env.sos_dict["rng"]
    assert isinstance(var, int)
    assert isinstance(env.sos_dict["step_rng"], list)
    assert env.sos_dict["step_rng"][-1] == var


def test_task_tags():
    """Test option tags of tasks"""
    import random

    tag = "tag{}".format(random.randint(1, 100000))
    with open("test_tags.sos", "w") as tt:
        tt.write("""
[10]
input: for_each={{'i': range(10)}}
task: tags='{}', trunk_size=2
sh: expand=True
echo {} {{i}}
""".format(tag, tag))
    wf = SoS_Script(filename="test_tags.sos").workflow()
    Base_Executor(
        wf,
        config={
            "sig_mode": "force",
            "script": "test_trunkworker.sos",
            "max_running_jobs": 10,
            "default_queue": "localhost",
            "worker_procs": ["4"],
        },
    ).run()
    ret = subprocess.check_output(
        "sos status -t {}".format(tag), shell=True).decode()
    assert len(ret.splitlines()) == 5, "Obtained {}".format(ret)
    # test multiple tags
    tag1 = "tag{}".format(random.randint(1, 100000))
    tag2 = "tag{}".format(random.randint(1, 100000))
    with open("test_tags.sos", "w") as tt:
        tt.write("""
[10]
input: for_each={{'i': range(2)}}
task: tags=['{}', '{}']
sh: expand=True
echo {} {{i}}
""".format(tag1, tag2, tag1))
    wf = SoS_Script(filename="test_tags.sos").workflow()
    Base_Executor(
        wf,
        config={
            "sig_mode": "force",
            "script": "test_trunkworker.sos",
            "max_running_jobs": 10,
            "worker_procs": ["4"],
            "default_queue": "localhost",
            "workflow": "default",
        },
    ).run()
    ret = subprocess.check_output(
        "sos status -t {}".format(tag2), shell=True).decode()
    assert len(ret.splitlines()) == 2, "Obtained {}".format(ret)


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_max_mem():
    """Test server restriction max_mem"""
    with pytest.raises(Exception):
        execute_workflow(
            """
            [10]
            task: mem='2G'
            print('test_max_mem')
            """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "docker_limited",
                "sig_mode": "force",
            },
        )


def test_local_runtime_max_walltime():
    """Test server max_walltime option"""
    with pytest.raises(Exception):
        execute_workflow(
            """
            [10]
            task:
            import time
            time.sleep(15)
            """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "local_limited",
                "sig_mode": "force",
            },
        )


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_max_cores():
    """Test server restriction max_cores"""
    with pytest.raises(Exception):
        execute_workflow(
            """
            [10]
            task: cores=8
            print('test_max_cores')
            """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "docker_limited",
                "sig_mode": "force",
            },
        )


@pytest.mark.skipIf(not has_docker, reason="Docker container not usable")
def test_override_max_cores():
    """Test use queue_args to override server restriction max_cores"""
    execute_workflow(
        """
        [10]
        task: cores=8
        print('test_override_max_cores')
        """,
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker_limited",
            "sig_mode": "force",
            "queue_args": {
                "cores": 1
            },
        },
    )


def test_list_hosts():
    """test list hosts using sos status -q"""
    for v in ["0", "1", "3", "4"]:
        output = subprocess.check_output(
            ["sos", "remote", "list", "-c", "~/docker.yml", "-v", v]).decode()
        assert "local_limited" in output, f"local_limited not in \n{output}\n for verbosity {v}"


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_max_walltime(purge_tasks):
    """Test server restriction max_walltime"""
    with pytest.raises(Exception):
        execute_workflow(
            """
            [10]
            task: walltime='1:00:00'
            print('test_max_walltime')
            """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "docker_limited",
                "sig_mode": "force",
            },
        )


def test_purge_all_with_option():
    """Test sos purge all with options such as -s completed"""
    with cd_new("temp_c"):
        with open("test.sos", "w") as tst:
            tst.write("""
input: for_each={'i': range(2)}
output: f'a{i}.txt'
task:
run: expand=True
echo temp_a
touch {_output}
""")
        subprocess.call("sos run test -s force -q localhost", shell=True)
        tasks = get_tasks()
        subprocess.call("sos purge -s failed", shell=True)
    # check tasks
    taskstatus = [
        x.split()[0] for x in subprocess.check_output(
            "sos status -v1", shell=True).decode().splitlines()
    ]
    assert all(x in taskstatus for x in tasks)
    # purge one of them
    subprocess.call(f"sos purge {tasks[0]}", shell=True)
    taskstatus = [
        x.split()[0] for x in subprocess.check_output(
            "sos status -v1", shell=True).decode().splitlines()
    ]
    assert tasks[0] not in taskstatus
    assert tasks[1] in taskstatus
    #
    subprocess.call("sos purge --all", shell=True)
    taskstatus = [
        x.split()[0] for x in subprocess.check_output(
            "sos status -v1", shell=True).decode().splitlines()
    ]
    assert tasks[1] not in taskstatus


def test_resubmit_task_with_different_walltime():
    """Test resubmission of tasks with different walltime #1019"""
    with cd_new("temp_walltime"):
        with open("test.sos", "w") as tst:
            tst.write("""
task: walltime='1m'
sh:
echo 0.1
""")
        subprocess.call("sos run test -s force -q localhost", shell=True)
        tasks = get_tasks()
        out = subprocess.check_output(f"sos status {tasks[0]} -v4", shell=True)
        assert "00:01:00" in out.decode()
        with open("test1.sos", "w") as tst:
            tst.write("""
task: walltime='2m'
sh:
echo 0.1
""")
        subprocess.call("sos run test1 -s force -q localhost", shell=True)
        new_tasks = get_tasks()
        assert tasks == new_tasks
        #
        out = subprocess.check_output(f"sos status {tasks[0]} -v4", shell=True)
        assert "00:02:00" in out.decode()


def test_task_no_signature(purge_tasks):
    """Test re-execution of tasks"""
    with cd_new("temp_signature"):
        with open("test.sos", "w") as tst:
            tst.write("""
task:
sh:
    sleep 2
""")
        subprocess.call("sos run test -s force -q localhost", shell=True)
        tasks = get_tasks()
        tf = TaskFile(tasks[0])
        assert not tf.has_signature()
        assert tf.status == "completed"
        #
        st = time.time()
        subprocess.call("sos run test -q localhost", shell=True)
        assert time.time() - st > 1


def test_task_with_signature(purge_tasks, clear_now_and_after):
    """Test re-execution of tasks"""
    # now with a real signature
    with cd_new("temp_signature"):
        with open("test.sos", "w") as tst:
            tst.write("""
output: 'a.txt'
task:
sh: expand=True
    sleep 3
    touch {_output}
""")
        subprocess.call("sos run test -s force -q localhost", shell=True)
        tasks = get_tasks()
        tf = TaskFile(tasks[0])
        assert tf.has_signature()
        assert tf.status == "completed"
        #
        st = time.time()
        subprocess.call("sos run test -q localhost", shell=True)
        assert time.time() - st < 2
        #
        clear_now_and_after('a.txt')
        st = time.time()
        subprocess.call("sos run test -q localhost", shell=True)
        assert time.time() - st > 2


def test_wrong_host():
    with pytest.raises(Exception):
        execute_workflow(
            """
            [10]
            task: walltime='1:00:00', queue='undefined'
            print('test_wrong_host')
            """,
            options={"default_queue": "localhost"})


def test_output_in_task():
    """Test passing _output to task #1136"""
    execute_workflow(
        """
        chunks  = [1,2]
        [1]
        input: for_each = 'chunks'
        output: f'{_chunks}.txt'
        _output.touch()

        [2]
        input: group_with = 'chunks'
        output: summary_stats = f'{_input}.summary', ld_matrix = f'{_input}.result'
        task:

        python3: expand="${ }"
            open("${_output['summary_stats']}", 'w').close()
            open("${_output['ld_matrix']}", 'w').close()
        """,
        options={"default_queue": "localhost"})


def test_repeated_tasks():
    """Test statement before task #1142 """
    for i in range(5):
        execute_workflow(
            """
            [1]
            input: for_each=dict(i=range(5))

            print(f'This is for {i}')
            task:  walltime='10m'
            print(f'this is task {i}')

            """,
            options={
                "sig_mode": "force",
                "default_queue": "localhost"
            })


def test_passing_parameters_to_task():
    """Test passing of parameters in global section to tasks #1155"""
    execute_workflow(
        """\
        [global]
        parameter: par=5
        def a():
            print(par)

        [default]
        task:
        a()
        """,
        options={
            "sig_mode": "force",
            "default_queue": "localhost"
        })


def test_trunk_size_with_stop_if(clear_now_and_after):
    """Test a case when some tasks are not submitted due to holes in slots #1159"""
    clear_now_and_after([f"{i+1}.txt" for i in range(5)])
    execute_workflow(
        """\
        [1]
        output: [f'{x+1}.txt' for x in range(5)]
        for i in range(5):
            name = f'{i+1}.txt'
            if i not in [0,1,2]:
                path(name).touch()
            else:
                with open(name, 'w') as f:
                    f.write('test it')

        [2]
        input: group_by = 1
        output: f'{_input:n}.out'
        stop_if(_input.stat().st_size==0, no_output=True)

        task: trunk_size = 80
        _output.touch()
        """,
        options={
            "sig_mode": "force",
            "default_queue": "localhost"
        })


def test_output_from_master_task():
    """Test splitting the output from master task #1203"""
    execute_workflow(
        """\
        l=[x for x in range(1,13)]

        [2]
        input: for_each = 'l'
        output: f'{_l}.out'

        task: trunk_size = 4
        _output.touch()

        [3]
        assert _input == f'{_index+1}.out'
        """,
        options={"default_queue": "localhost"})


@pytest.mark.skipIf(not has_docker, reason="Docker container not usable")
def test_remote_input_target(clear_now_and_after):
    """Test the use of remote target"""
    clear_now_and_after("vars.sh", "vars1.sh")
    execute_workflow(
        """
        [10]
        input: remote('/lib/init/vars.sh')
        output: f'vars1.sh'

        task:

        with open(_input, 'r') as inf, open(_output, 'w') as outf:
            outf.write(inf.read())
        """,
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )
    assert not os.path.isfile("vars.sh")
    assert os.path.isfile("vars1.sh")


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_delayed_interpolation(clear_now_and_after):
    """Test delayed interpolation with expression involving remote objects"""
    # purge all previous tasks
    clear_now_and_after('test.py', 'test.py.bak')
    execute_workflow(
        """
        [10]
        output: remote('test.py')
        task:
        run:
            touch test.py

        [20]
        output: remote(f"{_input:R}.bak")
        task:
        run: expand=True
            cp {_input} {_output}
        """,
        options={
            "config_file": "~/docker.yml",
            # do not wait for jobs
            "wait_for_task": True,
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )
    # this file is remote only
    assert not os.path.isfile("test.py")
    assert not os.path.isfile("test.py.bak")


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_remote_output_target(clear_now_and_after):
    """Test the use of remote target"""
    clear_now_and_after("vars.sh", "vars1.sh")
    execute_workflow(
        """
        [10]
        input: remote('/lib/init/vars.sh')
        output: remote(f'vars1.sh')

        task:

        with open(_input, 'r') as inf, open(_output, 'w') as outf:
            outf.write(inf.read())
        """,
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )
    assert not os.path.isfile("vars.sh")
    assert not os.path.isfile("vars1.sh")


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_remote_output_target_with_trunksize(clear_now_and_after):
    clear_now_and_after("vars.sh", "vars1.sh")

    execute_workflow(
        """\
        [10]
        import os
        input: remote('/lib/init/vars.sh'), remote('/lib/init/init-d-script'), group_by=1
        output: remote(os.path.basename(str(_input)))

        task: trunk_size=2

        with open(_input, 'r') as inf, open(_output, 'w') as outf:
            outf.write(inf.read())""",
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )
    assert not os.path.isfile("vars.sh")
    assert not os.path.isfile("init-d-script")


def test_runtime_max_walltime():
    """Test server max_walltime option"""
    with pytest.raises(Exception):
        execute_workflow(
            """
        [10]
        task:
        import time
        time.sleep(25)
        """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "docker_limited",
                "sig_mode": "force",
            },
        )


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_sync_master_task(clear_now_and_after):
    """Test sync input and output with remote host with trunksize"""
    clear_now_and_after([f"test_{i}.txt" for i in range(4)],
                        [f"test_{i}.bak" for i in range(4)])
    import random

    val = random.randint(1, 10000)
    execute_workflow(
        r"""
        parameter: g = 100

        [10]
        input: for_each=dict(i=range(4))
        output: f'test_{i}.txt'

        with open(f'test_{i}.txt', 'w') as tst:
            tst.write(f'test_{i}_{g}')

        [20]
        output: _input.with_suffix('.bak')

        task: trunk_size=2

        with open(_input, 'r') as inf, open(_output, 'w') as outf:
            outf.write(inf.read() + '.bak')
        """,
        args=["--g", str(val)],
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )
    # now check if
    for i in range(4):
        assert os.path.isfile(f"test_{i}.txt")
        with open(f"test_{i}.bak") as outf:
            assert outf.read() == f"test_{i}_{val}.bak"
        assert os.path.isfile(f"test_{i}.bak")
        with open(f"test_{i}.bak") as outf:
            assert outf.read() == f"test_{i}_{val}.bak"


def test_local_max_cores():
    """Test server restriction max_cores"""
    with pytest.raises(Exception):
        execute_workflow(
            """
        [10]
        task: cores=8
        print('test_local_max_cores')
        """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "local_limited",
                "sig_mode": "force",
            },
        )


def test_local_max_mem():
    """Test server restriction max_mem"""
    with pytest.raises(Exception):
        execute_workflow(
            """
        [10]
        task: mem='2G'
        print('test_local_max_mem')
        """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "local_limited",
                "sig_mode": "force",
            },
        )


def test_trunk_size_option(clear_now_and_after, purge_tasks):
    """Test option trunk_size"""
    clear_now_and_after([f"{i}.txt" for i in range(10)])
    execute_workflow(
        """
[10]
input: for_each={'I': range(10)}
task: trunk_size=5, walltime='10m'
run: expand=True
  echo {I} > {I}.txt
  sleep 0.1
""",
        options={
            "sig_mode": "force",
            "default_queue": "localhost",
        },
    )

    for i in range(10):
        assert os.path.isfile(f"{i}.txt")


def test_trunk_size_none_option(clear_now_and_after, purge_tasks):
    """Test option trunk_size"""
    clear_now_and_after([f"{i}.txt" for i in range(10)])
    # trunk size is None or 0, -1, intepreted as all tasks
    execute_workflow(
        """
[10]
input: for_each={'I': range(10)}
task: trunk_size=None, cores=1, walltime='10m'
run: expand=True
  echo {I} > {I}.txt
  sleep 0.1
""",
        options={
            "sig_mode": "force",
            "default_queue": "localhost",
        },
    )
    for i in range(10):
        assert os.path.isfile(f"{i}.txt")


def test_trunk_workers_option(clear_now_and_after, purge_tasks):
    """Test option trunk_workers"""
    clear_now_and_after([f"{i}.txt" for i in range(12)])

    execute_workflow(
        """
[10]
input: for_each={'I': range(12)}
task: trunk_size=6, trunk_workers=3, walltime='10m'
run: expand=True
  echo {I} > {I}.txt
  sleep 1
""",
        options={
            "max_running_jobs": 10,
            "worker_procs": ["4"],
            "default_queue": "localhost"
        },
    )
    for i in range(12):
        assert os.path.isfile(f"{i}.txt")


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_sync_input_output_and_rerun(clear_now_and_after, purge_tasks):
    """Test sync input and output with remote host"""
    clear_now_and_after([f'test_{i}.txt' for i in range(4)])
    clear_now_and_after([f'test_{i}.bak' for i in range(4)])

    import random
    val = random.randint(1, 10000)
    wf = """
        parameter: g = 100

        [10]
        input: for_each=dict(i=range(4))
        output: f'test_{i}.txt'

        with open(f'test_{i}.txt', 'w') as tst:
            tst.write(f'test_{i}_{g}')

        [20]
        output: _input.with_suffix('.bak')

        task:

        with open(_input, 'r') as inf, open(_output, 'w') as outf:
            outf.write(inf.read() + '.bak')
        """
    execute_workflow(
        wf,
        args=["--g", str(val)],
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )

    # now check if
    for i in range(4):
        assert os.path.isfile(f"test_{i}.txt")
        with open(f"test_{i}.bak") as outf:
            assert outf.read() == f"test_{i}_{val}.bak"

        assert os.path.isfile(f"test_{i}.bak")
        with open(f"test_{i}.bak") as outf:
            assert outf.read() == f"test_{i}_{val}.bak"
    #
    # test rerun the task file on local host
    clear_now_and_after([f'test_{i}.txt' for i in range(4)])
    clear_now_and_after([f'test_{i}.bak' for i in range(4)])

    execute_workflow(wf, args=["--g", str(val)], options={"sig_mode": "force"})

    for i in range(4):
        assert os.path.isfile(f"test_{i}.txt")
        with open(f"test_{i}.bak") as outf:
            assert outf.read() == f"test_{i}_{val}.bak"

        assert os.path.isfile(f"test_{i}.bak")
        with open(f"test_{i}.bak") as outf:
            assert outf.read() == f"test_{i}_{val}.bak"
