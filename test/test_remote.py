#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import subprocess

import pytest

from sos import execute_workflow
from sos.targets import file_target

has_docker = True
try:
    subprocess.check_output("docker ps | grep test_sos", shell=True).decode()
except subprocess.CalledProcessError:
    subprocess.call("sh build_test_docker.sh", shell=True)
    try:
        subprocess.check_output("docker ps | grep test_sos", shell=True).decode()
    except subprocess.CalledProcessError:
        print("Failed to set up a docker machine with sos")
        has_docker = False


def test_worker_procs():
    # test -j option
    execute_workflow(
        """
        [1]
        input: for_each=dict(i=range(10))

        bash: expand=True
        echo {i}
        """,
        options={
            "sig_mode": "force",
            "worker_proces": ["1", "localhost:2"]
        },
    )


def test_worker_procs_with_task():
    # test -j option
    execute_workflow(
        """
        [1]
        input: for_each=dict(i=range(10))

        task: trunk_size = 0

        bash: expand=True
        echo {i}
        """,
        options={
            "sig_mode": "force",
            "default_queue": "localhost",
            "worker_proces": ["1", "localhost:2"],
        },
    )


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_remote_execute(clear_now_and_after, script_factory):
    clear_now_and_after("result_remote.txt", "result_remote1.txt", "local.txt", "remote_exec.sos")

    test_remote_sos = script_factory("""
        [10]
        input: 'local.txt'
        output: 'result_remote.txt'
        task:

        run:
            cp local.txt result_remote.txt
            echo 'adf' >> 'result_remote.txt'
        """, filename='remote_exec.sos')
    with open("local.txt", "w") as w:
        w.write("something")

    assert 0 == subprocess.call(
        f"sos run {test_remote_sos} -c ~/docker.yml -r docker -s force",
        shell=True,
    )

    assert file_target("result_remote.txt").target_exists()
    # self.assertEqual(subprocess.call('sos preview result_remote.txt', shell=True), 0)
    with open("result_remote.txt") as w:
        content = w.read()
        assert "something" in content
        assert "adf" in content


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_remote_workflow_remote_queue(script_factory):
    test_r_q = script_factory('''
        input: for_each=dict(i=range(2))
        output: f'test_r_q_{i}.txt'

        task: walltime='10m', cores=1, mem='1G'
        sh: expand=True
            echo `pwd` > {_output}
            echo I am {i} >> {_output}
        ''')
    assert 0 == subprocess.call(f"sos run {test_r_q} -c ~/docker.yml -r ts -q ts", shell=True)


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_remote_exec(clear_now_and_after):
    clear_now_and_after("result_exec.txt")
    execute_workflow(
        """
        output: 'result_exec.txt'

        task:
        sh: expand=True
            echo Output: {_output} > {_output}
            echo PWD: `pwd`. >> {_output}
        """,
        options={
            "config_file": "~/docker.yml",
            "default_queue": "docker",
            "sig_mode": "force",
        },
    )
    assert file_target("result_exec.txt").target_exists()
    with open(file_target("result_exec.txt")) as res:
        result = res.read()
        assert "Output: result_exec.txt" in result
        assert f"PWD: {os.getcwd()}" in result


@pytest.mark.skipif(not has_docker, reason="Docker container not usable")
def test_remote_exec_workdir(clear_now_and_after):
    clear_now_and_after(file_target("result_workdir_wo_named.txt"))
    with pytest.raises(Exception):
        execute_workflow(
            """
        output: 'result_workdir_wo_named.txt'

        task: workdir='/other'
        sh: expand=True
            echo Output: {_output} > {_output}
            echo PWD: `pwd`. >> {_output}
        """,
            options={
                "config_file": "~/docker.yml",
                "default_queue": "docker",
                "sig_mode": "force",
            },
        )
