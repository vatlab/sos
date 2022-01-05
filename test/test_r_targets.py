#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import shutil
import subprocess

import pytest

from sos import execute_workflow


@pytest.mark.skipif(
    not shutil.which("Rscript") or "TRAVIS" in os.environ,
    reason="R not installed")
def test_r_library():
    """Test target R_Library"""
    execute_workflow("""
        [default]
        depends: R_library("dplyr", autoinstall=True)
        R:
            library('dplyr')
        """)


@pytest.mark.skipif(not shutil.which("Rscript"), reason="R not installed")
def test_depends_r_library():
    """Testing depending on R_library"""
    # first remove xtable package
    subprocess.call("R CMD REMOVE xtable", shell=True)
    execute_workflow("""
        [0]

        depends: R_library('xtable', autoinstall=True)
        R:
        library('xtable')
        ## Demonstrate data.frame
        tli.table <- xtable(cars)
        """)


@pytest.mark.skipif(not shutil.which("Rscript"), reason="R not installed")
def test_reexecution():
    """Test re-execution of steps with R_library"""
    subprocess.call("R CMD REMOVE xtable", shell=True)
    wf = """
    [1]
    depends: R_library("xtable", autoinstall=True)
    output: '1.txt'
    run: expand=True
        sleep 5
        touch {_output}
    """
    execute_workflow(wf)
    execute_workflow(wf)
