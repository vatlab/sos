#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action, SoS_ExecuteScript


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def julia(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified Julia script with command julia. This action accepts common
    action arguments such as input, active, workdir, docker_image and args. In
    particular, content of one or more files  specified by option input would be
    prepended before the specified script.
    """
    return SoS_ExecuteScript(script, interpreter or "julia", ".jl", args, entrypoint).run(**kwargs)
