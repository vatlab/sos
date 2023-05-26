#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action, SoS_ExecuteScript


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def bash(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified script using bash. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script."""
    return SoS_ExecuteScript(script, interpreter or ["/bin/bash", "bash"], ".sh", args, entrypoint).run(**kwargs)


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def csh(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified script using csh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script."""
    return SoS_ExecuteScript(script, interpreter or ["/bin/csh", "csh"], ".csh", args, entrypoint).run(**kwargs)


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def tcsh(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified script using tcsh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script."""
    return SoS_ExecuteScript(script, interpreter or ["/bin/tcsh", "tcsh"], ".sh", args, entrypoint).run(**kwargs)


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def zsh(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified script using zsh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script."""
    return SoS_ExecuteScript(script, interpreter or ["/bin/zsh", "zsh"], ".zsh", args, entrypoint).run(**kwargs)


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def sh(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified script using sh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script."""
    return SoS_ExecuteScript(script, interpreter or ["/bin/sh", "sh"], ".sh", args, entrypoint).run(**kwargs)
