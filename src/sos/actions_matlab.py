#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action, SoS_ExecuteScript


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def matlab(
        script,
        interpreter,
        args='''-nojvm -nodisplay -nosplash -nodesktop -r "try, run('{filename}'), catch me, fprintf('%s / %s\\n',me.identifier,me.message), exit(1), end, exit(0);"''',
        entrypoint="",
        **kwargs):
    """Execute specified script with command Matlab, with default options
    "-nojvm -nodisplay -nosplash -nodesktop -r". This action accepts common action arguments such as input,
    active, workdir, docker_image and args. In particular, content of one or more
    files  specified by option input would be
    prepended before the specified script.
    """
    return SoS_ExecuteScript(script, interpreter or "matlab", ".m", args, entrypoint).run(**kwargs)


@SoS_Action(acceptable_args=["script", "interpreter", "args", "entrypoint"])
def octave(script, interpreter="", args="", entrypoint="", **kwargs):
    """Execute specified script with command Matlab, with default options
    "-nodisplay -r". This action accepts common action arguments such as input,
    active, workdir, docker_image and args. In particular, content of one or more
    files  specified by option input would be
    prepended before the specified script.
    """
    return SoS_ExecuteScript(script, interpreter or "octave", ".m", args, entrypoint).run(**kwargs)
