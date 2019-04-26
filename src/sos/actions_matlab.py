#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action, SoS_ExecuteScript


@SoS_Action(acceptable_args=['script', 'args'])
def matlab(
        script,
        args='''-nojvm -nodisplay -nosplash -nodesktop -r "try, run('{filename}'), catch me, fprintf('%s / %s\\n',me.identifier,me.message), exit(1), end, exit(0);"''',
        **kwargs):
    '''Execute specified script with command Matlab, with default options
    "-nojvm -nodisplay -nosplash -nodesktop -r". This action accepts common action arguments such as input,
    active, workdir, docker_image and args. In particular, content of one or more
    files  specified by option input would be
    prepended before the specified script.
    '''
    return SoS_ExecuteScript(script, 'matlab', '.m', args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def octave(script, args='', **kwargs):
    '''Execute specified script with command Matlab, with default options
    "-nodisplay -r". This action accepts common action arguments such as input,
    active, workdir, docker_image and args. In particular, content of one or more
    files  specified by option input would be
    prepended before the specified script.
    '''
    return SoS_ExecuteScript(script, 'octave', '.m', args).run(**kwargs)
