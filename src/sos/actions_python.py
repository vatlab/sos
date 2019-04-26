#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action, SoS_ExecuteScript


@SoS_Action(acceptable_args=['script', 'args'])
def python(script, args='', **kwargs):
    '''Execute specified script using python (which can be python 2 or 3 depending on
    system configuration. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, 'python', '.py', args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def python2(script, args='', **kwargs):
    '''Execute specified script using python2, and python if python2 does
    not exist. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['python2', 'python2.7', 'python'], '.py',
                             args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def python3(script, args='', **kwargs):
    '''Execute specified script using python3, and python if python3 does
    not exist. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['python3', 'python'], '.py',
                             args).run(**kwargs)
