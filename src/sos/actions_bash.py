#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action, SoS_ExecuteScript


@SoS_Action(acceptable_args=['script', 'args'])
def bash(script, args='', **kwargs):
    '''Execute specified script using bash. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['/bin/bash', 'bash'], '.sh',
                             args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def csh(script, args='', **kwargs):
    '''Execute specified script using csh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['/bin/csh', 'csh'], '.csh',
                             args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def tcsh(script, args='', **kwargs):
    '''Execute specified script using tcsh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['/bin/tcsh', 'tcsh'], '.sh',
                             args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def zsh(script, args='', **kwargs):
    '''Execute specified script using zsh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['/bin/zsh', 'zsh'], '.zsh',
                             args).run(**kwargs)


@SoS_Action(acceptable_args=['script', 'args'])
def sh(script, args='', **kwargs):
    '''Execute specified script using sh. This action accepts common action arguments such as
    input, active, workdir, docker_image and args. In particular, content of one or more files
    specified by option input would be prepended before the specified script.'''
    return SoS_ExecuteScript(script, ['/bin/sh', 'sh'], '.sh',
                             args).run(**kwargs)
