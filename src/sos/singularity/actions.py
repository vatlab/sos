#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

from sos.actions import SoS_Action

from .client import SoS_SingularityClient


@SoS_Action(run_mode=['run', 'interactive'])
def singularity_build(script=None, src=None, dest=None, **kwargs):
    '''docker build command. By default a script is sent to the docker build command but
    you can also specify different parameters defined inu//docker-py.readthedocs.org/en/stable/api/#build
    '''
    singularity = SoS_SingularityClient()
    singularity.build(script, src, dest, **kwargs)
    return 0
