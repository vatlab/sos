#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from sos.actions import SoS_Action
from .client import DockerClient

@SoS_Action(run_mode=['run', 'interactive'])
def docker_build(dockerfile=None, **kwargs):
    '''docker build command. By default a script is sent to the docker build command but
    you can also specify different parameters defined inu//docker-py.readthedocs.org/en/stable/api/#build
    '''
    docker = DockerClient()
    docker.build(dockerfile, **kwargs)
    return 0

def docker_commit(**kwargs):
    docker = DockerClient()
    docker.commit(**kwargs)
    return 0

