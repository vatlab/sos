#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
#


###############################################################################
#                                                                             #
# To define new actions for your particular language                          #
#                                                                             #
# 1. Import SoS_Action and other needed functions from sos.actions            #
# 2. Define a function and decorate it with SoS_Action and appropriate        #
#    parameters. Please specify acceptable_args if the action accepts         #
#    a fixed set of arguments.                                                #
# 3. Modify setup.py and add the action to the [sos_actions] section          #
#    of entry_points.                                                         #
#                                                                             #
# The following is a very simple example that uses SoS_ExecuteScript to       #
# define action bash. Please refer to implementation of existing actions      #
# and the 'Extending_SoS' section of the SoS documentation                    #
# (vatlab.github.io/SOS/) for details.                                        #
#                                                                             #
###############################################################################

# from sos.actions import SoS_Action, SoS_ExecuteScript

# @SoS_Action(run_mode=['run', 'interactive'], acceptable_args=['script', 'args'])
# def bash(script, args='', **kwargs):
#    '''Execute specified script using bash.'''
#    return SoS_ExecuteScript(script, '/bin/bash', '.sh', args).run(**kwargs)

