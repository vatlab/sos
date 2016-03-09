#!/usr/bin/env python
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
#
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

import sys
from .sos_script import SoS_Script
from .utils import env

def sos_show(args, argv):
    try:
        script = SoS_Script(args.script, argv)
        workflow = script.workflow(args.workflow)
        print(workflow)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
    
def sos_run(args, argv):
    try:
        script = SoS_Script(args.script, argv)
        workflow = script.workflow(args.workflow)
        workflow.run()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
