#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)

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

ver = sys.version_info
if ver.major == 2 or (ver.major == 3 and (ver.minor, ver.micro) < (3, 0)):
    raise SystemError('SOS requires Python 3.3 or higher. Please upgrade your Python {}.{}.{}.'
        .format(ver.major, ver.minor, ver.micro))


__version__ = '0.5.1'
__py_version__ = '{}.{}.{}'.format(ver.major, ver.minor, ver.micro)
#
SOS_FULL_VERSION='{} for Python {}.{}.{}'.format(__version__, ver.major, ver.minor, ver.micro)
SOS_COPYRIGHT = '''SoS {} : Copyright (c) 2016 Bo Peng'''.format(__version__)
SOS_CONTACT = '''Please visit http://github.com/bpeng2000/SOS for more information.'''

# this import might fail when SoS is not properly installed
from .sos_script import *
from .utils import *
