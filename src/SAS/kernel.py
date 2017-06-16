#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/vatlab/SOS for more information.
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
import pandas as pd
from sos.utils import env

class sos_SAS:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'sas_kernel'
        self.background_color = 'teal'
        self.init_statements = ''

    def get_vars(self, names):
        #
        # get variables with names from env.sos_dict and create
        # them in the subkernel. The current kernel should be SAS
        for name in names:
            if not isinstance(env.sos_dict[name], pd.DataFrame):
                self.sos_kernel.warn('Cannot transfer a non DataFrame object {} of type {} to SAS'.format(name, env.sos_dict[name].__class__.__name__))
                continue
            # write to csv and read from SAS

    def put_vars(self, items, to_kernel=None):
        #
        for item in items:
            pass
        #
        if to_kernel == 'R':
            # write to xport and let R read it
        else:
            # others, write and read in Python
            pd.read_sas()

    def sessioninfo(self):
        # return information of the kernel
        pass



