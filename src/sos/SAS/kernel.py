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
import os
import shutil
import pandas as pd
import argparse
import tempfile
from sos.utils import env

from saspy.sasiostdio import SASsessionSTDIO

from io import BytesIO

class sos_SAS(SASsessionSTDIO):
    supported_kernels = {'SAS': ['sas']}
    background_color = '#9CD4F9'
    options = {}

    def __init__(self, sos_kernel, kernel_name='sas'):
        self.sos_kernel = sos_kernel
        self.kernel_name = kernel_name
        self.init_statements = ''
        #
        # we intentionally do not call SASsessionSTDIO's constructor, which needs to read
        # saspy configurations etc and will try to start SAS
        #
        #self.pid    = None
        #self.stdin  = None
        #self.stderr = None
        #self.stdout = None

        #self.sascfg   = SASconfigSTDIO(**kwargs)
        #self._log_cnt = 0
        #self._log     = ""
        #self._sb      = kwargs.get('sb', None)

        #self._startsas()
        self.stdin = BytesIO()
        self.pid = None
        # mimick a sascfg object with argparse Namespace
        self.sascfg = argparse.Namespace()
        self.sascfg.encoding = 'utf-8'
        self.sascfg.output = 'html'
        self._sb = argparse.Namespace()
        self._sb._dsopts = lambda x: ''

    def get_vars(self, names):
        #
        # get variables with names from env.sos_dict and create
        # them in the subkernel. The current kernel should be SAS
        for name in names:
            if not isinstance(env.sos_dict[name], pd.DataFrame):
                if self.sos_kernel._debug_mode:
                    self.sos_kernel.warn('Cannot transfer a non DataFrame object {} of type {} to SAS'.format(
                        name, env.sos_dict[name].__class__.__name__))
                continue
            # sas cannot handle columns with non-string header
            data = env.sos_dict[name].rename(columns={x:str(x) for x in env.sos_dict[name].columns})
            # convert dataframe to SAS
            self.dataframe2sasdata(data, name, "")
            sas_code = self.stdin.getvalue().decode('utf-8')
            self.stdin.seek(0)
            self.stdin.truncate()
            if self.sos_kernel._debug_mode:
                self.sos_kernel.warn('Executing\n{}'.format(sas_code))
            self.sos_kernel.run_cell(sas_code, True, False, on_error='Failed to put variable {} to SAS'.format(name))

    def put_vars(self, items, to_kernel=None):
        # put SAS dataset to Python as dataframe
        temp_dir = tempfile.mkdtemp()
        res = {}
        try:
            for idx, item in enumerate(items):
                try:
                    code = '''\
libname TEMP '{}';
Data TEMP.dat_{};
   set {};
run;
'''.format(temp_dir, idx, item)
                    # run the code to save file
                    self.sos_kernel.run_cell(code, True, False, on_error="Failed to get data set {} from SAS".format(item))
                    # check if file exists
                    saved_file = os.path.join(temp_dir, 'dat_{}.sas7bdat'.format(idx))
                    if not os.path.isfile(saved_file):
                        self.sos_kernel.warn('Failed to save dataset {} to {}'.format(item, saved_file))
                        continue
                    # now try to read it with Python
                    df = pd.read_sas(saved_file, encoding='utf-8')
                    res[item] = df
                except Exception as e:
                    self.sos_kernel.warn('Failed to get dataset {} from SAS: {}'.format(item, e))
        finally:
            # clear up the temp dir
            shutil.rmtree(temp_dir)
        return res

    def parse_response(self, html):
        # separate response into LOG (with class err) and LST (with class s)
        LOG = ''
        LST = ''
        for line in html.split('</span>'):
            if 'class="err"' in line:
                LOG += line.replace('<br>', '\n').replace('<span class="err">', '') + ' '
            elif 'class="s"' in line:
                LST += line.replace('<br>', '\n').replace('<span class="s">', '') + ' '
        return {'LOG':LOG, 'LST':LST}
        
    def sessioninfo(self):
        # return information of the kernel
        sas_code = '''\
PROC PRODUCT_STATUS;
run;
'''
        response = self.sos_kernel.get_response(sas_code, ('execute_result',))[0][1]['data']['text/html']
        return self.parse_response(response)['LOG']


