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
import argparse
from sos.utils import env

from saspy.sasiostdio import SASsessionSTDIO

from io import BytesIO

class sos_SAS(SASsessionSTDIO):
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'sas'
        self.background_color = '#dcb9b9'
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

    def get_vars(self, names):
        #
        # get variables with names from env.sos_dict and create
        # them in the subkernel. The current kernel should be SAS
        for name in names:
            if not isinstance(env.sos_dict[name], pd.DataFrame):
                if self.sos_kernel._debug_mode:
                    self.sos_kernel.warn('Cannot transfer a non DataFrame object {} of type {} to SAS'.format(name, env.sos_dict[name].__class__.__name__))
                continue
            # sas cannot handle columns with non-string header
            data = env.sos_dict[name].rename(columns={x:str(x) for x in env.sos_dict[name].columns})
            # convert dataframe to SAS
            self.dataframe2sasdata(data, name, "")
            sas_code = self.stdin.getvalue().decode('utf-8')
            self._reset_stdin()
            if self.sos_kernel._debug_mode:
                self.sos_kernel.warn('Executing\n{}'.format(sas_code))
            self.sos_kernel.run_cell(sas_code, True, False, on_error='Failed to put variable {} to SAS'.format(name))

    def put_vars(self, items, to_kernel=None):
        # put SAS dataset to Python as dataframe
        for item in items:
            # separate function from SASsessionSTDIO
            self.sasdata2dataframe_code(item)
            code = self.stdin.getvalue().decode('utf-8')
            self.sos_kernel.warn(code)
        return {}

    def sessioninfo(self):
        # return information of the kernel
        pass

    def _reset_stdin(self):
       self.stdin.seek(0)
       self.stdin.truncate()

    def sasdata2dataframe_code(self, table: str, libref: str ='', dsopts: dict ={}, **kwargs) -> '<Pandas Data Frame object>':
        '''
        This method exports the SAS Data Set to a Pandas Data Frame, returning the Data Frame object.
        table    - the name of the SAS Data Set you want to export to a Pandas Data Frame
        libref  - the libref for the SAS Data Set.
        port     - port to use for socket. Defaults to 0 which uses a random available ephemeral port
        '''

        if libref:
            tabname = libref+"."+table
        else:
            tabname = table

        code  = "proc sql; create view sasdata2dataframe as select * from "+tabname+self._sb._dsopts(dsopts)+";quit;\n"
        code += "data _null_; file STDERR;d = open('sasdata2dataframe');\n"
        code += "lrecl = attrn(d, 'LRECL'); nvars = attrn(d, 'NVARS');\n"
        code += "lr='LRECL='; vn='VARNUMS='; vl='VARLIST='; vt='VARTYPE='; vf='VARFMT=';\n"
        code += "put lr lrecl; put vn nvars; put vl;\n"
        code += "do i = 1 to nvars; var = varname(d, i); put var; end;\n"
        code += "put vt;\n"
        code += "do i = 1 to nvars; var = vartype(d, i); put var; end;\n"
        code += "run;"

        ll = self.submit(code, "text")

        l2 = ll['LOG'].rpartition("LRECL= ")
        l2 = l2[2].partition("\n")
        #lrecl = int(l2[0])

        l2 = l2[2].partition("VARNUMS= ")
        l2 = l2[2].partition("\n")
        nvars = int(l2[0])

        l2 = l2[2].partition("\n")
        varlist = l2[2].split("\n", nvars)
        del varlist[nvars]

        l2 = l2[2].partition("VARTYPE=")
        l2 = l2[2].partition("\n")
        vartype = l2[2].split("\n", nvars)
        del vartype[nvars]

        topts                 = dict(dsopts)
        topts['obs']        = 1
        topts['firstobs'] = ''

        code  = "data _null_; set "+tabname+self._sb._dsopts(topts)+";put 'FMT_CATS=';\n"
        for i in range(nvars):
            code += "_tom = vformatn('"+varlist[i]+"'n);put _tom;\n"
        code += "run;\n"
        return code

    def test(self, code, nvars, socks, port, tabname, dsopts, varlist, sas_date_fmts, sas_datetime_fmts, sas_time_fmts, vartype):
        ll = self.submit(code, "text")

        l2 = ll['LOG'].rpartition("FMT_CATS=")
        l2 = l2[2].partition("\n")
        varcat = l2[2].split("\n", nvars)
        del varcat[nvars]

        try:
            sock = socks.socket()
            sock.bind(("",port))
            port = sock.getsockname()[1]
        except OSError:
            print('Error try to open a socket in the sasdata2dataframe method. Call failed.')
            return None

        if self.sascfg.ssh:
            host = socks.gethostname()
        else:
            host = ''

        code  = ""
        code += "filename sock socket '"+host+":"+str(port)+"' lrecl=32767 recfm=v termstr=LF;\n"
        code += " data _null_; set "+tabname+self._sb._dsopts(dsopts)+";\n file sock; put "
        for i in range(nvars):
            code += "'"+varlist[i]+"'n "
            if vartype[i] == 'N':
                if varcat[i] in sas_date_fmts:
                    code += 'E8601DA10. '
                else:
                    if varcat[i] in sas_time_fmts:
                        code += 'E8601TM15.6 '
                    else:
                        if varcat[i] in sas_datetime_fmts:
                            code += 'E8601DT26.6 '
                        else:
                            code += 'best32. '
            if i < (len(varlist)-1):
                code += "'09'x "
        code += "; run;\n"

        sock.listen(0)
        self._asubmit(code, 'text')
        newsock = sock.accept()

        datar = ""
        while True:
            data = newsock[0].recv(4096)
            if len(data):
                datar += data.decode(self.sascfg.encoding)
            else:
                break

        newsock[0].shutdown(socks.SHUT_RDWR)
        newsock[0].close()
        sock.close()

        r = []
        for i in datar.splitlines():
            r.append(tuple(i.split(sep='\t')))

        df = pd.DataFrame.from_records(r, columns=varlist)

        for i in range(nvars):
            if vartype[i] == 'N':
                if varcat[i] not in sas_date_fmts + sas_time_fmts + sas_datetime_fmts:
                    df[varlist[i]] = pd.to_numeric(df[varlist[i]], errors='coerce')
                else:
                    df[varlist[i]] = pd.to_datetime(df[varlist[i]], errors='ignore')

        return df

