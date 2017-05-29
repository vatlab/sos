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


import argparse
import pandas
import numpy
from IPython.core.display import HTML

class Visualizer:
    def __init__(self, kernel, style):
        self.kernel = kernel
        if style is None:
            self.style = 'table'
            self.options = []
        else:
            self.style = style['style']
            self.options = style['options']

    def _parse_error(self, msg):
        self.kernel.warn(msg)

    def preview(self, df):
        if self.style == 'table':
            return self._handle_table(df)
        elif self.style == 'scatterplot':
            return self._handle_scatterplot(df)

    def get_tid(self):
        if not hasattr(self.kernel, '_tid'):
            self.kernel._tid = 1
        else:
            self.kernel._tid += 1
        return self.kernel._tid

    #
    # TABLE
    #
    def _get_table_parser(self):
        parser = argparse.ArgumentParser(prog='%preview -s table',
            description='''Preview Pandas DataFrame (or csv files) in tablular
            format with search and sort capacity''')
        parser.add_argument('-l', '--limit', type=int, default=200,
            help='''Limit the number of displayed records''')
        parser.error = self._parse_error
        return parser

    def _handle_table(self, df):
        parser = self._get_table_parser()
        try:
            args = parser.parse_args(self.options)
        except SystemExit:
            return

        if not isinstance(df, pandas.core.frame.DataFrame):
            raise ValuError('Not of DataFrame type')

        tid = self.get_tid()

        if df.shape[0] > args.limit:
            self.kernel.warn("Only the first {} of the {} records are previewed. Use option --limit to set a new limit.".format(args.limit, df.shape[0]))
        code = df.head(args.limit).to_html(index=True).replace('class=', 'id="dataframe_{}" class='.format(tid), 1)

        hr, rest = code.split('</tr>', 1)
        index_type = 'numeric' if isinstance(df.index, pandas.indexes.range.RangeIndex) else 'alphabetic'
        col_type = ['numeric' if numpy.issubdtype(x, numpy.number) else 'alphabetic' for x in df.dtypes]
        code = ''.join('''{} &nbsp; <i class="fa fa-sort" style="color:lightgray" onclick="sortDataFrame('{}', {}, '{}')"></th>'''.format(x,
            tid, idx,
            index_type if idx == 0 else col_type[idx-1]) if '<th' in x else x for idx,x in enumerate(hr.split('</th>')  )) + '</tr>' + rest

        code = """
    <div class='dataframe_container'>
    <input type="text" class='dataframe_input' id="search_{}" """.format(tid) + \
    """onkeyup="filterDataFrame('{}""".format(tid) + """')" placeholder="Search for names..">
    """ + code + '''</div>'''
        return {'text/html': HTML(code).data}, {}

    #
    # SCATTERPLOT
    #
    def get_scatterplot_parser(self):
        parser = argparse.ArgumentParser(prog='%preview -s scatterplot')

        parser.error = self._parse_error
        return parser


