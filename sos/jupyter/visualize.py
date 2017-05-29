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

from collections import defaultdict
from IPython.core.display import HTML
import json

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
    def _get_scatterplot_parser(self):
        parser = argparse.ArgumentParser(prog='%preview -s scatterplot')
        parser.add_argument('cols', nargs='*', help='''Columns to plot, which should all be numeric. If one
            column is specified, it is assumed to be a x-y plot with x being 1, 2, 3, .... If two or
            more columns (n) are specified, n-1 series will be plotted with the first column being the
            x axis, in which case an "_index" name can be used to specify 1, 2, 3, .... This option can be
            igured if the dataframe has only one or two columns.''')
        parser.add_argument('--ylim', nargs=2, help='''Range of y-axis''')
        parser.add_argument('--xlim', nargs=2, help='''Range of x-axis''')
        parser.add_argument('--width', default='800px', help='''Width of the plot.''')
        parser.add_argument('--height', default='600px', help='''Height of the plot.''')
        parser.add_argument('--show', nargs='+', help='''What to show in the plot,
            can be 'lines', 'points'. Default to both.''')
        parser.add_argument('-l', '--limit', default=2000, help='''Maximum number
            of records to plot.''')
        parser.error = self._parse_error
        return parser

    def _to_list(self, arr):
        if 'int' in arr.dtype.name:
            return [int(x) for x in arr]
        else:
            return [float(x) for x in arr]

    def _handle_scatterplot(self, df):
        parser = self._get_scatterplot_parser()
        try:
            args = parser.parse_args(self.options)
        except SystemExit:
            return

        if not isinstance(df, pandas.core.frame.DataFrame):
            raise ValuError('Not of DataFrame type')

        tid = self.get_tid()

        if df.shape[0] > args.limit:
            self.kernel.warn("Only the first {} of the {} records are plotted. Use option --limit to set a new limit.".format(args.limit, df.shape[0]))

        if not args.cols:
            if df.shape[1] == 1:
                args.cols = ['_index', df.columns[0]]
            elif df.shape[1] == 2:
                args.cols = df.columns
            else:
                raise ValueError('Please specify columns for plot. Available columns are {}'.format(
                    ' '.join(df.columns)))
        if len(args.cols) == 1:
            args.cols = ['_index', args.cols[0]]

        indexes = [str(x) for x in df.index]

        data = df.head(args.limit)
        nrow = data.shape[0]
        series = []

        # check datatype
        for col in args.cols:
            if col == '_index':
                continue
            if col not in data.columns:
                raise ValueError("Invalid column name {}".format(col))
            if not numpy.issubdtype(data[col].dtype, numpy.number):
                raise ValueError("Column {} is not of numeric type".format(col))

        all_series = []
        for col in args.cols[1:]:
            series = {}
            series['label'] = col
            if args.cols[0] == '_index':
                x = list(range(1, nrow + 1))
            else:
                x = self._to_list(data[args.cols[0]])
            if col == '_index':
                y = list(range(1, nrow + 1))
            else:
                y = self._to_list(data[col])
            series['data'] = [(_x,_y,_z) for _x,_y,_z in zip(x,y,indexes)]
            series['clickable'] = True
            series['hoverable'] = True
            all_series.append(series)

        options = defaultdict(dict)
        options['series']['lines'] = {'show': True if not args.show or 'lines' in args.show else False }
        options['series']['points'] = {'show': True if not args.show or 'points' in args.show else False }
        options['grid']['hoverable'] = True
        options['grid']['clickable'] = True

        if args.xlim:
            options['xaxis']['min'] = args.xlim[0]
            options['xaxis']['max'] = args.xlim[1]
        if args.ylim:
            options['yaxis']['min'] = args.ylim[0]
            options['yaxis']['max'] = args.ylim[1]

        code = """
    <div class='dataframe_container'>
    <div id="dataframe_scatterplot_{0}" width="{1}" height="{2}"></div>
    <script language="javascript" type="text/javascript" src="http://www.flotcharts.org/flot/jquery.flot.js"></script>
    <script>
    var plot = $.plot('dataframe_scatterplot_{0}', """.format(tid, args.width, args.height) + \
        json.dumps(all_series) + """, """ + json.dumps(options) + """)

        $("#dataframe_scatter_plot_{0}")""".format(tid) + """.bind("plothover", function (event, pos, item) {
        $("#tooltip").remove();
        if (item) {
            var tooltip = item.series.data[item.dataIndex][2];

            $('<div id="tooltip">' + tooltip + '</div>')
                .css({
                    position: 'absolute',
                    display: 'none',
                    top: item.pageY + 5,
                    left: item.pageX + 5,
                    border: '1px solid #fdd',
                    padding: '2px',
                    'background-color': '#fee',
                    opacity: 0.80 })
                .appendTo("body").fadeIn(200);


            showTooltip(item.pageX, item.pageY, tooltip);
        }
    });

        </script>
        </div>"""
        self.kernel.warn(code)
        return {'text/html': HTML(code).data}, {}
