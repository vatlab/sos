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

from itertools import tee
import operator

def is_sorted(iterable, compare=operator.le):
  a, b = tee(iterable)
  next(b, None)
  return all(map(compare, a, b))

class Visualizer:
    def __init__(self, kernel, style):
        self.kernel = kernel
        if style is None:
            self.style = 'table'
            self.options = []
        else:
            self.style = 'table' if style['style'] is None else style['style']
            self.options = style['options']

    def _parse_error(self, msg):
        self.kernel.warn(msg)

    def preview(self, df):
        if self.style == 'table':
            return self._handle_table(df)
        elif self.style == 'scatterplot':
            return self._handle_scatterplot(df)
        else:
            raise ValueError('Unknown style {}'.format(self.style))

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
            column is specified, it is assumed to be a x-y plot with x being 0, 1, 2, 3, .... If two or
            more columns (n) are specified, n-1 series will be plotted with the first column being the
            x axis, in which case an "_index" name can be used to specify 0, 1, 2, 3, .... This option can be
            igured if the dataframe has only one or two columns.''')
        parser.add_argument('--ylim', nargs=2, help='''Range of y-axis''')
        parser.add_argument('--xlim', nargs=2, help='''Range of x-axis''')
        parser.add_argument('--width', default='90%', help='''Width of the plot.''')
        parser.add_argument('--height', default='400px', help='''Height of the plot.''')
        parser.add_argument('--show', nargs='+', help='''What to show in the plot,
            can be 'lines', 'points' or both. Default to points, and lines if x-axis is
            sorted.''')
        parser.add_argument('-t', '--tooltip', nargs='*', help='''Fields to be shown in tooltip, in addition to
            the row index and point values that would be shown by default.''')
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
                args.cols = list(df.columns)
            else:
                raise ValueError('Please specify columns for plot. Available columns are {}'.format(
                    ' '.join(df.columns)))
        if len(args.cols) == 1:
            args.cols = ['_index', args.cols[0]]

        indexes = [str(x) for x in df.index]

        data = df.head(args.limit)
        nrow = data.shape[0]
        series = []

        if args.tooltip:
            args.tooltip = ['_index'] + args.tooltip
        else:
            args.tooltip = ['_index']

        # check datatype
        for col in args.cols + args.tooltip:
            if col == '_index':
                continue
            if col not in data.columns:
                raise ValueError("Invalid column name {}".format(col))
            if not numpy.issubdtype(data[col].dtype, numpy.number):
                raise ValueError("Column {} is not of numeric type".format(col))

        all_series = []
        if args.cols[0] == '_index':
            val_x = list(range(0, nrow))
        else:
            val_x = self._to_list(data[args.cols[0]])
        for col in args.cols[1:]:
            series = {}
            series['label'] = col
            if col == '_index':
                val_y = list(range(1, nrow + 1))
            else:
                val_y = self._to_list(data[col])

            tooltip = ['<br>'.join(['{}: {}'.format('index' if t == '_index' else t,
                idxvalue if t == '_index' else df[t][idx] ) for t in args.tooltip])
                for idx,idxvalue in enumerate(indexes)]
            series['data'] = [(x, y, z) for x, y, z in zip(val_x, val_y, tooltip)]
            series['clickable'] = True
            series['hoverable'] = True
            all_series.append(series)

        options = defaultdict(dict)
        options['series']['lines'] = {'show': is_sorted(val_x) if not args.show or 'lines' in args.show else False }
        options['series']['points'] = {'show': True if not args.show or 'points' in args.show else False }
        options['grid']['hoverable'] = True
        options['grid']['clickable'] = True

        # if there are actual indexes... and plot by x
        optional_style = ''
        if args.cols[0] == '_index' and not isinstance(df.index, pandas.indexes.range.RangeIndex):
            options['xaxis']['ticks'] = [[x,str(y)] for x,y in enumerate(df.index)]
            optional_style = r'''
var css = document.createElement("style");
css.type = "text/css";
css.innerHTML = '' +
   '#dataframe_scatterplot_{0} div.xAxis div.tickLabel {{\n' +
   '  transform: translateY(15px) translateX(15px) rotate(45deg);\n' +
   '  -ms-transform: translateY(15px) translateX(15px) rotate(45deg);\n' +
   '  -moz-transform: translateY(15px) translateX(15px) rotate(45deg);\n' +
   '  -webkit-transform: translateY(15px) translateX(15px) rotate(45deg);\n' +
   '  -o-transform: translateY(15px) translateX(15px) rotate(45deg);\n' +
   '  /*rotation-point:50% 50%;*/\n' +
   '  /*rotation:270deg;*/\n' +
   '}}\n'
document.body.appendChild(css);
'''.format(tid)

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
    function plotScatterPlot{0}() {{
        plot = $.plot('#dataframe_scatterplot_{0}', """.format(tid, args.width, args.height) + \
        json.dumps(all_series) + """, """ + json.dumps(options) + """)

    if ($('#tooltip').length == 0) {{
        $('#dataframe_scatterplot_{0}')""".format(tid) + """.append($("<div id='tooltip'></div>").css({
        position: "absolute",
        display: "none",
        border: "1px solid #fdd",
        padding: "2px",
        "background-color": "#fee",
        opacity: 0.80
        }));
    }
    """ + """
    $("#dataframe_scatterplot_{0}")""".format(tid) + """.bind("plothover", function (event, pos, item) {
            if (item) {
                $("#tooltip").html((item.series.label + ": (" + 
                    item.series.data[item.dataIndex][0].toString() + ", " +
                    item.series.data[item.dataIndex][1].toString() + ")<br>" +
                    item.series.data[item.dataIndex][2]).trim()).css({top: item.pageY+5, left: item.pageX+5})
                    .css('display', 'inline').fadeIn(200);
            } else {
                $("#tooltip").hide();
            }
    });
    }
""" + """
// we will wait for the div to be displayed on the HTML/Jupyter side before we plot
// the figure. This might not be necessary but at least this is a chance for us
// to resize the div and avoid some flot error

var dt = 100;
// the frontend might be notified before the table is inserted as results.
function showFigure{0}() {{
    if ( $('#dataframe_scatterplot_{0}').length === 0 || $.plot === undefined ) {{
          dt = dt * 1.5; // slow-down checks for datatable as time goes on;
          setTimeout("showFigure{0}", dt);
          return;
    }} else {{
        $('#dataframe_scatterplot_{0}').css('width', "{1}").css('height', "{2}");
        plotScatterPlot{0}();
    }}
}}
{3}
showFigure{0}()
</script>
</div>""".format(tid, args.width, args.height, optional_style)
        return {'text/html': HTML(code).data}, {}
