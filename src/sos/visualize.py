#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.
import argparse
import json
import operator
from collections import defaultdict
from itertools import tee

import numpy
import pandas


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
        if self.kernel:
            self.kernel.warn(msg)

    def preview(self, df):
        if self.style == 'table':
            return self._handle_table(df)
        elif self.style == 'scatterplot':
            return self._handle_scatterplot(df)
        else:
            raise ValueError(f'Unknown style {self.style}')

    def get_tid(self, vis_type):
        if not self.kernel:
            import random
            return random.randint(0, 1000000)
        if not hasattr(self.kernel, '_tid'):
            self.kernel._tid = {vis_type: 1}
        elif vis_type not in self.kernel._tid:
            self.kernel._tid[vis_type] = 1
        else:
            self.kernel._tid[vis_type] += 1
        return self.kernel._tid[vis_type]

    #
    # TABLE
    #
    def _get_table_parser(self):
        parser = argparse.ArgumentParser(
            prog='%preview -s table',
            description='''Preview Pandas DataFrame (or csv files) in tablular
            format with search and sort capacity''')
        parser.add_argument(
            '-l',
            '--limit',
            type=int,
            default=200,
            help='''Limit the number of displayed records. A negative value indicates
            showing all records.''')
        parser.error = self._parse_error
        return parser

    def _is_numeric_type(self, x):
        try:
            return numpy.issubdtype(x, numpy.number)
        except Exception:
            # this includes Pandas category type
            return False

    def _handle_table(self, df):
        parser = self._get_table_parser()
        try:
            args = parser.parse_args(self.options)
        except SystemExit:
            return

        if not isinstance(df, pandas.core.frame.DataFrame):
            raise ValueError('Not of DataFrame type')

        tid = self.get_tid('table')

        hint = ''
        if args.limit >= 0 and df.shape[0] > args.limit and args.limit == 200:
            hint += f'<div class="sos_hint">Only the first {args.limit} of the {df.shape[0]} records are previewed. Use option --limit to set a new limit.</div><br>'
        if args.limit >= 0:
            code = df.head(args.limit).to_html(index=True).replace(
                'class="', f'id="dataframe_{tid}" class="sos_dataframe ', 1)
        else:
            code = df.to_html(index=True).replace(
                'class="', f'id="dataframe_{tid}" class="sos_dataframe ', 1)

        hr, rest = code.split('</tr>', 1)
        index_type = 'numeric' if isinstance(
            df.index, pandas.RangeIndex) else 'alphabetic'
        col_type = [
            'numeric' if self._is_numeric_type(x) else 'alphabetic'
            for x in df.dtypes
        ]
        code = ''.join(
            '''{} &nbsp; <i class="fa fa-sort" style="color:lightgray" onclick="sortDataFrame('{}', {}, '{}')"></th>'''
            .format(x, tid, idx, index_type if idx ==
                    0 else col_type[idx - 1]) if '<th' in x else x
            for idx, x in enumerate(hr.split('</th>'))) + '</tr>' + rest

        # we put max-height 400px here because the notebook could be exported without using sos template
        # and associated css, resulting in very long table.
        code = f"""
    <div class='dataframe_container' style="max-height:400px">
    <input type="text" class='dataframe_input' id="search_{tid}" """ + \
            f"""onkeyup="filterDataFrame('{tid}""" + """')" placeholder="Search for names..">
    """ + code + '''</div>'''
        return {'text/html': hint + code}

    #
    # SCATTERPLOT
    #
    def _get_scatterplot_parser(self):
        parser = argparse.ArgumentParser(prog='%preview -s scatterplot')
        parser.add_argument(
            'cols',
            nargs='*',
            help='''Columns to plot, which should all be numeric. If one
            column is specified, it is assumed to be a x-y plot with x being 0, 1, 2, 3, .... If two or
            more columns (n) are specified, n-1 series will be plotted with the first column being the
            x axis, in which case an "_index" name can be used to specify 0, 1, 2, 3, .... This option can be
            igured if the dataframe has only one or two columns.''')
        parser.add_argument('--ylim', nargs=2, help='''Range of y-axis''')
        parser.add_argument('--xlim', nargs=2, help='''Range of x-axis''')
        parser.add_argument(
            '--log',
            choices=['x', 'y', 'xy', 'yx'],
            help='''Make x-axis, y-axis, or both to logarithmic''')
        parser.add_argument(
            '--width', default='50vw', help='''Width of the plot.''')
        parser.add_argument(
            '--height', default='38vw', help='''Height of the plot.''')
        parser.add_argument(
            '-b',
            '--by',
            nargs='+',
            help='''columns by which the data points are stratified.''')
        parser.add_argument(
            '--show',
            nargs='+',
            help='''What to show in the plot,
            can be 'lines', 'points' or both. Default to points, and lines if x-axis is
            sorted.''')
        parser.add_argument(
            '-t',
            '--tooltip',
            nargs='*',
            help='''Fields to be shown in tooltip, in addition to
            the row index and point values that would be shown by default.''')
        parser.add_argument(
            '-l',
            '--limit',
            default=2000,
            help='''Maximum number
            of records to plot.''')
        parser.error = self._parse_error
        return parser

    def _to_list(self, arr):
        if 'int' in arr.dtype.name:
            return [int(x) for x in arr]
        else:
            return [float(x) for x in arr]

    def _natural_ticks(self, rg):
        # given a range, get natural ticks such as 0.1, 1, 10, 100, ...
        import math
        # get integer by 1 or 2 (1, 10) or (1, 100)
        logl = math.floor(math.log10(rg[0]))
        logh = math.ceil(math.log10(rg[1]))
        # small scale, let flop decide
        if logh - logl < 3:
            return None
        return list(10**x for x in range(logl, logh + 1))

    def _handle_scatterplot(self, df):
        parser = self._get_scatterplot_parser()
        try:
            args = parser.parse_args(self.options)
        except SystemExit:
            return

        if not isinstance(df, pandas.core.frame.DataFrame):
            raise ValueError('Not of DataFrame type')

        tid = str(self.get_tid('scatterplot'))

        hint: str = ''
        if df.shape[0] > args.limit:
            hint += f'<div class="sos_hint">Only the first {args.limit} of the {df.shape[0]} records are plotted. Use option --limit to set a new limit.</div><br>'

        # replacing ' ' with &nbsp and '-' with unicode hyphen will disallow webpage to separate words
        # into lines
        indexes = [
            str(x).replace(' ', '&nbsp;').replace('-', '&#8209;')
            for x in df.index
        ]

        data = df.head(args.limit)
        nrow = data.shape[0]

        if not args.cols:
            if df.shape[1] == 1:
                args.cols = ['_index', df.columns[0]]
            elif df.shape[1] == 2:
                args.cols = list(df.columns)
            else:
                raise ValueError(
                    f'Please specify columns for plot. Available columns are {" ".join(df.columns)}'
                )

        if len(args.cols) == 1:
            args.cols = ['_index', args.cols[0]]

        if args.cols[0] == '_index':
            args.tooltip = [args.cols[0]] + \
                (args.tooltip if args.tooltip else [])
        else:
            args.tooltip = ['_index', args.cols[0]] + \
                (args.tooltip if args.tooltip else [])

        # check datatype
        for col in args.cols + args.tooltip + (args.by if args.by else []):
            if col == '_index':
                continue
            if col not in data.columns:
                raise ValueError(f"Invalid column name {col}")

        # tooltip and --by columns does not have to be numeric
        for col in args.cols:
            if col == '_index':
                continue
            if not self._is_numeric_type(data[col].dtype):
                raise ValueError(f"Column {col} is not of numeric type")

        if args.cols[0] == '_index':
            val_x = list(range(0, nrow))
        else:
            val_x = self._to_list(data[args.cols[0]])

        all_series = []

        if args.by:
            # create seris with _by
            vals = []
            for by_col in args.by:
                vals.append(sorted(list(set(data[by_col]))))
            # outer product
            import itertools
            categories = list(itertools.product(*vals))

        for col in args.cols[1:]:
            series = {}
            series['clickable'] = True
            series['hoverable'] = True

            if col == '_index':
                val_y = list(range(1, nrow + 1))
            else:
                val_y = self._to_list(data[col])

            tooltip = [
                '<br>'.join([
                    f'{"index" if t == "_index" else t}: {idxvalue if t == "_index" else df[t][idx]}'
                    for t in args.tooltip
                ])
                for idx, idxvalue in enumerate(indexes)
            ]

            all_data = [(x, y, z) for x, y, z in zip(val_x, val_y, tooltip)]

            if args.by:
                for cat in categories:
                    series = {}
                    series['label'] = col + \
                        ' (' + ' '.join(f'{x}={y}' for x,
                                        y in zip(args.by, cat)) + ')'
                    # find index of values that falls into the category
                    series['data'] = [
                        all_data[i]
                        for i in range(len(all_data))
                        if all(data[b][i] == v for b, v in zip(args.by, cat))
                    ]
                    if len(series['data']) > 0:
                        all_series.append(series)
            else:
                series['label'] = col
                series['data'] = all_data
                all_series.append(series)

        options = defaultdict(dict)
        options['xaxis'] = {}
        options['yaxis'] = {}
        options['series']['lines'] = {
            'show':
                is_sorted(val_x) and not args.by
                if not args.show or 'lines' in args.show else False
        }
        options['series']['points'] = {
            'show': True if not args.show or 'points' in args.show else False
        }
        options['grid']['hoverable'] = True
        options['grid']['clickable'] = True

        # if there are actual indexes... and plot by x
        class_name = 'scatterplot'
        if args.cols[0] == '_index' and not isinstance(df.index,
                                                       pandas.RangeIndex):
            options['xaxis']['ticks'] = [
                [x, str(y)] for x, y in enumerate(indexes)
            ]
            class_name = 'scatterplot_by_rowname'

        if args.xlim:
            options['xaxis']['min'] = args.xlim[0]
            options['xaxis']['max'] = args.xlim[1]
        if args.ylim:
            options['yaxis']['min'] = args.ylim[0]
            options['yaxis']['max'] = args.ylim[1]

        optfunc = ''
        if args.log and 'x' in args.log:
            range_x = [min(val_x), min(val_x)]
            optfunc = '''
                options['xaxis']['transform'] = function(v) { return Math.log(v); }
                options['xaxis']['inverseTransform'] = function(v) { return Math.exp(v); }
            '''
            ticks = self._natural_ticks(range_x)
            if ticks:
                optfunc += f'''
                    options['xaxis']['ticks'] = {ticks!r};
                    '''
            if not args.xlim:
                options['xaxis']['min'] = range_x[0]
                options['xaxis']['max'] = range_x[1]
        if args.log and 'y' in args.log:
            range_y = [
                min([
                    min([x[1] for x in series['data']]) for series in all_series
                ]),
                max([
                    max([x[1] for x in series['data']]) for series in all_series
                ])
            ]
            optfunc += '''
                options['yaxis']['transform'] = function(v) { return Math.log(v); };
                options['yaxis']['inverseTransform'] = function(v) { return Math.exp(v); };
            '''
            ticks = self._natural_ticks(range_y)
            if ticks:
                optfunc += f'''
                    options['yaxis']['ticks'] = {ticks!r};
                    '''
            # flot does not seems to scale correctly without min/max
            if not args.ylim:
                options['yaxis']['min'] = range_y[0]
                options['yaxis']['max'] = range_y[1]
        code = """
<div class='scatterplot_container'>
<div class='""" + class_name + """' id='dataframe_scatterplot_""" + tid + """' width='""" + args.width + """' height='""" + args.height + """'></div>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script language="javascript" type="text/javascript" src="https://www.flotcharts.org/flot/jquery.flot.js"></script>
<script>
    var options = """ + json.dumps(options) + """;""" + optfunc + """
    function plotScatterPlot""" + tid + """() {
        plot = $.plot('#dataframe_scatterplot_""" + tid + """', """ + json.dumps(
            all_series) + """, options)

    if ($('#dftooltip').length == 0) {
        $("<div id='dftooltip'></div>").css({
            position: "absolute",
            display: "none",
            border: "1px solid #fdd",
            padding: "2px",
            "background-color": "#fee",
            "z-index": 1000,
            opacity: 0.80
            }).appendTo("body");
    }
    $('#dataframe_scatterplot_""" + tid + """').bind("plothover", function (event, pos, item) {
            if (item) {
                $("#dftooltip").html((item.series.label + ": " +
                    item.series.data[item.dataIndex][1].toString() + "<br>" +
                    item.series.data[item.dataIndex][2]).trim()).css({top: item.pageY+5, left: item.pageX+5})
                    .fadeIn(200);
            } else {
                $("#dftooltip").hide();
            }
    });
    }
""" + """
// we will wait for the div to be displayed on the HTML/Jupyter side before we plot
// the figure. This might not be necessary but at least this is a chance for us
// to resize the div and avoid some flot error

var dt = 100;
// the frontend might be notified before the table is inserted as results.
function showFigure""" + tid + """() {
    if ( $('#dataframe_scatterplot_""" + tid + """').length === 0 || $.plot === undefined ) {
          dt = dt * 1.5; // slow-down checks for datatable as time goes on;
          setTimeout(showFigure""" + tid + """, dt);
          return;
   } else {
        $('#dataframe_scatterplot_""" + tid + """').css('width', '""" + args.width + """').css('height', '""" + args.height + """');
        plotScatterPlot""" + tid + """();
    }
}
showFigure""" + tid + """()
</script>
</div>"""
        return {'text/html': hint + code}
