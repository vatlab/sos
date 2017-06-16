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

from collections import Sequence
import tempfile
from sos.utils import short_repr, env
from IPython.core.error import UsageError

def homogeneous_type(seq):
    iseq = iter(seq)
    first_type = type(next(iseq))
    if first_type in (int, float):
        return True if all( (type(x) in (int, float)) for x in iseq ) else False
    else:
        return True if all( (type(x) is first_type) for x in iseq ) else False

#
#  support for %get
#
#  Converting a Python object to a R expression that will be executed
#  by the R kernel.
#
#
def _R_repr(obj):
    if isinstance(obj, bool):
        return 'TRUE' if obj else 'FALSE'
    elif isinstance(obj, (int, float, str)):
        return repr(obj)
    elif isinstance(obj, Sequence):
        if len(obj) == 0:
            return 'c()'
        # if the data is of homogeneous type, let us use c()
        # otherwise use list()
        # this can be confusion but list can be difficult to handle
        if homogeneous_type(obj):
            return 'c(' + ','.join(_R_repr(x) for x in obj) + ')'
        else:
            return 'list(' + ','.join(_R_repr(x) for x in obj) + ')'
    elif obj is None:
        return 'NULL'
    elif isinstance(obj, dict):
        return 'list(' + ','.join('{}={}'.format(x, _R_repr(y)) for x,y in obj.items()) + ')'
    elif isinstance(obj, set):
        return 'list(' + ','.join(_R_repr(x) for x in obj) + ')'
    else:
        import numpy
        import pandas
        if isinstance(obj, (numpy.intc, numpy.intp, numpy.int8, numpy.int16, numpy.int32, numpy.int64,\
                numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.float16, numpy.float32, \
                numpy.float64)):
            return repr(obj)
        elif isinstance(obj, numpy.matrixlib.defmatrix.matrix):
            try:
                import feather
            except ImportError:
                raise UsageError('The feather-format module is required to pass numpy matrix as R matrix'
                    'See https://github.com/wesm/feather/tree/master/python for details.')
            feather_tmp_ = tempfile.NamedTemporaryFile(suffix='.feather', delete=False).name
            feather.write_dataframe(pandas.DataFrame(obj).copy(), feather_tmp_)
            return 'data.matrix(..read.feather({!r}))'.format(feather_tmp_)
        elif isinstance(obj, numpy.ndarray):
            return 'c(' + ','.join(_R_repr(x) for x in obj) + ')'
        elif isinstance(obj, pandas.DataFrame):
            try:
                import feather
            except ImportError:
                raise UsageError('The feather-format module is required to pass pandas DataFrame as R data.frame'
                    'See https://github.com/wesm/feather/tree/master/python for details.')
            feather_tmp_ = tempfile.NamedTemporaryFile(suffix='.feather', delete=False).name
            try:
                data = obj.copy()
                # if the dataframe has index, it would not be transferred due to limitations
                # of feather. We will have to do something to save the index separately and 
                # recreate it. (#397)
                if isinstance(data.index, pandas.Index):
                    df_index = list(data.index)
                elif not isinstance(data.index, pandas.RangeIndex):
                    # we should give a warning here
                    df_index=None
                feather.write_dataframe(data, feather_tmp_)
            except:
                # if data cannot be written, we try to manipulate data
                # frame to have consistent types and try again
                for c in data.columns:
                    if not homogeneous_type(data[c]):
                        data[c] = [str(x) for x in data[c]]
                feather.write_dataframe(data, feather_tmp_)
                # use {!r} for path because the string might contain c:\ which needs to be
                # double quoted.
            return '..read.feather({!r}, index={})'.format(feather_tmp_, _R_repr(df_index))
        else:
            return repr('Unsupported datatype {}'.format(short_repr(obj)))



# R    length (n)    Python
# NULL        None
# logical    1    boolean
# integer    1    integer
# numeric    1    double
# character    1    unicode
# logical    n > 1    array
# integer    n > 1    array
# numeric    n > 1    list
# character    n > 1    list
# list without names    n > 0    list
# list with names    n > 0    dict
# matrix    n > 0    array
# data.frame    n > 0    DataFrame

R_init_statements = r'''
..py.repr.logical.1 <- function(obj) {
    if(obj)
        'True'
    else
        'False'
}
..py.repr.integer.1 <- function(obj) {
    as.character(obj)
}
..py.repr.double.1 <- function(obj) {
    as.character(obj)
}
..py.repr.character.1 <- function(obj) {
    paste0('r"""', obj, '"""')
}
..has.row.names <- function(df) {
  !all(row.names(df)==seq(1, nrow(df)))
}
..py.repr.dataframe <- function(obj) {
    if (!require("feather")) {
        install.packages('feather', repos='http://cran.stat.ucla.edu/')
        }
    library(feather)
    tf = tempfile('feather')
    write_feather(obj, tf)
    if (..has.row.names(obj)) {
        paste0("read_dataframe(r'", tf, "').set_index([", ..py.repr(row.names(obj)),"])")
    } else {
        paste0("read_dataframe(r'", tf, "')")
    }
}
..py.repr.matrix <- function(obj) {
    if (!require("feather")) {
        install.packages('feather', repos='http://cran.stat.ucla.edu/')
        }
    library(feather)
    tf = tempfile('feather')
    write_feather(as.data.frame(obj), tf)
    if (..has.row.names(obj)) {
       paste0("read_dataframe(r'", tf, "').set_index([", ..py.repr(row.names(obj)),"]).as_matrix()")
    } else {
       paste0("read_dataframe(r'", tf, "').as_matrix()")
    }
}
..py.repr.n <- function(obj) {
    paste("[",
        paste(sapply(obj, ..py.repr), collapse=','),
        "]")
}
..py.repr <- function(obj) {
    if (is.matrix(obj)) {
        ..py.repr.matrix(obj)
    } else if (is.data.frame(obj)) {
        ..py.repr.dataframe(obj)
    } else if (is.null(obj)) {
        'None'
    } else if (is.integer(obj)) {
        if (length(obj) == 1)
            ..py.repr.integer.1(obj)
        else
            paste("[", paste(obj, collapse=','), "]")
    } else if (is.double(obj)){
        if (length(obj) == 1)
            ..py.repr.double.1(obj)
        else
            paste("[", paste(obj, collapse=','), "]")
    } else if (is.character(obj)) {
        if (length(obj) == 1)
            ..py.repr.character.1(obj)
        else {
            paste("[", paste(sapply(obj, ..py.repr.character.1), collapse=','), "]")
        }
    } else if (is.logical(obj)) {
        if (length(obj) == 1)
            ..py.repr.logical.1(obj)
        else
            ..py.repr.n(obj)
    } else if (is.list(obj)) {
        # if the list has no name
        if (is.null(names(obj)))
            ..py.repr.n(obj)
        else {
            paste("{",
                  paste(sapply(names(obj), function (x)
                      paste(shQuote(gsub("\\.", "_", as.character(x))), ":", ..py.repr(obj[[x]]))),
                      collapse=','),
                  "}")
        }
    } else {
        "'Untransferrable variable'"
    }
}
..read.feather <- function(filename, index=NULL) {
    if (! suppressMessages(suppressWarnings(require("feather", quietly = TRUE)))) {
        try(install.packages('feather', repos='http://cran.stat.ucla.edu/'), silent=TRUE)
        if (!suppressMessages(suppressWarnings(require("feather"))))
            stop('Failed to install feather library')
    }
    suppressPackageStartupMessages(library(feather, quietly = TRUE))
    data = as.data.frame(read_feather(filename))
    if (!is.null(index))
        rownames(data) <- index
    return(data)
}
'''


class sos_R:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'ir'
        self.background_color = '#FDEDEC'
        self.init_statements = R_init_statements

    def get_vars(self, names):
        for name in names:
            if name.startswith('_'):
                self.sos_kernel.warn('Variable {} is passed from SoS to kernel {} as {}'.format(name, self.kernel_name, '.' + name[1:]))
                newname = '.' + name[1:]
            else:
                newname = name
            r_repr = _R_repr(env.sos_dict[name])
            self.sos_kernel.run_cell('{} <- {}'.format(newname, r_repr), True, False, on_error='Failed to get variable {} to R'.format(name))

    def put_vars(self, items, to_kernel=None):
        # first let us get all variables with names starting with sos
        response = self.sos_kernel.get_response('cat(..py.repr(ls()))', ('stream',), name=('stdout',))[0][1]
        all_vars = eval(response['text'])
        all_vars = [all_vars] if isinstance(all_vars, str) else all_vars

        items += [x for x in all_vars if x.startswith('sos')]

        for item in items:
            if '.' in item:
                self.sos_kernel.warn('Variable {} is put to SoS as {}'.format(item, item.replace('.', '_')))

        if not items:
            return {}

        py_repr = 'cat(..py.repr(list({})))'.format(','.join('{0}={0}'.format(x) for x in items))
        response = self.sos_kernel.get_response(py_repr, ('stream',), name=('stdout',))[0][1]
        expr = response['text']

        if to_kernel in ('Python2', 'Python3'):
            # directly to python3
            return '{}\nglobals().update({})'.format('from feather import read_dataframe\n' if 'read_dataframe' in expr else '', expr)
        # to sos or any other kernel
        else:
            # irkernel (since the new version) does not produce execute_result, only
            # display_data
            try:
                if 'read_dataframe' in expr:
                    # imported to be used by eval
                    from feather import read_dataframe
                    # suppress flakes warning
                    read_dataframe
                # evaluate as raw string to correctly handle \\ etc
                return eval(expr)
            except Exception as e:
                self.sos_kernel.warn('Failed to evaluate {!r}: {}'.format(expr, e))
                return None

    def sessioninfo(self):
        response = self.sos_kernel.get_response(r'cat(paste(capture.output(sessionInfo()), collapse="\n"))', ('stream',), name=('stdout',))[0]
        return response[1]['text']
