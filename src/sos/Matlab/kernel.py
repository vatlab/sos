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

from collections import Sequence
import tempfile
from sos.utils import short_repr, env
from IPython.core.error import UsageError

def homogeneous_type(seq):
    iseq = iter(seq)
    first_type = type(next(iseq))
    if first_type in (int, float):
        return True if all(isinstance(x, (int, float)) for x in iseq) else False
    else:
        return True if all(isinstance(x, first_type) for x in iseq) else False

#
#  support for %get
#
#  Converting a Python object to a Matlab expression that will be executed
#  by the Matlab kernel.
#
#
def _Matlab_repr(obj):
    if isinstance(obj, bool):
        return 'true' if obj else 'false'
    elif isinstance(obj, (int, float, str)):
        return repr(obj)
    elif isinstance(obj, Sequence):
        if len(obj) == 0:
            return '[]'
        # if the data is of homogeneous type, let us use []
        if homogeneous_type(obj):
            return '[' + ','.join(_Matlab_repr(x) for x in obj) + ']'
        else:
            return '{' + ','.join(_Matlab_repr(x) for x in obj) + '}'
    elif obj is None:
        return 'NaN'
    elif isinstance(obj, dict):
        return 'struct(' + ','.join('{},{}'.format(_Matlab_repr(x), _Matlab_repr(y)) for x,y in obj.items()) + ')'
    elif isinstance(obj, set):
        return '{' + ','.join(_Matlab_repr(x) for x in obj) + '}'
    else:
        import numpy
        import pandas
        if isinstance(obj, (numpy.intc, numpy.intp, numpy.int8, numpy.int16, numpy.int32, numpy.int64,\
                numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.float16, numpy.float32, \
                numpy.float64)):
            return repr(obj)
        elif isinstance(obj, numpy.matrixlib.defmatrix.matrix):
            try:
                import mat
            except ImportError:
                raise UsageError('The mat-format module is required to pass numpy matrix as R matrix'
                    'See https://github.com/wesm/mat/tree/master/python for details.')
            mat_tmp_ = tempfile.NamedTemporaryFile(suffix='.mat', delete=False).name
            mat.write_dataframe(pandas.DataFrame(obj).copy(), mat_tmp_)
            return 'data.matrix(..read.mat({!r}))'.format(mat_tmp_)
        elif isinstance(obj, numpy.ndarray):
            return 'c(' + ','.join(_Matlab_repr(x) for x in obj) + ')'
        elif isinstance(obj, pandas.DataFrame):
            try:
                import mat
            except ImportError:
                raise UsageError('The mat-format module is required to pass pandas DataFrame as R data.frame'
                    'See https://github.com/wesm/mat/tree/master/python for details.')
            mat_tmp_ = tempfile.NamedTemporaryFile(suffix='.mat', delete=False).name
            try:
                data = obj.copy()
                # if the dataframe has index, it would not be transferred due to limitations
                # of mat. We will have to do something to save the index separately and
                # recreate it. (#397)
                if isinstance(data.index, pandas.Index):
                    df_index = list(data.index)
                elif not isinstance(data.index, pandas.RangeIndex):
                    # we should give a warning here
                    df_index=None
                mat.write_dataframe(data, mat_tmp_)
            except Exception:
                # if data cannot be written, we try to manipulate data
                # frame to have consistent types and try again
                for c in data.columns:
                    if not homogeneous_type(data[c]):
                        data[c] = [str(x) for x in data[c]]
                mat.write_dataframe(data, mat_tmp_)
                # use {!r} for path because the string might contain c:\ which needs to be
                # double quoted.
            return '..read.mat({!r}, index={})'.format(mat_tmp_, _Matlab_repr(df_index))
        else:
            return repr('Unsupported datatype {}'.format(short_repr(obj)))



# Matlab    length (n)    Python
# NULL        NaN
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

Matlab_init_statements = r'''
__py_repr_logical_1 <- function(obj) {
    if(obj)
        'True'
    else
        'False'
}
__py_repr_integer_1 <- function(obj) {
    as_character(obj)
}
__py_repr_double_1 <- function(obj) {
    as_character(obj)
}
__py_repr_character_1 <- function(obj) {
    paste0('r"""', obj, '"""')
}
__has_row_names <- function(df) {
  !all(row_names(df)==seq(1, nrow(df)))
}
__py_repr_dataframe <- function(obj) {
    if (!require("mat")) {
        install_packages('mat', repos='http://cran_stat_ucla_edu/')
        }
    library(mat)
    tf = tempfile('mat')
    write_mat(obj, tf)
    if (__has_row_names(obj)) {
        paste0("read_dataframe(r'", tf, "')_set_index([", __py_repr(row_names(obj)),"])")
    } else {
        paste0("read_dataframe(r'", tf, "')")
    }
}
__py_repr_matrix <- function(obj) {
    if (!require("mat")) {
        install_packages('mat', repos='http://cran_stat_ucla_edu/')
        }
    library(mat)
    tf = tempfile('mat')
    write_mat(as_data_frame(obj), tf)
    if (__has_row_names(obj)) {
       paste0("read_dataframe(r'", tf, "')_set_index([", __py_repr(row_names(obj)),"])_as_matrix()")
    } else {
       paste0("read_dataframe(r'", tf, "')_as_matrix()")
    }
}
__py_repr_n <- function(obj) {
    paste("[",
        paste(sapply(obj, __py_repr), collapse=','),
        "]")
}
__py_repr <- function(obj) {
    if (is_matrix(obj)) {
        __py_repr_matrix(obj)
    } else if (is_data_frame(obj)) {
        __py_repr_dataframe(obj)
    } else if (is_null(obj)) {
        'None'
    } else if (is_integer(obj)) {
        if (length(obj) == 1)
            __py_repr_integer_1(obj)
        else
            paste("[", paste(obj, collapse=','), "]")
    } else if (is_double(obj)){
        if (length(obj) == 1)
            __py_repr_double_1(obj)
        else
            paste("[", paste(obj, collapse=','), "]")
    } else if (is_character(obj)) {
        if (length(obj) == 1)
            __py_repr_character_1(obj)
        else {
            paste("[", paste(sapply(obj, __py_repr_character_1), collapse=','), "]")
        }
    } else if (is_logical(obj)) {
        if (length(obj) == 1)
            __py_repr_logical_1(obj)
        else
            __py_repr_n(obj)
    } else if (is_list(obj)) {
        # if the list has no name
        if (is_null(names(obj)))
            __py_repr_n(obj)
        else {
            paste("{",
                  paste(sapply(names(obj), function (x)
                      paste(shQuote(gsub("\\_", "_", as_character(x))), ":", __py_repr(obj[[x]]))),
                      collapse=','),
                  "}")
        }
    } else {
        "'Untransferrable variable'"
    }
}
__read_mat <- function(filename, index=NULL) {
    if (! suppressMessages(suppressWarnings(require("mat", quietly = TRUE)))) {
        try(install_packages('mat', repos='http://cran_stat_ucla_edu/'), silent=TRUE)
        if (!suppressMessages(suppressWarnings(require("mat"))))
            stop('Failed to install mat library')
    }
    suppressPackageStartupMessages(library(mat, quietly = TRUE))
    data = as_data_frame(read_mat(filename))
    if (!is_null(index))
        rownames(data) <- index
    return(data)
}
'''


class sos_Matlab:
    def __init__(self, sos_kernel):
        self.sos_kernel = sos_kernel
        self.kernel_name = 'matlab_kernel'
        self.background_color = '#FA8072'
        self.init_statements = Matlab_init_statements

    def get_vars(self, names):
        for name in names:
            if name.startswith('_'):
                self.sos_kernel.warn('Variable {} is passed from SoS to kernel {} as {}'.format(name, self.kernel_name, '.' + name[1:]))
                newname = '.' + name[1:]
            else:
                newname = name
            matlab_repr = _Matlab_repr(env.sos_dict[name])
            self.sos_kernel.run_cell('{} = {}'.format(newname, matlab_repr), True, False, on_error='Failed to get variable {} to R'.format(name))

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

        try:
            if 'save_mat' in expr:
                # imported to be used by eval
                from scipy.io import read_mat
                # suppress flakes warning
                read_mat
            # evaluate as raw string to correctly handle \\ etc
            return eval(expr)
        except Exception as e:
            self.sos_kernel.warn('Failed to evaluate {!r}: {}'.format(expr, e))
            return None

    def sessioninfo(self):
        return self.sos_kernel.get_response(r'cat(paste(capture.output(sessionInfo()), collapse="\n"))', ('stream',), name=('stdout',))[0]
