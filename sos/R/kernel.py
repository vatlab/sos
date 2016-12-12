#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
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

from collections.abc import Sequence
import tempfile
from sos.utils import short_repr
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
            return 'data.matrix(read_feather("{}"))'.format(feather_tmp_)
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
                feather.write_dataframe(data, feather_tmp_)
            except:
                # if data cannot be written, we try to manipulate data
                # frame to have consistent types and try again
                for c in data.columns:
                    if not homogeneous_type(data[c]):
                        data[c] = [str(x) for x in data[c]]
                feather.write_dataframe(data, feather_tmp_)
            return 'read_feather("{}")'.format(feather_tmp_)
        else:
            return repr('Unsupported datatype {}'.format(short_repr(obj)))

def R_repr_of_py_obj(name, obj):
    new_name = '.' + name[1:] if name.startswith('_') else name
    r_repr = _R_repr(obj)
    #
    if 'read_feather' in r_repr:
        return new_name, '''if (!require("feather")) {{
            install.packages('feather', repos='http://cran.stat.ucla.edu/')
            }}
            library(feather)
            {} <- {}'''.format('.' + name[1:] if name.startswith('_') else name, r_repr)
    else:
        return new_name, '{} <- {}'.format(new_name, r_repr)

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
    options(useFancyQuotes=FALSE)
    dQuote(obj)
}
..py.repr.dataframe <- function(obj) {
    if (!require("feather")) {
        install.packages('feather', repos='http://cran.stat.ucla.edu/')
        }
    library(feather)
    tf = tempfile('feather')
    write_feather(obj, tf)
    paste0("read_dataframe('", tf, "')")
}
..py.repr.matrix <- function(obj) {
    if (!require("feather")) {
        install.packages('feather', repos='http://cran.stat.ucla.edu/')
        }
    library(feather)
    tf = tempfile('feather')
    write_feather(as.data.frame(obj), tf)
    paste0("read_dataframe('", tf, "').as_matrix()")
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
            ..py.repr.n(obj)
    } else if (is.double(obj)){
        if (length(obj) == 1)
            ..py.repr.double.1(obj)
        else
            ..py.repr.n(obj)
    } else if (is.character(obj)) {
        if (length(obj) == 1)
            ..py.repr.character.1(obj)
        else
            ..py.repr.n(obj)
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
            options(useFancyQuotes=FALSE)
            paste("{",
                  paste(sapply(names(obj), function (x)
                      paste(dQuote(gsub("\\.", "_", as.character(x))), ":", ..py.repr(obj[[x]]))),
                      collapse=','),
                  "}")
        }
    } else {
        "'Untransferrable variable'"
    }
}
'''


def py_repr_of_R_obj(items):
    return [x.replace('.', '_') for x in items], \
        '..py.repr(list({}))'.format(','.join('{0}={0}'.format(x) for x in items))

def py_from_R_repr(expr):
    '''
    Convert expression returned from R to python
    '''
    try:
        if 'read_dataframe' in expr:
            # imported to be used by eval
            from feather import read_dataframe
            # suppress flakes warning
            read_dataframe
        # the result is something like
        # [1] "{'a': 1}"
        return eval(eval(expr.split(' ', 1)[-1]))
    except Exception as e:
        raise UsageError('Failed to convert {} to Python object: {}'.format(expr, e))

class sos_R:
    def __init__(self):
        self.kernel_name = 'ir'
        self.init_statements = R_init_statements
        self.py_to_lan = R_repr_of_py_obj
        self.lan_to_py = py_repr_of_R_obj
        self.py_to_dict = py_from_R_repr
