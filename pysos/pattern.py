#!/usr/bin/env python3
#
# This file is part of Script of Scripts (SoS), a workflow system
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
import os
import re
import copy
import collections
from itertools import chain
from .utils import env

__all__ = ['expand_pattern']

#
# The following is adapted from snakemake so that the syntax of the pattern parameter
# would be the same as snakemake.
#
# https://bitbucket.org/snakemake/snakemake/src/22eff35627401d3bb243e068c5dce97107e7090b/snakemake/io.py?at=master&fileviewer=file-view-default
#
# I will re-implement it if there is any license issue with the code
#

_wildcard_regex = re.compile(
    "\{\s*(?P<name>\w+?)(\s*,\s*(?P<constraint>([^\{\}]+|\{\d+(,\d+)?\})*))?\s*\}")


def regex(filepattern):
    f = []
    last = 0
    wildcards = set()
    for match in _wildcard_regex.finditer(filepattern):
        f.append(re.escape(filepattern[last:match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "If multiple wildcards of the same name "
                    "appear in a string, eventual constraints have to be defined "
                    "at the first occurence and will be inherited by the others.")
            f.append("(?P={})".format(wildcard))
        else:
            wildcards.add(wildcard)
            f.append("(?P<{}>{})".format(wildcard, match.group("constraint") if
                                         match.group("constraint") else ".+"))
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)


def glob_wildcards(pattern, files=None):
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.
    """
    pattern = os.path.normpath(pattern)
    first_wildcard = re.search("{[^{]", pattern)
    dirname = os.path.dirname(pattern[:first_wildcard.start(
    )]) if first_wildcard else os.path.dirname(pattern)
    if not dirname:
        dirname = "."

    names = [match.group('name')
             for match in _wildcard_regex.finditer(pattern)]
    res = {x: [] for x in names}
    pattern = re.compile(regex(pattern))

    if files is None:
        files = ((os.path.join(dirpath, f) if dirpath != "." else f)
                 for dirpath, dirnames, filenames in os.walk(dirname)
                 for f in chain(filenames, dirnames))

    for f in files:
        match = re.match(pattern, f)
        if match:
            for name, value in match.groupdict().items():
                res[name].append(value)
    return res

def apply_wildcards(pattern,
                    wildcards,
                    fill_missing=False,
                    fail_dynamic=False,
                    dynamic_fill=None,
                    keep_dynamic=False):
    def format_match(match):
        name = match.group("name")
        try:
            value = wildcards[name]
            if fail_dynamic and value == dynamic_fill:
                raise RuntimeError(name)
            return str(value)  # convert anything into a str
        except KeyError as ex:
            if keep_dynamic:
                return "{{{}}}".format(name)
            elif fill_missing:
                return dynamic_fill
            else:
                raise RuntimeError('Wildcard apply error: {} ({})'.format(ex, wildcards))

    return re.sub(_wildcard_regex, format_match, pattern)

def extract_pattern(pattern, ifiles):
    '''This function match pattern to a list of input files, extract and return
    pieces of filenames as a list of variables with keys defined by pattern.'''
    res = glob_wildcards(pattern, [])
    for ifile in ifiles:
        matched = glob_wildcards(pattern, [ifile])
        for key in matched.keys():
            if not matched[key]:
                env.logger.warning('Filename {} does not match pattern {}. None returned.'.format(ifile, pattern))
                res[key].append(None)
            else:
                res[key].extend(matched[key])
    return res

def expand_pattern(pattern):
    '''This function expand patterns against the current namespace
    and return a list of filenames'''
    ofiles = []
    sz = None
    res = glob_wildcards(pattern, [])
    sz = None
    wildcard = [{}]
    for key in res.keys():
        if key not in env.sos_dict:
            raise ValueError('Undefined variable {} in pattern {}'.format(key, pattern))
        if not isinstance(env.sos_dict[key], str) and isinstance(env.sos_dict[key], collections.Sequence):
            if sz is None:
                sz = len(env.sos_dict[key])
                wildcard = [copy.deepcopy(wildcard[0]) for x in range(sz)]
            elif sz != len(env.sos_dict[key]):
                raise ValueError('Variables in output pattern should have the same length (other={}, len({})={})'
                    .format(sz, key, len(env.sos_dict[key])))
            for idx, value in enumerate(env.sos_dict[key]):
                wildcard[idx][key] = value
        else:
            for v in wildcard:
                v[key] = env.sos_dict[key]
    #
    for card in wildcard:
        ofiles.append(apply_wildcards(pattern, card, fill_missing=False,
           fail_dynamic=False, dynamic_fill=None, keep_dynamic=False))
    return ofiles
