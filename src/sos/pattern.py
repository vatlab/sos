#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import collections
import copy
import os
import re
import sys
from itertools import chain
from typing import Any, Dict, List, Optional, Union

from .syntax import SOS_WILDCARD
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


def regex(filepattern: str) -> str:
    f = []
    last = 0
    wildcards = set()
    for match in SOS_WILDCARD.finditer(filepattern):
        f.append(re.escape(filepattern[last:match.start()]))
        wildcard = match.group("name")
        if wildcard in wildcards:
            if match.group("constraint"):
                raise ValueError(
                    "If multiple wildcards of the same name "
                    "appear in a string, eventual constraints have to be defined "
                    "at the first occurence and will be inherited by the others."
                )
            f.append(f"(?P={wildcard})")
        else:
            wildcards.add(wildcard)
            f.append("(?P<{}>{})".format(
                wildcard,
                match.group("constraint")
                if match.group("constraint") else ".+"))
        last = match.end()
    f.append(re.escape(filepattern[last:]))
    f.append("$")  # ensure that the match spans the whole file
    return "".join(f)


def glob_wildcards(pattern: str, files: Optional[List[str]] = None
                  ) -> Dict[str, Union[List[Any], List[str]]]:
    """
    Glob the values of the wildcards by matching the given pattern to the filesystem.
    Returns a named tuple with a list of values for each wildcard.
    """
    pattern = os.path.normpath(pattern)
    if sys.platform == 'win32':
        # we perform path matching with / slash only
        pattern = pattern.replace('\\', '/')
    first_wildcard = re.search("{[^{]", pattern)
    dirname = os.path.dirname(pattern[:first_wildcard.start()]
                             ) if first_wildcard else os.path.dirname(pattern)
    if not dirname:
        dirname = "."

    names = [match.group('name') for match in SOS_WILDCARD.finditer(pattern)]
    res = {x: [] for x in names}
    pattern = re.compile(regex(pattern))

    if files is None:
        files = ((os.path.join(dirpath, f) if dirpath != "." else f)
                 for dirpath, dirnames, filenames in os.walk(dirname)
                 for f in chain(filenames, dirnames))

    for f in files:
        # we perform path matching with only / slash
        match = re.match(pattern, str(f).replace('\\', '/'))
        if match:
            for name, value in match.groupdict().items():
                res[name].append(value)
    return res


def apply_wildcards(pattern: str,
                    wildcards: Dict[str, Union[int, str]],
                    fill_missing: bool = False,
                    fail_dynamic: bool = False,
                    dynamic_fill: None = None,
                    keep_dynamic: bool = False) -> str:

    def format_match(match):
        name = match.group("name")
        try:
            value = wildcards[name]
            if fail_dynamic and value == dynamic_fill:
                raise RuntimeError(name)
            return str(value)  # convert anything into a str
        except KeyError as ex:
            if keep_dynamic:
                return f"{{{name}}}"
            elif fill_missing:
                return dynamic_fill
            else:
                raise RuntimeError(f'Wildcard apply error: {ex} ({wildcards})')

    return SOS_WILDCARD.sub(format_match, pattern)


def extract_pattern(pattern: str, ifiles: List[str]) -> Dict[str, any]:
    '''This function match pattern to a list of input files, extract and return
    pieces of filenames as a list of variables with keys defined by pattern.'''
    res = glob_wildcards(pattern, [])
    for ifile in ifiles:
        matched = glob_wildcards(pattern, [ifile])
        for key in matched.keys():
            if not matched[key]:
                #env.logger.warning('Filename {} does not match pattern {}. None returned.'.format(ifile, pattern))
                res[key].append(None)
            else:
                res[key].extend(matched[key])
    return res


def expand_pattern(pattern: str) -> List[str]:
    '''This function expand patterns against the current namespace
    and return a list of filenames'''
    ofiles = []
    sz = None
    res = glob_wildcards(pattern, [])
    sz = None
    wildcard = [{}]
    for key in res.keys():
        if key not in env.sos_dict:
            raise ValueError(f'Undefined variable {key} in pattern {pattern}')
        if not isinstance(env.sos_dict[key], str) and isinstance(
                env.sos_dict[key], collections.Sequence):
            if sz is None:
                sz = len(env.sos_dict[key])
                wildcard = [copy.deepcopy(wildcard[0]) for x in range(sz)]
            elif sz != len(env.sos_dict[key]):
                raise ValueError(
                    f'Variables in output pattern should have the same length (other={sz}, len({key})={len(env.sos_dict[key])})'
                )
            for idx, value in enumerate(env.sos_dict[key]):
                wildcard[idx][key] = value
        else:
            for v in wildcard:
                v[key] = env.sos_dict[key]
    #
    for card in wildcard:
        ofiles.append(
            apply_wildcards(
                pattern,
                card,
                fill_missing=False,
                fail_dynamic=False,
                dynamic_fill=None,
                keep_dynamic=False))
    return ofiles
