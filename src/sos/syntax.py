#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import keyword
import re

from typing import Callable, List

SOS_TARGETS_OPTIONS = [
    'group_by', 'paired_with', 'pattern', 'group_with', 'for_each',
    'remove_empty_groups'
]
SOS_INPUT_OPTIONS = ['concurrent']
SOS_OUTPUT_OPTIONS = ['group_with']
SOS_DEPENDS_OPTIONS: List = []
SOS_RUNTIME_OPTIONS = [
    'workdir', 'concurrent', 'active', 'walltime', 'nodes', 'cores', 'mem',
    'shared', 'env', 'prepend_path', 'queue', 'to_host', 'from_host',
    'map_vars', 'name', 'trunk_size', 'trunk_workers', 'tags'
]
SOS_ACTION_OPTIONS = [
    'workdir', 'container', 'engine', 'docker_image', 'docker_file', 'active',
    'input', 'output', 'allow_error', 'tracked', 'stdout', 'stderr',
    'default_env', 'env'
]

SOS_DIRECTIVES = ['input', 'output', 'depends', 'task', 'parameter']
SOS_SECTION_OPTIONS = ['provides', 'shared', 'workdir']

SOS_KEYWORDS = SOS_INPUT_OPTIONS + SOS_OUTPUT_OPTIONS + SOS_DEPENDS_OPTIONS + SOS_RUNTIME_OPTIONS \
    + SOS_ACTION_OPTIONS + SOS_DIRECTIVES + SOS_SECTION_OPTIONS

SOS_USAGES = {
    'input':
        '''
input: filename, filename, ... [group_by=GROUP] [filetype=FILETYPE]
          [paired_with=PAIRS] [for_each=VARS] [pattern=PATTEN]

Specify input targets of a SoS step.

See online documentation for details of variables.
''',
    'output':
        '''
output: target, target, ...

Specify output targets of a SoS step.

See online documentation for details of variables.
''',
    'depends':
        '''
depends: target, target, ...

Specify dependent targets of a SoS step.

See online documentation for details of variables.
'''
}

#
# This code is copied from Bazzar to compile regular expressions only when they
# are needed, to avoid compiling all regular expressions up front even when they
# are not used
#


class LazyRegex(object):
    """A proxy around a real regex, which won't be compiled until accessed."""

    # These are the parameters on a real _sre.SRE_Pattern object, which we
    # will map to local members so that we don't have the proxy overhead.
    _regex_attributes_to_copy = [
        '__copy__', '__deepcopy__', 'findall', 'finditer', 'match', 'scanner',
        'search', 'split', 'sub', 'subn'
    ]

    # We use slots to keep the overhead low. But we need a slot entry for
    # all of the attributes we will copy
    __slots__ = [
        '_real_regex',
        '_regex_args',
        '_regex_kwargs',
    ] + _regex_attributes_to_copy

    def __init__(self, *args, **kwargs) -> None:
        """Create a new proxy object, passing in the args to pass to re.compile
        :param args: The `*args` to pass to re.compile
        :param kwargs: The `**kwargs` to pass to re.compile
        """
        self._real_regex = None
        self._regex_args = args
        self._regex_kwargs = kwargs

    def _compile_and_collapse(self) -> None:
        """Actually compile the requested regex"""
        self._real_regex = self._real_re_compile(*self._regex_args,
                                                 **self._regex_kwargs)
        for attr in self._regex_attributes_to_copy:
            setattr(self, attr, getattr(self._real_regex, attr))

    def _real_re_compile(self, *args, **kwargs):
        """Thunk over to the original re.compile"""
        try:
            return re.compile(*args, **kwargs)
        except re.error as e:
            # raise ValueError instead of re.error as this gives a
            # cleaner message to the user.
            raise ValueError('"' + args[0] + '" ' + str(e))

    def __getstate__(self):
        """Return the state to use when pickling."""
        return {
            "args": self._regex_args,
            "kwargs": self._regex_kwargs,
        }

    def __setstate__(self, sdict):
        """Restore from a pickled state."""
        self._real_regex = None
        setattr(self, "_regex_args", sdict["args"])
        setattr(self, "_regex_kwargs", sdict["kwargs"])

    def __getattr__(self, attr: str) -> Callable:
        """Return a member from the proxied regex object.
        If the regex hasn't been compiled yet, compile it
        """
        if self._real_regex is None:
            self._compile_and_collapse()
        # Once we have compiled, the only time we should come here
        # is actually if the attribute is missing.
        return getattr(self._real_regex, attr)


# Regular expressions for parsing section headers and options
_SECTION_HEADER_TMPL = r'''
    ^\[\s*                             # [
    (?P<section_name>
    [a-zA-Z*0-9_-]+                    # name,
    \s*(\([^)]*\))?                    # optional alias
    (\s*,\s*                           # ,
    [a-zA-Z*0-9_]+                     # another name
    \s*(\([^)]*\))?                    # alias
    )*                                 # optional second or more
    )                                  # section names
    (:\s*                              # :
    (?P<section_option>.*)             # section options
    )?                                 # optional
    \]\s*$                             # ]
    '''

_GLOBAL_SECTION_HEADER_TMPL = r'''
    ^\[\s*                             # [
    global\s*                          # global
    \]\s*$                             # ]
    '''

_SECTION_NAME_TMPL = r'''
    ^\s*                               # start
    (?P<name>                          # optional name
    [a-zA-Z*]                          # alphabet or '*'
    ([-\w\d_*]*?                       # followed by alpha numeric, '-' or '*'
    [-a-zA-Z\d*])??                    # but last character cannot be _
    )?                                 # name is optional
    (?(name)                           # if there is name
    (_(?P<index>\d+))?                 #   optional _index
    |(?P<default_index>\d+))           # no name, then index
    \s*(\(\s*
    (?P<alias>[^)]+)                   # optional alias
    \s*\))?
    \s*$
    '''

_SUBWORKFLOW_TMPL = r'''
    ^\s*                               # leading space
    (?P<name>                          # name
    ([a-zA-Z]+\.)*                     # optional name
    ([a-zA-Z*]                         # cannot start with _ etc
    ([-\w\d_]*?))?)                    # can - have _ and digit
    (:(?P<steps>                       # index start from :
    [\d\s-]+))?                        # with - and digit
    \s*$                               # end
    '''

_SECTION_OPTION_TMPL = r'''
    ^\s*                               # start
    (?P<name>{})                       # one of the option names
    (\s*=\s*                           # =
    (?P<value>.+)                      # value
    )?                                 # value is optional
    \s*$
    '''.format('|'.join(SOS_SECTION_OPTIONS))

_FORMAT_LINE_TMPL = r'''
    ^                                  # from first column
    \#fileformat\s*=\s*                # starts with #fileformat=SOS
    (?P<format_name>\S*)               # format name
    \s*$                               # till end of line
    '''

_FORMAT_VERSION_TMPL = r'''
    ^                                  # from first column
    (?P<format_name>[a-zA-Z]+)         # format name
    (?P<format_version>[\d\.]+)        # any number and .
    \s*$                               # till end of line
    '''

_DIRECTIVE_TMPL = r'''
    ^                                  # from start of line
    (?P<directive_name>                #
    (?!({})\s*:)                       # not a python keyword followed by : (can be input)
    ({}                                # name of directive
    |[a-zA-Z][\w\d_]*))             #    or action
    \s*:\s*                            # followed by :
    (?P<directive_value>               # and values that can be
    ([^:.|=&@$^<>].*)?$)               # name followed by and values, which can be
                                       # constant 'a', "b", variable or function call
                                       # a(), or arbitrary expression (['a'...], dictionary, set
                                       # etc) which is difficult to match, so we use negative
                                       # pattern to exclude expressions starting with :, | etc
    '''.format('|'.join(keyword.kwlist), '|'.join(SOS_DIRECTIVES))

_INDENTED_ACTION_TMPL = r'''
    ^                                  # from start of line but allow space
    (?P<action_name>                   #
    (?!\s+({}|{})\s*:)                 # not a python keyword or SoS directive followed by :
    (\s+[a-zA-Z][\w\d_]*))             #    or action
    \s*:\s*                            # followed by :
    (?P<action_value>                  # and values that can be
    ([^:.|=&@$^<>].*)?$)               # name followed by and values, which can be
                                       # constant 'a', "b", variable or function call
                                       # a(), or arbitrary expression (['a'...], dictionary, set
                                       # etc) which is difficult to match, so we use negative
                                       # pattern to exclude expressions starting with :, | etc
    '''.format('|'.join(keyword.kwlist), '|'.join(SOS_DIRECTIVES))

_ASSIGNMENT_TMPL = r'''
    ^                                  # from start of line
    (?P<var_name>[\w_][\d\w_]*)        # variable name
    \s*=(?!=)\s*                            # assignment
    (?P<var_value>.*)             # variable content
    '''

_CONFIG_NAME = r'''
   ^
   [a-zA-Z]                            # first letter must be
   [a-zA-Z0-9_]*
   (\.[a-zA-Z]                         # followed by more names separated with '.'
   [a-zA-Z0-9_]*)*
   $
   '''

_SOS_MAGIC_TMPL = r'''                  # SOS magic
    ^%(dict                             # %dict
    |run                                # %run
    |paste                              # %paste
    |set                                # %set
    |get                                # %get
    |preview                            # %preview
    |with                               # %with
    |use                                # %use
    |shutdown                           # %shutdown
    )(\s+.*)?
    $
    '''

_SOS_CELL_TMPL = r'''                   # %cell
    ^%cell(\s+.*)?
    '''

_SOS_CELL_LINE_TMPL = r'''
    ^
    %cell(\s+                           # %cell
    (?P<cell_type>
    (markdown|code)                     # markdown or code
    )?
    (\s+
    (?P<cell_count>
    [0-9]+
    ))?
    (?P<metainfo>
    .*                                  # arbitrary stuff
    )?)?
    $
    '''

_INDENTED_TMPL = r'''
    ^                                   # start from beginning of string
    (
    \s*\n                               # empty lines are ignored
    )*
    (\s*)\S                             # match a line with a non-space character
    '''

_SOS_WILDCARD_TMPL = r'''
    \{
        \s*
        (?P<name>\w+?)
        (\s*,\s*
        (?P<constraint>
            ([^\{\}]+|\{\d+(,\d+)?\})
        *)
        )?
        \s*
    \}
    '''

_SOS_TAG_TMPL = r'''
    ^[\w_.-]+$                          # - is allowed in between
    '''

_SOS_LOGLINE = r'''
    ^                                   # 2017-11-01 14:26:16,145:
    \d\d\d\d\-\d\d\-\d\d\s+\d\d:\d\d:\d\d,\d\d\d:
'''

SOS_SECTION_HEADER = LazyRegex(_SECTION_HEADER_TMPL, re.VERBOSE)
SOS_GLOBAL_SECTION_HEADER = LazyRegex(_GLOBAL_SECTION_HEADER_TMPL, re.VERBOSE)
SOS_SECTION_NAME = LazyRegex(_SECTION_NAME_TMPL, re.VERBOSE)
SOS_SUBWORKFLOW = LazyRegex(_SUBWORKFLOW_TMPL, re.VERBOSE)
SOS_SECTION_OPTION = LazyRegex(_SECTION_OPTION_TMPL, re.VERBOSE)
SOS_FORMAT_LINE = LazyRegex(_FORMAT_LINE_TMPL, re.VERBOSE)
SOS_FORMAT_VERSION = LazyRegex(_FORMAT_VERSION_TMPL, re.VERBOSE)
SOS_DIRECTIVE = LazyRegex(_DIRECTIVE_TMPL, re.VERBOSE)
SOS_INDENTED_ACTION = LazyRegex(_INDENTED_ACTION_TMPL, re.VERBOSE)
SOS_ASSIGNMENT = LazyRegex(_ASSIGNMENT_TMPL, re.VERBOSE)
CONFIG_NAME = LazyRegex(_CONFIG_NAME, re.VERBOSE)
SOS_MAGIC = LazyRegex(_SOS_MAGIC_TMPL, re.VERBOSE)
SOS_CELL = LazyRegex(_SOS_CELL_TMPL, re.VERBOSE)
SOS_CELL_LINE = LazyRegex(_SOS_CELL_LINE_TMPL, re.VERBOSE)
INDENTED = LazyRegex(_INDENTED_TMPL, re.VERBOSE)
SOS_WILDCARD = LazyRegex(_SOS_WILDCARD_TMPL, re.VERBOSE)
SOS_TAG = LazyRegex(_SOS_TAG_TMPL, re.VERBOSE)
SOS_LOGLINE = LazyRegex(_SOS_LOGLINE, re.VERBOSE)
