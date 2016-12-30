#!/usr/bin/env python3
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
#
# Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
##
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

import re
import keyword

SOS_INPUT_OPTIONS = ['group_by', 'filetype', 'paired_with', 'for_each', 'pattern']
SOS_OUTPUT_OPTIONS = []
SOS_DEPENDS_OPTIONS = []
SOS_RUNTIME_OPTIONS = ['workdir', 'concurrent', 'active', 'walltime', 'env', 'prepend_path']
SOS_ACTION_OPTIONS = ['workdir', 'docker_image', 'docker_file', 'active', 'input', 'output']

SOS_DIRECTIVES = ['input', 'output', 'depends', 'task', 'parameter']
SOS_SECTION_OPTIONS = ['skip', 'sigil', 'provides', 'shared', 'workdir']

SOS_KEYWORDS = SOS_INPUT_OPTIONS + SOS_OUTPUT_OPTIONS + SOS_DEPENDS_OPTIONS + SOS_RUNTIME_OPTIONS \
    + SOS_ACTION_OPTIONS + SOS_DIRECTIVES + SOS_SECTION_OPTIONS

SOS_USAGES = {
    'input': '''
input: filename, filename, ... [group_by=GROUP] [filetype=FILETYPE] 
          [paired_with=PAIRS] [for_each=VARS] [pattern=PATTEN]

Specify input targets of a SoS step.

See online documentation for details of variables.
''',
    'output': '''
output: target, target, ...

Specify output targets of a SoS step.

See online documentation for details of variables.
''',
    'depends': '''
depends: target, target, ...

Specify dependent targets of a SoS step.

See online documentation for details of variables.
'''}

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
                 '__copy__', '__deepcopy__', 'findall', 'finditer', 'match',
                 'scanner', 'search', 'split', 'sub', 'subn'
                 ]

    # We use slots to keep the overhead low. But we need a slot entry for
    # all of the attributes we will copy
    __slots__ = ['_real_regex', '_regex_args', '_regex_kwargs',
                ] + _regex_attributes_to_copy

    def __init__(self, *args, **kwargs):
        """Create a new proxy object, passing in the args to pass to re.compile
        :param args: The `*args` to pass to re.compile
        :param kwargs: The `**kwargs` to pass to re.compile
        """
        self._real_regex = None
        self._regex_args = args
        self._regex_kwargs = kwargs

    def _compile_and_collapse(self):
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
            raise ValueError('"' + args[0] + '" ' +str(e))

    def __getstate__(self):
        """Return the state to use when pickling."""
        return {
            "args": self._regex_args,
            "kwargs": self._regex_kwargs,
            }

    def __setstate__(self, dict):
        """Restore from a pickled state."""
        self._real_regex = None
        setattr(self, "_regex_args", dict["args"])
        setattr(self, "_regex_kwargs", dict["kwargs"])

    def __getattr__(self, attr):
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

_SECTION_NAME_TMPL = '''
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

_SUBWORKFLOW_TMPL = '''
    ^\s*                               # leading space
    (?P<name>                          # name
    ([a-zA-Z]+\.)*                     # optional name
    [a-zA-Z*]                          # cannot start with _ etc
    ([-\w\d_]*?))                      # can - have _ and digit
    (:(?P<steps>                       # index start from :
    [\d\s-]+))?                        # with - and digit
    \s*$                               # end
    '''

_SECTION_OPTION_TMPL = '''
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
    |[a-zA-Z][\w\d_]*))                #    or action
    \s*:\s*                            # followed by :
    (?P<directive_value>               # and values that can be
    ([^:.|=&@$^<>].*)?$)               # name followed by and values, which can be
                                       # constant 'a', "b", variable or function call
                                       # a(), or arbitrary expression (['a'...], dictionary, set
                                       # etc) which is difficult to match, so we use negative
                                       # pattern to exclude expressions starting with :, | etc
    '''.format('|'.join(keyword.kwlist), '|'.join(SOS_DIRECTIVES))

_ASSIGNMENT_TMPL = r'''
    ^                                  # from start of line
    (?P<var_name>[\w_][\d\w_]*)        # variable name
    \s*=\s*                            # assignment
    (?P<var_value>.*)                  # variable content
    '''
_CONFIG_NAME = r'''
   ^
   [a-zA-Z]                            # first letter must be
   [a-zA-Z0-9_]*
   (\.[a-zA-Z]                         # followed by more names separated with '.'
   [a-zA-Z0-9_]*)*
   $
   '''


_SOS_STRU_TMPL = r'''                   # structural directive
    ^%(if                               # %if
    |elif                               # %elif
    |else                               # %else
    |endif                              # %endif
    |cell                               # %cell
    |set_options                        # %set_options
    |include                            # %include
    |from                               # %from
    )(\s+.*)?
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
    |restart                            # %restart
    )(\s+.*)?
    $
    '''

_SOS_IF_TMPL = r'''                     # %if
    ^%if\s+
    (?P<condition>.+)
    $
    '''

_SOS_ELIF_TMPL = r'''                   # %if
    ^%elif\s+
    (?P<condition>.+)
    $
    '''

_SOS_ELSE_TMPL = r'''                   # %else
    ^%else\s*
    $
    '''

_SOS_ENDIF_TMPL = r'''                  # %endif
    ^%endif\s*
    $
    '''

_SOS_CELL_TMPL = r'''                   # %cell
    ^%cell(\s+.*)?
    '''

_SOS_OPTIONS_TMPL = r'''                # %set_options
    ^%set_options(\s+
    (?P<options>.*)
    )?
    '''

_SOS_INCLUDE_TMPL = r'''
    ^%                                  # from start of line
    include\s+                          # keyword
    (?P<sos_files>
    [a-zA-Z][a-zA-Z0-9_.]*              # filename
    (\s+as\s+
    [a-zA-Z][a-zA-Z0-9_]*)?             # optional as ...
    (\s*,\s*
    [a-zA-Z][a-zA-Z0-9_.]*
    (\s+as\s+
    [a-zA-Z][a-zA-Z0-9_]*)?             # optional as ...
    )*)                                 # or more filename
    \s*$
    '''

_SOS_FROM_INCLUDE_TMPL = r'''
    ^%                                  # from start of line
    from
    \s+
    (?P<sos_file>
    [a-zA-Z][a-zA-Z0-9_.]*)\s+          # filename
    include\s+                          # include
    (?P<names>
    (\*|
    [a-zA-Z][a-zA-Z0-9_]*               # workflow as aslias
    (\s+as\s+
    [a-zA-Z][a-zA-Z0-9_]*)?
    (\s*,\s*                            # parate by ,
    [a-zA-Z][a-zA-Z0-9_]*               # additional workflow as aslias
    (\s+as\s+
    [a-zA-Z][a-zA-Z0-9_]*)?
    )*)) 
    \s*$
    '''

_SOS_AS_TMPL = r'''
    (?P<name>[a-zA-Z][a-zA-Z0-9_]*)     # name of item
    (\s+as\s+                           #  as 
    (?P<alias>
    [a-zA-Z][a-zA-Z0-9_]*))?            # optional alias
    '''

_SOS_CELL_LINE_TMPL = r'''
    ^
    %cell(\s+                           # %cell
    (?P<cell_type>
    (markdown|code)                     # markdown or code
    )
    (\s+\d+\s+)?)?                      # arbitrary stuff
    $
    '''

_INDENTED_TMPL = r'''
    ^                                   # start from beginning of string
    (
    \s*\n                               # empty lines are ignored
    )*
    (\s*)\S                             # match a line with a non-space character
    '''

# Format specifier that can be used at the end of the string to convert
# result (!s|r|q) and control output format (:.2f etc). The pattern
# is constructed according to Python format mini language.
_FORMAT_SPECIFIER_TMPL = r'''
    ^                                   # start of expression
    (?P<expr>.*?)                       # any expression
    (?P<conversion>!\s*                 # conversion starting with !
    [srqabden,]+                        # conversion, q, a, b, n, and , are added by SoS
    )?
    (?P<format_spec>:\s*                # format_spec starting with :
    (?P<fill>.?[<>=^])?                 # optional fill|align
    (?P<sign>[-+ ])?                    # optional sign
    \#?                                 #
    0?                                  #
    (?P<width>\d+)?                     # optional width
    ,?                                  # optional ,
    (?P<precision>\.\d+)?               # optional precision
    (?P<type>[bcdeEfFgGnosxX%])?        # optional type
    )?                                  # optional format_spec
    \s*$                                # end of tring
    '''


# we handle simple cases in an easier way to avoid linear search each time.
# simple case means ${ } as sigil, and there is nothing but variable name
# within it.
#
_SIMPLE_SUB_TMPL = r'''
    (?<!                                # if not preceded by
    \\                                  # a back slash
    )
    \$\{                                # left sigil
    (                                   # capture variable name
    [_a-zA-Z]\w*                        # alpha numeric with no leading numeric
    )
    \}                                  # right sigil
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


SOS_SECTION_HEADER = LazyRegex(_SECTION_HEADER_TMPL, re.VERBOSE)
SOS_SECTION_NAME = LazyRegex(_SECTION_NAME_TMPL, re.VERBOSE)
SOS_SUBWORKFLOW = LazyRegex(_SUBWORKFLOW_TMPL, re.VERBOSE)
SOS_SECTION_OPTION = LazyRegex(_SECTION_OPTION_TMPL, re.VERBOSE)
SOS_FORMAT_LINE = LazyRegex(_FORMAT_LINE_TMPL, re.VERBOSE)
SOS_FORMAT_VERSION = LazyRegex(_FORMAT_VERSION_TMPL, re.VERBOSE)
SOS_DIRECTIVE = LazyRegex(_DIRECTIVE_TMPL, re.VERBOSE)
SOS_ASSIGNMENT = LazyRegex(_ASSIGNMENT_TMPL, re.VERBOSE)
CONFIG_NAME = LazyRegex(_CONFIG_NAME, re.VERBOSE)
SOS_AS = LazyRegex(_SOS_AS_TMPL, re.VERBOSE)
SOS_STRU = LazyRegex(_SOS_STRU_TMPL, re.VERBOSE)
SOS_MAGIC = LazyRegex(_SOS_MAGIC_TMPL, re.VERBOSE)
SOS_IF = LazyRegex(_SOS_IF_TMPL, re.VERBOSE)
SOS_ELIF = LazyRegex(_SOS_ELIF_TMPL, re.VERBOSE)
SOS_ELSE = LazyRegex(_SOS_ELSE_TMPL, re.VERBOSE)
SOS_ENDIF = LazyRegex(_SOS_ENDIF_TMPL, re.VERBOSE)
SOS_OPTIONS = LazyRegex(_SOS_OPTIONS_TMPL, re.VERBOSE)
SOS_CELL = LazyRegex(_SOS_CELL_TMPL, re.VERBOSE)
SOS_INCLUDE = LazyRegex(_SOS_INCLUDE_TMPL, re.VERBOSE)
SOS_FROM_INCLUDE = LazyRegex(_SOS_FROM_INCLUDE_TMPL, re.VERBOSE)
SOS_CELL_LINE = LazyRegex(_SOS_CELL_LINE_TMPL, re.VERBOSE)
INDENTED = LazyRegex(_INDENTED_TMPL, re.VERBOSE)
FORMAT_SPECIFIER = LazyRegex(_FORMAT_SPECIFIER_TMPL, re.VERBOSE | re.DOTALL)
SIMPLE_SUB = LazyRegex(_SIMPLE_SUB_TMPL, re.VERBOSE | re.DOTALL)
SOS_WILDCARD = LazyRegex(_SOS_WILDCARD_TMPL, re.VERBOSE)
