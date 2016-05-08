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

SOS_INPUT_OPTIONS = ['group_by', 'skip', 'filetype', 'paired_with', 'for_each', 'pattern', 'dynamic']
SOS_OUTPUT_OPTIONS = ['pattern', 'dynamic']
SOS_DEPENDS_OPTIONS = ['pattern', 'dynamic']
SOS_RUNTIME_OPTIONS = ['workdir', 'concurrent', 'active']
SOS_ACTION_OPTIONS = ['workdir', 'docker_image', 'docker_file', 'active']

SOS_DIRECTIVES = ['input', 'output', 'depends', 'task']
SOS_SECTION_OPTIONS = ['alias', 'skip', 'sigil', 'target']
SOS_PARAMETERS_SECTION_NAME = 'parameters'

SOS_REPORT_PREFIX = '!'

# Regular expressions for parsing section headers and options
_SECTION_HEADER_TMPL = r'''
    ^\[\s*                             # [
    (?P<section_name>[\d\w_,*\s]+)     # digit, alphabet, _ and ,
    (:\s*                              # :
    (?P<section_option>.*)             # section options
    )?                                 # optional
    \]\s*$                             # ]
    '''

_SECTION_NAME_TMPL = '''
    ^\s*                               # start
    (?P<name>                          # optional name
    [a-zA-Z*]                          # alphabet or '*'
    ([\w\d_*]*?                        # followed by alpha numeric or '*'
    [a-zA-Z\d*])??                     # but last character cannot be _
    )?                                 # name is optional
    (?(name)                           # if there is name
    (_(?P<index>\d+))?                 #   optional _index
    |(?P<default_index>\d+))           # no name, then index
    \s*$
    '''

_SUBWORKFLOW_TMPL = '''
    ^\s*                               # leading space
    (?P<name>                          # name
    ([a-zA-Z]+\.)*                     # optional name
    [a-zA-Z*]                          # cannot start with _ etc
    ([\w\d_]*?))                       # can have _ and digit
    (_(?P<steps>                       # index start from _
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
    (?P<format_name>.*)                # format name
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
    (?P<directive_value>.*)            # and values
    '''.format('|'.join(keyword.kwlist), '|'.join(SOS_DIRECTIVES))

_ASSIGNMENT_TMPL = r'''
    ^                                  # from start of line
    (?P<var_name>[\w_][\d\w_]*)        # variable name
    \s*=\s*                            # assignment
    (?P<var_value>.*)                  # variable content
    '''

SOS_SECTION_HEADER = re.compile(_SECTION_HEADER_TMPL, re.VERBOSE)
SOS_SECTION_NAME = re.compile(_SECTION_NAME_TMPL, re.VERBOSE)
SOS_SUBWORKFLOW = re.compile(_SUBWORKFLOW_TMPL, re.VERBOSE)
SOS_SECTION_OPTION = re.compile(_SECTION_OPTION_TMPL, re.VERBOSE)
SOS_FORMAT_LINE = re.compile(_FORMAT_LINE_TMPL, re.VERBOSE)
SOS_FORMAT_VERSION = re.compile(_FORMAT_VERSION_TMPL, re.VERBOSE)
SOS_DIRECTIVE = re.compile(_DIRECTIVE_TMPL, re.VERBOSE)
SOS_ASSIGNMENT = re.compile(_ASSIGNMENT_TMPL, re.VERBOSE)

