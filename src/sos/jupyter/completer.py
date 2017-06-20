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

import os
import glob

import rlcompleter
from sos.utils import env

def last_valid(line):
    text = line
    for char in (' ', '\t', '"', "'", '=', '('):
        if text.endswith(char):
            text = ''
        elif char in text:
            text = text.rsplit(char, 1)[-1]
    return text

class SoS_MagicsCompleter:
    def __init__(self, kernel):
        self.kernel = kernel

    def get_completions(self, line):
        text = last_valid(line)

        if not text.strip():
            if line.startswith('%get'):
                return text, [x for x in env.sos_dict.keys() if x not in \
                    self.kernel.original_keys and not x.startswith('_')]
            elif any(line.startswith(x) for x in ('%use', '%with', '%restart')):
                return text, ['sos'] + list(self.kernel.supported_languages.keys())
            else:
                return None
        elif text.startswith('%') and line.startswith(text):
            return text, ['%' + x + ' ' for x in self.kernel.ALL_MAGICS if x.startswith(text[1:])]
        elif any(line.startswith(x) for x in ('%use', '%with', '%restart')):
            return text, [x for x in self.kernel.supported_languages.keys() if x.startswith(text)]
        elif line.startswith('%get'):
            return text, [x for x in env.sos_dict.keys() if x.startswith(text) \
                and x not in self.kernel.original_keys and not x.startswith('_')]
        else:
            return None

class SoS_PathCompleter:
    '''PathCompleter.. The problem with ptpython's path completor is that
    it only matched 'text_before_cursor', which would not match cases such
    as %cd ~/, which we will need.'''
    def __init__(self):
        pass

    def get_completions(self, line):
        text = last_valid(line)

        if not text.strip():
            return text, glob.glob('*')
        else:
            matches = glob.glob(os.path.expanduser(text) + '*')
            if len(matches) == 1 and matches[0] == os.path.expanduser(text) \
                and os.path.isdir(os.path.expanduser(text)):
                return text, glob.glob(os.path.expanduser(text) + '/*')
            else:
                return text, matches

class PythonCompleter:
    def __init__(self):
        pass

    def get_completions(self, line):
        text = last_valid(line)
        
        completer = rlcompleter.Completer(env.sos_dict._dict)
        return text, completer.global_matches(text)

class SoS_Completer(object):
    def __init__(self, kernel):
        self.completers = [
            SoS_MagicsCompleter(kernel),
            SoS_PathCompleter(),
            PythonCompleter(),
        ]

    def complete_text(self, code, cursor_pos = None):
        if cursor_pos is None:
            cursor_pos = len(code)

        # get current line before cursor
        doc = code[:cursor_pos].rpartition('\n')[2]

        for c in self.completers:
            matched = c.get_completions(doc)
            if matched is None:
                continue
            elif isinstance(matched, tuple):
                if matched[1]:
                    return matched
            else:
                # iterator ...
                matched = list(matched)
                if matched:
                    return code[matched[0].start_position:], [x.text for x in matched]
        # No match
        return '', []

