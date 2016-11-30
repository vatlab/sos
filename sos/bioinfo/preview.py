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


from sos.utils import short_repr

def preview_bam(filename):
    import pysam
    res = ''
    with pysam.AlignmentFile(filename, 'rb') as bam:
        headers = bam.header
        for record_type in ('RG', 'PG', 'SQ'):
            if record_type not in headers:
                continue
            else:
                records = headers[record_type]
            res += record_type + ':\n'
            for i, record in enumerate(records):
                if type(record) == str:
                    res += '  ' + short_repr(record) + '\n'
                elif type(record) == dict:
                    res += '  '
                    for idx, (k, v) in enumerate(record.items()):
                        if idx < 4:
                            res += '{}: {}    '.format(k, short_repr(v))
                        elif idx == 4:
                            res += '...'
                            break
                if i > 4:
                    res += '\n  ...\n'
                    break
                else:
                    res += '\n'
    return res

