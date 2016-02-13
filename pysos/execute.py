#!/usr/bin/env python
#
# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS
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

## import os
## import sys
## import subprocess
## from multiprocessing import Process
## import threading
## import pipes
## import pprint
## #
## import glob
## import hashlib
## import shlex
## import argparse
## import logging
## import shutil
## import tempfile
## import tarfile
## import gzip
## import bz2
## import zipfile
## import time
## import re
## import csv
## import platform
## import logging
## import random
## import traceback
## from multiprocessing import Process
## from collections import namedtuple, MutableMapping
## from itertools import tee, izip, combinations
## 
## from .utils import *
## 
## from .project import PipelineDescription, Project
## 
## if sys.version_info.major == 2:
##     # for parallel execution of steps
##     import Queue
## else:
##     import queue as Queue
##     from ucsctools_py3 import showTrack
## 
## try:
##     import pysam
##     hasPySam = True
## except (ImportError, ValueError) as e:
##     hasPySam = False
## 
## 
## 
## class NamedList:
##     '''This class implements a named list to assist users in inputting a large
##     number of items to a pipeline. Because these lists are often stored in a text 
##     or excel file and are associated with a name or some other meta information,
##     users can input all them in the format of three ':' seprated parts
## 
##     1. name (optional): name of the list.
##     2. values: values can be either a space or comma separated list such as A1 A2 ... 
##        or A1,A2,A3,... or a field in a file in the format of fieldname@filename?query.
##        Fieldname should be one or more fields from the file joined with a non-alphabetical
##        character (e.g. '+'). Filename should be a file in .csv, text (tab delimited) or 
##        EXCEL file. Query is a query that can limit the items from the retrieved list.
##        A query can consist other fields in the input data file with the entire file treated
##        as a retional database.
##     3. meta (optional): meta information for the list, which can be extra comment, weight (relative
##        to other list), and others, must be in the format of comma separated name=value
## 
##     Examples:
##         age@phenotype.xls
##             A column called phenotype from "phenotype.xls". Age will be treated as 
##             name of the list.
## 
##         deceased:dead@phenotype.xls
##             A list with name "deceased" retrieved from column "dead" of "phenotype.xls".
## 
##         id1,id2,id3,id4
##             A list of IDs.
## 
##         IDs@phenotype.xls?recurrance=="1"
##            get a list of IDs (column IDs) from "phenotype.xls" where the recurrance column has
##            value "1"
## 
##         age@phenotype.xls:weight==1000
##             A named list with meta information 1000
## 
##     The named list has attribute
##         name
##             Name of the list, default to "" unless default_name is given.
## 
##         items
##             A list of items.
## 
##         meta
##             Optional meta information as a dictionary. Default to {}.
## 
##     '''
##     def __init__(self, value_string, default_name=""):
##         self.name = default_name
##         self.items = []
##         self.meta = {}
##         # if it is None etc
##         if value_string is None:
##             return
##         if not isinstance(value_string, str):
##             if len(value_string) == 1:
##                 value_string = value_string[0]
##             else:
##                 # if multiple items are passed, treat directly
##                 # as list of strings.
##                 self.items = value_string
##                 return
##         #
##         self._parse(value_string, default_name)
## 
##     def _parse(self, value_string, default_name):
##         if not value_string:
##             return
##         # single space separated string
##         if re.match('([\w\d-]+\s+)+[\w\d-]+', value_string):
##             self.items = value_string.split()
##             return
##         #
##         # comma separated named list
##         matched = re.match('^([\w\d-]+:)*((([\w\d-]+,)+[\w\d-]+)|([^:]*)@([^:?]*)(\?([^:]*))*)(:([^:]+))*$', value_string)
##         if matched is None:
##             raise ValueError(('"{}" is not a valid named list / query string, which should be name (optional) + comma separated list or '
##                 'colname@filename with optional query string (?), with optional meta. Three parts should be separated by :.')
##                 .format(value_string))
##         name = matched.group(1)
##         comma_list = matched.group(3)
##         colname = matched.group(5)
##         filename = matched.group(6)
##         query = matched.group(8)
##         meta = matched.group(10)
##         if name is not None:
##             self.name = name.rstrip(':')
##         if meta is not None:
##             matched = re.match('[\d\w_]+\s*=\s*[^=,]*(,[\d\w_]+\s*=\s*[^=,]*)*$', meta)
##             if not matched:
##                 raise ValueError('Meta information is not in the format of key=value: {}'.format(meta))
##             self.meta = {x.split('=')[0].strip(): x.split('=')[1].strip() for x in meta.split(',')}
##         #
##         if comma_list is not None:
##             self.items = comma_list.split(',')
##         else:
##             import pandas as pd
##             filename = os.path.expanduser(filename)
##             if not os.path.isfile(filename):
##                 raise ValueError('File does not exist: {}'.format(filename))
##             # pandas can be slow to import
##             if filename.endswith('.csv'):
##                 data = pd.read_csv(filename)
##             else:
##                 data = pd.read_excel(filename)
##             # convert everything to str, however, there are cases where 1 is converted to 1.0 when missing value
##             # is present (na is considered float, thus the behavior)
##             data = data.applymap(str)
##             env.logger.trace("{} records are loaded from {}".format(data.shape[0], filename))
##             #env.logger.trace(data)
##             #
##             # if query?
##             if query is not None:
##                 if re.match('.*[\d\w_]+\s*=\s*[\d\w_]+.*', query):
##                     raise ValueError('Syntax "a=b" is not allowed. Please use "a==b" instead: {}'.format(query))
##                 try:
##                     pre_filter = data.shape[0]
##                     data = data.query(query)
##                     if pre_filter != data.shape[0]:
##                         env.logger.info('{} out of {} records are removed by filter {}'.format(pre_filter - data.shape[0], 
##                             pre_filter, query))
##                 except Exception as e:
##                     raise ValueError('Failed to apply query "{}" to data file {}: {}'
##                         .format(query, filename, e))
##             #
##             values = None
##             for col in re.split('([^\w\d_])', colname):
##                 if re.match('[^\w\d_]', col):
##                     if values is None:
##                         raise ValueError('Leading non-ascii word is not allowed. {}'.format(colname))
##                     else:
##                         values = [x+col for x in values]
##                     continue
##                 if col not in data.columns:
##                     raise ValueError('File {} does not have column {}. Available columns are {}'
##                         .format(filename, col, ', '.join(list(data.columns))))
##                 if values is None:
##                     values = list(data[col].fillna(''))
##                 else:
##                     values = [x+y for x,y in zip(values, data[col].fillna(''))]
##             #
##             self.items = values
##             if self.name == default_name:
##                 self.name = colname
## 
## def rvec(vec):
##     return 'c(' + ','.join([repr(x) for x in vec]) + ')'
## 
## class EmitInput:
##     '''An input emitter that emits input files individually, in pairs, or 
##     altogether.'''
##     def __init__(self, group_by='all', select=True, skip=False, pass_unselected=True):
##         '''Select input files of certain types, group them, and send input files
##         to action. Selection criteria can be True (all input file types, default),
##         'False' (select no input file, but an empty list will still be passed to
##         pipeline action), 'fastq' (check content of files), or one or
##         more file extensions (e.g. ['.sam', '.bam']).  Eligible files are by default
##         sent altogether (group_by='all') to action (${INPUT} equals ${INPUT#} where
##         # is the index of step, but can also be sent individually (group_by='single',
##         ${INPUT} equals to a list of a single file) or in pairs 
##         (group_by='paired', e.g. filename_1.txt and filename_2.txt), pairwise for
##         (a0, a1), (a1, a2), (a2, a3) ..., or combinations for all combination of
##         different input files. Unselected files are by default passed directly as
##         output of a step. If skip is set to True, the whole step will be skipped'''
##         self.group_by = group_by
##         if type(select) == str:
##             if select not in ['fastq', 'bam', 'sam'] and not str(select).startswith('.'):
##                 raise ValueError("Value to option select can only be True/False, "
##                     "'fastq', or a file extension with leading '.': '{}' provided."
##                     .format(select))
##             self.select = [select]
##         elif select in [True, False]:
##             self.select = select
##         else:
##             for s in select:
##                 if s not in ['fastq', 'bam', 'sam'] and not str(s).startswith('.'):
##                     raise ValueError("Value to option select can only be True/False, "
##                         "'fastq', or a file extension with leading '.': '{}' provided."
##                         .format(s))
##             self.select = select
##         self.skip = skip
##         self.pass_unselected = pass_unselected
## 
##     def _isFastq(self, filename):
##         try:
##             if not os.path.isfile(filename) and not os.path.isfile(filename + '.file_info'):
##                 raise RuntimeError('File not found')
##             fl = FileInfo(filename).firstline()
##             if fl is None:
##                 env.logger.info('Cannot detect the type of file because the {} has been removed.'
##                     .format(filename))
##                 return filename.lower().split('.')[-1] not in ['bam', 'sam', 'gz', 'zip']
##             if not fl.startswith('@'):
##                 return False
##             if filename.endswith('.gz'):
##                 env.logger.warning('{}: compressed fastq file might not be '
##                     'acceptable to downstream analysis.'.format(filename))
##         except Exception as e:
##             env.logger.debug('Input file {} is not in fastq format: {}'.format(filename, e))
##             return False
##         return True
## 
##     def _is_paired(self, f1, f2, at=None):
##         if len(f1) != len(f2):
##             return False
##         if f1 >= f2:
##             return False
##         diffs = [x != y for x,y in zip(f1, f2)]
##         if sum(diffs) != 1:
##             return False
##         diff_at = diffs.index(True)
##         if sorted([f1[diff_at], f2[diff_at]]) != ['1', '2']:
##             return False
##         if at is not None and diff_at not in at:
##             return False
##         return True
## 
##     def _pairByReadNames(self, selected, unselected):
##         # we should pair files by actual read names
##         read_map = {}
##         for filename in selected:
##             read = FileInfo(filename).firstline().strip()
##             if read[:-1] in read_map:
##                 if read.endswith('1'):
##                     if read_map[read[:-1]][0] is not None:
##                         raise RuntimeError('Fastq file {} has the same first read as {}'
##                             .format(filename, read_map[read[:-1]][0]))
##                     else:
##                         read_map[read[:-1]][0] = filename
##                 elif read.endswith('2'):
##                     if read_map[read[:-1]][1] is not None:
##                         raise RuntimeError('Fastq file {} has the same first read as {}'
##                             .format(filename, read_map[read[:-1]][1]))
##                     else:
##                         read_map[read[:-1]][1] = filename
##                 else:
##                     raise RuntimeError('Fastq file {} is not paired because its read name does '
##                         'not end with 1 or 2'.format(filename))
##             else:
##                 if read.endswith('1'):
##                     read_map[read[:-1]] = [filename, None]
##                 elif read.endswith('2'):
##                     read_map[read[:-1]] = [None, filename]
##                 else:
##                     raise RuntimeError('Fastq file {} is not paired because its read name does '
##                         'not end with 1 or 2'.format(filename))
##         # now, let us go through files
##         pairs = []
##         for read, filenames in read_map.items():
##             if filenames[0] is None:
##                 raise RuntimeError('Fastq file {} is not paired (no matching read is found)'
##                     .format(filenames[0]))
##             elif filenames[1] is None:
##                 raise RuntimeError('Fastq file {} is not paired (no matching read is found)'
##                     .format(filenames[1]))
##             else:
##                 if not self._is_paired(filenames[0], filenames[1]):
##                     env.logger.warning('{} and {} contain paired reads but the filenames '
##                         'do not follow illumina filename convention'
##                         .format(filenames[0], filenames[1]))
##                 pairs.append(filenames)
##         return sorted(pairs), unselected     
## 
##     def _pairByFileName(self, selected, unselected):
##             #
##             # there is a possibility that one name differ at multiple parts
##             # with another name. e.g
##             #
##             #      A1_TAGCTT_L007_R1_001.fastq.gz
##             #
##             # differ with the following two names by a number
##             #
##             #      A1_TAGCTT_L007_R1_002.fastq.gz
##             #      A1_TAGCTT_L007_R2_001.fastq.gz
##             # 
##             # the code below tries to find good pairs first, then use matched
##             # locations to pair others
##             if not selected:
##                 env.logger.warning('No file matching type "{}" is selected for pairing.'
##                     .format(self.select))
##                 return [], unselected
##             all_pairs = [[x,y] for x in selected for y in selected if self._is_paired(x,y)]
##             unpaired = [x for x in selected if not any([x in y for y in all_pairs])]
##             if unpaired:
##                 raise ValueError('Failed to pair input filenames: {} is not paired'
##                     'with any other names.'.format(', '.join(unpaired)))
##             uniquely_paired = [x for x in selected if sum([x in y for y in all_pairs]) == 1]
##             # if some filenames are uniquely paired, we can use them to identify
##             # index locations.
##             if uniquely_paired:
##                 pairs = [x for x in all_pairs if x[0] in uniquely_paired or x[1] in uniquely_paired]
##                 if len(pairs) != all_pairs:
##                     # find the differentiating index of existing pairs
##                     diff_at = set([[i != j for i,j in zip(x[0],x[1])].index(True) for x in pairs])
##                     # use the diff_at locations to screen the rest of the pairs
##                     pairs.extend([x for x in all_pairs if x not in pairs and self._is_paired(x[0], x[1], diff_at)])
##                     #
##                     if len(pairs) * 2 != len(selected):
##                         unpaired = [x for x in selected if not any([x in y for y in pairs])]
##                         raise ValueError('Failed to pair input files because they '
##                             'match multiple filenames: {}'.format(', '.join(unpaired)))
##                 return sorted(pairs), unselected
##             else:
##                 # all filenames match to multiple names, so we try to get all
##                 # differentiating indexes and see if one of them can pair filenames
##                 # perfectly. We start from the end because we assume that _1 _2
##                 # are close to the end of filenames.
##                 #
##                 diff_at = set([[i != j for i,j in zip(x[0],x[1])].index(True) for x in all_pairs])
##                 acceptable_diff_at = []
##                 for d in diff_at:
##                     # try to pair all names at this location.
##                     pairs = [x for x in all_pairs if self._is_paired(x[0], x[1], [d])]
##                     if len(pairs) * 2 != len(selected):
##                         continue
##                     # all filename should appear once and only once
##                     if not all([sum([x in y for y in pairs]) == 1 for x in selected]):
##                         continue
##                     acceptable_diff_at.append(d)
##                 # fortunately, only one perfect pairing is found
##                 if len(acceptable_diff_at) == 1:
##                     pairs = [x for x in all_pairs if self._is_paired(x[0], x[1], acceptable_diff_at)]
##                     return sorted(pairs), unselected
##                 elif len(acceptable_diff_at) > 1:
##                     env.logger.warning('There are {} ways to match all filenames '
##                         'perfectly. The one using a latter differentiating index '
##                         'is used.'.format(len(acceptable_diff_at)))
##                     diff_at = sorted(list(acceptable_diff_at))[-1]
##                     pairs = [x for x in all_pairs if self._is_paired(x[0], x[1], [diff_at])]
##                     return sorted(pairs), unselected
##                 else:
##                     raise ValueError('All filenames match multiple names but no differentiating '
##                         'index can pair filenames perfectly.') 
## 
##     def __call__(self, ifiles, pipeline=None):
##         if self.skip:
##             return [], ifiles
##         selected = []
##         unselected = []
##         for filename in ifiles:
##             match = False
##             if self.select is True:
##                 match = True
##             elif self.select is False:
##                 pass
##             else:   # list of types
##                 for t in self.select:
##                     if t == 'fastq':
##                         if self._isFastq(filename):
##                             match = True
##                             break
##                     if filename.lower().endswith('.' + t.lstrip('.').lower()):
##                         match = True
##                         break
##             #
##             if match:
##                 selected.append(filename)
##             elif self.pass_unselected:
##                 unselected.append(filename)
##         #
##         # for this special case, the step is skipped
##         if self.group_by == 'single':
##             return [[x] for x in selected], unselected
##         elif self.group_by == 'all':
##             return [selected], unselected
##         elif self.group_by == 'paired':
##             if 'fastq' in self.select:
##                 try:
##                     return self._pairByReadNames(selected, unselected)
##                 except Exception as e:
##                     # if failed to pair by read name, pair by filenames
##                     env.logger.warning('Failed to pair fastq files by read names. '
##                         'Trying to pair files by filenames: {}'.format(e))
##             else:
##                 # this should not happen becase we do not need to pair non-fastq files 
##                 # at this point, but I will leave the code here anyway.
##                 env.logger.warning('It is unsafe to pair input files by names instead of '
##                     'their content. Please add option select="fastq" if you need to '
##                     'pair input fastq files')
##             return self._pairByFileName(selected, unselected)
##         elif self.group_by == 'pairwise':
##             f1, f2 = tee(selected)
##             next(f2, None)
##             return [list(x) for x in izip(f1, f2)], unselected
##         elif self.group_by == 'combinations':
##             return [list(x) for x in combinations(selected, 2)], unselected
## 
## 
## class PipelineAction:
##     '''Base class for all pipeline actions. If one or more output files
##     are specified, the pipeline will record the runtime signature of
##     this action in a file ``${cache_dir}/path/to/$OUTPUT[0].exe_info``, which 
##     consists of the MD5 signature of input and output files, command used, 
##     and additional information such as start and end time of execution,
##     standard and error outputs. An action will be skipped if the action
##     is re-run with the same input, output and command.
## 
##     NOTE: The ``__call__`` function or a pipeline action implements the 
##     runtime signature feature and calls function ``execute`` for actual work.
##     User-defined actions should either override the ``__call__`` function
##     (without the runtime signature feature) or function 
##     ``_execute(self, ifiles, pipeline)`` (with runtime signature feature).
##     User can also define function ``_bypass(self, ifiles, pipeline)`` if the
##     step is bypassed due to identical execution signatures.
##     '''
##     def __init__(self, cmd='', output=[]):
##         '''
##         Parameters:
##             cmd (string or list of strings):
##                 one or more commands to be executed. It should capture
##                 the name and all options used by the command.
## 
##             output (string or list of strings):
##                 Output files. If at least one output file is specified,
##                 the runtime signature of this action will be saved to
##                 $output[0].exe_info. The output directory will be created
##                 if it does not exist.
##         '''
##         # multiple command is not allowed.
##         if not cmd:
##             self.cmd = []
##         elif isinstance(cmd, str):
##             self.cmd = [' '.join(cmd.split('\n'))]
##         else:
##             self.cmd = [' '.join(x.split('\n')) for x in cmd]
##         #
##         if not output:
##             self.output = []
##         elif isinstance(output, str):
##             self.output = [output]
##         else:
##             self.output = output
##         #
##         self.runtime = RuntimeFiles(self.output)
## 
##     def _bypass(self, ifiles, pipeline=None):
##         '''Function called by ``__call__`` if the step is bypassed due to identical
##         execution signature. This function can be used, for example, to set pipeline
##         variable even when the step is not executed.'''
##         return True
## 
##     def _execute(self, ifiles, pipeline=None):
##         '''Function called by ``__call__`` for actual action performed on ifiles. A user-defined
##         action should re-define __call__ or this function. This funciton should return ``True`` if
##         the action is completed successfully, ``False`` for pending (signature will be written later,
##         and raise an exception for errors. '''
##         raise RuntimeError('Please define your own execute function in an derived class of PipelineAction.')
## 
##     def _write_info(self, pipeline=None):
##         if not self.output:
##             return
##         with open(self.runtime.proc_info, 'a') as exe_info:
##             exe_info.write('#End: {}\n'.format(time.asctime(time.localtime())))
##             for f in self.output:
##                 if not os.path.isfile(f):
##                     raise RuntimeError('Output file {} does not exist after completion of the job.'.format(f))
##                 # for performance considerations, use partial MD5
##                 exe_info.write('{}\t{}\t{}\n'.format(self._useVars(f, pipeline.VARS), os.path.getsize(f),
##                     calculateMD5(f, partial=True)))
##             # write standard output to exe_info
##             exe_info.write('\n\nSTDOUT\n\n')
##             if os.path.isfile(self.runtime.proc_out):
##                 with open(self.runtime.proc_out) as stdout:
##                     for line in stdout:
##                         exe_info.write(line)
##             # write standard error to exe_info
##             exe_info.write('\n\nSTDERR\n\n')
##             if os.path.isfile(self.runtime.proc_err):
##                 with open(self.runtime.proc_err) as stderr:
##                     for line in stderr:
##                         exe_info.write(line)
##         # if command succeed, remove all out_ and err_ files, 
##         self.runtime.clear()
## 
##     def _useVars(self, text, pipeline_vars={}):
##         # reduce extra newline, space etc
##         text = ' '.join(text.split())
##         #
##         # if there are pipeline vars, try to use pipeline vars to replace the 
##         # text 
##         if pipeline_vars:
##             files_and_dirs = []
##             for key, item in pipeline_vars.items():
##                 if isinstance(item, str) and (os.path.isfile(item) or os.path.isdir(item)):
##                     files_and_dirs.append([key, item])
##             # sort by length
##             files_and_dirs = sorted(files_and_dirs, key=lambda x: -len(x[1]))
##             #
##             # try to subsitute
##             pieces = text.split()
##             for key, item in files_and_dirs:
##                 for idx,p in enumerate(pieces):
##                     if p.startswith(item) and (len(p) == len(item) or not (p[len(item)].isalpha() or p[len(item)].isdigit())):
##                         pieces[idx] = '${{{}}}{}'.format(key, pieces[idx][len(item):])
##             text = ' '.join(pieces)
##         return text                   
##             
##     def __call__(self, ifiles, pipeline=None):
##         '''Execute action with input files ``ifiles`` with runtime information
##         stored in ``pipeline``. This function is called by the pipeline and calls
##         user-defined ``execute`` function.
## 
##         Parameters:
## 
##             ifiles (string or list of strings):
##                 input file names
## 
##             pipeline (an pipeline object):
##                 An Pipeline object for which the action is executed. The action
##                 can set or retrieve runtime information from a dictionary 
##                 ``pipeline.VARS``.
## 
##         Result:
##             An action returns output files (parameter ``output`` of the action)
##             if any output is given. Otherwise input files (``ifiles``) are passed
##             through and returned.
##         '''
##         if self.output:
##             if os.path.isfile(self.runtime.proc_info):
##                 with open(self.runtime.proc_info) as exe_info:
##                     cmd = exe_info.readline().strip()
##                 if self._useVars(cmd, pipeline.VARS) == self._useVars('; '.join(self.cmd), pipeline.VARS) \
##                     and existAndNewerThan(self.output, ifiles + pipeline.step_dependent_files,
##                     md5file=self.runtime.proc_info, pipeline=pipeline):
##                     env.logger.info('Reuse existing {}'.format(', '.join(self.output)))
##                     self._bypass(ifiles, pipeline)
##                     if self.output:
##                         return self.output
##                     else:
##                         return ifiles
##             # create directory if output directory does not exist
##             for d in [os.path.split(os.path.abspath(x))[0] for x in self.output]:
##                 if not os.path.isdir(d):
##                     try:
##                         os.makedirs(d)
##                     except Exception as e:
##                         raise RuntimeError('Failed to create directory {} for output file: {}'.format(d, e))
##         # We cannot ignore this step, but do we have all the input files?
##         # If not, we will have to rewind the execution
##         for ifile in ifiles:
##             if not os.path.isfile(ifile):
##                 env.logger.warning('Rewind execution because input file {} does not exist.'.format(ifile))
##                 raise RewindExecution(ifile)
##         #
##         if self.output:
##             with open(self.runtime.proc_info, 'w') as exe_info:
##                 exe_info.write('{}\n'.format(self._useVars('; '.join(self.cmd), pipeline.VARS)))
##                 exe_info.write('#Start: {}\n'.format(time.asctime(time.localtime())))
##                 for f in ifiles + pipeline.step_dependent_files:
##                     # for performance considerations, use partial MD5
##                     exe_info.write('{}\t{}\t{}\n'.format(self._useVars(f, pipeline.VARS), os.path.getsize(f),
##                         calculateMD5(f, partial=True)))
##         # now, run the job, write info if it is successfully finished.
##         # Otherwise the job might be forked and it will record the signature by itself.
##         ret = self._execute(ifiles, pipeline)
##         if ret not in [True, False]:
##             env.logger.warning('User defined execute function of a PipelineAction should return True or False')
##         if ret:
##             self._write_info(pipeline)#
##         #
##         if self.output:
##             return self.output
##         else:
##             return ifiles
## 
## # for backward compatibility
## SkiptableAction=PipelineAction
## 
## try:
##     from simulation import *
##     hasSimuPOP = True
## except ImportError as e:
##     hasSimuPOP = False
## 
## class SequentialActions(PipelineAction):
##     '''Define an action that calls a list of actions, specified by Action1,
##     Action2 etc. This allows the specification of multiple small tasks in
##     a single pipeline step.
## 
##     NOTE: this action is automatically applied if a list or tuple of actions
##     are specified in the SPEC file (e.g. action=Action1(), Action2()).
## 
## 
##     Examples:
##         action=CheckCommands('bowtie'), CheckOutput('bowtie --version', '1.1.1')
## 
##     '''
##     def __init__(self, actions):
##         '''
##         Parameters:
##             actions (a tuple or list of actions):
##                 A list of actions that will be applied to 
##         '''
##         self.actions = []
##         for a in actions:
##             if hasattr(a, '__call__'):
##                 self.actions.append(a.__call__)
##             else:
##                 self.actions.append(a)
## 
##     def __call__(self, ifiles, pipeline=None):
##         '''
##         Pass ifiles to the first action, take its output and pass it to
##         the second action, and so on. Return the output from the last
##         action as the result of this ``SequentialAction``.
##         '''
##         for a in self.actions:
##             # the input of the next action is the output of the
##             # previous action. However, ${INPUT} is the same for all 
##             # options because it is substitute before ...
##             ifiles = a(ifiles, pipeline)
##         # return the output of the last action
##         return ifiles
## 
## 
## class CheckVariantToolsVersion(PipelineAction):
##     '''Check the version of variant tools and determine if it is
##     recent enough to execute the pipeline.
## 
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Raises:
##         Fail if the version of variant tools used to execute the
##         pipeline is older than the specified version.
## 
##     Examples:
##         action=CheckVariantToolsVersion('2.5.0')
## 
##     '''
##     def __init__(self, version=''):
##         '''
##         Parameters:
##             version (string):
##                 Oldest version of variant tools that can be used
##                 to execute this pipeline
##         '''
##         self.min_version = version
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         vtools_version = [int(x) for x in re.sub('\D', ' ', pipeline.VARS['vtools_version']).split()]
##         # e.g. minimal 2.2.0, vtools 2.1.1
##         if [int(x) for x in re.sub('\D', ' ', self.min_version).split()] > vtools_version:
##             raise RuntimeError('Version {} is required to execute this pipeline. '
##                 'Please upgrade your installation of variant tools (version {})'
##                 .format(self.min_version, pipeline.VARS['vtools_version']))
##         return ifiles
## 
## 
## class ImportModules(PipelineAction):
##     '''Import functions and action from a Python module. This action passed input
##     files to output and does not change the pipeline.
##     
##     File Flow: Input passthrough, but import symbols to pipeline.
## 
##                      Pipeline
##                         ^
##             INPUT =============> INPUT
## 
##     Raises:
##         Raise a RuntimeError if one or more modules can not
##         be imported.
## 
##     Examples:
##         action=ImportModules('DNASeq_tools.py')
##         action=ImportModules(['DNASeq_tools.py', 'simuPOP.demography'])
##     '''
##     def __init__(self, modules=[], script=''):
##         '''Import one or more modules to be used by the existing pipeline. 
##         
##         Parameters:
##             module (string or list of strings):
##                 One or more module, which can be either the name of a system module
##                 or a .py file. In the latter case, Variant Tools will try to locate
##                 the file directly (a full path can be given), look for the module in
##                 the path of the pipeline (if a local pipeline is used), or download
##                 from the Variant Tools Repository under directory pipeline.
## 
##             script (string or list of strings):
##                 One or more in-line script that defines Python functions or customized
##                 actions that will be used in this pipeline. This allows users to define
##                 actions and utility functions that do not need to be shared with other
##                 pipelines but might be used repeatedly in this pipeline. Otherwise a
##                 ExecutePipelineFunction action can be used.
##         '''
##         if isinstance(modules, str):
##             self.modules = [modules]
##         else:
##             self.modules = modules
##         #
##         if isinstance(script, str):
##             self.script = script
##         else:
##             self.script = '\n'.join(script)
## 
##     def __call__(self, ifiles, pipeline=None):
##         for module in self.modules:
##             # this is a path to a .py file
##             if module.endswith('.py'):
##                 if os.path.isfile(module):
##                     pyfile = module
##                 # if the .py file locates in the same directory as the pipeline file
##                 elif pipeline is not None \
##                     and os.path.isfile(os.path.join(os.path.split(pipeline.spec_file)[0], module)):
##                     pyfile = os.path.join(os.path.split(pipeline.spec_file)[0], module)
##                 else:
##                     # try to download it from online
##                     try:
##                         pyfile = downloadFile('simulation/{}'.format(module))
##                     except Exception as e:
##                         try:
##                             pyfile = downloadFile('pipeline/{}'.format(module))
##                         except Exception as e:
##                             raise ValueError('Failed to download required python module {}: {}'.format(module, e))
##                 try:
##                     p,f = os.path.split(os.path.abspath(os.path.expanduser(pyfile)))
##                     sys.path.append(p)
##                     local_dict = __import__(f[:-3] if f.endswith('.py') else f, globals(), locals(), module.split('.', 1)[-1:])
##                     env.logger.info('{} symbols are imported form module {}'.format(len(local_dict.__dict__), module))
##                     pipeline.GLOBALS.update(local_dict.__dict__)
##                 except Exception as e:
##                     raise RuntimeError('Failed to import module {}: {}'.format(module, e))
##             # now a system module
##             else:
##                 try:
##                     # allow loading from current directory
##                     sys.path.append(os.getcwd())
##                     local_dict = __import__(module, globals(), locals(), module.split('.', 1)[-1:])
##                     env.logger.info('{} symbols are imported form module {}'.format(len(local_dict.__dict__), module))
##                     pipeline.GLOBALS.update(local_dict.__dict__)
##                 except ImportError as e:
##                     raise RuntimeError('Failed to import module {}: {}'.format(module, e))
##         # script
##         if self.script:
##             try:
##                 local_dict = {}
##                 exec(self.script, globals(), local_dict)
##                 env.logger.info('{} symbols are imported form inline script'.format(len(local_dict)))
##                 pipeline.GLOBALS.update(local_dict)
##             except Exception as e:
##                 raise RuntimeError('Failed to execute script "{}": {}'.format(self.script, e))
## 
##         return ifiles
## 
## 
## class CheckCommands(PipelineAction):
##     '''Check the existence of specified commands and raise an error if one of
##     the commands does not exist.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT 
## 
##     Raises:
##         A RuntimeError will be raised if a command is not found.
## 
##     Examples:
##         action=CheckCommands('java')
##         action=CheckCommands(['java', 'tophat2'])
##     '''
##     def __init__(self, commands):
##         '''
##         Parameters:
##             commands (string or list of strings):
##                 Name of one of more commands to be checked. No option is allowed.
##         '''
##         PipelineAction.__init__(self)
##         if type(commands) == type(''):
##             self.commands= [commands]
##         else:
##             self.commands= commands
## 
##     def __call__(self, ifiles, pipeline=None):
##         for cmd in self.commands:
##             if which(cmd) is None:
##                 raise RuntimeError('Command {} does not exist. Please install it and try again.'
##                     .format(cmd))
##             else:
##                 env.logger.info('Command {} is located.'.format(cmd))
##         return ifiles
## 
## 
## class CheckOutput(PipelineAction):
##     '''Run a command and check if its output matches at least one of specified
##     patterns. The pipeline will be terminated if failIfMismatch is set to True
##     (default). Otherwise a warning message will be printed.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Raises:
##         Raise a RuntimeError if the output of command does not match
##         any of the patterns, if ``failIfMismatch`` is set to ``True``.
## 
##     Examples:
##         action=CheckOutput('tophat2 --version', ['v2.0.13', 'v2.0.14'])
## 
##         # if strict_version is a command line parameter
##         action=CheckOutput('samtools', '0.1.19', %(strict_version)s)
##     '''
##     def __init__(self, command, patterns, failIfMismatch=True):
##         '''
##         Parameters:
##             command (string):
##                 A command (with or without options)
## 
##             patterns (string or list of strings):
##                 One or more patterns (usually a piece of version string)
##                 that will be compared to the output of ``command``
## 
##             failIfMismatch (boolean):
##                 If set to ``True`` (default), the action will terminate the
##                 pipeline if the output of command does not match any of the
##                 patterns. Otherwise a warning message will be printed when 
##                 the output of command does not match any of the patterns.
##         '''
##         self.command = command
##         if isinstance(patterns, str):
##             self.patterns = [patterns]
##         else:
##             self.patterns = patterns
##         self.fail = failIfMismatch
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         try:
##             # do not use subprocess.check_output because I need to get
##             # output even when the command returns non-zero return code
##             p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
##                 stderr=subprocess.PIPE, shell=True)
##             odata, edata = p.communicate()
##             output = odata.decode() + edata.decode()
##             env.logger.trace('Output of command "{}" is "{}"'
##                 .format(self.command, output))
##         except Exception as e:
##             raise RuntimeError('Failed to execute command "{}": {}'
##                 .format(self.cmd, e))
##         #
##         if all([re.search(x, output, re.MULTILINE) is None for x in self.patterns]):
##             msg = ('Output of command "{}" ("{}") does not ' + 
##                     'match specified regular expression {}.').format(self.command,
##                         ' '.join(output[:40].split()), ' or '.join(self.patterns))
##             if self.fail:
##                 raise RuntimeError(msg)
##             else:
##                 env.logger.warning(msg)
##         return ifiles
## 
## class CheckFiles(PipelineAction):
##     '''Check the existence of specified files and raise an
##     error if one of the files does not exist.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Raises:
##         Raise a RuntimeError if any of the files is not found.
## 
##     Example:
##         # assume gatk_path is a command line argument
##         action=CheckFile('%(gatk_path)s/GenomeAnalysisTK.jar',
##             'Please point --gatk_path to a directory with GenomeAnalysisTK.jar')
##     '''
##     def __init__(self, files, message=''):
##         '''
##         Parameters:
##             files (string or list of strings):
##                 One or more files to check.
## 
##             message (string):
##                 A message when one of the files cannot be found.
##         '''
##         if type(files) == str:
##             self.files = [files]
##         else:
##             self.files = files
##         self.message = message
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         for f in self.files:
##             if os.path.isfile(os.path.expanduser(f)):
##                 env.logger.info('{} is located.'.format(f))
##             else:
##                 raise RuntimeError('Cannot locate {}: {}'.format(f, self.message))
##         return ifiles
## 
## 
## class CheckDirs(PipelineAction):
##     '''Check the existence of specified directories and raise an
##     error if one of the directories does not exist.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Raises:
##         Raise a RuntimeError if any of the directories is not found.
## 
##     Example:
##         action=CheckDirs('${cmd_output}',
##             'Value of parameter --output need to be an existing directory')
##     '''
##     def __init__(self, dirs, message=''):
##         '''
##         Parameters:
##             files (string or list of strings):
##                 One or more directories to check.
## 
##             message (string):
##                 A message when one of the directories cannot be found.
##         '''
## 
##         if type(dirs) == str:
##             self.dirs = [dirs]
##         else:
##             self.dirs = dirs
##         self.message = message
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         for d in self.dirs:
##             if os.path.isdir(d):
##                 env.logger.info('Directory {} is located.'.format(d))
##             else:
##                 raise RuntimeError('Cannot locate directory {}. {}'.format(d, self.message))
##         return ifiles
## 
## 
## class TerminateIf(PipelineAction):
##     '''Terminate a pipeline if a condition is not met.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Raises:
##         A RuntimeError will be raised to terminate the pipeline if
##         the condition is met.
## 
##     Examples:
##         action=TerminateIf(not '${cmd_output}', 'No --output is specified.')
##     '''
##     def __init__(self, cond, message=''):
##         '''
##         Parameters:
##             cond (boolean):
##                 True or False. In practice, ``cond`` is usually
##                 a lambda function that checks the existence of a file or value
##                 of a pipeline variable.
## 
##             message (string):
##                 A message to be outputted when the condition is met.
##         '''
##         self.cond = cond
##         self.message = message
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         if self.cond:
##             raise RuntimeError(self.message)
##         return ifiles
## 
## 
## class WarnIf(PipelineAction):
##     '''Send a warning message if a condition is not met.
## 
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Examples:
##         action=WarnIf('%(LGD)' == 'NA', 'Default value of parameter --LGD is used.')
##     '''
##     def __init__(self, cond, message):
##         '''
##         Parameters:
##             cond (boolean):
##                 True or False. In practice, ``cond`` is usually
##                 a lambda function that checks the existence of a file or value
##                 of a pipeline variable.
## 
##             message (string):
##                 A message to be outputted when the condition is met.
##         '''
##         self.cond = cond
##         self.message = message
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         if self.cond:
##             env.logger.warning(self.message)
##         return ifiles
## 
## 
## class OutputText(PipelineAction):
##     '''Write specified text to standard output, or a file if a filename is
##     specified. The text can be a list of strings. A new line is added 
##     automatically to each line of the text.
## 
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Examples:
##         action=OutputText('Hey, the biggest part is done.')
##     
##     '''
##     def __init__(self, text='', output=None, mode='a'):
##         '''
##         Parameters:
##             text (string or list of strings):
##                 Text to be written to output.
## 
##             output (a file name or None):
##                 Output files. The text will be written to standard output
##                 if no output is specified.
## 
##             mode (string):
##                 Mode to open file. 'a' for append and 'w' for overwrite.
## 
##         '''
##         if not isinstance(text, str):
##             self.text = ''.join([str(x) + '\n' for x in text])
##         else:
##             self.text = text + '\n'
##         self.filename = filename
##         self.mode = mode
##         PipelineAction.__init__(self, 'OutputText', filename if filename is not None else '')
## 
##     def __call__(self, ifiles, pipeline=None):
##         if self.filename is not None:
##             with open(self.filename, self.mode) as output:
##                 output.write(self.text)
##         else:
##             sys.stdout.write(self.text)        
##         return ifiles
## 
## 
## class FieldsFromTextFile(PipelineAction):
##     '''Read a text file, guess its delimeter, field name (from header)
##     and create field descriptions. If a vcf file is encountered, all
##     fields will be exported.
## 
##     File Flow: extract format of input and output format.
## 
##         INPUT ==> Get Format ==> OUTPUT
## 
##     Raises:
##         Raise a RuntimeError if this action failed to guess format (fields)
##         from the input file.
## 
##     Examples:
##         action=FieldsFromTextFile('format.txt')
## 
##     '''
##     def __init__(self, output):
##         '''
##         Parameters:
##             output:
##                 Output file that records the format of the input files.
##         '''
##         PipelineAction.__init__(self, 'FieldsFromTextFile', output)
## 
##     def _execute(self, ifiles, pipeline=None):
##         if len(ifiles) > 1:
##             env.logger.warning('Only the format of the first input file would be outputted.')
##         try:
##             if ifiles[0].endswith('.vcf') or ifiles[0].endswith('.vcf.gz'):
##                 showTrack(ifiles[0], self.output[0])
##             else:
##                 with open(self.output[0], 'w') as fo:
##                     csv_dialect = csv.Sniffer().sniff(open(ifiles[0], 'rU').read(8192))
##                     fo.write('delimiter="{}"\n\n'.format(csv_dialect.delimiter.replace('\t', r'\t')))
##                     values = []
##                     with open(ifiles[0], 'rU') as fi:
##                         reader = csv.reader(fi, dialect=csv_dialect)
##                         headers = reader.next()
##                         values = [[] for x in headers]
##                         for line in reader:
##                             for idx in range(len(headers)):
##                                 values[idx].append(line[idx])
##                             if len(values[0]) > 100:
##                                 break
##                     #
##                     for idx, header in enumerate(headers):
##                         fo.write('[{}]\n'.format(validFieldName(header)))
##                         fo.write('index={}\n'.format(idx+1))
##                         fo.write('type={}\n\n'.format(typeOfValues(values[idx])))
##         except Exception as e:
##             raise RuntimeError('Failed to guess fields from {}: {}'.format(ifiles[0], e))
##         #
##         return True
##        
## class RewindExecution(Exception):
##     pass
## 
## class NullAction(PipelineAction):
##     '''A pipeline action that does nothing. This is usually used when the goal
##     of the step is to change input, output, or assign variables to pipelines.
##     The action will be assumed if an empty action line is given.
## 
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Example:
##         action=
##         action=NullAction()
##     '''
##     def __init__(self, output=[], *args, **kwargs):
##         '''A null action that does nothing.'''
##         PipelineAction.__init__(self, cmd='NullAction', output=output)
## 
##     def _execute(self, ifiles, pipeline=None):
##         return True
##         
## class MonitorThread(threading.Thread):
##     def __init__(self, group=None, target=None, name=None,
##                  args=(), kwargs={}, Verbose=None):
##         threading.Thread.__init__(self, group, target, name, args, kwargs, Verbose)
##         self._return = None
## 
##     def run(self):
##         if self._Thread__target is not None:
##             self._return = self._Thread__target(*self._Thread__args,
##                                                 **self._Thread__kwargs)
##     def join(self):
##         threading.Thread.join(self)
##         return self._return
## 
## class SharedProcess:
##     def __init__(self, runtime):
##         self.runtime = runtime
## 
##     def update_prog(self):
##         while True:
##             with open(self.runtime.proc_prog, 'w') as prog:
##                 prog.write(' ')
##             time.sleep(30)
## 
##     def __enter__(self):
##         # wait for the availability of lock
##         start_time = time.time()
##         prog_time = None
##         while True:
##             # there is a progress file
##             if os.path.isfile(self.runtime.proc_prog):
##                 if prog_time is None:
##                     env.logger.trace('Job started with progress file {}'.format(self.runtime.proc_prog))
##                 os.remove(self.runtime.proc_prog)
##                 prog_time = time.time()
##             #
##             if prog_time is None or time.time() - prog_time > 60:
##                 # no progress by another thread
##                 break
##             #
##             if os.path.isfile(self.runtime.proc_done):
##                 break
##         #
##         # start monitoring process
##         self.update_process = Process(target=self.update_prog)
##         self.update_process.start()
## 
##     def __exit__(self, type, value, traceback):
##         with open(self.runtime.proc_done, 'w') as prog:
##             prog.write('1')
##         self.runtime.clear(['prog'])
##         self.update_process.terminate()
## 
## class RunCommand(PipelineAction):
##     """This action execute specified commands. If the pipeline is running
##     in parallel mode and a submitter is specified, it will use the submitter
##     command to execute the commands in a separate job. 
##     
##     File Flow:
## 
##         Input passthrough if no output file is specified.
##             INPUT ====> INPUT
##         Generate output if one or more output files are specified.
##             INPUT ==> CMD ==> OUTPUT
##         
##     Raises:
##         Raises an error if an command fails to execute.
## 
##     Examples:
##         # simple commands without checking output
##         action=RunCommand(cmd='vtools init myproj -f')
##         
##         action=RunCommand(cmd=[
##             '[ -d ${DIR1} ] || mkdir -p ${DIR1}',
##             '[ -d ${DIR2} ] || mkdir -p ${DIR2}'
##             ])
##         
##         # multiple commands, change working directory
##         # check output
##         action=RunCommand([
##         	'update_blastdb.pl human_genomic --decompress || true',
##         	'update_blastdb.pl nt --decompress || true',
##         	],
##             working_dir='${NCBI_RESOURCE_DIR}/blast/db',
##         	output=['${NCBI_RESOURCE_DIR}/blast/db/human_genomic.nal',
## 		        '${NCBI_RESOURCE_DIR}/blast/db/nt.nal']
##             )
##         
##         # run command in background, with pipes
##         action=RunCommand('''samtools view -h ${INPUT}
##            | awk '$6 ~/N/' | awk '{ if ($9 ~ /^-/) {print $1"\t-"} else print $1"\t+" }'
##            | sort -T ${TEMP_DIR} -u | wc -l > ${ALIGNMENT_OUT}/junction.count''',
##            output='${ALIGNMENT_OUT}/junction.count',
##            submitter='sh {} &')
## 
##     """
##     def __init__(self, cmd='', output=[], working_dir=None, submitter=None, wait=True, max_jobs=None):
##         '''This action accepts one (a string) or more command (a list of strings)
##         and executes them in a shell environment, possibly as a separate job. 
##         
##         Parameters:
##             cmd (string or list of strings):
##                 One or more commands to execute. The commands will be executed in
##                 shell mode so pipes are allowed.
## 
##             output (string or list of strings):
##                 Expected output files of the action. If specified, the execution
##                 signature will be created to record the input, output and command
##                 of the action, and ignore the action if the signature matches of 
##                 a previous run.
## 
##             working_dir (None or string):
##                 Working directory of the command. Variant Tools will change to
##                 this directory before executing the commands if a valid directory
##                 is passed.
## 
##             submitter (None or string):
##                 If a submitter is specified and the pipeline is executed in multi-job
##                 mode (e.g. --jobs 2), a shell script will be written with the commands
##                 to be executed. The submitter command will be executed with ``{}`` in
##                 parameter ``submitter`` replaced by the name of shell script. For
##                 example, submitter='sh {} &' will run the job as a background job,
##                 and submitter='qsub -q long < {}' will submit the shell script to the
##                 long queue of a cluster system. Because the pipeline will be terminated
##                 if the submitter command fails, `qsub new_job ... && false` can be used
##                 to replace the running process by start a new job and terminate the
##                 existing process intentionally.
## 
##             wait (True, False, or number of seconds):
##                 If a job is submitted, whether or not wait it to be completed. The default
##                 is True, meaning that the master thread will continue to execute until
##                 is will be waiting for the outcome of this command. If you set this parameter
##                 to False, the pipeline execute will be stopped and you can re-run the
##                 pipeline till the subcommand is completed. You can also set a number
##                 to let the master thread wait for a pre-determined period of time. This 
##                 option is useful if the subprocess might die.
## 
##             max_jobs: (deprecated)
##         '''
##         self.submitter = submitter
##         self.working_dir = working_dir
##         self.wait = wait
##         if type(output) == str:
##             self.output = [os.path.expanduser(output)]
##         else:
##             self.output = [os.path.expanduser(x) for x in output]
##         #
##         self.runtime = RuntimeFiles(self.output)
##         self.start_time = time.time()
##         if not cmd:
##             cmd = ['echo "None command executed."']
##         PipelineAction.__init__(self, cmd=cmd, output=output)
## 
##     def _elapsed_time(self):
##         '''Return the elapsed time in human readable format since start time'''
##         second_elapsed = int(time.time() - self.start_time)
##         days_elapsed = second_elapsed // 86400
##         return ('{} days '.format(days_elapsed) if days_elapsed else '') + \
##             time.strftime('%H:%M:%S', time.gmtime(second_elapsed % 86400))
##       
##     def _run_command(self):
##         '''Call a list of external command cmd, raise an error if any of them
##         fails. '''
##         if self.runtime.proc_lck:
##             env.lock(self.runtime.proc_lck, str(os.getpid()))
##         for cur_cmd in self.cmd:
##             if self.working_dir is not None and not os.path.isdir(self.working_dir):
##                 os.makedirs(self.working_dir)
##             env.logger.info('Running ``{}``{}'.format(cur_cmd,
##                 ' under {}'.format(self.working_dir) if self.working_dir else ''))
##             ret = subprocess.call(cur_cmd, shell=True, 
##                 stdout=None if self.runtime.proc_out is None else open(self.runtime.proc_out, 'w'),
##                 stderr=None if self.runtime.proc_err is None else open(self.runtime.proc_err, 'w'),
##                 cwd=self.working_dir)
##             if ret < 0:
##                 if self.output:
##                     try:
##                         env.unlock(self.runtime.proc_lck, str(os.getpid()))
##                     except:
##                         pass
##                 raise RuntimeError("Command '{}' was terminated by signal {} after executing {}"
##                     .format(cur_cmd, -ret, self._elapsed_time()))
##             elif ret > 0:
##                 if self.output:
##                     with open(self.runtime.proc_err) as err:
##                         for line in err.read().split('\n')[-50:]:
##                             env.logger.error(line)
##                     try:
##                         env.unlock(self.runtime.proc_lck, str(os.getpid()))
##                     except:
##                         pass
##                 raise RuntimeError("Execution of command '{}' failed after {} (return code {})."
##                     .format(cur_cmd, self._elapsed_time(), ret))
## 
##     def _monitor(self):
##         start_time = time.time()
##         prog_time = None
##         while True:
##             if os.path.isfile(self.runtime.proc_prog):
##                 if prog_time is None:
##                     env.logger.trace('Job started with progress file {}'.format(self.runtime.proc_prog))
##                 os.remove(self.runtime.proc_prog)
##                 prog_time = time.time()
##             #
##             if prog_time is None:
##                 # if the job has not been started for 10 minutes, quite
##                 if time.time() - start_time > 600:
##                     return('Background job has not been started after 10 minutes.')
##             else:
##                 if time.time() - prog_time > 60:
##                     return('Background job has not updated it progress for 1 minutes.')
##             #
##             if os.path.isfile(self.runtime.proc_done):
##                 break
##             else:
##                 if self.wait is False:
##                     return('Do not wait for the completion of submitted job (wait=False).')
##                 if self.wait is not True and isinstance(self.wait, int) and prog_time is not None and time.time() - prog_time > self.wait:
##                     return('Quitted after waiting {} seconds.'.format(self.wait))
##                 time.sleep(10)
##         try:
##             env.unlock(self.runtime.proc_lck, str(os.getpid()))
##         except:
##             env.logger.warning('Failed to remove lock for file {}'.format(self.output[0]))
##             pass
##         try:
##             with open(self.runtime.proc_done) as done:
##                 ret = int(done.read().strip())
##         except Exception as e:
##             return('Failed to retrive return information for forked process from {}. {}'
##                 .format(self.runtime.proc_done, e))
##         #
##         if ret < 0:
##             return("Command '{}' was terminated by signal {} after executing {}"
##                 .format('; '.join(self.cmd), -ret, self._elapsed_time()))
##         elif ret > 0:
##             if self.output:
##                 with open(self.runtime.proc_err) as err:
##                     for line in err.read().split('\n')[-50:]:
##                         env.logger.error(line)
##             return("Execution of command '{}' failed after {} (return code {})."
##                 .format('; '.join(self.cmd), self._elapsed_time(), ret))
##         # remove the .done file
##         if not self.output[0] in self.pipeline.THREADS:
##             return('Output is not waited by any threads')
##         # DO NOT POP FROM ANOTHER THREAD, this will cause race condition
##         # (unless we use thread safe dictionry). In this case, we only need
##         # to monitor the status of threads from the master threads.
##         #    self.pipeline.THREADS.pop(self.output[0])
##         #
##         # the thread will end here
##         env.logger.trace('Thread for output {} ends.'.format(self.output[0]))
##         self.runtime.clear(['done'])
##         return('')
## 
##     def _submit_command(self):
##         '''Submit a job and wait for its completion.'''
##         # use full path because the command might be submitted to a remote machine
##         if os.path.isfile(self.runtime.proc_done):
##             os.remove(self.runtime.proc_done)
##         if self.runtime.proc_lck:
##             env.lock(self.runtime.proc_lck, str(os.getpid()))
##         #
##         if os.path.isfile(self.runtime.proc_cmd):
##             with open(self.runtime.proc_cmd) as old_cmd:
##                 old_script = old_cmd.read()
##         else:
##             old_script = None
##         # create a batch file for execution
##         with open(self.runtime.proc_cmd, 'w') as sh_file:
##             sh_file.write('#PBS -o {}\n'.format(os.path.abspath(self.runtime.proc_out)))
##             sh_file.write('#PBS -e {}\n'.format(os.path.abspath(self.runtime.proc_err)))
##             sh_file.write('#PBS -N {}\n'.format(os.path.basename(self.output[0])))
##             #sh_file.write('#PBS -N {}.{}_{}\n'.format(self.runtime.proc_err))
##             sh_file.write('#PBS -V\n')
##             # we try to reproduce the environment as much as possible becaus ehte
##             # script might be executed in a different environment
##             for k, v in os.environ.items():
##                 if any([k.startswith(x) for x in ('SSH', 'PBS', '_')]) or not k.replace('_', '').isalpha():
##                     continue
##                 sh_file.write('export {}="{}"\n'.format(k, v.replace('\n', '\\n')))
##             #
##             sh_file.write('\ncd {}\n'.format(os.path.abspath(os.getcwd())))
##             if self.working_dir is not None:
##                 sh_file.write('[ -d {0} ] || mkdir -p {0}\ncd {0}\n'.format(os.path.abspath(self.working_dir)))
##             #
##             sh_file.write('''
## progress() {{
##   while true
##   do
##     touch {}
##     sleep 30
##   done
## }}
## 
## progress &
## MYSELF=$!
## '''.format(self.runtime.proc_prog))
##             sh_file.write('\n'.join(self.cmd))
##             #
##             sh_file.write('\n\nCMD_RET=$?\nif [ $CMD_RET == 0 ]; then vtools admin --record_exe_info {} {}; fi\n'
##                 .format(os.getpid(), ' '.join(self.output)))
##             # a signal to show the successful completion of the job
##             sh_file.write('\nrm -f {}\nkill $MYSELF >/dev/null 2>&1\necho $CMD_RET > {}\n'
##                 .format(self.runtime.proc_prog, self.runtime.proc_done))
##         #
##         if old_script is not None:
##             #with open(self.runtime.proc_cmd) as new_cmd:
##             #     if old_script == new_cmd.read():
##             #         env.logger.debug('Identical script {}'.format(self.runtime.proc_cmd))
##                      # if there is no change in command
##                      other_prog = glob.glob(os.path.abspath(self.output[0]) + '.working_*')
##                      if other_prog:
##                          for op in other_prog:
##                              # if the working file is less than 2 minutes old, ...
##                              if time.time() - os.path.getmtime(op) < 120:
##                                  env.logger.info('Another process appears to be working on {}, checking ...'.format(self.output[0]))
##                                  last_time = os.path.getmtime(op)
##                                  time.sleep(60)
##                                  # if the working file does not change after 60 seconds
##                                  if os.path.getmtime(op) != last_time:
##                                      raise RuntimeError('Failed to submit job because a job is currently running or has been failed within 2 minutes. Status file is {} (pid is {})'.format(op, os.getpid()))
##         # try to submit command
##         if '{}' in self.submitter:
##             submit_cmd = self.submitter.replace('{}', self.runtime.proc_cmd)
##         else:
##             submit_cmd = self.submitter
##         #
##         env.logger.info('Running job {} with command "{}" from directory {}'.format(
##             self.runtime.proc_cmd, submit_cmd, os.getcwd()))
##         ret = subprocess.call(submit_cmd, shell=True,
##             stdout=open(self.runtime.proc_out, 'w'), stderr=open(self.runtime.proc_err, 'w'),
##             cwd=self.working_dir)
##         if ret != 0:
##             try:
##                 env.unlock(self.runtime.proc_out, str(os.getpid()))
##             except:
##                 pass
##         # 
##         if ret < 0:
##             raise RuntimeError("Failed to submit job {} due to signal {} (submitter='{}')" .format(self.runtime.proc_cmd, -ret, self.submitter))
##         elif ret > 0:
##             if os.path.isfile(self.runtime.proc_err):
##                  with open(self.runtime.proc_err) as err:
##                      msg = err.read()
##             else:
##                 msg = ''
##             raise RuntimeError("Failed to submit job {} using submiter '{}': {}".format(self.runtime.proc_cmd, self.submitter, msg))
##         else:
##             t = MonitorThread(target=self._monitor)
##             t.daemon = True
##             t.start()
##             if self.output[0] in self.pipeline.THREADS:
##                 raise RuntimeError('Two spawned jobs have the same self.output[0] file {}'.format(self.output[0]))
##             self.pipeline.THREADS[self.output[0]] = t
## 
## 
##     def _execute(self, ifiles, pipeline=None):
##         # substitute cmd by input_files and output_files
##         if pipeline.jobs > 1 and self.submitter is not None and not self.output:
##             env.logger.warning('Fail to execute in parallel because no output is specified.')
##         #
##         self.pipeline = pipeline
##         # Submit the job on a cluster system
##         # 1. if there is output, otherwise we cannot track the status of the job
##         # 2. if a submit command is specified
##         # 3. if --jobs with a value greater than 1 is used.
##         if self.output and pipeline.jobs > 1 and self.submitter is not None:
##             self._submit_command()
##             return False
##         else:
##             self._run_command()
##             return True
## 
## class ExecutePythonCode(PipelineAction):
##     '''This action execute a piece of Python code under the pipeline namespace, which
##     means all pipeline variables will be available to the code. This action provides a
##     way to implement pipeline actions on the fly. Arbitary parameters can be passed
##     and be made available to the script. Pipeline variables are also available to the 
##     script as a variabe "pvars".
##    
##     File Flow:
## 
##         Input passthrough if no output file is specified.
##             INPUT ====> INPUT
##         Generate output if one or more output files are specified.
##             INPUT ==> CMD ==> OUTPUT
##         
##     Raises:
##         Raises an error if the python code fails to execute.
## 
##     '''
##     def __init__(self, script='', output=[], modules=[], export=None, **kwargs):
##         '''This action accepts one or a list of strings and execute it as a piece of Python
##         code. Pipeline variables are made available as a dictionary "pvars".
## 
##         Parameters:
##             script (string or list of strings):
##                 One or more strings to execute. List of strings will be concatenated by new
##                 lines.
## 
##             modules (string or list of strings):
##                 Modules to import for this step. It is similar to action ImportModules but
##                 the imported symbols are available to this action only.
## 
##             export (None or filename):
##                 A filename to which the execute code will be exported. This option is useful
##                 for debugging because pipeline variables will be prepended to the script so
##                 that the exported script could be executed with no or minimal modification.
## 
##             kwargs (additional parameters):
##                 Any additional kwargs will be passed to the function executed.
##         '''
##         if not script:
##             env.logger.warning('No valid script is specified.')
##             script = ''
##         if isinstance(script, str):
##             self.script = script
##         else:
##             self.script = '\n'.join(script)
##         #
##         m = hashlib.md5()
##         m.update(self.script.encode('utf-8'))
##         #
##         self.kwargs = kwargs
##         self.modules = modules
##         self.export = export
##         #
##         PipelineAction.__init__(self, cmd='python -e {} {}'.format(m.hexdigest(), kwargs), output=output)
## 
##     def _execute(self, ifiles, pipeline=None):
##         if self.export is not None:
##             with open(self.export, 'w') as exported_script:
##                 exported_script.write('#!/usr/env python\n')
##                 exported_script.write('#\n#Script exported by action ExecutePythonCode\n')
##                 # modules?
##                 exported_script.write(''.join(['import {}\n'.format(x) for x in ('sys', 'os', 're', 'glob')]))
##                 exported_script.write('\nfrom variant_tools.pipeline import *\n')
##                 for module in self.modules:
##                     exported_script.write('import {}\n'.format(os.path.basename(module)[:-3] if module.endswith('.py') else module))
##                 #
##                 # pipeline variables
##                 exported_script.write('pvars=')
##                 pprint.pprint(pipeline.VARS.dict(), stream=exported_script)
##                 exported_script.write('\n')
##                 #
##                 exported_script.write('# Pipeline variables are case insensitive\n')
##                 exported_script.write('pvars.update({x.upper():y for x,y in pvars.items()})\n')
##                 exported_script.write('pvars.update({x.lower():y for x,y in pvars.items()})\n')
##                 #
##                 # script
##                 exported_script.write(self.script)
##             env.logger.info('Python code exported to ``{}``'.format(self.export))
##         for module in self.modules:
##             # this is a path to a .py file
##             if module.endswith('.py'):
##                 if os.path.isfile(module):
##                     pyfile = module
##                 # if the .py file locates in the same directory as the pipeline file
##                 elif pipeline is not None \
##                     and os.path.isfile(os.path.join(os.path.split(pipeline.spec_file)[0], module)):
##                     pyfile = os.path.join(os.path.split(pipeline.spec_file)[0], module)
##                 else:
##                     # try to download it from online
##                     try:
##                         pyfile = downloadFile('simulation/{}'.format(module))
##                     except Exception as e:
##                         try:
##                             pyfile = downloadFile('pipeline/{}'.format(module))
##                         except Exception as e:
##                             raise ValueError('Failed to download required python module {}: {}'.format(module, e))
##                 try:
##                     p,f = os.path.split(os.path.abspath(os.path.expanduser(pyfile)))
##                     sys.path.append(p)
##                     local_dict = __import__(f[:-3] if f.endswith('.py') else f, globals(), locals(), module.split('.', 1)[-1:])
##                     env.logger.info('{} symbols are imported form module {}'.format(len(local_dict.__dict__), module))
##                     pipeline.GLOBALS.update(local_dict.__dict__)
##                 except Exception as e:
##                     raise RuntimeError('Failed to import module {}: {}'.format(module, e))
##             # now a system module
##             else:
##                 try:
##                     # allow loading from current directory
##                     sys.path.append(os.getcwd())
##                     local_dict = __import__(module, globals(), locals(), module.split('.', 1)[-1:])
##                     env.logger.info('{} symbols are imported form module {}'.format(len(local_dict.__dict__), module))
##                     pipeline.GLOBALS.update(local_dict.__dict__)
##                 except ImportError as e:
##                     raise RuntimeError('Failed to import module {}: {}'.format(module, e))
##         env.logger.info('Executing Python script:\n{}'.format(self.script))
##         try:
##             globals().update(pipeline.GLOBALS)
##             local_dict = self.kwargs
##             local_dict['pvars'] = pipeline.VARS
##             exec(self.script, globals(), local_dict)
##         except Exception as e:
##             ex_type, ex, tb = sys.exc_info()
##             traceback.print_tb(tb)
##             raise RuntimeError('Failed to execute script: {}'.format(e))
##         return True
## 
## 
## class ExecuteScript(PipelineAction):
##     """This action execute specified in-line script with specified command (bash, python, perl
##     etc). If the pipeline is running in parallel mode and a submitter is specified, it will use
##     the submitter command to execute the commands in a separate job. This action is the base
##     action for ExecuteRScript, ExecutePerlScript and ExecutePythonScript and is usually not
##     used directly.
##     
##     
##     File Flow:
## 
##         Input passthrough if no output file is specified.
##             INPUT ====> INPUT
##         Generate output if one or more output files are specified.
##             INPUT ==> CMD ==> OUTPUT
##         
##     Raises:
##         Raises an error if an command fails to execute.
## 
##     """
##     def __init__(self, script='', interpreter='', args='', output=[], working_dir=None, 
##          export=None, submitter=None, suffix=None, wait=True):
##         '''This action accepts one or a list of strings, write them to a temporary file
##         and executes them by a interpreter, possibly as a separate job. 
##         
##         Parameters:
##             script (string or list of strings):
##                 One or more strings to execute. List of strings will be concatenated by new
##                 lines. The complete script will be written to a temporary file to be executed
##                 by an interpreter.
## 
##             interpreter (string):
##                 An interpreter that will be used to execute the script. It is usually just a command
##                 but more complex command line is allowed with '{}' replaced by the path to the
##                 temporary script.
## 
##             args (string or list of strings):
##                 Command line arguments which can be a single string or a list of strings.
##                 Filenames will be properly quoted if needed.
## 
##             output (string or list of strings):
##                 Expected output files of the action. If specified, the execution
##                 signature will be created to record the input, output and command
##                 of the action, and ignore the action if the signature matches of 
##                 a previous run.
## 
##             working_dir (None or string):
##                 Working directory of the command. Variant Tools will change to
##                 this directory before executing the commands if a valid directory
##                 is passed.
## 
##             export (None or string):
##                 A filename to which the script will be exported before execution. 
##                 This option makes it easier to debug the script because the script in the
##                 spec file might contain pipeline variables.
## 
##             submitter (None or string):
##                 If a submitter is specified and the pipeline is executed in multi-job
##                 mode (e.g. --jobs 2), the script will be submitted by the submitter.
##                 The submitter command will be executed with ``{}`` in
##                 parameter ``submitter`` replaced by the name of shell script. For
##                 example, submitter='sh {} &' will run the job as a background job,
##                 and submitter='qsub -q long < {}' will submit the shell script to the
##                 long queue of a cluster system. Because the pipeline will be terminated
##                 if the submitter command fails, `qsub new_job ... && false` can be used
##                 to replace the running process by start a new job and terminate the
##                 existing process intentionally.
## 
##             suffix (None or string):
##                 An optional suffix (file extension) to the temporary script.
## 
##             wait (True, False, or number of seconds):
##                 If a job is submitted, whether or not wait it to be completed. The default
##                 is True, meaning that the master thread will continue to execute until
##                 is will be waiting for the outcome of this command. If you set this parameter
##                 to False, the pipeline execute will be stopped and you can re-run the
##                 pipeline till the subcommand is completed. You can also set a number
##                 to let the master thread wait for a pre-determined period of time. This 
##                 option is useful if the subprocess might die.
##         '''
##         self.interpreter = interpreter
##         self.submitter = submitter
##         self.working_dir = working_dir
##         self.wait = wait
##         if type(output) == str:
##             self.output = [os.path.expanduser(output)]
##         else:
##             self.output = [os.path.expanduser(x) for x in output]
##         #
##         self.runtime = RuntimeFiles(self.output)
##         self.start_time = time.time()
##         if not script:
##             raise ValueError('No valid script is specified.')
##         if isinstance(script, str):
##             self.script = script
##         else:
##             self.script = '\n'.join(script)
##         #
##         self.args = args
##         env.logger.info('Executing\n{}'.format(self.script))
##         #
##         m = hashlib.md5()
##         m.update(self.script.encode('utf-8'))
##         #
##         self.script_file = tempfile.NamedTemporaryFile(mode='w+t', suffix=suffix, delete=False).name
##         with open(self.script_file, 'w') as script_file:
##             script_file.write(self.script)
##         #
##         if export is not None:
##             with open(export, 'w') as exported_script:
##                 exported_script.write(self.script)
##                 env.logger.info('Script exported to ``{}``'.format(export))
##         PipelineAction.__init__(self, cmd='{} {}'.format(interpreter, m.hexdigest()), output=output)
## 
##     def __del__(self):
##         try:
##             os.remove(self.script_file)
##         except Exception as e:
##             env.logger.debug('Failed to remove temporary script file {}: {}'.format(self.script_file, e))
## 
##     def _elapsed_time(self):
##         '''Return the elapsed time in human readable format since start time'''
##         second_elapsed = int(time.time() - self.start_time)
##         days_elapsed = second_elapsed // 86400
##         return ('{} days '.format(days_elapsed) if days_elapsed else '') + \
##             time.strftime('%H:%M:%S', time.gmtime(second_elapsed % 86400))
##       
##     def _run_command(self):
##         '''Call a list of external command cmd, raise an error if any of them
##         fails. '''
##         if self.runtime.proc_lck:
##             env.lock(self.runtime.proc_lck, str(os.getpid()))
##         if '{}' in self.interpreter:
##             cmd = self.interpreter.replace('{}', pipes.quote(self.script_file))
##         else:
##             cmd = self.interpreter + ' ' + pipes.quote(self.script_file) + \
##                 (self.args if isinstance(self.args, str) else ' '.join(pipes.quote(x) for x in self.args))
##         env.logger.info('Running ``{}``'.format(cmd))
##         ret = subprocess.call(cmd, shell=True, 
##             stdout=None if self.runtime.proc_out is None else open(self.runtime.proc_out, 'w'),
##             stderr=None if self.runtime.proc_err is None else open(self.runtime.proc_err, 'w'),
##             cwd=self.working_dir)
##         if ret < 0:
##             if self.output:
##                 try:
##                     env.unlock(self.runtime.proc_lck, str(os.getpid()))
##                 except:
##                     pass
##             raise RuntimeError("Command '{}' was terminated by signal {} after executing {}"
##                 .format(cmd, -ret, self._elapsed_time()))
##         elif ret > 0:
##             if self.output:
##                 with open(self.runtime.proc_err) as err:
##                     for line in err.read().split('\n')[-50:]:
##                         env.logger.error(line)
##                 try:
##                     env.unlock(self.runtime.proc_lck, str(os.getpid()))
##                 except:
##                     pass
##             raise RuntimeError("Execution of command '{}' failed after {} (return code {})."
##                 .format(cmd, self._elapsed_time(), ret))
##         else:
##             # write standard out to terminal
##             if self.runtime.proc_out:
##                 with open(self.runtime.proc_out) as proc_out:
##                     for line in proc_out:
##                         env.logger.info(line.rstrip())
##             if self.runtime.proc_err:
##                 with open(self.runtime.proc_err) as proc_err:
##                     for line in proc_err:
##                         env.logger.warning(line.rstrip())
## 
##     def _monitor(self):
##         start_time = time.time()
##         prog_time = None
##         while True:
##             if os.path.isfile(self.runtime.proc_prog):
##                 if prog_time is None:
##                     env.logger.trace('Job started with progress file {}'.format(self.runtime.proc_prog))
##                 os.remove(self.runtime.proc_prog)
##                 prog_time = time.time()
##             #
##             if prog_time is None:
##                 # if the job has not been started for 10 minutes, quite
##                 if time.time() - start_time > 600:
##                     return('Background job has not been started after 10 minutes.')
##             else:
##                 if time.time() - prog_time > 60:
##                     return('Background job has not updated it progress for 1 minutes.')
##             if os.path.isfile(self.runtime.proc_done):
##                 break
##             else:
##                 if self.wait is False:
##                     return
##                 if self.wait is not True and isinstance(self.wait, int) and prog_time is not None and time.time() - prog_time > self.wait:
##                     return('Quitted after waiting {} seconds.'.format(self.wait))
##                 time.sleep(10)
##         try:
##             env.unlock(self.runtime.proc_lck, str(os.getpid()))
##         except:
##             env.logger.warning('Failed to remove lock for file {}'.format(self.output[0]))
##             pass
##         with open(self.runtime.proc_done) as done:
##             ret = int(done.read().strip())
##         #
##         if ret < 0:
##             return("Command '{}' was terminated by signal {} after executing {}"
##                 .format('; '.join(self.cmd), -ret, self._elapsed_time()))
##         elif ret > 0:
##             if self.output:
##                 with open(self.runtime.proc_err) as err:
##                     for line in err.read().split('\n')[-50:]:
##                         env.logger.error(line)
##             return("Execution of command '{}' failed after {} (return code {})."
##                 .format('; '.join(self.cmd), self._elapsed_time(), ret))
##         # remove the .done file
##         if not self.output[0] in self.pipeline.THREADS:
##             return('Output is not waited by any threads')
##         # DO NOT POP FROM ANOTHER THREAD, this will cause race condition
##         # (unless we use thread safe dictionry). In this case, we only need
##         # to monitor the status of threads from the master threads.
##         #    self.pipeline.THREADS.pop(self.output[0])
##         #
##         # the thread will end here
##         env.logger.info('{} has been successfully generated.'.format(self.output[0]))
##         self.runtime.clear(['done'])
##         return('')
## 
##     def _submit_command(self):
##         '''Submit a job and wait for its completion.'''
##         # use full path because the command might be submitted to a remote machine
##         #
##         if os.path.isfile(self.runtime.proc_done):
##             os.remove(self.runtime.proc_done)
##         if self.runtime.proc_lck:
##             env.lock(self.runtime.proc_lck, str(os.getpid()))
##         #
##         if os.path.isfile(self.runtime.proc_cmd):
##             with open(self.runtime.proc_cmd) as old_cmd:
##                 old_script = old_cmd.read()
##         else:
##             old_script = None
##         # create a batch file for execution
##         with open(self.runtime.proc_cmd, 'w') as sh_file:
##             sh_file.write('#PBS -o {}\n'.format(os.path.abspath(self.runtime.proc_out)))
##             sh_file.write('#PBS -e {}\n'.format(os.path.abspath(self.runtime.proc_err)))
##             sh_file.write('#PBS -N {}\n'.format(os.path.basename(self.output[0])))
##             #sh_file.write('#PBS -N {}.{}_{}\n'.format(self.runtime.proc_err))
##             sh_file.write('#PBS -V\n')
##             # we try to reproduce the environment as much as possible becaus ehte
##             # script might be executed in a different environment
##             for k, v in os.environ.items():
##                 if any([k.startswith(x) for x in ('SSH', 'PBS', '_')]) or not k.replace('_', '').isalpha():
##                     continue
##                 sh_file.write('export {}="{}"\n'.format(k, v.replace('\n', '\\n')))
##             #
##             sh_file.write('\ncd {}\n'.format(os.path.abspath(os.getcwd())))
##             if self.working_dir is not None:
##                 sh_file.write('[ -d {0} ] || mkdir -p {0}\ncd {0}\n'.format(os.path.abspath(self.working_dir)))
## #
##             sh_file.write('''
## progress() {{
##   while true
##   do
##     touch {}
##     sleep 30
##   done
## }}
## 
## progress &
## MYSELF=$!
## '''.format(self.runtime.proc_prog))
## 
##             # interpreter
##             if '{}' in self.interpreter:
##                 sh_file.write(self.interpreter.replace('{}', pipes.quote(self.script_file)) +
##                     (self.args if isinstance(self.args, str) else ' '.join(pipes.quote(x) for x in self.args)) + '\n')
##             else:
##                 sh_file.write(self.interpreter + ' ' + pipes.quote(self.script_file) + \
##                     (self.args if isinstance(self.args, str) else ' '.join(pipes.quote(x) for x in self.args)) + '\n')
##             #
##             sh_file.write('\n\nCMD_RET=$?\nif [ $CMD_RET == 0 ]; then vtools admin --record_exe_info {} {}; fi\n'
##                 .format(os.getpid(), ' '.join(self.output)))
##             # a signal to show the successful completion of the job
##             sh_file.write('\nrm -f {}\nkill $MYSELF >/dev/null 2>&1\necho $CMD_RET > {}\n'
##                 .format(self.runtime.proc_prog, self.runtime.proc_done))
##         #
##         # try to submit command
##         if '{}' in self.submitter:
##             submit_cmd = self.submitter.replace('{}', self.runtime.proc_cmd)
##         else:
##             submit_cmd = self.submitter
##         #
##         if old_script is not None:
##             #with open(self.runtime.proc_cmd) as new_cmd:
##             #     if old_script == new_cmd.read():
##                      # if there is no change in command
##                      other_prog = glob.glob(os.path.abspath(self.output[0]) + '.working_*')
##                      if other_prog:
##                          for op in other_prog:
##                              if time.time() - os.path.getmtime(op) < 120:
##                                  env.logger.info('Another process appears to be working on {}, checking ...'.format(self.output[0]))
##                                  last_time = os.path.getmtime(op)
##                                  time.sleep(60)
##                                  # if the working file does not change after 60 seconds
##                                  if os.path.getmtime(op) != last_time:
##                                      raise RuntimeError('Failed to submit job because a job is currently running or has been failed within 2 minutes. Status file is {} (pid is {})'.format(op, os.getpid()))
##         env.logger.info('Running job {} with command "{}" from directory {}'.format(
##             self.runtime.proc_cmd, submit_cmd, os.getcwd()))
##         ret = subprocess.call(submit_cmd, shell=True,
##             stdout=open(self.runtime.proc_out, 'w'), stderr=open(self.runtime.proc_err, 'w'),
##             cwd=self.working_dir)
##         if ret != 0:
##             try:
##                 env.unlock(self.runtime.proc_out, str(os.getpid()))
##             except:
##                 pass
##         # 
##         if ret < 0:
##             raise RuntimeError("Failed to submit job {} due to signal {} (submitter='{}')" .format(self.runtime.proc_cmd, -ret, self.submitter))
##         elif ret > 0:
##             if os.path.isfile(self.runtime.proc_err):
##                  with open(self.runtime.proc_err) as err:
##                      msg = err.read()
##             else:
##                 msg = ''
##             raise RuntimeError("Failed to submit job {} using submiter '{}': {}".format(self.runtime.proc_cmd, self.submitter, msg))
##         else:
##             t = MonitorThread(target=self._monitor)
##             t.daemon = True
##             t.start()
##             if self.output[0] in self.pipeline.THREADS:
##                 raise RuntimeError('Two spawned jobs have the same self.output[0] file {}'.format(self.output[0]))
##             self.pipeline.THREADS[self.output[0]] = t
## 
## 
##     def _execute(self, ifiles, pipeline=None):
##         # substitute cmd by input_files and output_files
##         if pipeline.jobs > 1 and self.submitter is not None and not self.output:
##             env.logger.warning('Fail to execute in parallel because no output is specified.')
##         #
##         self.pipeline = pipeline
##         # Submit the job on a cluster system
##         # 1. if there is output, otherwise we cannot track the status of the job
##         # 2. if a submit command is specified
##         # 3. if --jobs with a value greater than 1 is used.
##         if self.output and pipeline.jobs > 1 and self.submitter is not None:
##             self._submit_command()
##             return False
##         else:
##             self._run_command()
##             return True
## 
## 
## class ExecuteRScript(ExecuteScript):
##     '''Execute in-line R script using Rscript as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='', output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='Rscript', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.R', wait=wait)
## 
## class ExecuteShellScript(ExecuteScript):
##     '''Execute in-line shell script using bash as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='', output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='bash', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.sh', wait=wait)
## 
## class ExecuteCShellScript(ExecuteScript):
##     '''Execute in-line shell script using bash as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='', output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='tcsh', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.csh', wait=wait)
## 
## class ExecutePythonScript(ExecuteScript):
##     '''Execute in-line python script using python as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='', output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='python', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.py', wait=wait)
## 
## class ExecutePython3Script(ExecuteScript):
##     '''Execute in-line python script using python3 as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='', output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='python3', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.py', wait=wait)
## 
## class ExecutePerlScript(ExecuteScript):
##     '''Execute in-line perl script using perl as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='',  output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='perl', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.perl', wait=wait)
## 
## 
## class ExecuteRubyScript(ExecuteScript):
##     '''Execute in-line perl script using perl as interpreter. Please
##     check action ExecuteScript for more details.
##     '''
##     def __init__(self, script='', args='',  output=[], export=None, working_dir=None, submitter=None, wait=True):
##         ExecuteScript.__init__(self, script=script, interpreter='ruby', args=args,
##             output=output, export=export, working_dir=working_dir, submitter=submitter,
##             suffix='.rb', wait=wait)
## 
## 
## class CheckRLibraries(ExecuteRScript):
##     '''Check the existence of specified R libraries. If an package is not available,
##     it will try to install it using "install.package" and from bioconductor. The pipeline
##     will raise an error if one of the library is not available and cannot be installed.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT 
## 
##     Raises:
##         A RuntimeError will be raised if a R library is not installed.
## 
##     Examples:
##         action=CheckRLibraries('edgeR')
##         action=CheckRLibraries(['edgeR', 'AIMS'])
##     '''
##     def __init__(self, libraries):
##         '''
##         Parameters:
##             libraries (string or list of strings):
##                 Name of one of more R libraries to be checked.
##         '''
##         PipelineAction.__init__(self)
##         if type(libraries) == type(''):
##             self.libraries= [libraries]
##         else:
##             self.libraries= libraries
##         # script
##         # get temp filename
##         self.output_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False).name
##         script = r'''
##         for (package in c({0})) {{
##             if (require(package, character.only=TRUE, quietly=TRUE)) {{
##                 write(paste(package, "AVAILABLE"), file="{1}", append=TRUE)
##                 next
##             }} else {{
##                 install.packages(package, repos="http://cran.us.r-project.org", 
##                     quiet=TRUE)
##             }}
##             # if the package still does not exist
##             if (!require(package, character.only=TRUE, quietly=TRUE)) {{
##                 source("http://bioconductor.org/biocLite.R")
##                 biocLite(package, ask=FALSE)
##             }}
##             # if it still does not exist, write the package name to output
##             if (require(package, character.only=TRUE, quietly=TRUE)) {{
##                 write(paste(package, "INSTALLED"), file="{1}", append=TRUE)
##             }} else {{
##                 write(paste(package, "MISSING"), file="{1}", append=TRUE)
##             }}
##         }}
##         '''.format(', '.join(['"{}"'.format(lib) for lib in self.libraries]), 
##             self.output_file)
##         ExecuteRScript.__init__(self, script=script, output=self.output_file)
## 
##     def __call__(self, ifiles, pipeline=None):
##         ExecuteRScript.__call__(self, ifiles, pipeline)
##         with open(self.output_file) as tmp:
##             count = 0
##             for line in tmp:
##                 lib, status = line.split()
##                 if status.strip() == "MISSING":
##                     env.logger.error('R Library {} is not available and cannot be installed.'.format(lib))
##                     count += 1
##                 elif status.strip() == 'AVAILABLE':
##                     env.logger.info('R library {} is available'.format(lib))
##                 elif status.strip() == 'INSTALLED':
##                     env.logger.info('R library {} has been installed'.format(lib))
##                 else:
##                     raise RuntimeError('This should not happen: {}'.format(line))
##         #
##         try:
##             os.remove(self.output_file)
##         except:
##             pass
##         if count > 0:
##             raise RuntimeError("One or more R libraries are not available.")
##         return ifiles
## 
## 
## class DecompressFiles(PipelineAction):
##     '''This action gets a list of input files from input file, decompressing
##     input files (.tar.gz, .zip, etc) if necessary. The decompressed files
##     are returned as output. One particular feature of this action is that
##     it records content of large tar or tar.gz files to a manifest file and
##     ignores the step if the manifest file exists.
## 
##     File Flow: Decompress input files
## 
##         INPUT ==> Decompress ==> OUTPUT
## 
##     Examples:
##         action=DecompressFiles()
##     '''
##     def __init__(self, dest_dir=None):
##         '''
##         Parameters:
##             dest_dir (None or string):
##                 Destination directory, default to current directory.
##         '''
##         self.dest_dir = dest_dir if dest_dir else '.'
##         PipelineAction.__init__(self)
## 
##     def _decompress(self, filename):
##         '''If the file ends in .tar.gz, .tar.bz2, .bz2, .gz, .tgz, .tbz2, decompress
##         it to dest_dir (current directory if unspecified), and return a list of files.
##         Uncompressed files will be returned untouched. If the destination files exist
##         and newer, this function will return immediately.'''
##         mode = None
##         if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tar.bz2'):
##             mode = 'r:gz'
##         elif filename.lower().endswith('.tbz2') or filename.lower().endswith('.tgz'):
##             mode = 'r:bz2'
##         elif filename.lower().endswith('.tar'):
##             mode = 'r'
##         elif filename.lower().endswith('.gz'):
##             dest_file = os.path.join(self.dest_dir, os.path.basename(filename)[:-3])
##             if existAndNewerThan(ofiles=dest_file, ifiles=filename):
##                 env.logger.info('Using existing decompressed file {}'.format(dest_file))
##             else:
##                 env.logger.info('Decompressing {} to {}'.format(filename, dest_file))
##                 with gzip.open(filename, 'rb') as gzinput, open(TEMP(dest_file), 'wb') as output:
##                     content = gzinput.read(10000000)
##                     while content:
##                         output.write(content)
##                         content = gzinput.read(10000000)
##                 # only rename the temporary file to the right one after finishing everything
##                 # this avoids corrupted files
##                 os.rename(TEMP(dest_file), dest_file)
##             return [dest_file]
##         elif filename.lower().endswith('.bz2'):
##             dest_file = os.path.join(self.dest_dir, os.path.basename(filename)[:-4])
##             if existAndNewerThan(ofiles=dest_file, ifiles=filename):
##                 env.logger.warning('Using existing decompressed file {}'.format(dest_file))
##             else:
##                 env.logger.info('Decompressing {} to {}'.format(filename, dest_file))
##                 with bz2.BZ2File(filename, 'rb') as bzinput, open(TEMP(dest_file), 'wb') as output:
##                     content = bzinput.read(10000000)
##                     while content:
##                         output.write(content)
##                         content = bzinput.read(10000000)
##                 # only rename the temporary file to the right one after finishing everything
##                 # this avoids corrupted files
##                 os.rename(TEMP(dest_file), dest_file)
##             return [dest_file]
##         elif filename.lower().endswith('.zip'):
##             bundle = zipfile.ZipFile(filename)
##             bundle.extractall(self.dest_dir)
##             env.logger.info('Decompressing {} to {}'.format(filename, self.dest_dir))
##             return [os.path.join(self.dest_dir, name) for name in bundle.namelist()]
##         #
##         # if it is a tar file
##         if mode is not None:
##             env.logger.info('Extracting fastq sequences from tar file {}'
##                 .format(filename))
##             #
##             # MOTE: open a compressed tar file can take a long time because it needs to scan
##             # the whole file to determine its content. I am therefore creating a manifest
##             # file for the tar file in the dest_dir, and avoid re-opening when the tar file
##             # is processed again.
##             manifest = RuntimeFiles(filename).manifest
##             all_extracted = False
##             dest_files = []
##             if existAndNewerThan(ofiles=manifest, ifiles=filename):
##                 all_extracted = True
##                 for f in [x.strip() for x in open(manifest).readlines()]:
##                     dest_file = os.path.join(self.dest_dir, os.path.basename(f))
##                     if existAndNewerThan(ofiles=dest_file, ifiles=filename):
##                         dest_files.append(dest_file)
##                         env.logger.info('Using existing extracted file {}'.format(dest_file))
##                     else:
##                         all_extracted = False
##             #
##             if all_extracted:
##                 return dest_files
##             #
##             # create a temporary directory to avoid corrupted file due to interrupted decompress
##             try:
##                 os.mkdir(os.path.join(self.dest_dir, 'tmp'))
##             except:
##                 # directory might already exist
##                 pass
##             #
##             dest_files = []
##             with tarfile.open(filename, mode) as tar:
##                 # only extract files
##                 files = [x.name for x in tar.getmembers() if x.isfile()]
##                 # save content to a manifest
##                 with open(manifest, 'w') as manifest:
##                     for f in files:
##                         manifest.write(f + '\n')
##                 for f in files:
##                     # if there is directory structure within tar file, decompress all to the current directory
##                     dest_file = os.path.join(self.dest_dir, os.path.basename(f))
##                     dest_files.append(dest_file)
##                     if existAndNewerThan(ofiles=dest_file, ifiles=filename):
##                         env.logger.info('Using existing extracted file {}'.format(dest_file))
##                     else:
##                         env.logger.info('Extracting {} to {}'.format(f, dest_file))
##                         tar.extract(f, os.path.join(self.dest_dir, 'tmp'))
##                         # move to the top directory with the right name only after the file has been properly extracted
##                         shutil.move(os.path.join(self.dest_dir, 'tmp', f), dest_file)
##                 # set dest_files to the same modification time. This is used to
##                 # mark the right time when the files are created and avoid the use
##                 # of archieved but should-not-be-used files that might be generated later
##                 [os.utime(x, None) for x in dest_files]
##             return dest_files
##         # return source file if nothing needs to be decompressed
##         return [filename]
##         
##     def __call__(self, ifiles, pipeline=None):
##         # decompress input files and return a list of output files
##         filenames = []
##         for filename in ifiles:
##             filenames.extend(self._decompress(filename))
##         filenames.sort()
##         return filenames
## 
## 
## class RemoveIntermediateFiles(PipelineAction):
##     '''This action removes specified files (not the step input files) and replaces
##     them with their signature (file size, md5 signature etc). A pipeline can bypass
##     completed steps with these files as input or output by checking the signatures.
##     In contrast, the steps would have to be re-run if the files are removed from the
##     file system. 
##     
##     File Flow: Input passthrough. Specified files are replaced by their signature.
## 
##         INPUT ====> INPUT
## 
##     Examples:
##         action=RemoveIntermediateFiles('${OUTPUT200}')
##         action=RemoveIntermediateFiles('${OUTPUT200} ${OUTPUT330}')
##         action=RemoveIntermediateFiles(['${OUTPUT200}', '${OUTPUT330}'])
##     '''
##     def __init__(self, files):
##         '''Replace ``files`` with their signatures. This pipeline passes its 
##         input to output and does not change the flow of pipeline.
##  
##         Parameters:
##             files (string or list of strings)
##                 One or more files to be removed. Multiple files can be specified
##                 in the same string if they are separated by spaces.
## 
##         '''
##         if isinstance(files, str):
##             self.files_to_remove = [files]
##         else:
##             self.files_to_remove = files
##         PipelineAction.__init__(self)
## 
##     def _getFiles(self):
##         for name in self.files_to_remove:
##             files = shlex.split(name)
##             for f in files:
##                 yield f
## 
##     def __call__(self, ifiles, pipeline=None):
##         env.logger.trace('Remove intermediate files {}'.format(' '.join(self.files_to_remove)))
##         for f in self._getFiles():
##             if not os.path.isfile(f):
##                 if os.path.isfile(f + '.file_info'):
##                     env.logger.info('Keeping existing {}.file_info.'.format(f))
##                 else:
##                     raise RuntimeError('Failed to create {}.file_info: Missing input file.'
##                         .format(f))
##             else:
##                 FileInfo(f).save()
##                 env.logger.info('Replace {0} with {0}.file_info'.format(f))
##                 try:
##                     os.remove(f)
##                 except e:
##                     env.logger.warning('Failed to remove intermediate file {}'.format(f))
##         return ifiles
## 
## 
## class LinkToDir(PipelineAction):
##     '''Create hard links of input files to a specified directory. This is 
##     usually used to link input files to a common cache directory so that 
##     all operations can be performed on that directory.
## 
##     File Flow: Link input files to specified destination directory.
## 
##         INPUT == LINK ==> DEST_DIR/INPUT
## 
##     Examples:
##         action=LinkToDir('cache')
## 
##     '''
##     def __init__(self, dest_dir):
##         '''
##         Parameters:
##             dest_dir (string):
##                 A directory to which input files will be linked to.
##                 The directory will be created if it does not exist.
##         '''
##         self.dest = dest_dir
##         if not os.path.isdir(self.dest):
##             env.logger.info('Creating directory {}'.format(self.dest))
##             try:
##                 os.makedirs(self.dest)
##             except Exception as e:
##                 raise RuntimeError('Failed to create directory {}: {}'.format(self.dest, e))
##             if not os.path.isdir(self.dest):
##                 raise RuntimeError('Failed to create directory {}: {}'.format(self.dest, e))
##         PipelineAction.__init__(self)
## 
##     def __call__(self, ifiles, pipeline=None):
##         ofiles = []
##         for filename in ifiles:
##             path, basename = os.path.split(filename)
##             if not os.path.isfile(filename):
##                 if os.path.isfile(filename + '.file_info'):
##                     dest_file = os.path.join(self.dest, basename) + '.file_info'
##                     if os.path.isfile(dest_file):
##                         if not os.path.samefile(filename + '.file_info', dest_file):
##                             os.remove(dest_file)
##                             env.logger.info('Linking {} to {}'.format(filename, self.dest))
##                             os.link(filename + '.file_info', os.path.join(self.dest, basename) + '.file_info')
##                         else:
##                             env.logger.trace('Reusing existing linked file_info file: {}'
##                                 .format(os.path.join(self.dest, basename) + '.file_info'))
##                     else:
##                         env.logger.info('Linking {} to {}'.format(filename, self.dest))
##                         os.link(filename + '.file_info', os.path.join(self.dest, basename) + '.file_info')
##                 else:
##                     raise RuntimeError('Failed to link {} to directory {}: file does not exist'
##                         .format(filename, self.dest))
##             else:
##                 dest_file = os.path.join(self.dest, basename)
##                 if os.path.isfile(dest_file):
##                     if not os.path.samefile(filename, dest_file):
##                         os.remove(dest_file)
##                         env.logger.info('Linking {} to {}'.format(filename, self.dest))
##                         os.link(filename, dest_file)
##                     else:
##                         env.logger.trace('Reusing existing linked file: {}'.format(dest_file))
##                 else:
##                     env.logger.info('Linking {} to {}'.format(filename, self.dest))
##                     os.link(filename, dest_file)
##             ofiles.append(os.path.join(self.dest, basename))
##         return ofiles
## 
## 
## class DownloadResource(PipelineAction):
##     '''Download resources to specified destination directory. dest_dir can
##     be a full path name or a directory relative to 
##     $local_resource/pipeline_resource where $local_resource is the local
##     resource directory of the project (default to ~/.variant_tools,
##     see runtime option local_resource for details). The default pipeline 
##     resource directory is $local_resource/pipeline_resource/NAME where NAME
##     is the name of the pipeline.
##     
##     File Flow: Input passthrough. 
## 
##         INPUT ====> INPUT
## 
##     Examples:
##         action=DownloadResource(resource='ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz',
##              dest_dir="${LOCAL_RESOURCE}/iGenomes")
## 
##         action=DownloadResource(resource='ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz
##             ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz.md5',
##             dest_dir='${LOCAL_RESOURCE/GATK')
##     
##     NOTE:
##         1. If FILE.md5 file is downloaded, it will be used to validate FILE.
##         2. The resources will be automatically decompressed if decompress=True (default). You would get both 
##             FILE and FILE.gz if you downloaded FILE.gz.
##     '''
##     def __init__(self, resource, dest_dir, output=[], decompress=True):
##         '''Download resources from specified URLs in ``resource``. 
## 
##         Parameters:
##             dest_dir:
##                 Directory where the downloaded resources will be placed.
##         '''
##         self.resource = [x for x in resource.split() if x]
##         if not dest_dir or type(dest_dir) != str:
##             raise ValueError('Invalid resource directory {}'.format(dest_dir))
##         else:
##             self.pipeline_resource = os.path.expanduser(dest_dir)
##         try:
##             if not os.path.isdir(self.pipeline_resource):
##                 os.makedirs(self.pipeline_resource)
##         except:
##             raise RuntimeError('Failed to create pipeline resource directory '
##                 .format(self.pipeline_resource))
##         self.decompress=decompress
##         PipelineAction.__init__(self, cmd='Download Resource {} to {}'.format(resource, dest_dir),
##             output=output)
## 
##     def __call__(self, ifiles, pipeline=None):
##         saved_dir = os.getcwd()
##         os.chdir(self.pipeline_resource)
##         ofiles, md5files = self._downloadFiles(ifiles)
##         self._validate(md5files)
##         os.chdir(saved_dir)
##         return ofiles
## 
##     def _validate(self, md5_files):
##         if md5_files:
##             prog = ProgressBar('Validating md5 signature', sum([x[1] for x in md5_files]))
##             mismatched_files = []
##             for filename, s in md5_files:
##                 try:
##                     downloaded_md5 = open(filename + '.md5').readline().split()[0]
##                     calculated_md5 = calculateMD5(filename, partial=False)
##                     if downloaded_md5 != calculated_md5:
##                         mismatched_files.append(filename)
##                 except Exception as e:
##                     env.logger.warning('Failed to verify md5 signature of {}: {}'
##                         .format(filename[:-4], e))
##                 prog.update(prog.count + s)
##             prog.done()
##             if mismatched_files:
##                 env.logger.warning('md5 signature of {} mismatch. '
##                       'Please remove {} and try again.'
##                       .format(', '.join(mismatched_files),
##                       'this file' if len(mismatched_files) == 1 else 'these files'))
## 
##     def _downloadFiles(self, ifiles):
##         '''Download resource'''
##         # decompress all .gz files
##         skipped = []
##         md5_files = []
##         for cnt, URL in enumerate(sorted(self.resource)):
##             filename = URL.rsplit('/', 1)[-1]
##             dest_file = os.path.join(self.pipeline_resource, filename)
##             try:
##                 if os.path.isfile(dest_file):
##                     skipped.append(filename)
##                 else:
##                     downloadURL(URL, dest_file, False,
##                         message='{}/{} {}'.format(cnt+1, len(self.resource), filename))
##             except KeyboardInterrupt as e:
##                 raise e
##             except Exception as e:
##                 raise RuntimeError('Failed to download {}: {} {}'
##                     .format(filename, type(e).__name__, e))
##             #
##             if filename.endswith('.tar.gz'):
##                 manifest_file = RuntimeFiles(filename).manifest
##                 env.logger.trace('Checking manifest {}'.format(manifest_file))
##                 decompress = not os.path.isfile(manifest_file)
##                 if not decompress:
##                     with open(manifest_file) as mf:
##                         for item in mf:
##                             if not os.path.isfile(item.strip()):
##                                 decompress = True
##                                 break
##                 if decompress:
##                     with tarfile.open(filename, 'r:gz') as tar: 
##                         s = delayedAction(env.logger.info, 'Extracting {}'.format(filename))
##                         tar.extractall(self.pipeline_resource)
##                         del s
##                         # only extract files
##                         files = [x.name for x in tar.getmembers() if x.isfile()]
##                         # save content to a manifest
##                         with open(manifest_file, 'w') as manifest:
##                             for f in files:
##                                 manifest.write(f + '\n')
##             elif filename.endswith('.gz'):
##                 if not existAndNewerThan(ofiles=filename[:-3], ifiles=filename):
##                     s = delayedAction(env.logger.info,
##                         'Decompressing {}'.format(filename))
##                     decompressGzFile(filename, inplace=False, force=True)
##                     del s
##             elif filename.endswith('.zip'):
##                 manifest_file = RuntimeFiles(filename).manifest
##                 env.logger.trace('Checking manifest {}'.format(manifest_file))
##                 decompress = not os.path.isfile(manifest_file)
##                 if not decompress:
##                     with open(manifest_file) as mf:
##                         for item in mf:
##                             if not os.path.isfile(item.strip()):
##                                 decompress = True
##                                 break
##                 if decompress:
##                     env.logger.trace('Decompressing {}'.format(filename))
##                     s = delayedAction(env.logger.info,
##                         'Decompressing {}'.format(filename))
##                     bundle = zipfile.ZipFile(filename)
##                     bundle.extractall(os.path.dirname(filename))
##                     with open(manifest_file, 'w') as manifest:
##                         for f in bundle.namelist():
##                             manifest.write(f + '\n')
##             #
##             if filename.endswith('.md5') and os.path.isfile(filename[:-4]):
##                 md5_files.append([filename[:-4], os.path.getsize(filename[:-4])])
##         if skipped:
##             env.logger.info('Using {} existing resource files under {}.'
##                 .format(', '.join(skipped), self.pipeline_resource))
##         return ifiles, md5_files
##  
## 
## class _CaseInsensitiveDict(MutableMapping):
##     """A case-insensitive ``dict``-like object that
##     1. limits the type of items to string or list of strings.
##     2. returns '' is the key does not exist (yield a warning)
##     3. allows attribute-like access
##     """
##     def __init__(self, data=None, **kwargs):
##         self._store = dict()
##         if data is None:
##             data = {}
##         self.update(data, **kwargs)
## 
##     def __setitem__(self, key, value):
##         # Use the uppercased key for lookups, but store the actual
##         # key alongside the value.
##         reset = 'reset' if key.upper() in self._store else 'set' # and value != self._store[key.upper()][1] else 'set'
##         if not isinstance(value, (str, list, tuple)):
##             value = str(value)
##             env.logger.warning('Pipeline variable {} is converted to "{}"'.format(key, value))
##         if isinstance(value, (list, tuple)) and not all([isinstance(x, str) for x in value]):
##             raise ValueError('Only string or list of strings are allowed for pipeline variables: {} for key {}'.format(value, key))
##         self._store[key.upper()] = (key, value)
##         if isinstance(value, str) or len(value) <= 2 or len(str(value)) < 50:
##             env.logger.debug('Pipeline variable ``{}`` is {} to ``{}``'.format(key, reset, str(value)))
##         else: # should be a list or tuple
##             val = str(value).split(' ')[0] + ' ...] ({} items)'.format(len(value))
##             env.logger.debug('Pipeline variable ``{}`` is {} to ``{}``'.format(key, reset, val))
## 
##     def __contains__(self, key):
##         return key.upper() in self._store
## 
##     def dict(self):
##         return {x:y for x,y in self._store.values()}
## 
##     def __setattr__(self, key, value):
##         if key == '_store':
##             self.__dict__[key] = value
##         else:
##             self.__setitem__(key, value)
## 
##     def __getattr__(self, key):
##         return self.__getitem__(key)
## 
##     def __getitem__(self, key):
##         try:
##             return self._store[key.upper()][1]
##         except:
##             env.logger.warning('Pipeline variable "{}" does not exist. A blank string is returned.'.format(key))
##             return ''
## 
##     def __delitem__(self, key):
##         del self._store[key.upper()]
## 
##     def __iter__(self):
##         return (casedkey for casedkey, mappedvalue in self._store.values())
## 
##     def __len__(self):
##         return len(self._store)
## 
##     def upper_items(self):
##         """Like iteritems(), but with all uppercase keys."""
##         return (
##             (upperkey, keyval[1])
##             for (upperkey, keyval)
##             in self._store.items()
##         )
## 
##     def __eq__(self, other):
##         if isinstance(other, collections.Mapping):
##             other = _CaseInsensitiveDict(other)
##         else:
##             return NotImplemented
##         # Compare insensitively
##         return dict(self.upper_items()) == dict(other.upper_items())
## 
##     # Copy is required
##     def copy(self):
##          return _CaseInsensitiveDict(self._store.values())
## 
##     def __repr__(self):
##         return '%s(%r)' % (self.__class__.__name__, dict(self.items()))
## 
## 
## class Pipeline:
##     '''The Variant Tools pipeline class. Its instance will be passed to each action
##     to provide runtime information. An action should not change any attribute of
##     the pipeline, except for setting additional variables through its ``VARS``
##     dictionary. Note that VARS is a case-insensitive dictionary but it is generally
##     recommended to use CAPTICAL names for pipeline variables. '''
##     def __init__(self, name, extra_args=[], pipeline_type='pipeline', verbosity=None, jobs=1):
##         self.pipeline = PipelineDescription(name, extra_args, pipeline_type)
##         self.spec_file = self.pipeline.spec_file
##         self.verbosity = verbosity
##         self.jobs = jobs
## 
##     def limit_steps(self, psteps, allowed_steps):
##         '''Restrict steps of a pipeline using allowed_steps'''
##         all_steps = {int(x.index):False for x in psteps}
##         #
##         for item in allowed_steps.split(','):
##             # remove space
##             item = ''.join([x for x in item if x != ' '])
##             if item.isdigit():
##                 # pipeline:100
##                 all_steps[int(item)] = True
##             elif '-' in item and item.count('-') == 1:
##                 l, u = item.split('-')
##                 if (l and not l.isdigit()) or (u and not u.isdigit()) or \
##                     (l and u and int(l) > int(u)):
##                     raise ValueError('Invalid pipeline step item {}'.format(item))
##                 # pipeline:-100, pipeline:100+ or pipeline:10-100
##                 if not l:
##                     l = min(all_steps.keys())
##                 if not u:
##                     u = max(all_steps.keys())
##                 #
##                 for key in all_steps.keys():
##                     if key >= int(l) and key <= int(u):
##                         all_steps[key] = True
##             else:
##                 raise ValueError('Invalid pipeline step item {}'.format(item))
##         # 
##         # disable limited steps
##         for idx in range(len(psteps)):
##             if not all_steps[int(psteps[idx].index)]:
##                 psteps[idx].options.append('skip')
##         env.logger.warning('Steps {} are skipped due to restriction {}'
##             .format(','.join([str(x) for x in all_steps.keys() if not all_steps[x]]), allowed_steps))
## 
##     def execute(self, pname, **kwargs):
##         allowed_steps = None
##         if not pname:
##             pname = ''
##         else:
##             if ':' in pname:
##                 pname, allowed_steps = pname.split(':', 1)
##         if not pname:
##             if len(self.pipeline.pipelines) == 1:
##                 pname = self.pipeline.pipelines.keys()[0]
##             elif 'default' in self.pipeline.pipelines:
##                 pname = 'default'
##             else:
##                 raise ValueError('Name of pipeline should be specified because '
##                     '{}.pipeline defines more than one pipelines without a default one. '
##                     'Available pipelines are: {}.'.format(self.pipeline.name,
##                     ', '.join(self.pipeline.pipelines.keys())))
##         elif pname not in self.pipeline.pipelines.keys():
##             raise ValueError('Pipeline {} is undefined in configuraiton file '
##                 '{}. Available pipelines are: {}'.format(pname,
##                 self.pipeline.name, ', '.join(self.pipeline.pipelines.keys())))
##         #
##         psteps = self.pipeline.pipelines[pname]
##         if allowed_steps is not None:
##             self.limit_steps(psteps, allowed_steps)
##         #
##         # the project will be opened when needed
##         with Project(mode=['ALLOW_NO_PROJ', 'READ_ONLY'], verbosity=self.verbosity) as proj:
##             self.VARS = _CaseInsensitiveDict(
##                 home=os.path.expanduser('~'),
##                 temp_dir=env.temp_dir,
##                 cache_dir=env.cache_dir,
##                 local_resource=env.local_resource,
##                 ref_genome_build=proj.build if proj.build is not None else '',
##                 pipeline_name=pname,
##                 spec_file=self.spec_file,
##                 model_name=pname,
##                 vtools_version=proj.version,
##                 working_dir=os.getcwd(),
##                 pipeline_format=self.pipeline.pipeline_format)
##         
##         if not os.path.isdir(env.cache_dir):
##             os.makedirs(env.cache_dir)
##         # these are command line options
##         if float(self.pipeline.pipeline_format) <= 1.0:
##             if 'cmd_input' in self.pipeline.commandline_opts:
##                 if not self.pipeline.commandline_opts['cmd_input']:
##                     self.pipeline.commandline_opts['cmd_input'] = []
##                 else:
##                     self.pipeline.commandline_opts['cmd_input'] = self.pipeline.commandline_opts['cmd_input'].split(',')
##             if 'cmd_output' in self.pipeline.commandline_opts:
##                 if not self.pipeline.commandline_opts['cmd_output']:
##                     self.pipeline.commandline_opts['cmd_output'] = []
##                 else:
##                     self.pipeline.commandline_opts['cmd_output'] = self.pipeline.commandline_opts['cmd_output'].split(',')
##         self.VARS.update(self.pipeline.commandline_opts)
##         self.VARS.update({k:str(v) for k,v in kwargs.items()})
##         if 'cmd_input' not in self.VARS:
##             self.VARS['cmd_input'] = []
##         if 'cmd_output' not in self.VARS:
##             self.VARS['cmd_output'] = []
##         # if there is a output file, write log to .log
##         if self.VARS['cmd_output'] and 'logfile' not in self.VARS:
##             self.VARS['logfile'] = self.VARS['cmd_output'][0] + '.log'
##         #
##         self.GLOBALS = {}
##         self.GLOBALS.update(globals())
##         self.THREADS = {}
##         # we need to put self.pipeline.pipeline_vars in self.VARS because
##         # they might refer to each other
##         self.VARS.update(self.pipeline.pipeline_vars)
##         for key, val in self.pipeline.pipeline_vars.items():
##             #if key in ('vtools_version', 'spec_file', 'home', 'pipeline_name', 'model_name'):
##             #    raise ValueError('Cannot reset read-only pipeline variable {}'.format(key))
##             self.VARS[key] = substituteVars(val, self.VARS, self.GLOBALS, asString=False)
##         for key, val in self.VARS.items():
##             if key == 'working_dir' and val != os.getcwd():
##                 env.logger.warning('Changing working directory to {}'.format(val))
##                 os.chdir(val)
##             if key == 'cache_dir' and val != env.cache_dir:
##                 env.logger.warning('Changing cache directory to {}'.format(val))
##                 env.cache_dir = val
##         #
##         if 'logfile' in self.VARS:
##             env.logger.info('Logging information is saved to {}'.format(self.VARS['logfile']))
##             if '/' in self.VARS['logfile']:
##                 d = os.path.split(self.VARS['logfile'])[0]
##                 if not os.path.isdir(d):
##                     env.logger.info('Making directory {} for output file'.format(d))
##                     os.makedirs(d)
##             ch = logging.FileHandler(self.VARS['logfile'].lstrip('>'), mode = 'a')
##             ch.setLevel(logging.DEBUG)
##             ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
##             env.logger.addHandler(ch)
##         #
##         ifiles = self.VARS['cmd_input']
##         step_index = 0
##         rewind_count = 0
##         while True:
##             # step_index can jump back and forth depending on the 
##             # execution status of each step
##             command = psteps[step_index]
##             if 'skip' in command.options:
##                 step_index += 1
##                 env.logger.info('Step {}.{}_{} is skipped'
##                     .format(self.pipeline.name, pname, command.index))
##                 if step_index == len(psteps):
##                     break
##                 step_output = []
##                 continue
##             self.VARS['pipeline_step'] = command.index
##             env.logger.info('Executing ``{}.{}_{}``: {}'
##                 .format(self.pipeline.name, pname, command.index, 
##                     ' '.join(command.comment.split())))
##             # init
##             for key, val in command.init_action_vars:
##                 self.VARS[key] = substituteVars(val, self.VARS, self.GLOBALS, asString=False)
##             # substitute ${} variables
##             emitter = None
##             if 'no_input' in command.options or 'independent' in command.options:
##                 step_input = []
##                 step_named_input = []
##             elif command.input is None or not command.input.strip():
##                 step_input = ifiles
##                 step_named_input = [['', ifiles]]
##             else:
##                 command_input_line = substituteVars(command.input, self.VARS, self.GLOBALS)
##                 if ':' in command_input_line:
##                     input_line, emitter_part = command_input_line.split(':', 1)
##                 else:
##                     input_line = command_input_line
##                     emitter_part = ''
##                 #
##                 if not input_line.strip():
##                     step_input = ifiles
##                     step_named_input = [['', ifiles]]
##                 else:
##                     # look for pattern of name=filenames
##                     pieces = re.split('([\w\d_]+\s*=)', input_line)
##                     step_named_input = []
##                     for piece in pieces:
##                         if not piece.strip():
##                             continue
##                         if piece.endswith('='):
##                             step_named_input.append([piece[:-1].strip(), []])
##                         else:
##                             expanded_files = sum([glob.glob(os.path.expanduser(x)) for x in shlex.split(piece)], [])
##                             if not expanded_files:
##                                 raise ValueError('{} does not expand to any valid file.'.format(piece))
##                             if not step_named_input:
##                                 step_named_input.append(['', expanded_files]) 
##                             else:
##                                 step_named_input[-1][1].extend(expanded_files)
##                     #
##                     step_input = sum([x[1] for x in step_named_input if x[0] == ''], [])
##                     if not step_input:
##                         step_input = ifiles
##                 if emitter_part:
##                     try:
##                         # remove ${INPUT} because it is determined by the emitter
##                         if 'input' in self.VARS:
##                             self.VARS.pop('input')
##                         emitter = eval('EmitInput({})'.format(
##                             substituteVars(emitter_part, self.VARS, self.GLOBALS)
##                             ), globals(), self.GLOBALS)
##                     except Exception as e:
##                         raise ValueError('Failed to interpret input emit options "{}"'.format(e))
##             #
##             #
##             self.VARS['input{}'.format(command.index)] = step_input
##             self.VARS['input'] = step_input
##             self.step_dependent_files = []
##             for n, f in step_named_input:
##                 if n:
##                     self.VARS['input{}_{}'.format(command.index, n)] = f
##                     self.VARS['input_{}'.format(n)] = f
##                     self.step_dependent_files.extend(f)
##             if self.step_dependent_files:
##                 env.logger.debug('Step dependent files are {}'.format(', '.join(self.step_dependent_files)))
##             #
##             saved_dir = os.getcwd()
##             for opt in command.options:
##                 matched = re.match('^input_alias\s*=\s*([\w\d_]+)$', opt)
##                 if matched:
##                     self.VARS[matched.group(1)] = step_input
##                     for n, f in step_named_input:
##                         if n:
##                             self.VARS['{}_{}'.format(matched.group(1), n)] = f
##                 matched = re.match('^working_dir\s*=\s*(\S+)$', opt)
##                 if matched:
##                     working_dir = os.path.expanduser(matched.group(1))
##                     if not os.path.isdir(working_dir):
##                         raise ValueError('Invalid working directory: {}'.format(working_dir))
##                     env.logger.info('Use working directory ``{}`` for {}_{}'.format(working_dir, pname, command.index))
##                     os.chdir(working_dir)
##             #
##             env.logger.trace('INPUT of step {}_{}: {}'
##                     .format(pname, command.index, step_input))
##             # 
##             # now, group input files
##             if not command.input_emitter:
##                 if emitter is None:  # if not defined in input:
##                     emitter = EmitInput()
##             else:
##                 if emitter is not None:
##                     raise ValueError('Cannot define input emitter in both input and input_emitter')
##                 try:
##                     # ${CMD_INPUT} etc can be used.
##                     emitter = eval(substituteVars(command.input_emitter, self.VARS, self.GLOBALS), globals(), self.GLOBALS)
##                 except Exception as e:
##                     raise RuntimeError('Failed to group input files: {}'
##                         .format(e))
##             # pass Pipeline itself to emitter
##             igroups, step_output = emitter(step_input, self)
##             try:
##                 for ig in igroups:
##                     if ig != self.VARS['input']:
##                         self.VARS['input'] = ig
##                     if not ig and float(self.pipeline.pipeline_format) <= 1.0:
##                         env.logger.trace('Step skipped due to no input file (for pipeline format < 1.0 only)')
##                         continue
##                     # pre action variables are evaluated for each ig because they might involve
##                     # changing ${input}
##                     for key, val in command.pre_action_vars:
##                         self.VARS[key] = substituteVars(val, self.VARS, self.GLOBALS, asString=False)
##                     action = substituteVars(command.action, self.VARS, self.GLOBALS)
##                     env.logger.trace('Emitted input of step {}_{}: {}'
##                         .format(pname, command.index, ig))
##                     env.logger.trace('Action of step {}_{}: {}'
##                         .format(pname, command.index, action))
##                     # check if the input file is ready. This is used for
##                     # parallel execution of the pipeline while the input file
##                     # might be worked on by another job
##                     for ifile in ig:
##                         # is ifile in any of the output files?
##                         if ifile in self.THREADS:
##                             # wait for the thread to complete
##                             env.logger.info('Waiting for the input file {} to be available.'
##                                 .format(ifile))
##                             #while self.THREADS[ifile].isAlive():
##                             ret = self.THREADS[ifile].join()
##                             # thread closed, remove from self.THREADS
##                             self.THREADS.pop(ifile)
##                             if ret:
##                                 raise RuntimeError('Failed to generate {}: {}'.format(ifile, ret))
##                         if not (os.path.isfile(ifile) or os.path.isfile(ifile + '.file_info')):
##                             #raise RewindExecution(ifile)
##                             raise RuntimeError('Non-existent input file {} due to ongoing or failed background job'.format(ifile))
##                     #
##                     if not action.strip():
##                         action = 'NullAction()'
##                     action = eval('(' + action + ')', globals(), self.GLOBALS)
##                     if isinstance(action, (tuple, list)):
##                         action = SequentialActions(action)
##                     if not issubclass(action.__class__, PipelineAction):
##                         env.logger.warning('Pipeline action {} is not a subclass of PipelineAction'.format(action.__class__))
##                     # pass the Pipeline object itself to action
##                     # this allows the action to have access to pipeline variables
##                     # and other options
##                     if 'blocking' in command.options:
##                         self.runtime = RuntimeFiles('{}_{}'.format(pname, command.index))
##                         with SharedProcess(self.runtime) as protection:
##                             ofiles = action(ig, self)
##                     else:
##                         ofiles = action(ig, self)
##                     if type(ofiles) == str:
##                         step_output.append(ofiles)
##                     else:
##                         step_output.extend(ofiles)
##                 # wait for all pending jobs to finish
##                 self.VARS['output{}'.format(command.index)] = step_output
##                 for opt in command.options:
##                     matched = re.match('^output_alias\s*=\s*([\w\d_]+)$', opt)
##                     if matched:
##                         env.logger.debug('Setting variable {} to {}'.format(matched.group(1), step_output))
##                         self.VARS[matched.group(1)] = step_output
##                 env.logger.trace('OUTPUT of step {}_{}: {}'
##                     .format(pname, command.index, step_output))
##                 for f in step_output:
##                     if not (os.path.isfile(f) or os.path.isfile(f + '.file_info') or f in self.THREADS):
##                         raise RuntimeError('Output file {} does not exist after '
##                             'completion of step {}_{} (working directory: {})'
##                             .format(f, pname, command.index, os.getcwd()))
##                 for key, val in command.post_action_vars:
##                     self.VARS[key] = substituteVars(val, self.VARS, self.GLOBALS, asString=False)
##                 #
##                 # In case of passthrough, the original input files will be passed to 
##                 # the next step regardless what has been produced during the step.
##                 if 'independent' not in command.options:
##                     ifiles = step_output
##                 # this step is successful, go to next
##                 os.chdir(saved_dir)
##                 step_index += 1
##                 env.logger.debug('Step {}.{}_{} is executed successfully.'
##                     .format(self.pipeline.name, pname, command.index))
##                 if step_index == len(psteps):
##                     break
##             except RewindExecution:
##                 rewind_count += 1
##                 if rewind_count >= 3:
##                     raise RuntimeError('Failed to execute pipeline {}.{}: excessive '
##                         'rewind during execution.'.format(self.pipeline.name, pname))
##                 # unfortunately, an input file has been removed (replaced by .file_info) but
##                 # a later steps need it. We will have to figure out how to create this 
##                 # file by looking backward ...
##                 to_be_regenerated = [x for x in step_input if not os.path.isfile(x)]
##                 # we need to check if this file is actually generated at all before
##                 # otherwise a misspecified input file would cause the whole pipline 
##                 # to start from step 1 again and again
##                 all_input_and_output_files = []
##                 for k,v in self.VARS.items():
##                     if ((k.startswith('INPUT') and k[5:].isdigit() and int(k[5:]) < int(command.index)) or \
##                         (k.startswith('OUTPUT') and k[6:].isdigit() and int(k[6:]) < int(command.index))) and \
##                          isinstance(v, (tuple, list)):
##                         all_input_and_output_files.extend(v)
##                 #
##                 for x in to_be_regenerated:
##                     if x not in all_input_and_output_files:
##                         raise RuntimeError('Specified input file "{}" does not exist and is not '
##                             'generated by any previous step.'.format(x))
##                     # remove all fony files so that they will be re-generated
##                     if os.path.isfile(x + '.file_info'):
##                         os.remove(x + '.file_info')
##                 remaining = [x for x in to_be_regenerated]
##                 env.logger.debug('Missing input file {}'.format(', '.join(remaining)))
##                 while step_index > 0:
##                     step_index -= 1
##                     command = psteps[step_index]
##                     # remove all fony files so that they will be re-generated
##                     for x in self.VARS['input{}'.format(command.index)]:
##                         if os.path.isfile(x + '.file_info'):
##                             env.logger.debug('Remove file info {}'.format(x + '.file_info'))
##                             os.remove(x + '.file_info')
##                     # if any of the input file does not exist, go back further
##                     if not all([os.path.isfile(x) for x in self.VARS['input{}'.format(command.index)]]):
##                         env.logger.debug('Not all input files are available: {}'.format(', '.join(self.VARS['input{}'.format(command.index)])))
##                         continue
##                     # check if a real file can be generated at this step
##                     remaining = [x for x in remaining if x not in self.VARS['output{}'.format(command.index)]]
##                     if not remaining:
##                         break
##                 if step_index > 1:
##                     ifiles = self.VARS['output{}'.format(psteps[step_index - 1].index)]
##                 else:
##                     ifiles = self.VARS['cmd_input']
##                 env.logger.warning('Rewinding to ``{}.{}_{}``: input files {} need to be re-generated.'
##                     .format(self.pipeline.name, pname, command.index, ', '.join(to_be_regenerated)))
##                 os.chdir(saved_dir)
##             except Exception as e:
##                 env.logger.debug('Failed to execute step {}.{}_{}.'
##                     .format(self.pipeline.name, pname, command.index))
##                 raise RuntimeError('Failed to execute step {}_{}: {}'
##                     .format(pname, command.index, e))
##             #
##             # clear variables that are local to step
##             for n, f in step_named_input:
##                 if n:
##                     self.VARS.pop('input_{}'.format(n))
##         #
##         # at the end of pipeline wait for all threads to complete
##         if self.THREADS:
##             for k, v in self.THREADS.items():
##                 env.logger.trace('Waiting for {} to be completed.'.format(k))
##                 while v.isAlive():
##                     v.join(5)
##                 # thread closed, remove from self.THREADS
##                 self.THREADS.pop(ifile)
##         env.logger.info('Execution of pipeline {}.{} is successful with output {}'
##             .format(self.pipeline.name, pname, ', '.join(step_output)))
## 
## 
## def executeArguments(parser):
##     parser.add_argument('specfile', metavar='SPECFILE',
##         help='''Name of a pipeline configuration file, which can be a
##             path to a .pipeline file (with or without extension) or one
##             of the online pipelines listed by command "vtools show pipelines".
##             For backward compatibility, if no input and output files are
##             specified (options --input and --output), values of this option 
##             is treated as a SQL query that will be executed against the project
##             database, with project genotype database attached as "genotype" and
##             annotation databases attached by their names.''')
##     parser.add_argument('pipelines', nargs='*', metavar='PIPELINES',
##         help='''Name of one or more pipelines defined in SPECFILE, which can be
##             ignored if the SPECFILE only defines one pipeline. One or more steps
##             can be specified in the form of 'pipeline:5' (step_5 only),
##             'pipeline:-5' (up to step 5), 'pipeline:5-' (from step 5),
##             'pipeline:2,5' (step 2 and 5), 'pipeline:2-5' (step 2 to 5). This
##             essentially adds an option "skip" to the unselected pipeline steps and
##             it is up to the user to ensure that the pipeline is executable
##             with only a subset of steps. Please use command "vtools show pipeline
##             SPECFILE" for details of available pipelines in SPECFILE, including
##             pipeline-specific parameters that could be used to change the default
##             behavior of the pipelines.''')
##     parser.add_argument('-j', '--jobs', default=1, type=int,
##         help='''Execute the pipeline in parallel model if a number other than
##             1 is specified. In this mode, the RunCommand action will create
##             a shell script and submit the job using a command specified by
##             option ``submitter``,  if this parameter is defined.''')
##     parser.add_argument('-d', '--delimiter', default='\t',
##         help='''Delimiter used to output results of a SQL query.''')
## 
## def execute(args):
##     # to keep backward compatibility, the vtools execute command
##     # can execute a SQL query and a pipeline
##     def executeQuery():
##         with Project(verbosity=args.verbosity) as proj:
##             # if there is no output, 
##             proj.db.attach('{}_genotype'.format(proj.name), 'genotype')
##             # for backward compatibility
##             proj.db.attach('{}_genotype'.format(proj.name))
##             cur = proj.db.cursor()
##             # 
##             query = args.specfile + ' ' + ' '.join(args.pipelines)
##             if query.upper().startswith('SELECT'):
##                 env.logger.trace('Analyze statement: "{}"'.format(query))
##                 cur.execute('EXPLAIN QUERY PLAN ' + query)
##                 for rec in cur:
##                     env.logger.trace('\t'.join([str(x) for x in rec]))
##             # really execute the query
##             try:
##                 cur.execute(query)
##             except Exception as e:
##                 raise RuntimeError('Failed to execute SQL query "{}": {}'
##                     .format(query, e))
##             proj.db.commit()
##             sep = args.delimiter
##             for rec in cur:
##                 print(sep.join(['{}'.format(x) for x in rec]))
##     #
##     def executePipeline():                
##         pipeline = Pipeline(args.specfile, extra_args=args.unknown_args, verbosity=args.verbosity,
##             jobs=args.jobs)
##         # unspecified
##         if not args.pipelines:
##             pipeline.execute(None)
##         else:
##             for name in args.pipelines:
##                 pipeline.execute(name)
##     # 
##     try:
##         env.verbosity = args.verbosity
##         env.logger = None
##         # definitely a pipeline
##         if args.specfile.endswith('.pipeline') or args.unknown_args:
##             executePipeline()
##         # definitely a sql query
##         elif args.delimiter != '\t':
##             executeQuery()
##         else:
##             try:
##                 # try to execute as a SQL query
##                 executeQuery()
##             except (RuntimeError, ValueError) as e:
##                 env.logger.debug('Failed to execute {} as SQL query: {}'
##                     .format(' '.join(args.pipelines), e))
##                 executePipeline()
##     except Exception as e:
##         env.unlock_all()
##         env.logger.error(e)
##         sys.exit(1)
## 
## if __name__ == '__main__':
##     pass
