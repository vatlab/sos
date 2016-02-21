#!/usr/bin/env python
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

import os
import sys
import glob
import logging
import getpass
import textwrap
import tempfile
import shutil
import argparse
import threading
import Queue
import signal
import time
import re
import tarfile
import urllib
import inspect
import pydoc

from .utils import *


# define a field type
Field = namedtuple('Field', ['name', 'index', 'adj', 'fmt', 'type', 'comment'])
Column = namedtuple('Column', ['index', 'field', 'adj', 'comment'])
#
# see http://varianttools.sourceforge.net/Calling/New for details
#
PipelineCommand = namedtuple('PipelineCommand', ['index', 'options', 'input',
    'input_emitter', 'action', 'init_action_vars', 'pre_action_vars', 'post_action_vars', 'comment'])
#
# How field will be use in a query. For example, for field sift, it is
# connection clause will be:
#   field = dbNSFP.sift
#   table = dbNSFP.dbNSFP  # annotation
#   link = chr=dbNSFP.chr AND pos=dbNSFP.h18pos
#
FieldConnection = namedtuple('FieldConnection', ['field', 'table', 'link'])

class MyConfigParser(RawConfigParser):
    def __init__(self, *args, **kwargs):
        RawConfigParser.__init__(self, *args, **kwargs)

    def optionxform(self, x):
        return str(x)

    def items(self, section, raw=False, vars={}):
        res = RawConfigParser.items(self, section)
        return res

    def get(self, section, item, raw=False, vars={}):
        res = RawConfigParser.get(self, section, item)
        if re.search(r'%\(\w+\)s', res):
            env.logger.warning('The use of %(VAR)s variable is deprecated. Please use ${{}} instead: {} ..'
                .format(' '.join(res.split())[:40]))
            new_res = re.sub(r'%\((\w+)\)s', r'${\1}', res)  
            env.logger.debug('Replacing "{}" with "{}"'.format(res, new_res))
            return new_res
        return res

class PipelineDescription:
    def __init__(self, name, extra_args=[], pipeline_type='pipeline'):
        '''Pipeline configuration file'''
        self.description = None
        self.pipeline_format = '1.0'
        self.pipeline_vars = {}
        self.pipelines = {}
        self.pipeline_descriptions = {}
        self.pipeline_type = pipeline_type
        self.commandline_opts = {}
        # a strange piece of text and replaces newlines in triple
        # quoted text so that the text can be read as a single string.
        # Newlines in the processed values will be translated back.
        self.newline_PH = chr(5) + chr(6) + chr(7)
        self.semicolon_PH = chr(16) + chr(17)
        #
        if os.path.isfile(name + '.pipeline'):
            self.name = os.path.split(name)[-1]
            self.commandline_opts = self.parseArgs(name + '.pipeline', extra_args)
            self.parsePipeline(name + '.pipeline', defaults=self.commandline_opts) 
        elif name.endswith('.pipeline') and os.path.isfile(name):
            self.name = os.path.split(name)[-1][:-9]
            self.commandline_opts = self.parseArgs(name, extra_args)
            self.parsePipeline(name, defaults=self.commandline_opts) 
        else:
            # not found, try online
            if name.endswith('.pipeline'):
                url = '{}/{}'.format(pipeline_type, name)
            else:
                url = '{}/{}.pipeline'.format(pipeline_type, name)
            try:
                pipeline = downloadFile(url, quiet=True)
            except Exception as e:
                raise ValueError('Failed to download pipeline specification '
                    'file {}.pipeline: {}'.format(name, e))
            self.name = name
            self.commandline_opts = self.parseArgs(pipeline, extra_args)
            self.parsePipeline(pipeline, defaults=self.commandline_opts)

    def _translateConfigText(self, filename):
        # We would like to keep everything between triple quotes literal. There is no easy way 
        # to achive that so we have to convert the text and then convert back after
        # the data is read.
        if not hasattr(self, 'config_text'):
            with open(filename, 'r') as inputfile:
                self.config_text = inputfile.read().replace(';', self.semicolon_PH)
            #
            for quote in ('"""', "'''"):
                pieces = re.split(quote, self.config_text)
                for i in range(1, len(pieces), 2):
                    # if ''' starts from a new line (this should not happen, automatically pad it
                    if pieces[i-1].endswith('\n'):
                        pieces[i-1] = pieces[i-1] + ' '
                    # automatically add r to ''' ''' quotes
                    if quote == "'''" and not pieces[i-1].endswith('r'):
                        pieces[i-1] = pieces[i-1] + 'r'
                    # replace string with an unlikely character
                    pieces[i] = pieces[i].replace('\\\n', '').replace('\n', self.newline_PH)
                self.config_text = quote.join(pieces)
            # handling comments
            #
            # This will allow the use of 
            #
            # [header]
            # # comment_text
            #
            # instead of
            #
            # [header]
            # comment: comment_text
            pieces = re.split('(\n\[.+\]\s*)', self.config_text)
            for idx, piece in enumerate(pieces):
                if re.match('^\s*\[[*\w\d,\s:=-]+\]\s*$', piece) and idx + 1 != len(pieces):
                    comment = []
                    non_comment = []
                    has_comment = False
                    for line in pieces[idx + 1].split('\n'):
                        if line.startswith('#') and not non_comment:
                            comment.append(line.lstrip('#'))
                        else:
                            if line.startswith('comment:') or line.startswith('comment='):
                                # if there is existing ...
                                has_comment = True
                                break
                            non_comment.append(line)
                    if not has_comment and comment:
                        pieces[idx + 1] = '\n'.join(non_comment) + '\n' + 'comment=' + ' '.join(comment) + '\n'
            self.config_text = '\n'.join(pieces)
            # now, the [DEFAULT] section
            #
            # we can convert
            #
            # par=default
            #     comment 
            #
            # automatically to
            #
            # par=default
            # par_comment=comment
            #
            #
            # for the regular section
            #
            # if we have
            #
            # [section]
            # input:
            # somethingelse(
            #
            # we change it to
            #
            # [section]
            # input:
            # action: somethingelse(
            #
            # automatically.
            #
            # 
            has_pipeline_description_section = False
            has_pipeline_description = False
            pieces = []
            par = None
            in_default = False
            in_comment = False
            comment_count = 0
            second_comment = []
            for line in self.config_text.split('\n'):
                pieces.append(line)
                #
                if re.match('^\[\s*pipeline desciption\]\s*$', line):
                    has_pipeline_description_section = True
                if re.match('^description\s*[=\:]', line):
                    has_pipeline_description = True
                if line.startswith('#'):
                    if not in_comment:
                        in_comment = True
                        comment_count += 1
                else:
                    in_comment = False
                #
                if not in_comment and not re.match('^\s+', line) and not has_pipeline_description_section:
                    pieces.insert(len(pieces)-1, '[pipeline description]')
                    has_pipeline_description_section = True
                #
                if in_comment and comment_count == 2:
                    second_comment.append('    ' + line.lstrip('#'))
                #
                if re.match('^\[\s*DEFAULT\s*\]\s*$', line):
                    in_default = True
                    continue
                if not in_default:
                    # not in default section, we expand action
                    # automatically
                    if (not line.startswith('#')) and (re.match('^[\w\d_]+\s*\(', line) or re.match('^\${', line)):
                        pieces[-1] = 'action: ' + line
                    continue
                elif re.match('^\[', line):
                    in_default = False
                    continue
                #
                matched = re.match('^(\w+[\w\d_]*)\s*[=:]', line)
                if matched:
                    par = matched.group(1)
                elif par is not None:
                    # if not matched, but par is True, must be the next line
                    if line.startswith(' ') or line.startswith('\t'):
                        pieces[-1] = '{}_comment : {}'.format(par, line.lstrip())
                        par = None
            self.config_text = '\n'.join(pieces)
            if not has_pipeline_description:
                self.config_text = self.config_text.replace('[pipeline description]', '[pipeline description]\ndescription:\n{}'.format('\n'.join(second_comment)))
            #with open(os.path.join(env.temp_dir, 'pipeline_executed.tmp'), 'w') as tmp:
            #    tmp.write(self.config_text)
        return self.config_text
        
    def parseArgs(self, filename, fmt_args):
        with open(filename) as pp:
            for line in pp:
                if not line.startswith('#'):
                    break
                if line.startswith('##fileformat='):
                    m = re.match('##fileformat=\D*([\d.]+)', line)
                    if m is None:
                        raise ValueError('Pipeline format string should have format ##fileformat=PIPELINEx.xx: {} detected'
                            .format(line))
                    self.pipeline_format = m.group(1)
        #
        env.logger.debug('Pipeline version {}'.format(self.pipeline_format))
        # We used format interpolation in older version 
        try:
            if float(self.pipeline_format) <= 1.0:
                if sys.version_info.major == 2:
                    fmt_parser = SafeConfigParser()
                else:
                    fmt_parser = ConfigParser(strict=False)
                fmt_parser.read(filename)
            else:
                # and now we only use pipeline variables.
                fmt_parser = MyConfigParser()
                fmt_parser.readfp(StringIO(self._translateConfigText(filename)))
        except Exception as e:
            msg = repr(e).split('\n')
            if msg[-1].strip().startswith('[line'):
                line_no = int(msg[-1].strip()[6:].split(']')[0])
                lines = self._translateConfigText(filename).split('\n')
                if line_no > 2:
                    env.logger.error('{}: {}'.format(line_no-1, lines[line_no - 2]))
                env.logger.error('{}: {}'.format(line_no, lines[line_no-1]))
                if line_no < len(lines):
                    env.logger.error('{}: {}'.format(line_no + 1, lines[line_no]))
            raise
        parameters = fmt_parser.items('DEFAULT')
        parser = argparse.ArgumentParser(prog='vtools CMD --pipeline {}'
            .format(os.path.split(filename)[-1]),
            description='Parameters to override parameters of existing steps.')
        self.parameters = []
        if 'input' not in [x[0] for x in parameters]:
            parser.add_argument('-i', '--input', help='Input of pipeline as variable ${cmd_input}', nargs='*', default=[])
        else:
            par_help = [x[1] for x in parameters if x[0] == 'input_comment']
            parser.add_argument('-i', '--input', help=par_help[0] if par_help else '', nargs='*', default=[])
        if 'output' not in [x[0] for x in parameters]:
            parser.add_argument('-o', '--output', help='Output of pipeline as variable ${cmd_ontput}', nargs='*', default=[])
        else:
            par_help = [x[1] for x in parameters if x[0] == 'output_comment']
            parser.add_argument('-o', '--output', help=par_help[0] if par_help else '', nargs='*', default=[])
        for par in parameters:
            # $NAME_comment is used for documentation only
            if par[0].endswith('_comment') or par[0] in ('input', 'output'):
                continue
            if par[0].lower() in ('home', 'cwd', 'cmd_input', 'cmd_output', 'temp_dir', 'cache_dir', 'local_resource',
                    'ref_genome_build', 'pipeline_name', 'spec_file', 'model_name', 'vtools_version', 'pipeline_format'):
                raise ValueError('Command option {} is reserved and cannot be specified from command line.'.format(par[0]))
            par_help = [x[1] for x in parameters if x[0] == par[0] + '_comment']
            self.parameters.append((par[0], par[1], par_help[0] if par_help else ''))
            parser.add_argument('--{}'.format(par[0]), help=self.parameters[-1][2],
                nargs='*', default=par[1])
        args = vars(parser.parse_args(fmt_args))
        if  float(self.pipeline_format) <= 1.0:
            for key,value in args.items():
                if not isinstance(value, str):
                    args[key] = ','.join(value)
        if 'input' in args:
            args['cmd_input'] = args['input']
            args.pop('input')
        if 'output' in args:
            args['cmd_output'] = args['output']
            args.pop('output')
        return args

    def parsePipeline(self, filename, defaults):
        self.spec_file = filename
        if float(self.pipeline_format) <= 1.0:
            if sys.version_info.major == 2:
                parser = SafeConfigParser()
            else:
                parser = ConfigParser(strict=False)
            parser.optionxform = str
        else:
            # and now we only use pipeline variables.
            parser = MyConfigParser()
        # this allows python3 to read .pipeline file with non-ascii characters,
        # but there is no simple way to make it python2 compatible.
        #with open(filename, 'r', encoding='UTF-8') as inputfile:
        #    parser.readfp(inputfile)
        try:
            parser.readfp(StringIO(self._translateConfigText(filename)))
        except:
            env.logger.error(self._translateConfigText(filename))
            raise
        # sections?
        sections = parser.sections()
        if 'pipeline description' not in sections:
            raise ValueError("Missing section 'pipeline description' in "
                "configuration file {}".format(filename))
        #
        for section in sections:
            if section.lower() == 'pipeline description':
                for item in parser.items(section, vars=defaults):
                    if item[0] == 'description':
                        self.description = item[1].strip()
                        if (self.description.startswith("r'''") or self.description.startswith("'''")) and self.description.endswith("'''"):
                            self.description = self.description[(4 if self.description.startswith('r') else 3):-3].replace(self.newline_PH, '<br>')
                        elif (self.description.startswith('r"""') or self.description.startswith('"""')) and self.description.endswith('"""'):
                            self.description = self.description[(4 if self.description.startswith('r') else 3):-3].replace(self.newline_PH, '\n')

                    elif item[0].endswith('_description'):
                        self.pipeline_descriptions[item[0].strip().rsplit('_', 1)[0]] = item[1]
                    elif item[0] in defaults or item[0].endswith('_comment'):
                        pass
                    else:
                        self.pipeline_vars[item[0]] = item[1]
            else:
                #
                # section header can contain multiple steps
                # [A_1,B_1,*_3]
                try:
                    section_headers = [x.strip() for x in section.split(':', 1)[0].split(',')]
                    for header in section_headers:
                        if not re.match('^([\w*_][\w\d*_]*_)?[\d]+$', header) and not re.match('^[\w][\w\d]*$', header):
                            raise ValueError('Invalid section header "{}"'.format(section))
                    #
                    pnames = [x.strip().rsplit('_', 1)[0] if '_' in x and x.rsplit('_',1)[-1].isdigit() else ('default' if x.isdigit() else x) for x in section_headers]
                    pidxs = [x.strip().rsplit('_', 1)[1] if '_' in x and x.rsplit('_',1)[-1].isdigit() else (x if x.isdigit() else '0') for x in section_headers]
                    #
                    if ':' in section:
                        options = [x.strip() for x in section.split(':', 1)[-1].split(',')]
                        for opt in options:
                            if opt not in ['no_input', 'independent', 'skip', 'blocking'] and not re.match('^(output_alias|input_alias|action)\s*=\s*([\w\d_]+)$', opt) \
                                and not re.match('^working_dir\s*=\s*(\S+)$', opt):
                                env.logger.warning('Unrecognized section option: {}'.format(opt))
                    else:
                        options = []
                except Exception as e:
                    raise ValueError('Invalid section name {} in pipeline description file {}: {}'
                        .format(section, filename, e))
                if not all([x.isdigit() for x in pidxs]):
                    raise ValueError('Index of a pipeline step should be an integer: {} provided'
                        .format(', '.join(pidxs)))
                for pname in pnames:
                    if pname not in self.pipelines:
                        self.pipelines[pname] = []
                try:
                    items = [x[0] for x in parser.items(section, raw=True)]
                    #if 'action' not in items:
                    #    raise ValueError('Missing item "action" in section {}.'.format(section))
                    has_input = False
                    step_init_vars = []
                    step_pre_vars = []
                    step_post_vars = []
                    before_input_action = True
                    before_action = True
                    for item in items:
                        if item.endswith('_comment'):
                            continue
                        if item not in ['input_emitter', 'comment'] + defaults.keys():
                            #env.logger.warning('ITEM {}'.format(item))
                            if item == 'input':
                                before_input_action = False
                                has_input = True
                                continue
                            elif item == 'action':
                                before_action = False
                                before_input_action = False
                                continue
                            if before_input_action:
                                step_init_vars.append([item, parser.get(section, item, vars=defaults)])
                            elif before_action:
                                step_pre_vars.append([item, parser.get(section, item, vars=defaults)])
                            else:
                                step_post_vars.append([item, parser.get(section, item, vars=defaults)])
                    #env.logger.warning('INIT VAR {}'.format(step_init_vars))
                    #env.logger.warning('PRE ACTION VAR {}'.format(step_pre_vars))
                    #env.logger.warning('POST ACTION VAR {}'.format(step_post_vars))
                    # if no input, assume post_input, pre-action
                    if not has_input:
                        step_pre_vars = step_init_vars
                        step_init_vars = []
                    for pname,pidx in zip(pnames, pidxs):
                        command = PipelineCommand(index=pidx,
                            options=options,
                            input=parser.get(section, 'input', vars=defaults).replace(self.newline_PH, '\n').replace(self.semicolon_PH, ';') if 'input' in items else None,
                            input_emitter=parser.get(section, 'input_emitter', vars=defaults).replace(self.newline_PH, '\n').replace(self.semicolon_PH, ';') if 'input_emitter' in items else '',
                            action=parser.get(section, 'action', vars=defaults).replace(self.newline_PH, '\n').replace(self.semicolon_PH, ';') if 'action' in items else '',
                            init_action_vars=step_init_vars,
                            pre_action_vars=step_pre_vars,
                            post_action_vars=step_post_vars,
                            comment=parser.get(section, 'comment', raw=True).replace(self.newline_PH, '\n').replace(self.semicolon_PH, ';') if 'comment' in items else '')
                        self.pipelines[pname].append(command)
                except Exception as e:
                    raise ValueError('Invalid section {}: {}'.format(section, e))
        # for pipelines with all * sections, look for a description or use default name
        not_wildname = [y for y in self.pipelines.keys() if '*' not in y and '?' not in y]
        # if all names are * ...
        if not not_wildname:
            if self.pipeline_descriptions:
                self.pipelines.update({x:[] for x in self.pipeline_descriptions})
            else:
                self.pipelines.update({'default':[]})
        # process wild cast pipelines
        for wildname in [x for x in self.pipelines.keys() if '*' in x or '?' in x]:
            for pname in [y for y in self.pipelines.keys() if '*' not in y and '?' not in y]:
                if matchName(wildname, pname):
                    self.pipelines[pname].extend(self.pipelines[wildname])
        #
        self.pipelines = {x:y for x,y in self.pipelines.items() if '*' not in x and '?' not in x}
        # sort steps
        for pname in self.pipelines:
            self.pipelines[pname].sort(key=lambda x: int(x[0].strip().rsplit('_')[-1]))
        # 
        # validate
        for pname in self.pipeline_descriptions:
            if pname not in self.pipelines.keys():
                env.logger.warning('Invalid item {0}_description because pipeline '
                    '"{0}" is not defined in this file (available: {1}).'
                    .format(pname, ', '.join(self.pipelines.keys())))
        for pname, pipeline in self.pipelines.items():
            if pname not in self.pipeline_descriptions:
                #if pname != 'default':
                #    env.logger.warning('No description for {} {} is available.'.format(self.pipeline_type, pname))
                self.pipeline_descriptions[pname] = ''
            for idx, cmd in enumerate(pipeline):
                if cmd is None:
                    raise ValueError('Invalid pipeline {}. Step {} is left unspecified.'
                        .format(pname, idx+1))
                for opt in cmd.options:
                    matched = re.match('^action\s*=\s*([\w\d_]+)$', opt)
                    if matched:
                        header = matched.group(1)
                        if not re.match('^([\w*_][\w\d*_]*_)?[\d]+$', header) and not re.match('^[\w][\w\d]*$', header):
                            raise ValueError('Invalid section header for option {}'.format(opt))
                        #
                        pn = header.strip().rsplit('_', 1)[0] if '_' in header else ('default' if header.isdigit() else header)
                        pi = header.strip().rsplit('_', 1)[1] if '_' in header else (header if header.isdigit() else '0')
                        if pn not in self.pipelines:
                            raise ValueError('Cannot find pipeline {} for option {}'.format(pn, opt))
                        found = False
                        for step in self.pipelines[pn]:
                            if step.index == pi:
                                if cmd.action.strip() != '':
                                    raise ValueError('No action should be specified if option action is used for step {}_{}'.format(pname, cmd.index))
                                env.logger.info('Using action for step [[{}_{}]] for step [[{}_{}]]'.format( pn, pi, pname, cmd.index))
                                # have to re-create the whole object
                                pipeline[idx] = PipelineCommand(
                                    index=cmd.index,
                                    options=cmd.options,
                                    input=cmd.input,
                                    input_emitter=cmd.input_emitter,
                                    action=step.action,
                                    init_action_vars=cmd.init_action_vars,
                                    pre_action_vars=cmd.pre_action_vars,
                                    post_action_vars=cmd.post_action_vars,
                                    comment=cmd.comment, 
                                )
                                cmd = pipeline[idx]
                                found = True
                                break
                        if not found:
                            raise ValueError('Cannot find step {} for option {}'.format(pi, opt))
                if not cmd.action:
                    raise ValueError('Missing or empty action for step {} of pipeline {}'
                        .format(cmd.index, pname))
                # step.comment might have expression with pipeline_name and pipeline_step
                if '${' in cmd.comment:
                    pipeline[idx] = PipelineCommand(
                        index=cmd.index,
                        options=cmd.options,
                        input=cmd.input,
                        input_emitter=cmd.input_emitter,
                        action=cmd.action,
                        init_action_vars=cmd.init_action_vars,
                        pre_action_vars=cmd.pre_action_vars,
                        post_action_vars=cmd.post_action_vars,
                        comment=substituteVars(cmd.comment, 
                            {'pipeline_name': pname, 
                             'pipeline_step': cmd.index,
                             'pipeline_format': self.pipeline_format},
                            {})
                        )
     
    def describe(self):
        textWidth = max(60, getTermWidth())
        if self.description is not None:
            # separate \n\n 
            for paragraph in dehtml(self.description).split('\n\n'):
                print('\n'.join(textwrap.wrap(paragraph, width=textWidth)))
        #
        text = 'Available {}: {}'.format(
            'simulation models' if self.pipeline_type == 'simulation' else 'pipelines',
            ', '.join(sorted(self.pipelines.keys())))
        print('\n' + '\n'.join(textwrap.wrap(text, width=textWidth, subsequent_indent=' '*8)))
        for pname, pipeline in sorted(self.pipelines.items()):
            paragraphs = dehtml(self.pipeline_descriptions[pname]).split('\n\n')
            print('\n' + '\n'.join(textwrap.wrap('{} "{}":  {}'
                .format(
                'Model' if self.pipeline_type == 'simulation' else 'Pipeline',
                pname, paragraphs[0]), width=textWidth)))
            for paragraph in paragraphs[1:]:
                print('\n'.join(textwrap.wrap(paragraph, width=textWidth)))
            for idx, step in enumerate(pipeline):
                # hide a step if there is no comment
                if step.comment:
                    text = '{:<22}'.format('  {}_{}:'.format(pname, step.index)) + step.comment
                    print('\n'.join(textwrap.wrap(text, width=textWidth, subsequent_indent=' '*22)))
        #
        if self.parameters:
            print('\n{} parameters:'.format('Model' if self.pipeline_type == 'simulation' else 'Pipeline'))
            for item in self.parameters:
                #
                text = '  ' + item[0] + \
                    (' '*(22-len(item[0])-2) if len(item[0])<20 else ' ') + \
                    (item[2] + ' ' if item[2] else '') + \
                    ('(default: {})'.format(item[1]) if item[1] else '')
                print('\n'.join(textwrap.wrap(text, subsequent_indent=' '*22,
                    width=textWidth)))
        #else:
        #    print('\nNo configurable parameter is defined for this {}.\n'
        #        .format('model' if self.pipeline_type == 'simulation' else 'pipeline'))


if __name__ == '__main__':
    # for testing purposes only. The main interface is provided in vtools
    pass
