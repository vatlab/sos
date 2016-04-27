#!/usr/bin/env python3
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

import os
import sys
import yaml
import atexit
import fnmatch
import webbrowser

from pygments import highlight
from pygments.lexers import PythonLexer, get_lexer_by_name
from pygments.formatters import HtmlFormatter

from .utils import env, get_traceback
from .sos_script import SoS_Script, SoS_Workflow
from .sos_executor import Sequential_Executor
from io import StringIO

#
# subcommmand show
#
class ContinuousHtmlFormatter(HtmlFormatter):
    def _wrap_tablelinenos(self, inner):
        dummyoutfile = StringIO()
        lncount = 0
        for t, line in inner:
            if t:
                lncount += 1
            dummyoutfile.write(line)

        fl = self.linenostart
        mw = len(str(lncount + fl - 1))
        sp = self.linenospecial
        st = self.linenostep
        la = self.lineanchors
        aln = self.anchorlinenos
        nocls = self.noclasses
        if sp:
            lines = []

            for i in range(fl, fl+lncount):
                if i % st == 0:
                    if i % sp == 0:
                        if aln:
                            lines.append('<a href="#%s-%d" class="special">%*d</a>' %
                                         (la, i, mw, i))
                        else:
                            lines.append('<span class="special">%*d</span>' % (mw, i))
                    else:
                        if aln:
                            lines.append('<a href="#%s-%d">%*d</a>' % (la, i, mw, i))
                        else:
                            lines.append('%*d' % (mw, i))
                else:
                    lines.append('')
            ls = '\n'.join(lines)
        else:
            lines = []
            for i in range(fl, fl+lncount):
                if i % st == 0:
                    if aln:
                        lines.append('<a href="#%s-%d">%*d</a>' % (la, i, mw, i))
                    else:
                        lines.append('%*d' % (mw, i))
                else:
                    lines.append('')
            ls = '\n'.join(lines)

        # in case you wonder about the seemingly redundant <div> here: since the
        # content in the other cell also is wrapped in a div, some browsers in
        # some configurations seem to mess up the formatting...
        if nocls:
            yield 0, ('<tr><td><div class="linenodiv" '
                      'style="background-color: #f0f0f0; padding-right: 10px">'
                      '<pre style="line-height: 125%">' +
                      ls + '</pre></div></td><td class="code">')
        else:
            yield 0, ('<tr><td class="linenos"><div class="linenodiv"><pre>' +
                      ls + '</pre></div></td><td class="code">')
        yield 0, dummyoutfile.getvalue()
        yield 0, '</td></tr>'
        self.linenostart += lncount

def view_script(script, script_file, style):

    if isinstance(script, SoS_Workflow):
        # workflow has
        #   filename
        #   name
        #   description
        #   sections
        #   auxillary_sections
        title = script.name
    else:
        # script has
        #   filename
        #   overall description 
        #   workflow descriptions
        #   sections
        title = os.path.basename(script_file)
    #
    # we just get pieces of code, not the complete HTML file.
    #
    html_file = '.sos/{}.html'.format(os.path.basename(script_file))
    formatter = ContinuousHtmlFormatter(cssclass="source", full=False, linenos=True)
    with open(html_file, 'w') as html:
        html.write('''<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
 <title>{}</title>
<meta http-equiv="content-type" content="text/html; charset=None">
<style type="text/css">
pre {{ margin-top: 0px; margin-bottom: 0px; }}
{}
 </style>
</head>
<body>'''.format(title, formatter.get_style_defs()))
        #

        html.write('<table class="%stable">' % formatter.cssclass)
        if script.sections and script.sections[0].global_def:
            html.write('{}\n'.format(highlight(
                script.sections[0].global_def,
                PythonLexer(), formatter)))
        for section in script.sections:
            html.write(highlight('[{}]'.format(section.names), PythonLexer(), formatter) + '\n')
            if section.comment:
                html.write(highlight('# ' + section.comment, PythonLexer(), formatter) + '\n')
            for stmt in section.statements:
                if stmt[0] == ':':
                    html.write('<span style="directive">{}</span>: {}\n'.format(stmt[1], highlight(stmt[2], PythonLexer(),
                        HtmlFormatter(cssclass="source", full=False, linenos=True))))
                elif stmt[0] == '=':
                    html.write(highlight(stmt[2], PythonLexer(), formatter) + '\n')
                else:
                    html.write(highlight(stmt[1], PythonLexer(), formatter) + '\n')
            if section.task:
                html.write(highlight(section.task, PythonLexer(), formatter) + '\n')
        html.write('</table>')
        html.write('''\n</body></html>\n''')
    url = 'file://{}'.format(os.path.abspath(html_file))
    env.logger.info('Viewing {} in a browser'.format(url))
    webbrowser.open(url, new=2)


def sos_show(args, workflow_args):
    try:
        script = SoS_Script(filename=args.script)
        if args.workflow:
            workflow = script.workflow(args.workflow)
            if not args.html:
                workflow.show()
            else:
                view_script(workflow, args.script, args.style)
        elif not args.html:
            script.show()
        else:
            view_script(script, args.script, args.style)
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)

#
# subcommand dryrun
#
def sos_dryrun(args, workflow_args):
    args.__max_jobs__ = 1
    args.__dryrun__ = True
    args.__prepare__ = False
    args.__run__ = False
    args.__rerun__ = False
    args.__config__ = None
    sos_run(args, workflow_args)

#
# subcommand prepare
#
def sos_prepare(args, workflow_args):
    args.__max_jobs__ = 1
    args.__dryrun__ = True
    args.__prepare__ = True
    args.__run__ = False
    args.__rerun__ = False
    args.__config__ = None
    sos_run(args, workflow_args)

#
# subcommand run
#
def sos_run(args, workflow_args):
    env.max_jobs = args.__max_jobs__
    env.verbosity = args.verbosity
    # kill all remainging processes when the master process is killed.
    atexit.register(env.cleanup)
    # default mode: run in dryrun mode
    args.__run__ = not (args.__rerun__ or args.__prepare__ or args.__dryrun__)
    #
    if args.__run__ or args.__rerun__:
        args.__prepare__ = True
    #
    # always run in dryrun mode
    env.run_mode = 'dryrun'
    # if this is not the last step, use verbosity 1 (warning)
    #if args.__prepare__:
    #    env.verbosity = min(args.verbosity, 1)
    #else:
    #
    try:
        script = SoS_Script(filename=args.script)
        workflow = script.workflow(args.workflow)
        executor = Sequential_Executor(workflow)
        executor.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
    except Exception as e:
        if args.verbosity and args.verbosity > 2:
            sys.stderr.write(get_traceback())
        env.logger.error(e)
        sys.exit(1)
    # then prepare mode
    if args.__prepare__:
        # if this is not the last step, use verbosity 1 (warning)
        #if args.__run__ or args.__rerun__:
        #    env.verbosity = min(args.verbosity, 1)
        #else:
        #    env.verbosity = args.verbosity
        #
        env.run_mode = 'prepare'
        try:
            script = SoS_Script(filename=args.script)
            workflow = script.workflow(args.workflow)
            executor.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
        except Exception as e:
            if args.verbosity and args.verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error(e)
            sys.exit(1)
    # then run mode
    if args.__run__ or args.__rerun__:
        env.run_mode = 'run'
        # env.verbosity = args.verbosity
        if args.__rerun__:
            env.sig_mode = 'ignore'
        try:
            script = SoS_Script(filename=args.script)
            workflow = script.workflow(args.workflow)
            executor = Sequential_Executor(workflow)
            executor.run(workflow_args, cmd_name='{} {}'.format(args.script, args.workflow), config_file=args.__config__)
        except Exception as e:
            if args.verbosity and args.verbosity > 2:
                sys.stderr.write(get_traceback())
            env.logger.error(e)
            sys.exit(1)

#
# subcommand config
#
def sos_config(args, workflow_args):
    if workflow_args:
        raise RuntimeError('Unrecognized arguments {}'.format(' '.join(workflow_args)))
    #
    if args.__global_config__:
        config_file = os.path.expanduser('~/.sos/config.json')
    elif args.__config_file__:
        config_file = os.path.expanduser(args.__config_file__)
    else:
        config_file = os.path.expanduser('.sos/config.json')
    if args.__get_config__ is not None:
        if os.path.isfile(config_file):
            try:
                with open(config_file) as config:
                    cfg = yaml.safe_load(config)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to parse sos config file {}, is it in YAML/JSON format? ({}}'.format(config_file, e))
                sys.exit(1)
            for option in (args.__get_config__ if args.__get_config__ else ['*']):
                for k, v in cfg.items():
                    if fnmatch.fnmatch(k, option):
                        print('{}\t{!r}'.format(k, v))
    elif args.__unset_config__:
        if os.path.isfile(config_file):
            try:
                with open(config_file) as config:
                    cfg = yaml.safe_load(config)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to parse sos config file {}, is it in YAML/JSON format? ({})'.format(config_file, e))
                sys.exit(1)
        else:
            env.logger.error('Config file {} does not exist'.format(config_file))
        #
        unset = []
        for option in args.__unset_config__:
            for k in cfg.keys():
                if fnmatch.fnmatch(k, option):
                    unset.append(k)
                    print('Unset {}'.format(k))
        #
        if unset:
            for k in set(unset):
                cfg.pop(k)
            # 
            if unset:
                with open(config_file, 'w') as config:
                    config.write(yaml.safe_dump(cfg, default_flow_style=False))
    elif args.__set_config__:
        if os.path.isfile(config_file):
            try:
                with open(config_file) as config:
                    cfg = yaml.safe_load(config)
                if cfg is None:
                    cfg = {}
            except Exception as e:
                env.logger.error('Failed to sos config file {}, is it in YAML/JSON format? ({})'.format(config_file, e))
                sys.exit(1)
        else:
            cfg = {}
        #
        for option in args.__set_config__:
            k, v = option.split('=', 1)
            try:
                v = eval(v)
            except Exception as e:
                env.logger.error('Cannot interpret option {}. Please quote the string if it is a string option. ({})'.format(option, e))
                sys.exit(1)
            cfg[k] = v
            print('Set {} to {!r}'.format(k, v))
        #
        with open(config_file, 'w') as config:
            config.write(yaml.safe_dump(cfg, default_flow_style=False))
        


