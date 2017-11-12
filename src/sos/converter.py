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
import sys
import argparse
import webbrowser

from io import StringIO
from textwrap import dedent

from pygments import highlight
from pygments.token import Name, Keyword
from pygments.lexers import PythonLexer, TextLexer, get_lexer_by_name, guess_lexer
from pygments.formatters import HtmlFormatter, TerminalFormatter
from pygments.styles import get_all_styles
from pygments.util import shebang_matches

from .utils import env, pretty_size
from .actions import get_actions
from .parser import SoS_Script
from .syntax import SOS_INPUT_OPTIONS, SOS_OUTPUT_OPTIONS, SOS_DEPENDS_OPTIONS, \
    SOS_RUNTIME_OPTIONS, SOS_SECTION_OPTIONS

#
# subcommmand show
#
# CSS style for the output of html format, extracted from github :-)
inline_css = '''
<style type="text/css">
html {
  font-family: sans-serif;
  -ms-text-size-adjust: 100%;
  -webkit-text-size-adjust: 100%;
}

body {
  margin: 0;
}

a {
  background-color: transparent;
  -webkit-text-decoration-skip: objects;
}

a:active,
a:hover {
  outline-width: 0;
}

img {
  border-style: none;
}

pre {
  font-family: monospace, monospace;
  font-size: 1em;
}

button::-moz-focus-inner,
[type="button"]::-moz-focus-inner,
[type="reset"]::-moz-focus-inner,
[type="submit"]::-moz-focus-inner {
  border-style: none;
  padding: 0;
}

button:-moz-focusring,
[type="button"]:-moz-focusring,
[type="reset"]:-moz-focusring,
[type="submit"]:-moz-focusring {
  outline: 1px dotted ButtonText;
}

table {
  border-spacing: 0;
  border-collapse: collapse;
}

td {
  padding: 0;
}

* {
  box-sizing: border-box;
}

body {
  font: 13px / 1.4 Helvetica, arial, nimbussansl, liberationsans, freesans, clean, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";
  color: #333;
  background-color: #fff;
}

a {
  color: #4078c0;
  text-decoration: none;
}

a:hover,
a:active {
  text-decoration: underline;
}

pre {
  margin-top: 0;
  margin-bottom: 0;
  font: 12px Consolas, "Liberation Mono", Menlo, Courier, monospace;
}

.container {
  width: 980px;
  margin-right: auto;
  margin-left: auto;
}

.container::before {
  display: table;
  content: "";
}

.container::after {
  display: table;
  clear: both;
  content: "";
}

.right {
  float: right !important;
}

::-moz-placeholder {
  color: #999;
}

:-ms-input-placeholder {
  color: #999;
}

::placeholder {
  color: #999;
}

.form-select::-ms-expand {
  opacity: 0;
}

.avatar {
  display: inline-block;
  overflow: hidden;
  line-height: 1;
  vertical-align: middle;
  border-radius: 3px;
}

.commit-tease {
  position: relative;
  padding: 10px;
  margin-bottom: -1px;
  line-height: 20px;
  color: #68777d;
  background-color: #f2f9fc;
  border: 1px solid #c9e6f2;
  border-radius: 3px;
  border-bottom-right-radius: 0;
  border-bottom-left-radius: 0;
}

.commit-tease .message {
  color: inherit;
}

.commit-tease .avatar {
  margin-top: -1px;
}

.commit-tease-sha {
  display: inline-block;
  font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  color: #445055;
}

.full-commit .btn-outline:not(:disabled):hover {
  color: #4078c0;
  border: 1px solid #4078c0;
}

.blob-wrapper {
  overflow-x: auto;
  overflow-y: hidden;
  border-bottom-right-radius: 3px;
  border-bottom-left-radius: 3px;
}

.blob-num {
  width: 1%;
  min-width: 50px;
  padding-right: 10px;
  padding-left: 10px;
  font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  font-size: 12px;
  line-height: 18px;
  color: rgba(0,0,0,0.3);
  text-align: right;
  white-space: nowrap;
  vertical-align: bottom;
  cursor: pointer;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
  border: solid #eee;
  border-width: 0 1px 0 0;
}

.blob-num:hover {
  color: rgba(0,0,0,0.6);
}

.blob-num::before {
  content: attr(data-line-number);
}

.blob-code {
  position: relative;
  padding-right: 10px;
  padding-left: 10px;
  vertical-align: top;
  white-space: nowrap;
}

.sos-comment {
  /* background: #f7f7f7; */
}
.sos-report {
  background: #FAFAD2;
  padding-top: 5px;
}
.sos-header {
  margin-top: 5px;
  margin-bottom: 0px;
  padding-top: 5px;
  padding-bottom: 5px;
  background-color: #ddd;
  /* border-bottom-color: #dde4e6; */
  color: darkblue;
  border-top-style: solid;

}
.sos-directive {
  padding-top: 5px;
  /* background-color: #fff9ea; */
  border-color: #f1c0c0;
}
.sos-statement {
  padding-top: 5px;
  /* background-color: #f2f9fc; */
}
.sos-error {
  /* background: #ffdddd; */
  color: red;
}
.sos-script {
  position: relative;
  padding-right: 10px;
  padding-left: 30px;
  vertical-align: top;
  /* background-color: #eff0f1; */
}


.blob-code-inner {
  overflow: visible;
  font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  font-size: 12px;
  color: #333;
  word-wrap: normal;
  white-space: pre;
}

.blob-code-inner::before {
  content: "";
}

body {
  min-width: 1020px;
  word-wrap: break-word;
}

.user-mention {
  font-weight: bold;
  color: #333;
  white-space: nowrap;
}

.prose-diff>.markdown-body :not(li.moved).removed {
  color: #a33;
  text-decoration: line-through;
  background: #ffeaea;
}

.prose-diff>.markdown-body :not(.github-user-ins):not(li.moved).removed {
  text-decoration: line-through;
}

.prose-diff>.markdown-body :not(li.moved).added {
  background: #eaffea;
}

.prose-diff>.markdown-body :not(.github-user-del):not(li.moved).added li:not(.moved):not(.github-user-del).added {
  text-decoration: none;
}

:checked+.radio-label {
  position: relative;
  z-index: 1;
  border-color: #4078c0;
}

.select-menu-text-filter input::-moz-placeholder {
  color: #aaa;
}

.select-menu-text-filter input:-ms-input-placeholder {
  color: #aaa;
}

.select-menu-text-filter input::placeholder {
  color: #aaa;
}

.tab-size[data-tab-size="8"] {
  -moz-tab-size: 8;
  -o-tab-size: 8;
  tab-size: 8;
}

# override source
.source {
  background: auto;
}

.file {
  position: relative;
  margin-top: 20px;
  margin-bottom: 15px;
  border: 1px solid #ddd;
  border-radius: 3px;
}

.file-header {
  padding: 5px 10px;
  background-color: #f7f7f7;
  border-bottom: 1px solid #d8d8d8;
  border-top-left-radius: 2px;
  border-top-right-radius: 2px;
}

.file-header::before {
  display: table;
  content: "";
}

.file-header::after {
  display: table;
  clear: both;
  content: "";
}

.file-actions {
  float: right;
  padding-top: 3px;
}

.file-info {
  float: left;
  font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  font-size: 12px;
  line-height: 32px;
}

.file-info-divider {
  display: inline-block;
  width: 1px;
  height: 18px;
  margin-right: 3px;
  margin-left: 3px;
  vertical-align: middle;
  background-color: #ddd;
}

.header-search-input::-ms-clear {
  display: none;
}

.org-validate-group input::-ms-clear {
  display: none;
}


'''

template_pre_style = '''
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
    <meta charset='utf-8'>
 <title>%s</title>
<meta http-equiv="content-type" content="text/html; charset=None">
'''

template_pre_table = '''
 </style>
</head>
<body>


'''

template_header = '''
<div class="container new-discussion-timeline experiment-repo-nav">

 <div class="commit-tease">
      <span class="right">
      </span>
 <div>
        <a class="header-logo-invertocat" href="https://vatlab.github.io/SOS" data-hotkey="g d" aria-label="Homepage">
        <img src="http://vatlab.github.io/SOS/img/sos_icon.svg" alt="sos_icon">
        </a>

        <span class="user-mention">%s</span>
      </div>
      </div>

  <div class="repository-content">
  <div class="file">
  <div class="file-header">

  <div class="file-info">
      %s
      <span class="file-info-divider"></span>
    %s
  </div>
  <div class="file-actions">
  </div>
  </div>
 <div itemprop="text" class="blob-wrapper data type-python">

'''

template_post_table = '''
</div>
</div>
</div>
</body></html>
'''


class SoS_Lexer(PythonLexer):
    """
    A Python lexer with SOS keywords, used for SoS directive
    """

    name = 'Sript of Scripts'
    aliases = ['sos']

    # override the mimetypes to not inherit them from python
    mimetypes = []
    filenames = ['*.sos']
    mimetypes = ['text/x-sos', 'application/x-sos']

    PythonLexer.tokens['root'].insert(0, (r'(^\w+)\s*:', Keyword.Namespace))

    EXTRA_KEYWORDS = set( SOS_INPUT_OPTIONS + SOS_OUTPUT_OPTIONS +
        SOS_DEPENDS_OPTIONS + SOS_RUNTIME_OPTIONS + SOS_SECTION_OPTIONS +
        get_actions())

    def get_tokens_unprocessed(self, text):
        for index, token, value in \
                PythonLexer.get_tokens_unprocessed(self, text):
            if token is Name and value in self.EXTRA_KEYWORDS:
                yield index, Keyword.Pseudo, value
            else:
                yield index, token, value

    def analyse_text(text):
        return (shebang_matches(text, r'sos-runner') or \
            '#fileformat=SOS' in text[:1000])

class ContinuousHtmlFormatter(HtmlFormatter):
    def _wrap_tablelinenos(self, inner):
        # change wrap_table to do not output '<table' and '</table>' because
        # we will output several tables and use a single table tag around them.
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
            yield 0, ('<tr><td class="blob-num js-line-number"><div class="linenodiv"><pre>' +
                      ls + '</pre></div></td><td class="code">')
        yield 0, dummyoutfile.getvalue()
        yield 0, '</td></tr>'
        self.linenostart += lncount

def write_html_content(content_type, content, formatter, html):
    # dedent content but still keeps empty lines
    old_class = formatter.cssclass
    nlines = len(content)
    content = dedent('\n'.join(content))
    # ' ' to keep pygments from removing empty lines
    # split, merge by \n can introduce one additional line
    content = [' \n' if x == '' else x + '\n' for x in content.split('\n')][:nlines]
    #
    if content_type == 'COMMENT':
        formatter.cssclass = 'source blob-code sos-comment'
        html.write('\n'.join(highlight(x, SoS_Lexer(), formatter) for x in content))
    elif content_type in ('REPORT', 'report'):
        formatter.cssclass = 'source blob-code sos-report'
        html.write('{}\n'.format(highlight(''.join(content),
            TextLexer(), formatter)))
    elif content_type == 'SECTION':
        formatter.cssclass = 'source blob-code sos-header'
        html.write('{}\n'.format(highlight(''.join(content),
            SoS_Lexer(), formatter)))
    elif content_type == 'DIRECTIVE':
        formatter.cssclass = 'source blob-code sos-directive'
        html.write('{}\n'.format(highlight(''.join(content),
            SoS_Lexer(), formatter)))
    elif content_type == 'STATEMENT':
        formatter.cssclass = 'source blob-code sos-statement'
        html.write('{}\n'.format(highlight(''.join(content),
            SoS_Lexer(), formatter)))
    elif content_type == 'ERROR':
        formatter.cssclass = 'source blob-code sos-error '
        html.write('{}\n'.format(highlight(''.join(content),
            SoS_Lexer(), formatter)))
    else:
        formatter.cssclass = 'source blob-code sos-script '
        if content_type == 'run':
            content_type = 'bash'
        elif content_type == 'node':
            content_type = 'JavaScript'
        elif content_type == 'report':
            content_type = 'text'
        try:
            lexer = get_lexer_by_name(content_type)
        except Exception:
            try:
                lexer = guess_lexer(''.join(content))
            except Exception:
                lexer = TextLexer()
        html.write('{}\n'.format(highlight((''.join(content)),
            lexer, formatter)))
    formatter.cssclass = old_class

#
# utility function
#
def transcribe_script(filename=None, content=None):
    with StringIO() as transcript:
        try:
            SoS_Script(filename=filename, content=content, transcript=transcript)
        except Exception as e:
            raise ValueError('Failed to parse {}: {}'.format(filename, e))
        return transcript.getvalue().splitlines()

#
# Converter to HTML
#
#
def get_script_to_html_parser():
    parser = argparse.ArgumentParser('sos convert FILE.sos FILE.html (or --to html)',
        description='''Convert sos file to html format with syntax highlighting,
        and save the output either to a HTML file or view it in a broaser.''')
    parser.add_argument('--raw', help='''URL to the raw sos file, which will be linked
        to filenames in the HTML output''')
    parser.add_argument('--style', choices=list(get_all_styles()),
        help='''Pygments style for the HTML output.''',
        default='default')
    parser.add_argument('--linenos', action='store_true',
        help='''Display lineno to the left of the source code''')
    parser.add_argument('-v', '--view', action='store_true',
        help='''Open the output file in a broswer. In case no html file is specified,
        this option will display the HTML file in a browser, instead of writing its
        content to standard output.''')
    return parser

def script_to_html(script_file, html_file, args=None, unknown_args=None):
    '''
    Convert sos file to html format with syntax highlighting, and
    either save the output either to a HTML file or view it in a broaser.
    This converter accepts additional parameters --style or pygments
    styles, --linenos for displaying linenumbers, and a parameter --raw
    to embed a URL to the raw sos file.
    '''
    import tempfile

    if unknown_args:
        raise ValueError('Unrecognized parameter {}'.format(unknown_args))
    if args:
        formatter = ContinuousHtmlFormatter(cssclass="source", full=False,
            **{x:y for x,y in vars(args).items() if x != ('raw', 'view')})
    else:
        formatter = ContinuousHtmlFormatter(cssclass="source", full=False)

    with StringIO() as html:
        html.write(template_pre_style % os.path.basename(script_file))
        html.write(inline_css)
        # remove background definition so that we can use our own
        html.write('\n'.join(x for x in formatter.get_style_defs().split('\n') if 'background' not in x))
        html.write(template_pre_table)
        if args and args.raw:
            raw_link = f'<a href="{args.raw}" class="commit-tease-sha">{script_file}</a>'
            script_link = f'<a href="{args.raw}">{os.path.basename(script_file)}</a>'
        else:
            raw_link = script_file
            script_link = os.path.basename(script_file)
        html.write(template_header % (script_link, raw_link,
            pretty_size(os.path.getsize(script_file))))
        #
        html.write('<table class="highlight tab-size js-file-line-container">')
        content = []
        content_type = None
        content_number = None
        next_type = None
        for line in transcribe_script(script_file):
            line_type, line_no, script_line = line.split('\t', 2)
            # Does not follow section because it has to be one line
            if line_type == 'FOLLOW' and content_type in (None, 'SECTION'):
                line_type = 'COMMENT'
            if content_type == line_type or line_type == 'FOLLOW':
                if next_type is not None and not script_line.rstrip().endswith(','):
                    formatter.linenostart = content_number
                    write_html_content(content_type, content, formatter, html)
                    content = [script_line]
                    content_type = next_type
                    content_number = int(line_no)
                    next_type = None
                else:
                    content.append(script_line)
            else:
                if content:
                    formatter.linenostart = content_number
                    write_html_content(content_type, content, formatter, html)
                if line_type.startswith('SCRIPT_'):
                    content_type = 'DIRECTIVE'
                    next_type = line_type[7:]
                else:
                    content_type = line_type
                content_number = int(line_no)
                content = [script_line]
        if content:
            formatter.linenostart = content_number
            write_html_content(content_type, content, formatter, html)
        html.write('</table>')
        html.write(template_post_table)

        html_content = html.getvalue()

    # need to save to a temp file to view
    if html_file is None and args and args.view:
        html_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.html', delete=False).name

    if html_file is None:
        sys.stdout.write(html_content)
    else:
        with open(html_file, 'w') as out:
            out.write(html_content)
        env.logger.info(f'SoS script saved to {html_file}')
        #
    if args and args.view:
        url = f'file://{html_file}'
        env.logger.info(f'Viewing {url} in a browser')
        webbrowser.open(url, new=2)

#
# Converter to Markdown
#
#

def markdown_content(content_type, content, fh):
    # write content to a file
    import re
    if content_type == 'COMMENT':
        fh.write('{}\n'.format(''.join([x.lstrip('#').strip() + '  \n'
                                        for x in content if not x.startswith('#!')])))
    elif content_type == 'REPORT':
        fh.write(f'{"".join(content)}\n')
    elif content_type == 'SECTION':
        fh.write(f'## {"".join(content)}\n')
    elif content_type == 'DIRECTIVE':
        fh.write('{}\n'.format(''.join(['**{}**  \n'.format(re.sub(r'(\$|_)', r'`\1`', x).strip())
                                        for x in content])))
    elif content_type == 'STATEMENT':
        fh.write(f'```python\n{"".join(content)}\n```\n')
    elif content_type == 'ERROR':
        fh.write(f'{"".join(content)}\n')
    else:
        if content_type == 'run':
            content_type = 'bash'
        elif content_type == 'node':
            content_type = 'JavaScript'
        elif content_type == 'report':
            content_type = ''
        fh.write(f'```{content_type}\n{"".join(content)}```\n')

def get_script_to_markdown_parser():
    parser = argparse.ArgumentParser('sos convert FILE.sos FILE.md (or --to md)',
        description='''Convert SOS script to a markdown format with scripts
            quoted in markdown syntax.''')
    return parser

def script_to_markdown(script_file, markdown_file, style_args=None, unknown_args=None):
    '''
    Convert SOS scriot to a markdown file with syntax highlighting.
    '''
    if unknown_args:
        raise ValueError(f'Unrecognized parameter {unknown_args}')

    if not markdown_file:
        markdown = sys.stdout
    else:
        markdown = open(markdown_file, 'w')
    # remove background definition so that we can use our own
    content = []
    content_type = None
    # content_number = None
    next_type = None
    for line in transcribe_script(script_file):
        line_type, _, script_line = line.split('\t', 2)
        # Does not follow section because it has to be one line
        if line_type == 'FOLLOW' and content_type in (None, 'SECTION'):
            line_type = 'COMMENT'
        if content_type == line_type or line_type == 'FOLLOW':
            if next_type is not None and not script_line.rstrip().endswith(','):
                markdown_content(content_type, content, markdown)
                content = [script_line]
                content_type = next_type
                # content_number = int(line_no)
                next_type = None
            else:
                content.append(script_line)
        else:
            if content:
                markdown_content(content_type, content, markdown)
            if line_type.startswith('SCRIPT_'):
                content_type = 'DIRECTIVE'
                next_type = line_type[7:]
            else:
                content_type = line_type
            # content_number = int(line_no)
            content = [script_line]
    if content:
        markdown_content(content_type, content, markdown)
    if markdown != sys.stdout:
        markdown.close()
        env.logger.info(f'SoS script saved to {markdown_file}')

#
# Output to terminal
#
def write_content(content_type, content, formatter, fh=sys.stdout):
    #
    nlines = len(content)
    content = dedent(''.join(content))
    # ' ' to keep pygments from removing empty lines
    # split, merge by \n can introduce one additional line
    content = [' \n' if x == '' else x + '\n' for x in content.split('\n')][:nlines]
    #
    if content_type == 'COMMENT':
        fh.write(highlight(''.join(content), SoS_Lexer(), formatter))
    elif content_type in ('REPORT', 'report'):
        fh.write(highlight(''.join(content), TextLexer(), formatter))
    elif content_type == 'SECTION':
        fh.write(highlight(''.join(content), SoS_Lexer(), formatter))
    elif content_type == 'DIRECTIVE':
        fh.write(highlight(''.join(content), SoS_Lexer(), formatter))
    elif content_type == 'STATEMENT':
        fh.write(highlight(''.join(content), SoS_Lexer(), formatter))
    elif content_type == 'ERROR':
        fh.write(highlight(''.join(content), SoS_Lexer(), formatter))
    else:
        if content_type == 'run':
            content_type = 'bash'
        elif content_type == 'node':
            content_type = 'JavaScript'
        elif content_type == 'report':
            content_type = 'text'
        try:
            lexer = get_lexer_by_name(content_type)
        except Exception:
            try:
                lexer = guess_lexer(''.join(content))
            except Exception:
                lexer = TextLexer()
        fh.write(highlight((''.join(content)), lexer, formatter))

def get_script_to_term_parser():
    parser = argparse.ArgumentParser('sos convert FILE.sos --to term',
        description='Write script to terminal with syntax highlighting.')
    parser.add_argument('--bg', choices=['light', 'dark'],
        help='Color theme of the output',
        default='light')
    parser.add_argument('--linenos', action='store_true',
        help='Display lineno to the left of the script')
    return parser

def script_to_term(script_file, output_file, args, unknown_args=None):
    '''
    Write script to terminal. This converter accepts additional parameters
    --bg [light|dark] for light or dark theme, and --linenos for output
    lineno.
    '''

    if unknown_args:
        raise ValueError('Unrecognized parameter {}'.format(unknown_args))

    formatter = TerminalFormatter(**vars(args))
    # remove background definition so that we can use our own
    content = []
    content_type = None
    # content_number = None
    next_type = None
    for line in transcribe_script(script_file):
        line_type, _, script_line = line.split('\t', 2)
        # Does not follow section because it has to be one line
        if line_type == 'FOLLOW' and content_type in (None, 'SECTION'):
            line_type = 'COMMENT'
        if content_type == line_type or line_type == 'FOLLOW':
            if next_type is not None and not script_line.rstrip().endswith(','):
                write_content(content_type, content, formatter)
                content = [script_line]
                content_type = next_type
                # content_number = int(line_no)
                next_type = None
            else:
                content.append(script_line)
        else:
            if content:
                write_content(content_type, content, formatter)
            if line_type.startswith('SCRIPT_'):
                content_type = 'DIRECTIVE'
                next_type = line_type[7:]
            else:
                content_type = line_type
            # content_number = int(line_no)
            content = [script_line]
    if content:
        write_content(content_type, content, formatter)
