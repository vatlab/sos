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
import argparse
import webbrowser
import re

from io import StringIO
from textwrap import dedent

from pygments import highlight
from pygments.token import Name, Keyword
from pygments.lexers import PythonLexer, TextLexer, get_lexer_by_name, guess_lexer
from pygments.formatters import HtmlFormatter, TerminalFormatter
from pygments.styles import get_all_styles
from pygments.util import shebang_matches

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

from .utils import env
from .actions import get_actions
from .kernel import SoS_Exporter
from .sos_syntax import SOS_INPUT_OPTIONS, SOS_OUTPUT_OPTIONS, SOS_DEPENDS_OPTIONS, \
    SOS_RUNTIME_OPTIONS, SOS_SECTION_OPTIONS, SOS_SECTION_HEADER

__all__ = ['SoS_Lexer']

#
# subcommmand show
#
# CSS style for the output of html format, extracted from github :-)
inline_css = '''
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
  background: #f7f7f7;
}
.sos-report {
  background: #FAFAD2;
  padding-top: 5px;
}
.sos-header {
  margin-top: 20px;
  margin-bottom: 0px;
  padding-top: 5px;
  background-color: #C0C0C0;
  border-bottom-color: #dde4e6;

}
.sos-directive {
  padding-top: 5px;
  background-color: #fff9ea;
  border-color: #f1c0c0;
}
.sos-statement {
  padding-top: 5px;
  background-color: #f2f9fc;
}
.sos-error {
  background: #ffdddd;
}
.sos-script {
  position: relative;
  padding-right: 10px;
  padding-left: 30px;
  vertical-align: top;
  background-color: #eff0f1;
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
<style type="text/css">
'''

template_pre_table = '''

 </style>
</head>
<body>


<div class="container new-discussion-timeline experiment-repo-nav">

 <div class="commit-tease">
      <span class="right">
      </span>
 <div>
        <a class="header-logo-invertocat" href="https://github.com/BoPeng/SOS" data-hotkey="g d" aria-label="Homepage">
        <svg xmlns="http://www.w3.org/2000/svg" style="width: 25px; height: 25px" viewBox="0 0 64 64" enable-background="new 0 0 64 64"><path fill="#ec1c24" d="M63.791,56.913c0,3.802-3.076,6.882-6.876,6.882H6.874C3.078,63.795,0,60.715,0,56.913V6.877    C0,3.08,3.078,0,6.874,0h50.042c3.8,0,6.876,3.08,6.876,6.877V56.913z"/><path fill="#c32129" d="m0 48.574v8.339c0 3.802 3.078 6.882 6.874 6.882h50.04c3.8 0 6.876-3.08 6.876-6.882v-50.04c-10.894 39.569-50.05 42.08-63.792 41.697"/><g fill="#fff"><path d="m20.445 36.649c0 1.311-.337 2.49-1.015 3.536-.678 1.051-1.667 1.869-2.973 2.463-1.305.59-2.849.888-4.636.888-2.146 0-3.911-.405-5.306-1.215-.986-.584-1.789-1.357-2.404-2.335-.617-.973-.927-1.919-.927-2.836 0-.532.186-.989.553-1.371.374-.38.844-.568 1.418-.568.462 0 .856.148 1.178.444.319.296.592.733.821 1.319.278.689.574 1.268.896 1.729.321.467.776.849 1.356 1.148.582.304 1.35.455 2.298.455 1.305 0 2.361-.304 3.174-.911.815-.607 1.227-1.368 1.227-2.277 0-.719-.221-1.303-.66-1.754-.439-.448-1.01-.791-1.705-1.028-.697-.239-1.628-.492-2.792-.758-1.56-.365-2.867-.792-3.919-1.28-1.054-.49-1.885-1.156-2.502-2-.619-.844-.929-1.893-.929-3.148 0-1.195.325-2.255.979-3.186.65-.927 1.593-1.64 2.828-2.14 1.238-.498 2.688-.749 4.357-.749 1.334 0 2.485.167 3.457.496.974.333 1.783.772 2.426 1.319.643.548 1.111 1.126 1.405 1.726.3.602.447 1.191.447 1.764 0 .523-.186.993-.554 1.417-.375.417-.835.629-1.387.629-.502 0-.888-.127-1.148-.38-.265-.249-.545-.663-.853-1.234-.396-.823-.866-1.462-1.422-1.918-.555-.46-1.44-.692-2.667-.692-1.136 0-2.051.251-2.747.749-.693.5-1.046 1.102-1.046 1.799 0 .439.121.811.356 1.127.237.317.562.591.978.817.416.227.835.406 1.26.535.425.125 1.128.312 2.104.561 1.223.286 2.334.602 3.323.948.996.345 1.838.766 2.536 1.26.696.494 1.238 1.117 1.63 1.875.389.755.585 1.68.585 2.777"/><path d="m31.57 21.08c2.249 0 4.184.457 5.8 1.372 1.614.915 2.835 2.212 3.664 3.895.829 1.686 1.246 3.663 1.246 5.933 0 1.682-.228 3.208-.68 4.579-.457 1.372-1.137 2.563-2.048 3.572-.907 1.01-2.023 1.775-3.347 2.309-1.324.534-2.841.801-4.548.801-1.699 0-3.223-.274-4.565-.82-1.344-.549-2.465-1.323-3.36-2.322-.899-.995-1.579-2.194-2.038-3.6-.457-1.4-.69-2.915-.69-4.546 0-1.667.241-3.202.719-4.595.48-1.393 1.172-2.577 2.081-3.554.905-.98 2.02-1.727 3.319-2.246 1.304-.52 2.786-.778 4.447-.778m6.253 11.172c0-1.593-.257-2.968-.772-4.136-.516-1.164-1.246-2.045-2.202-2.643-.952-.598-2.044-.896-3.278-.896-.88 0-1.693.167-2.439.496-.743.331-1.385.811-1.924 1.444-.539.631-.964 1.44-1.275 2.424-.312.981-.466 2.086-.466 3.311 0 1.232.155 2.352.466 3.345.312.997.75 1.826 1.318 2.484.568.657 1.223 1.146 1.958 1.472.737.329 1.542.49 2.42.49 1.127 0 2.161-.28 3.105-.845.944-.562 1.695-1.431 2.252-2.608.555-1.173.837-2.619.837-4.338"/><path d="m60.38 36.649c0 1.311-.339 2.49-1.016 3.536-.676 1.051-1.665 1.869-2.97 2.463-1.304.59-2.85.888-4.637.888-2.146 0-3.91-.405-5.302-1.215-.989-.584-1.793-1.357-2.41-2.335-.615-.973-.925-1.919-.925-2.836 0-.532.185-.989.553-1.371.372-.38.845-.568 1.416-.568.465 0 .858.148 1.18.444.319.296.594.733.821 1.319.274.689.576 1.268.894 1.729.323.467.772.849 1.357 1.148.582.304 1.347.455 2.299.455 1.305 0 2.361-.304 3.174-.911.815-.607 1.225-1.368 1.225-2.277 0-.719-.219-1.303-.66-1.754-.438-.448-1.01-.791-1.7-1.028-.7-.239-1.629-.492-2.794-.758-1.561-.365-2.867-.792-3.921-1.28-1.053-.49-1.884-1.156-2.503-2-.619-.844-.927-1.893-.927-3.148 0-1.195.325-2.255.98-3.186.648-.927 1.591-1.64 2.829-2.14 1.232-.498 2.687-.749 4.354-.749 1.337 0 2.486.167 3.463.496.972.333 1.778.772 2.422 1.319.639.548 1.108 1.126 1.402 1.726.302.602.447 1.191.447 1.764 0 .523-.185.993-.553 1.417-.374.417-.835.629-1.388.629-.501 0-.888-.127-1.148-.38-.264-.249-.546-.663-.85-1.234-.398-.823-.866-1.462-1.427-1.918-.553-.46-1.44-.692-2.665-.692-1.137 0-2.054.251-2.751.749-.691.5-1.04 1.102-1.04 1.799 0 .439.115.811.356 1.127.235.317.561.591.974.817.419.227.837.406 1.262.535.424.125 1.129.312 2.101.561 1.225.286 2.336.602 3.327.948.991.345 1.838.766 2.535 1.26.696.494 1.239 1.117 1.627 1.875.394.755.589 1.68.589 2.777"/></g> </svg>
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
    content = dedent(''.join(content))
    # ' ' to keep pygments from removing empty lines
    # split, merge by \n can introduce one additional line
    content = [' \n' if x == '' else x + '\n' for x in content.split('\n')][:nlines]
    #
    if content_type == 'COMMENT':
        formatter.cssclass = 'source blob-code sos-comment'
        html.write('{}\n'.format(highlight(''.join(content),
            SoS_Lexer(), formatter)))
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
    elif content_type == 'ASSIGNMENT':
        formatter.cssclass = 'source blob-code sos-statement'
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
            content_type == 'text'
        try:
            lexer = get_lexer_by_name(content_type)
        except:
            try:
                lexer = guess_lexer(''.join(content))
            except:
                lexer = TextLexer()
        html.write('{}\n'.format(highlight((''.join(content)),
            lexer, formatter)))
    formatter.cssclass = old_class

#
# Converter to HTML
#
#
def parse_html_args(style_args):
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw')
    parser.add_argument('--style', choices=list(get_all_styles()), default='default')
    parser.add_argument('--linenos', action='store_true')
    try:
        args = parser.parse_args(style_args)
    except Exception as e:
        raise RuntimeError('Unrecognized style argument {}: {}'
            .format(' '.join(style_args), e))
    return vars(args)

def script_to_html(transcript, script_file, html_file, style_args):
    '''
    Write a HTML file with the transcript of a SOS file.
    '''
    if html_file == '__BROWSER__':
        html_file = '.sos/{}.html'.format(os.path.basename(script_file))
    sargs = parse_html_args(style_args)
    formatter = ContinuousHtmlFormatter(cssclass="source", full=False,
        **{x:y for x,y in sargs.items() if x != 'raw'})
    with open(html_file, 'w') as html:
        html.write(template_pre_style % os.path.basename(script_file))
        html.write(inline_css)
        # remove background definition so that we can use our own
        html.write('\n'.join(x for x in formatter.get_style_defs().split('\n') if 'background' not in x))
        if 'raw' in sargs and sargs['raw']:
            raw_link = '<a href="{}" class="commit-tease-sha">{}</a>'.format(sargs['raw'], script_file)
        else:
            raw_link = script_file
        html.write(template_pre_table % (os.path.basename(script_file), raw_link,
            '{:.1f} KB'.format(os.path.getsize(script_file) / 1024)))
        #
        html.write('<table class="highlight tab-size js-file-line-container">')
        with open(transcript) as script:
            content = []
            content_type = None
            content_number = None
            next_type = None
            for line in script:
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
    #
    if html_file.startswith('.sos'):
        url = 'file://{}'.format(os.path.abspath(html_file))
        env.logger.info('Viewing {} in a browser'.format(url))
        webbrowser.open(url, new=2)
    else:
        env.logger.info('SoS script saved to {}'.format(html_file))


def workflow_to_html(workflow, script_file, html_file, style_args):
    '''
    Write a HTML file with the transcript of a SOS file.
    '''
    if html_file == '__BROWSER__':
        html_file = '.sos/{}.html'.format(os.path.basename(script_file))
    sargs = parse_html_args(style_args)
    formatter = ContinuousHtmlFormatter(cssclass="source", full=False,
        **{x:y for x,y in sargs.items() if x != 'raw'})
    with open(html_file, 'w') as html:
        html.write(template_pre_style % ('{} - {}'.format(os.path.basename(script_file), workflow)))
        html.write(inline_css)
        # remove background definition so that we can use our own
        html.write('\n'.join(x for x in formatter.get_style_defs().split('\n') if 'background' not in x))
        if 'raw' in sargs and sargs['raw']:
            raw_link = '<a href="{}" class="commit-tease-sha">{}</a>'.format(sargs['raw'], script_file)
        else:
            raw_link = script_file
        html.write(template_pre_table % (os.path.basename(script_file), raw_link, ''))
        #
        html.write('<table class="highlight tab-size js-file-line-container">')
        if workflow.sections and workflow.sections[0].global_def:
            write_html_content('STATEMENT', workflow.sections[0].global_def, formatter, html)
        #
        for section in workflow.sections:
            write_html_content('SECTION', '[{}_{}]'.format(section.name, section.index), formatter, html)
            #
            if section.comment:
                write_html_content('COMMENT', '#' + section.comment, formatter, html)
            for stmt in section.statements:
                if stmt[0] == ':':
                    write_html_content('DIRECTIVE', '{} : {}'.format(stmt[1], stmt[2]), formatter, html)
                elif stmt[0] == '=':
                    write_html_content('STATEMENT', '{} = {}'.format(stmt[1], stmt[2]), formatter, html)
                else:
                    write_html_content('STATEMENT', stmt[1].strip(), formatter, html)
            if section.task:
                #write_html_content('DIRECTIVE', 'task:', formatter, html)
                write_html_content('STATEMENT', section.task.strip(), formatter, html)
        html.write('</table>')
        html.write(template_post_table)
    #
    if html_file.startswith('.sos'):
        url = 'file://{}'.format(os.path.abspath(html_file))
        env.logger.info('Viewing {} in a browser'.format(url))
        webbrowser.open(url, new=2)
    else:
        env.logger.info('Workflow saved to {}'.format(html_file))

#
# Converter to Markdown
#
#

def markdown_content(content_type, content, fh):
    # write content to a file
    if content_type == 'COMMENT':
        fh.write('{}\n'.format(''.join([x.lstrip('#').strip() + '  \n'
                                        for x in content if not x.startswith('#!')])))
    elif content_type == 'REPORT':
        fh.write('{}\n'.format(''.join(content)))
    elif content_type == 'SECTION':
        fh.write('## {}\n'.format(''.join(content)))
    elif content_type == 'DIRECTIVE':
        fh.write('{}\n'.format(''.join(['**{}**  \n'.format(re.sub(r'(\$|_)', r'`\1`', x).strip())
                                        for x in content])))
    elif content_type == 'ASSIGNMENT':
        fh.write('```python\n{}\n```\n'.format(''.join(content))),
    elif content_type == 'STATEMENT':
        fh.write('```python\n{}\n```\n'.format(''.join(content))),
    elif content_type == 'ERROR':
        fh.write('{}\n'.format(''.join(content))),
    else:
        if content_type == 'run':
            content_type = 'bash'
        elif content_type == 'node':
            content_type = 'JavaScript'
        elif content_type == 'report':
            content_type == ''
        fh.write('```{}\n{}```\n'.format(content_type, ''.join(content)))

def script_to_markdown(transcript, script_file, markdown_file):
    '''
    Write a markdown file with the transcript of a SOS file.
    '''
    with open(markdown_file, 'w') as markdown:
        # remove background definition so that we can use our own
        with open(transcript) as script:
            content = []
            content_type = None
            # content_number = None
            next_type = None
            for line in script:
                line_type, line_no, script_line = line.split('\t', 2)
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
    #
    env.logger.info('SoS script saved to {}'.format(markdown_file))


def workflow_to_markdown(workflow, script_file, markdown_file):
    '''
    Write a markdown file with the transcript of a SOS file.
    '''
    with open(markdown_file, 'w') as markdown:
        if workflow.sections and workflow.sections[0].global_def:
            markdown_content('STATEMENT', workflow.sections[0].global_def, markdown)
        #
        for section in workflow.sections:
            markdown_content('SECTION', '[{}_{}]'.format(section.name, section.index), markdown)
            #
            if section.comment:
                markdown_content('COMMENT', '#' + section.comment, markdown)
            for stmt in section.statements:
                if stmt[0] == ':':
                    markdown_content('DIRECTIVE', '{} : {}'.format(stmt[1], stmt[2]), markdown)
                elif stmt[0] == '=':
                    markdown_content('STATEMENT', '{} = {}'.format(stmt[1], stmt[2]), markdown)
                else:
                    markdown_content('STATEMENT', stmt[1].strip(), markdown)
            if section.task:
                markdown_content('STATEMENT', section.task.strip(), markdown)
            markdown_content('COMMENT', '\n', markdown)
    #
    env.logger.info('Workflow saved to {}'.format(markdown_file))


#
# Converter from Notebook
#

def notebook_to_script(notebook_file, sos_file, convert_args=[]):
    '''
    convert a ipython notebook.
    '''
    exporter = SoS_Exporter()
    notebook = nbformat.read(notebook_file, nbformat.NO_CONVERT)
    output, resource = exporter.from_notebook_node(notebook, {})
    if sos_file == '__STDOUT__':
        sys.stdout.write(output)
    else:
        with open(sos_file, 'w') as sos:
            sos.write(output)
        env.logger.info('SoS script saved to {}'.format(sos_file))

#
# Converter to Notebook
#
def add_cell(cells, cell_src_code, cell_count):
    # if a section consist of all report, report it as a markdown cell
    if not cell_src_code:
        return
    if (cell_src_code[0].startswith('!') or not cell_src_code[0].strip() or SOS_SECTION_HEADER.match(cell_src_code[0])) and \
        all(not x.strip() or x.startswith('!') for x in cell_src_code[1:]):
        md = ''
        for line in cell_src_code:
            if SOS_SECTION_HEADER.match(line):
                continue
            else:
                md += line[2:]
        cells.append(new_markdown_cell(source=md))
        return cell_count
    else:
        cells.append(
             new_code_cell(
                 source=''.join(cell_src_code),
                 execution_count=cell_count)
        )
        return cell_count

def script_to_notebook(transcript, script_file, notebook_file):
    '''
    Write a notebook file with the transcript of a SOS file.
    '''
    cells = []
    cell_count = 1
    CELL_LINE = re.compile('^#cell\s+(markdown|code)(\s+\d+\s+)?$')
    with open(transcript) as script:
        cell_src_code = []
        content = []
        content_type = None
        # content_number = None
        next_type = None
        for line in script:
            line_type, line_no, script_line = line.split('\t', 2)
            if line_type == 'COMMENT':
                if script_line.startswith('#!'):
                    continue
                if script_line.startswith('#fileformat='):
                    if not script_line[12:].startswith('SOS'):
                        raise RuntimeError('{} is not a SoS script according to #fileformat line.'.format(script_file))
                    continue
            if CELL_LINE.match(script_line):
                # eat a new line before the CELL_LINE
                if content and content[-1] == '\n':
                    content = content[:-1]
                if not any(x.strip() for x in content):
                    continue
                cell_src_code.extend(content)
                add_cell(cells, cell_src_code, cell_count)
                parts = script_line.split()
                if len(parts) > 2 and parts[2].isdigit():
                    cell_count = int(parts[2])
                if parts[1] == 'markdown':
                    content_type == 'REPORT'
                else:
                    content_type == 'SECTION'
                content = []
                cell_src_code = []
                continue
            # Does not follow section because it has to be one line
            if line_type == 'FOLLOW' and content_type in (None, 'SECTION'):
                line_type = 'COMMENT'
            if content_type == line_type or line_type == 'FOLLOW':
                if next_type is not None and not script_line.rstrip().endswith(','):
                    if content_type == 'SECTION' and cell_src_code:
                        # wrap up a cell
                        add_cell(cells, cell_src_code, cell_count)
                        cell_count += 1
                        cell_src_code = content
                    else:
                        cell_src_code.extend(content)
                    content = [script_line]
                    content_type = next_type
                    # content_number = int(line_no)
                    next_type = None
                else:
                    content.append(script_line)
            else:
                if content:
                    if content_type == 'SECTION' and cell_src_code:
                        # wrap up a cell
                        add_cell(cells, cell_src_code, cell_count)
                        cell_count += 1
                        cell_src_code = content
                    else:
                        cell_src_code.extend(content)
                if line_type.startswith('SCRIPT_'):
                    content_type = 'DIRECTIVE'
                    next_type = line_type[7:]
                else:
                    content_type = line_type
                # content_number = int(line_no)
                content = [script_line]
    if content:
        add_cell(cells, cell_src_code, cell_count)
    #
    nb = new_notebook(cells = cells,
        metadata = {
            'kernelspec' : {
                "display_name": "SoS",
                "language": "sos",
                "name": "sos"
            },
            "language_info": {
                "file_extension": ".sos",
                "mimetype": "text/x-sos",
                "name": "sos",
                "pygments_lexer": "python"
            }
        }
    )
    if notebook_file == '__STDOUT__':
        nbformat.write(nb, sys.stdout, 4)
    else:
        with open(notebook_file, 'w') as notebook:
            nbformat.write(nb, notebook, 4)
        env.logger.info('SoS script saved to {}'.format(notebook_file))


def workflow_to_notebook(workflow, script_file, notebook_file):
    '''
    Write a notebook file with the transcript of a SOS file.
    '''
    cells = []
    cell_count = 1
    if workflow.sections and workflow.sections[0].global_def:
        cells.append(
            new_code_cell(
            source = workflow.sections[0].global_def,
            execution_count=cell_count)
        )
        cell_count += 1
    #
    for section in workflow.sections:
        # FIXME: section option is ignored
        cell_src_code = '[{}_{}]\n'.format(section.name, section.index)
        #
        if section.comment:
            cell_src_code += '# ' + section.comment + '\n'
        for stmt in section.statements:
            if stmt[0] == ':':
                cell_src_code += '{}: {}\n'.format(stmt[1], stmt[2])
            elif stmt[0] == '=':
                cell_src_code += '{} = {}\n'.format(stmt[1], stmt[2])
            else:
                cell_src_code += stmt[1].strip() + '\n'
        if section.task:
            cell_src_code += 'task:\n' + section.task.strip() + '\n'
        cells.append(
            new_code_cell(
            source = cell_src_code,
            execution_count=cell_count)
        )
        cell_count += 1
    #
    nb = new_notebook(cells = cells,
        metadata = {
            'kernelspec' : {
                "display_name": "SoS",
                "language": "sos",
                "name": "sos"
            },
            "language_info": {
                "file_extension": ".sos",
                "mimetype": "text/x-sos",
                "name": "sos",
                "pygments_lexer": "python"
            }
        }
    )
    with open(notebook_file, 'w') as notebook:
        nbformat.write(nb, notebook, 4)
    #
    env.logger.info('Workflow saved to {}'.format(notebook_file))

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
    elif content_type == 'ASSIGNMENT':
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
            content_type == 'text'
        try:
            lexer = get_lexer_by_name(content_type)
        except:
            try:
                lexer = guess_lexer(''.join(content))
            except:
                lexer = TextLexer()
        fh.write(highlight((''.join(content)), lexer, formatter))

def parse_term_args(style_args):
    parser = argparse.ArgumentParser()
    parser.add_argument('--bg', choices=['light', 'dark'], default='light')
    parser.add_argument('--linenos', action='store_true')
    try:
        args = parser.parse_args(style_args)
    except Exception as e:
        raise RuntimeError('Unrecognized style argument {}: {}'
            .format(' '.join(style_args), e))
    return vars(args)

def script_to_term(transcript, script_file, style_args):
    '''
    Write script to terminal
    '''
    sargs = parse_term_args(style_args)
    env.logger.trace('Using style argument {}'.format(sargs))
    formatter = TerminalFormatter(**sargs)
    # remove background definition so that we can use our own
    with open(transcript) as script:
        content = []
        content_type = None
        # content_number = None
        next_type = None
        for line in script:
            line_type, line_no, script_line = line.split('\t', 2)
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

def workflow_to_term(workflow, script_file, style_args):
    '''
    Write a workflow to terminal
    '''
    sargs = parse_term_args(style_args)
    env.logger.trace('Using style argument {}'.format(sargs))
    formatter = TerminalFormatter(**sargs)
    if workflow.sections and workflow.sections[0].global_def:
        write_content('STATEMENT', workflow.sections[0].global_def, formatter)
    #
    for section in workflow.sections:
        write_content('SECTION', '[{}_{}]'.format(section.name, section.index), formatter)
        #
        if section.comment:
            write_content('COMMENT', '#' + section.comment, formatter)
        for stmt in section.statements:
            if stmt[0] == ':':
                write_content('DIRECTIVE', '{} : {}'.format(stmt[1], stmt[2]), formatter)
            elif stmt[0] == '=':
                write_content('STATEMENT', '{} = {}'.format(stmt[1], stmt[2]), formatter)
            else:
                write_content('STATEMENT', stmt[1].strip(), formatter)
        if section.task:
            write_content('STATEMENT', section.task.strip(), formatter)
        write_content('COMMENT', '\n', formatter)
