#!/usr/bin/env python3
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import argparse
import os
import sys
from io import StringIO
from textwrap import dedent
from ._version import __version__

import nbformat
from pygments import highlight
from pygments.formatters import HtmlFormatter, TerminalFormatter
from pygments.lexers import (PythonLexer, TextLexer, get_lexer_by_name,
                             guess_lexer)
from pygments.styles import get_all_styles
from pygments.token import Keyword, Name
from pygments.util import shebang_matches

from .actions import get_actions
from .parser import SoS_Script
from .syntax import (SOS_DEPENDS_OPTIONS, SOS_INPUT_OPTIONS,
                     SOS_OUTPUT_OPTIONS, SOS_RUNTIME_OPTIONS,
                     SOS_SECTION_HEADER, SOS_SECTION_OPTIONS)
from .utils import env, pretty_size


class SoS_Lexer(PythonLexer):
    """
    A Python lexer with SOS keywords, it is not used by sos right now but
    needs to be kept for compatibility with older notebooks.
    """
    name = 'Sript of Scripts'
    aliases = ['sos']

    # override the mimetypes to not inherit them from python
    mimetypes = []
    filenames = ['*.sos']
    mimetypes = ['text/x-sos', 'application/x-sos']

    PythonLexer.tokens['root'].insert(0, (r'(^\w+)\s*:', Keyword.Namespace))

    EXTRA_KEYWORDS = set(SOS_INPUT_OPTIONS + SOS_OUTPUT_OPTIONS +
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
        return (shebang_matches(text, r'sos-runner') or
                '#fileformat=SOS' in text[:1000])


#
# Converter to HTML
#
#
codemirror_themes = [
    'default',
    '3024-day',
    '3024-night',
    'abcdef',
    'ambiance',
    'base16-dark',
    'base16-light',
    'bespin',
    'blackboard',
    'cobalt',
    'colorforth',
    'darcula',
    'dracula',
    'duotone-dark',
    'duotone-light',
    'eclipse',
    'elegant',
    'erlang-dark',
    'gruvbox-dark',
    'hopscotch',
    'icecoder',
    'idea',
    'isotope',
    'lesser-dark',
    'liquibyte',
    'lucario',
    'material',
    'mbo',
    'mdn-like',
    'midnight',
    'monokai',
    'neat',
    'neo',
    'night',
    'oceanic-next',
    'panda-syntax',
    'paraiso-dark',
    'paraiso-light',
    'pastel-on-dark',
    'railscasts',
    'rubyblue',
    'seti',
    'shadowfox',
    'solarized dark',
    'solarized light',
    'the-matrix',
    'tomorrow-night-bright',
    'tomorrow-night-eighties',
    'ttcn',
    'twilight',
    'vibrant-ink',
    'xq-dark',
    'xq-light',
    'yeti',
    'zenburn'
]


def get_script_to_html_parser():
    parser = argparse.ArgumentParser('sos convert FILE.sos FILE.html (or --to html)',
                                     description='''Convert sos file to html format with syntax highlighting,
        and save the output either to a HTML file or view it in a broaser.''')
    parser.add_argument('--raw', help='''URL to the raw sos file, which will be linked
        to filenames in the HTML output''')
    parser.add_argument('--style', choices=codemirror_themes,
                        help="Code mirror themes for the output",  default='default')
    parser.add_argument('--linenos', action='store_true',
                        help=argparse.SUPPRESS)
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
    from jinja2 import Environment, PackageLoader, select_autoescape
    environment = Environment(
        loader=PackageLoader('sos', 'templates'),
        autoescape=select_autoescape(['html', 'xml'])
    )
    template = environment.get_template('html_report.tpl')

    with open(script_file) as script:
        content = script.read()
    context = {
        'filename': script_file,
        'basename': os.path.basename(script_file),
        'script': content,
        'sos_version': __version__,
        'linenos': args.linenos if args and hasattr(args, 'linenos') else True,
        'raw': args.raw if args and hasattr(args, 'raw') else '',
        'theme': args.style if args and hasattr(args, 'style') else 'default',
    }
    html_content = template.render(context)
    if html_file is None:
        sys.stdout.write(html_content)
    else:
        with open(html_file, 'w') as out:
            out.write(html_content)
        env.logger.info(f'SoS script saved to {html_file}')
        #
    if args and args.view:
        import webbrowser
        url = f'file://{os.path.abspath(html_file)}'
        env.logger.info(f'Viewing {url} in a browser')
        webbrowser.open(url, new=2)


def extract_workflow(notebook_file):
    nb = nbformat.read(notebook_file, nbformat.NO_CONVERT)
    cells = nb.cells
    content = '#!/usr/bin/env sos-runner\n#fileformat=SOS1.0\n\n'
    for cell in cells:
        if cell.cell_type != "code":
            continue
        # Non-sos code cells are also ignored
        if 'kernel' in cell.metadata and cell.metadata['kernel'] not in ('sos', 'SoS', None):
            continue
        lines = cell.source.split('\n')
        valid_cell = False
        for idx,line in enumerate(lines):
            if valid_cell or (line.startswith('%include') or line.startswith('%from')):
                content += line + '\n'
            elif SOS_SECTION_HEADER.match(line):
                valid_cell = True
                # look retrospectively for comments
                c = idx - 1
                comment = ''
                while c >= 0 and lines[c].startswith('#'):
                    comment = lines[c] + '\n' + comment
                    c -= 1
                content += comment + line + '\n'
        if valid_cell:
            content += '\n'
    return content
