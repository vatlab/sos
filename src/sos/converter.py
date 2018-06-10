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
    parser.add_argument('--style', help='''deprecated option''',
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
        'script': content,
        'sos_version': __version__,
        'linenos': args.linenos,
        'raw': args.raw,
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
    content = [' \n' if x == '' else x +
               '\n' for x in content.split('\n')][:nlines]
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
        for line in lines:
            if valid_cell or (line.startswith('%include') or line.startswith('%from')):
                content += line + '\n'
            elif SOS_SECTION_HEADER.match(line):
                valid_cell = True
                content += line + '\n'
        if valid_cell:
            content += '\n'
    return content
