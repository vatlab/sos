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

import sys
import argparse

from io import StringIO

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook
from nbconvert.exporters import Exporter

from sos.utils import env
from sos.sos_syntax import SOS_SECTION_HEADER, SOS_CELL_LINE


#
# Converter from Notebook
#

def get_notebook_to_script_parser():
    parser = argparse.ArgumentParser('sos convert FILE.ipynb FILE.sos (or --to sos)',
        description='''Export Jupyter notebook with a SoS kernel to a 
        .sos file. The cells are presented in the .sos file as
        cell structure lines, which will be ignored if executed
        in batch mode ''')
    parser.add_argument('--reorder', action='store_true',
        help='''Reorder cells according to execution count''')
    parser.add_argument('--reset-index', action='store_true',
        help='''Reset index to 1, 2, 3, ... regardless of existing
            execution counts''')
    parser.add_argument('--add-header', action='store_true',
        help='''Add default section headers to cells without
            headers, which would help running the converted
            script in batch mode''')
    parser.add_argument('--no-index', action='store_true',
        help='''Do not output any index''')
    parser.add_argument('--metadata', action='store_true',
        help='''Output metadata''')
    parser.add_argument('--remove-magic', action='store_true',
        help='''Remove magic lines from the output''')
    return parser


# This class cannot be defined in .kernel because it would cause some
# weird problem with unittesting not able to resolve __main__
class SoS_Exporter(Exporter):
    def __init__(self, config=None, reorder=False, reset_index=False, add_header=False,
            no_index=False, no_metadata=False, remove_magic=False, 
            **kwargs):
        self.reorder = reorder
        self.reset_index = reset_index
        self.add_header = add_header
        self.no_index = no_index
        self.no_metadata = no_metadata
        self.remove_magic = remove_magic
        self.output_extension = '.sos'
        self.output_mimetype = 'text/x-sos'
        Exporter.__init__(self, config, **kwargs)

    def from_notebook_cell(self, cell, fh, idx = 0):
        if self.no_metadata:
            meta = ''
        else:
            meta = ' '.join('{}={}'.format(x,y) for x,y in cell.metadata.items())
        if not hasattr(cell, 'execution_count') or cell.execution_count is None or self.no_index:
            fh.write('\n%cell {} {}\n'.format(cell.cell_type, meta))
        else:
            idx += 1
            fh.write('\n%cell {} {} {}\n'.format(cell.cell_type,
                                              idx if self.reset_index else cell.execution_count, meta))
        if cell.cell_type == 'code':
            if any(cell.source.startswith(x) for x in ('%run', '%restart', '%dict', '%get', '%use', '%with', '%set', '%paste', '%matplotlib', '%edit')):
                if self.remove_magic:
                    cell.source = '\n'.join(cell.source.split('\n')[1:])
            if self.add_header and not any([SOS_SECTION_HEADER.match(x) for x in cell.source.split('\n')]):
                cell.source = '[{}]\n'.format(idx if self.reset_index else cell.execution_count) + cell.source
            fh.write(cell.source.strip() + '\n')
        elif cell.cell_type == "markdown":
            fh.write('\n'.join('#! ' + x for x in cell.source.split('\n')) + '\n')
        return idx

    def from_notebook_node(self, nb, resources, **kwargs):
        #
        if self.reorder:
            unnumbered_cells = {x: y for x, y in enumerate(nb.cells)
                              if not hasattr(y, 'execution_count') or y.execution_count is None}
            numbered_cells = [y for y in nb.cells
                              if hasattr(y, 'execution_count') and y.execution_count is not None]
            numbered_cells = sorted(numbered_cells, key = lambda x: x.execution_count)
            cells = []
            for idx in range(len(nb.cells)):
                if idx in unnumbered_cells:
                    cells.append(unnumbered_cells[idx])
                else:
                    cells.append(numbered_cells.pop(0))
        else:
            cells = nb.cells
        with StringIO() as fh:
            fh.write('#!/usr/bin/env sos-runner\n')
            fh.write('#fileformat=SOS1.0\n')
            idx = 0
            for cell in cells:
                idx = self.from_notebook_cell(cell, fh, idx)
            content = fh.getvalue()
        resources['output_extension'] = '.sos'
        return content, resources


def notebook_to_script(notebook_file, sos_file, args=None, unknown_args=[]):
    '''
    Convert a ipython notebook to sos format. This converter accepts options
    --reorder to reorder cells according to executing order, --reset-index
    to reset executing counts, --add-header to add a SoS header, --no-index
    to ignore indexes, --remove-magic to remove ipynb-only magics, and 
    --md-to-report to convert markdown cells to sos report actions.
    '''
    if unknown_args:
        raise ValueError('Unrecognized parameter {}'.format(' '.join(unknown_args)))
    if args:
        exporter = SoS_Exporter(reorder=args.reorder, reset_index=args.reset_index,
                            add_header=args.add_header, no_index=args.no_index,
                            no_metadata=args.no_metadata, remove_magic=args.remove_magic)
    else:
        exporter = SoS_Exporter()
    notebook = nbformat.read(notebook_file, nbformat.NO_CONVERT)
    output, resource = exporter.from_notebook_node(notebook, {})
    if not sos_file:
        sys.stdout.write(output)
    else:
        with open(sos_file, 'w') as sos:
            sos.write(output)
        env.logger.info('SoS script saved to {}'.format(sos_file))

#
# Converter to Notebook
#
def get_script_to_notebook_parser():
    parser = argparse.ArgumentParser('sos convert FILE.sos FILE._ipynb (or --to ipynb)',
        description='''Convert a sos script to Jupyter notebook (.ipynb)
            so that it can be opened by Jupyter notebook.''')
    parser.add_argument('-e', '--execute', action='store_true',
        help='''Execute the notebook as if running "Cell -> Run all" from the
            Jupyter notebook interface''')
    return parser

def add_cell(cells, content, cell_type, cell_count, metainfo):
    # if a section consist of all report, report it as a markdown cell
    if not content:
        return
    if cell_type not in ('code', 'markdown'):
        env.logger.warning('Unrecognized cell type {}, code assumed.'.format(cell_type))
    if cell_type == 'markdown' and any(x.strip() and not x.startswith('#! ') for x in content):
        env.logger.warning('Markdown lines not starting with #!, code cell assumed.')
        cell_type = 'code'
    #
    if cell_type == 'markdown':
        cells.append(new_markdown_cell(source=''.join([x[3:] for x in content]).strip(),
            metadata=metainfo))
    else:
        cells.append(
             new_code_cell(
                 # remove any trailing blank lines...
                 source=''.join(content).strip(),
                 execution_count=cell_count,
                 metadata=metainfo)
        )

def script_to_notebook(script_file, notebook_file, args=None, unknown_args=[]):
    '''
    Convert a sos script to iPython notebook (.ipynb) so that it can be opened
    by Jupyter notebook.
    '''
    if unknown_args:
        raise ValueError('Unrecognized parameter {}'.format(' '.join(unknown_args)))
    cells = []
    cell_count = 1
    cell_type = 'code'
    metainfo = {}
    content = []

    with open(script_file) as script:
        first_block = True
        for line in script:
            if line.startswith('#') and first_block:
                if line.startswith('#!'):
                    continue
                if line.startswith('#fileformat='):
                    if not line[12:].startswith('SOS'):
                        raise RuntimeError('{} is not a SoS script according to #fileformat line.'.format(script_file))
                    continue

            first_block = False

            mo = SOS_CELL_LINE.match(line)
            if mo:
                # get ride of empty content
                if not any(x.strip() for x in content):
                    content = []

                if content:
                    add_cell(cells, content, cell_type, cell_count, metainfo)

                cell_type = mo.group('cell_type')
                if not cell_type:
                    cell_type = 'code'
                cc = mo.group('cell_count')
                if cc:
                    cell_count = int(cc)
                else:
                    cell_count += 1
                metainfo = mo.group('metainfo')
                if metainfo:
                    pieces = [piece.split('=',1) for piece in metainfo.split()]
                    for idx,piece in enumerate(pieces):
                        if len(piece) == 1:
                            env.logger.warning('Incorrect metadata {}'.format(piece))
                            pieces[idx].append('')
                        if piece[1] == 'True':
                            pieces[idx][1] = True
                        elif piece[1] == 'False':
                            pieces[idx][1] = False
                    metainfo = {x:y for x,y in pieces}
                else:
                    metainfo = {}
                content = []
                continue
            else:
                content.append(line)
    #
    if content and any(x.strip() for x in content):
        add_cell(cells, content, cell_type, cell_count, metainfo)
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
    err = None
    if args and args.execute:
        from nbconvert.preprocessors import ExecutePreprocessor, CellExecutionError
        ep = ExecutePreprocessor(timeout=600, kernel_name='sos')
        try:
            ep.preprocess(nb, {'metadata': {'path': '.'}})
        except CellExecutionError as e:
            err = e
    #
    if not notebook_file:
        nbformat.write(nb, sys.stdout, 4)
    else:
        with open(notebook_file, 'w') as notebook:
            nbformat.write(nb, notebook, 4)
        env.logger.info('Jupyter notebook saved to {}'.format(notebook_file))
    if err:
        raise err


#
# notebook to HTML
#

def export_notebook(exporter_class, to_format, notebook_file, output_file, unknown_args=[]):

    import os
    import subprocess
    if not os.path.isfile(notebook_file):
        raise RuntimeError('{} does not exist'.format(notebook_file))
    if not output_file:
        import tempfile
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.' + to_format).name
        tmp_stderr = tempfile.NamedTemporaryFile(delete=False, suffix='.' + to_format).name
        with open(tmp_stderr, 'w') as err:
            ret = subprocess.call(['jupyter', 'nbconvert', notebook_file, '--to', to_format,
                '--output', tmp] + unknown_args, stderr=err)
        with open(tmp_stderr) as err:
            err_msg = err.read()
        if ret != 0:
            env.logger.error(err_msg) 
            env.logger.error('Failed to convert {} to {} format'.format(notebook_file, to_format))
        else:
            # identify output files
            dest_file = err_msg.rsplit()[-1]
            if not os.path.isfile(dest_file):
                env.logger.error(err_msg)
                env.logger.error('Failed to get converted file.')
            else:
                with open(dest_file, 'rb') as tfile:
                    sys.stdout.buffer.write(tfile.read())
        try:
            os.remove(tmp)
        except:
            pass
    else:
        ret = subprocess.call(['jupyter', 'nbconvert', notebook_file, '--to', to_format,
            '--output', output_file] + unknown_args)
        if ret != 0:
            env.logger.error('Failed to convert {} to {} format'.format(notebook_file, to_format))
        else:
            env.logger.info('Output saved to {}'.format(output_file))

  
def get_notebook_to_html_parser():
    parser = argparse.ArgumentParser('sos convert FILE.ipynb FILE.html (or --to html)',
        description='''Export Jupyter notebook with a SoS kernel to a 
        .html file. Additional command line arguments are passed directly to 
        command "jupyter nbconvert --to html" so please refer to nbconvert manual for
        available options.''')
    return parser

def notebook_to_html(notebook_file, output_file, sargs=None, unknown_args=[]):
    from nbconvert.exporters.html import HTMLExporter
    export_notebook(HTMLExporter, 'html', notebook_file, output_file, unknown_args)

def get_notebook_to_pdf_parser():
    parser = argparse.ArgumentParser('sos convert FILE.ipynb FILE.pdf (or --to pdf)',
        description='''Export Jupyter notebook with a SoS kernel to a 
        .pdf file. Additional command line arguments are passed directly to 
        command "jupyter nbconvert --to pdf" so please refer to nbconvert manual for
        available options.''')
    return parser

def notebook_to_pdf(notebook_file, output_file, sargs=None, unknown_args=[]):
    from nbconvert.exporters.pdf import PDFExporter
    export_notebook(PDFExporter, 'pdf', notebook_file, output_file, unknown_args)

def get_notebook_to_md_parser():
    parser = argparse.ArgumentParser('sos convert FILE.ipynb FILE.md (or --to md)',
        description='''Export Jupyter notebook with a SoS kernel to a 
        markdown file. Additional command line arguments are passed directly to 
        command "jupyter nbconvert --to markdown" so please refer to nbconvert manual for
        available options.''')
    return parser

def notebook_to_md(notebook_file, output_file, sargs=None, unknown_args=[]):
    from nbconvert.exporters.markdown import MarkdownExporter
    export_notebook(MarkdownExporter, 'markdown', notebook_file, output_file, unknown_args)


