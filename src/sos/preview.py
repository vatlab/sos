#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import argparse
import base64
import io
import pkg_resources
from sos.utils import dehtml, env, dot_to_gif, linecount_of_file


def get_previewers():
    # Note: data is zest.releaser specific: we want to pass
    # something to the plugin
    group = 'sos_previewers'
    result = []
    for entrypoint in pkg_resources.iter_entry_points(group=group):
        # if ':' in entry point name, it should be a function
        try:
            name, priority = entrypoint.name.split(',', 1)
            priority = int(priority)
        except Exception as e:
            env.logger.warning(
                f'Ignore incorrect previewer entry point {entrypoint}: {e}')
            continue
        # If name points to a function in a module. Let us try to import the module
        if ':' in name:
            import importlib
            try:
                mod, func = name.split(':')
                imported = importlib.import_module(mod)
                result.append((getattr(imported, func), entrypoint, priority))
            except ImportError:
                env.logger.warning(f'Failed to load function {mod}:{func}')
        else:
            result.append((name, entrypoint, priority))
    #
    result.sort(key=lambda x: -x[2])
    # we put string infront of functions so that filename matching will be used before
    # content matching. For example, the xlsx file is actually a zip file so it could be
    # previewed as a zip file first.
    return [x for x in result if isinstance(x[0], str) and x[0] != '*'] + \
        [x for x in result if not isinstance(x[0], str)] + \
        [x for x in result if x[0] == '*']


def _preview_img_parser():
    parser = argparse.ArgumentParser(prog='%preview *.pdf')
    parser.add_argument(
        '--width',
        help='Width of the previewed image, can be in any HTML units such as px, em.'
    )
    parser.add_argument(
        '--height',
        help='Height of the previewed image, can be in any HTML units such as px, em.'
    )
    parser.error = lambda msg: env.logger.warning(msg)
    return parser


def preview_img(filename, kernel=None, style=None):
    with open(filename, 'rb') as f:
        image = f.read()
    import imghdr
    image_type = imghdr.what(None, image)
    image_data = base64.b64encode(image).decode('ascii')

    args = None
    meta = {}
    if style is not None and 'options' in style:
        parser = _preview_img_parser()
        try:
            args = parser.parse_args(style['options'])
            meta.update({
                'image/png':
                    dict(([['width', args.width]] if args.width else []) +
                         ([['height', args.height]] if args.height else []))
            })
        except SystemExit:
            return

    if image_type != 'png':
        try:
            if image_type == 'gif':
                return {'image/png': image_data}, meta
            else:
                from wand.image import Image
                img = Image(filename=filename)
                return {
                    'image/' + image_type:
                        image_data,
                    'image/png':
                        base64.b64encode(img._repr_png_()).decode('ascii')
                }, meta
        except Exception:
            return {'image/' + image_type: image_data}, meta
    else:
        return {'image/' + image_type: image_data}, meta


def preview_svg(filename, kernel=None, style=None):
    with open(filename, 'r') as f:
        image_data = f.read()
    return {'image/svg+xml': image_data}


def _preview_pdf_parser():
    parser = argparse.ArgumentParser(prog='%preview *.pdf')
    parser.add_argument(
        '--pages', nargs='+', type=int, help='Pages of the PDF to preview.')
    parser.add_argument(
        '--width',
        help='Width of the previewed image, can be in any HTML units such as px, em.'
    )
    parser.add_argument(
        '--height',
        help='Height of the previewed image, can be in any HTML units such as px, em.'
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='resolution of the converted png preview')
    parser.error = lambda msg: env.logger.warning(msg)
    return parser


def preview_pdf(filename, kernel=None, style=None):
    use_png = False
    warn = kernel.warn if kernel is not None else env.logger.warning
    if style is not None and 'style' in style:
        if style['style'] != 'png':
            if style['style'] is not None:
                warn(
                    f'Option --style of PDF preview only accept parameter png: {style["style"]} provided'
                )
        else:
            use_png = True
    args = None
    if style is not None and 'options' in style:
        parser = _preview_pdf_parser()
        try:
            args = parser.parse_args(style['options'])
        except SystemExit:
            return
    meta = {}
    embed_options = ''
    if args and (args.width or args.height):
        meta.update({
            'image/png':
                dict(([['width', args.width]] if args.width else []) +
                     ([['height', args.height]] if args.height else []))
        })
        embed_options += (f'width="{args.width}" ' if args.width else ' ') + \
            (f'height="{args.height}" ' if args.height else ' ')
    if use_png:
        try:
            # this import will fail even if wand is installed
            # if imagemagick is not installed properly.
            from wand.image import Image
            img = Image(filename=filename, resolution=args.dpi)

            if img.width == 0 or img.height == 0:
                raise ValueError('Image appears to have zero width or height')
            nPages = len(img.sequence)
            pages = list(range(nPages))

            if args and args.pages is not None:
                pages = [x - 1 for x in args.pages]
                for p in pages:
                    if p >= nPages:
                        warn(
                            f'Page {p} out of range of the pdf file ({nPages} pages)'
                        )
                        pages = list(range(nPages))
                        break
            # single page PDF
            if len(pages) == 1 and pages[0] == 0:
                return {
                    'image/png':
                        base64.b64encode(img._repr_png_()).decode('ascii')
                }, meta
            elif len(pages) == 1:
                # if only one page
                return {
                    'image/png':
                        base64.b64encode(
                            Image(img.sequence[pages[0]])._repr_png_()
                        ).decode('ascii')
                }, meta
            else:
                widths = [img.sequence[p].width for p in pages]
                heights = [img.sequence[p].height for p in pages]
                image = Image(width=max(widths), height=sum(heights))
                for i, p in enumerate(pages):
                    image.composite(
                        img.sequence[p], top=sum(heights[:i]), left=0)
                image.format = 'png'
                with io.BytesIO() as out:
                    image.save(file=out)
                    img_data = out.getvalue()
                return {
                    'image/png': base64.b64encode(img_data).decode('ascii')
                }, meta
        except Exception as e:
            warn(e)
            return {
                'text/html':
                    f'<embed src="{filename}" {embed_options} type="application/pdf" />'
            }
    else:
        # by default use iframe, because PDF figure can have multiple pages (#693)
        return {
            'text/html':
                f'<embed src="{filename}" {embed_options} type="application/pdf" />'
        }


def preview_html(filename, kernel=None, style=None):
    with open(filename) as html:
        content = html.read()
    return {'text/html': content, 'text/plain': dehtml(content)}


def _get_txt_parser():
    parser = argparse.ArgumentParser(prog='%preview *.txt')
    parser.add_argument(
        '-l',
        '--limit',
        type=int,
        default=5,
        help='''Maximum number
        of lines to preview.''')
    parser.error = lambda msg: env.logger.warning(msg)
    return parser


def preview_txt(filename, kernel=None, style=None):
    if style is not None and 'options' in style:
        parser = _get_txt_parser()
        try:
            args = parser.parse_args(style['options'])
        except SystemExit:
            return
        limit = args.limit
    else:
        limit = 5

    nlines = linecount_of_file(filename)
    content = 'HINT: {} line{}{}\n'.format(
        nlines, 's' if nlines > 1 else '',
        f' ({limit} displayed, see --limit)' if nlines > limit else '')

    with open(filename, 'r') as fin:
        if limit < 0:
            content += fin.read()
        else:
            for _ in range(limit):
                content += fin.readline()
    return content


def preview_csv(filename, kernel=None, style=None):
    import pandas
    from .visualize import Visualizer
    data = pandas.read_csv(filename)
    return Visualizer(kernel, style).preview(data)


def preview_xls(filename, kernel=None, style=None):
    import pandas
    from .visualize import Visualizer
    data = pandas.read_excel(filename)
    return Visualizer(kernel, style).preview(data)


def preview_zip(filename, kernel=None, style=None):
    import zipfile
    names = zipfile.ZipFile(filename).namelist()
    return f'{len(names)} files\n' + '\n'.join(
        names[:5]) + ('\n...' if len(names) > 5 else '')


def preview_tar(filename, kernel=None, style=None):
    import tarfile
    with tarfile.open(filename, 'r:*') as tar:
        # only extract files
        names = [x.name for x in tar.getmembers() if x.isfile()]
    return f'{len(names)} files\n' + '\n'.join(
        names[:5]) + ('\n...' if len(names) > 5 else '')


def preview_gz(filename, kernel=None, style=None):
    import tarfile
    # .tar.gz
    if tarfile.is_tarfile(filename):
        return preview_tar(filename, kernel, style)
    import gzip
    content = b''
    with gzip.open(filename, 'rb') as fin:
        for _ in range(5):
            content += fin.readline()
    try:
        return content.decode()[:2000]
    except Exception:
        return 'binary data'


def preview_md(filename, kernel=None, style=None):
    import markdown
    with open(filename) as fin:
        text = fin.read()
    html = markdown.markdown(text)
    return {'text/html': html, 'text/plain': text}


def preview_dot(filename, kernel=None, style=None):
    data = dot_to_gif(
        filename, warn=kernel.warn if kernel else env.logger.warning)
    # according to https://github.com/ipython/ipython/issues/10045
    # I have to use 'image/png' instead of 'image/gif' to get the gif displayed.
    return {'image/png': data}
