#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import argparse
import base64

import pkg_resources
from sos.utils import dehtml, env


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
            env.logger.warning(f'Ignore incorrect previewer entry point {entrypoint}: {e}')
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


def preview_img(filename, kernel=None, style=None):
    with open(filename, 'rb') as f:
        image = f.read()

    import imghdr
    image_type = imghdr.what(None, image)
    image_data = base64.b64encode(image).decode('ascii')
    if image_type != 'png':
        try:
            if image_type == 'gif':
                return {'image/png': image_data}
            else:
                from wand.image import Image
                img = Image(filename=filename)
                return {'image/' + image_type: image_data,
                        'image/png': base64.b64encode(img._repr_png_()).decode('ascii')}
        except Exception:
            return {'image/' + image_type: image_data}
    else:
        return {'image/' + image_type: image_data}


def preview_svg(filename, kernel=None, style=None):
    with open(filename, 'r') as f:
        image_data = f.read()
    return {'image/svg+xml': image_data}


def _preview_pdf_parser():
    parser = argparse.ArgumentParser(prog='%preview *.pdf')
    parser.add_argument('--pages', nargs='+', type=int,
                        help='Pages of the PDF to preview.')
    parser.error = lambda msg: env.logger.warning(msg)
    return parser


def preview_pdf(filename, kernel=None, style=None):
    use_png = False
    if style is not None and 'style' in style:
        if style['style'] != 'png':
            if kernel is not None and style['style'] is not None:
                kernel.warn(
                    f'Option --style of PDF preview only accept parameter png: {style["style"]} provided')
        else:
            use_png = True
    if use_png:
        try:
            # this import will fail even if wand is installed
            # if imagemagick is not installed properly.
            from wand.image import Image
            img = Image(filename=filename)
            nPages = len(img.sequence)
            pages = list(range(nPages))

            if style is not None and 'options' in style:
                parser = _preview_pdf_parser()
                try:
                    args = parser.parse_args(style['options'])
                except SystemExit:
                    return
                if args.pages is not None:
                    pages = [x - 1 for x in args.pages]
                    for p in pages:
                        if p >= nPages:
                            if kernel is not None:
                                kernel.warn(
                                    f'Page {p} out of range of the pdf file ({nPages} pages)')
                            pages = list(range(nPages))
                            break
            # single page PDF
            if len(pages) == 1 and pages[0] == 0:
                return {
                    'image/png': base64.b64encode(img._repr_png_()).decode('ascii')}
            elif len(pages) == 1:
                # if only one page
                return {
                    'image/png': base64.b64encode(Image(img.sequence[pages[0]])._repr_png_()).decode('ascii')}
            else:
                image = Image(width=img.width, height=img.height * len(pages))
                for i, p in enumerate(pages):
                    image.composite(
                        img.sequence[p],
                        top=img.height * i,
                        left=0
                    )
                return {
                    'image/png': base64.b64encode(image._repr_png_()).decode('ascii')}
        except Exception as e:
            if kernel is not None:
                kernel.warn(e)
            return {'text/html':
                    f'<iframe src={filename} width="100%"></iframe>'}
    else:
        # by default use iframe, because PDF figure can have multiple pages (#693)
        # try to get width and height
        try:
            from wand.image import Image
            img = Image(filename=filename)
            return {'text/html':
                    f'<iframe src={filename} width="800px" height="{img.height/img.width * 800}px"></iframe>'}
        except Exception as e:
            kernel.warn(e)
            return {'text/html':
                    f'<iframe src={filename} width="100%"></iframe>'}


def preview_html(filename, kernel=None, style=None):
    with open(filename) as html:
        content = html.read()
    return {'text/html': content,
            'text/plain': dehtml(content)}


def _get_txt_parser():
    parser = argparse.ArgumentParser(prog='%preview *.txt')
    parser.add_argument('-l', '--limit', type=int, default=5, help='''Maximum number
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

    content = ''
    with open(filename, 'r') as fin:
        if limit < 0:
            content = fin.read()
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
    return f'{len(names)} files\n' + '\n'.join(names[:5]) + ('\n...' if len(names) > 5 else '')


def preview_tar(filename, kernel=None, style=None):
    import tarfile
    with tarfile.open(filename, 'r:*') as tar:
        # only extract files
        names = [x.name for x in tar.getmembers() if x.isfile()]
    return f'{len(names)} files\n' + '\n'.join(names[:5]) + ('\n...' if len(names) > 5 else '')


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
    import os
    import glob
    import tempfile
    from graphviz import Source
    with open(filename) as dot, tempfile.TemporaryDirectory() as tempDirectory:
        src = Source(dot.read())
        src.format = 'png'
        outfile = src.render(filename='sosDot', directory=tempDirectory)
        # dot command can generate more than outfiles returned by the render function
        pngFiles = glob.glob(os.path.join(tempDirectory, f'sosDot*.png'))
        if len(pngFiles) == 1:
            with open(outfile, 'rb') as content:
                data = content.read()
            return {'image/png': base64.b64encode(data).decode('ascii')}
        else:
            import imageio
            # create a gif files from multiple png files
            pngFiles.sort(key=lambda x: int(os.path.basename(x)[:-3].split('.')[1] or 0))
            # getting the maximum size
            try:
                images = [imageio.imread(x) for x in pngFiles]
            except Exception as e:
                if kernel:
                    kernel.warn(f'Failed to read gng file: {e}')
                with open(pngFiles[-1], 'rb') as content:
                    data = content.read()
                return {'image/png': base64.b64encode(data).decode('ascii')}
            maxWidth = max([x.shape[0] for x in images])
            maxHeight = max([x.shape[1] for x in images])
            if images[0].shape[0] < maxWidth or images[0].shape[1] < maxHeight:
                try:
                    from PIL import Image, ImageOps
                    newFirstImg = ImageOps.expand(Image.open(pngFiles[0]), border=(
                        0, 0, (maxHeight - images[0].shape[1]), (maxWidth - images[0].shape[0])), fill=0xFFFFFF)
                    newFirstImg.save(pngFiles[0], directory=tempDirectory)
                    # replace the original small one to the expanded one
                    images[0] = imageio.imread(pngFiles[0])
                except Exception as e:
                    if kernel:
                        kernel.warn(f'Failed to resize gif file: {e}')
            # create a gif file from images
            gifFile = os.path.join('sosDot.gif')
            try:
                imageio.mimsave(gifFile, images, duration=0.5)
            except Exception as e:
                if kernel:
                    kernel.warn(f'Failed to generate gif animation: {e}')
                with open(pngFiles[-1], 'rb') as content:
                    data = content.read()
                return {'image/png': base64.b64encode(data).decode('ascii')}
            with open(gifFile, 'rb') as f:
                image = f.read()
            # according to https://github.com/ipython/ipython/issues/10045
            # I have to use 'image/png' instead of 'image/gif' to get the gif displayed.
            return {'image/png': base64.b64encode(image).decode('ascii')}
