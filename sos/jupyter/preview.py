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
import base64
from IPython.core.display import HTML
from sos.utils import env, dehtml

def preview_img(filename):
    with open(filename, 'rb') as f:
        image = f.read()

    import imghdr
    image_type = imghdr.what(None, image)
    image_data = base64.b64encode(image).decode('ascii')
    if image_type != 'png':
        try:
            from wand.image import Image
            img = Image(filename=filename)
            return { 'image/' + image_type: image_data,
                'image/png': base64.b64encode(img._repr_png_()).decode('ascii') }
        except Exception:
            return { 'image/' + image_type: image_data }
    else:
        return { 'image/' + image_type: image_data }

def preview_pdf(filename):
    try:
        # this import will fail even if wand is installed
        # if imagemagick is not installed properly.
        from wand.image import Image
        img = Image(filename=filename)
        return {
            'text/html': HTML('<iframe src={0} width="100%"></iframe>'.format(filename)).data,
            'image/png': base64.b64encode(img._repr_png_()).decode('ascii') }
    except Exception as e:
        env.logger.error(e)
        return { 'text/html':
            HTML('<iframe src={0} width="100%"></iframe>'.format(filename)).data}

def preview_html(filename):
    with open(filename) as html:
        content = html.read()
    return { 'text/html': content,
        'text/plain': dehtml(content) }

def preview_txt(filename):
    try:
        content = ''
        with open(filename, 'r') as fin:
            for line in range(5):
                content += fin.readline()
        return content
    except:
        return ''

def preview_csv(filename):
    try:
        import pandas
        data = pandas.read_csv(filename)
        html = data._repr_html_()
        return { 'text/html': HTML(html).data}
    except Exception as e:
        env.logger.warning(e)
        return ''

def preview_xls(filename):
    try:
        import pandas
        data = pandas.read_excel(filename)
        html = data._repr_html_()
        return { 'text/html': HTML(html).data}
    except Exception as e:
        env.logger.warning(e)
        return ''

def preview_zip(filename):
    import zipfile
    zip = zipfile.ZipFile(filename)
    names = zip.namelist()
    return '{} files\n'.format(len(names)) + '\n'.join(names[:5]) + ('\n...' if len(names) > 5 else '')

def preview_tar(filename):
    import tarfile
    with tarfile.open(filename, 'r:*') as tar:
        # only extract files
        names = [x.name for x in tar.getmembers() if x.isfile()]
    return '{} files\n'.format(len(names)) + '\n'.join(names[:5]) + ('\n...' if len(names) > 5 else '')

def preview_gz(filename):
    import gzip
    content = b''
    with gzip.open(filename, 'rb') as fin:
        for line in range(5):
            content += fin.readline()
    try:
        return content.decode()[:2000]
    except:
        return 'binary data'

def preview_md(filename):
    import markdown
    try:
        with open(filename) as fin:
            text = fin.read()
        html = markdown.markdown(text)
        return {'text/html': HTML(html).data, 'text/plain': text}
    except:
        return ''
    
def preview_dot(filename):
    from graphviz import Source
    with open(filename) as dot:
        src = Source(dot.read())
    src.format='png'
    outfile = src.render()
    with open(outfile, 'rb') as content:
        data = content.read()
    return {'image/png': base64.b64encode(data).decode('ascii') }
