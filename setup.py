#!/usr/bin/env python
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

import sys, os
import shutil
from setuptools import find_packages, setup
from distutils import log
from setuptools.command.install import install

# obtain version of SoS
with open('sos/_version.py') as version:
    for line in version:
        if line.startswith('__version__'):
            __version__ = eval(line.split('=')[1])
            break

kernel_json = {
    "argv":         ["python", "-m", "sos.jupyter.kernel", "-f", "{connection_file}"],
    "display_name": "SoS",
    "language":     "sos",
}

class InstallWithConfigurations(install):
    def run(self):
        # Regular installation
        install.do_egg_install(self)

        # copy sos.vim and sos-detect.vim to .vim
        vim_syntax_dir = os.path.expanduser('~/.vim/syntax')
        vim_syntax_file = os.path.join(vim_syntax_dir, 'sos.vim')
        if not os.path.isdir(vim_syntax_dir):
            os.makedirs(vim_syntax_dir)
        shutil.copy('misc/sos.vim', vim_syntax_file)
        #
        vim_ftdetect_dir = os.path.expanduser('~/.vim/ftdetect')
        vim_ftdetect_file = os.path.join(vim_ftdetect_dir, 'sos.vim')
        if not os.path.isdir(vim_ftdetect_dir):
            os.makedirs(vim_ftdetect_dir)
        shutil.copy('misc/sos-detect.vim', vim_ftdetect_file)
        # copy vim-ipython to .vim/ftplugin
        vim_plugin_dir = os.path.expanduser('~/.vim/ftplugin/sos')
        if not os.path.isdir(vim_plugin_dir):
            os.makedirs(vim_plugin_dir)
        shutil.copy('misc/vim-ipython/ipy.vim', os.path.join(vim_plugin_dir, 'ipy.vim'))
        shutil.copy('misc/vim-ipython/vim_ipython.py', os.path.join(vim_plugin_dir, 'vim_ipython.py'))
        #
        # at this point, jupyter and ipython should have been installed.
        import json
        try:
            from jupyter_client.kernelspec import KernelSpecManager as KS
        except ImportError:
            from ipykernel.kernelspec import KernelSpecManager as KS
        from IPython.utils.tempdir import TemporaryDirectory
        from IPython.paths import get_ipython_dir
        #
        # copy ipython magic to ~/.ipython/extensions
        ext_dir = os.path.join(get_ipython_dir(), 'extensions')
        ext_file = os.path.join(ext_dir, 'sos_magic.py')
        if not os.path.isdir(ext_dir):
            os.makedirs(ext_dir)
        prof_dir = os.path.join(get_ipython_dir(), 'profile_sos')
        prof_file = os.path.join(prof_dir, 'ipython_config.py')
        if not os.path.isdir(prof_dir):
            os.makedirs(prof_dir)
        #
        shutil.copy('sos/jupyter/sos_magic.py', ext_file)
        shutil.copy('sos/jupyter/sos_ipython_profile.py', prof_file)
        #
        log.info('\nSoS is installed and configured to use with vim, ipython, and Jupyter.')
        log.info('Use "set syntax=sos" to enable syntax highlighting.')
        log.info('Use "ipython --profile sos" to start ipython with sos magic.')
        #
        # Now write the kernelspec
        with TemporaryDirectory() as td:
            os.chmod(td, 0o755)  # Starts off as 700, not user readable
            shutil.copy('sos/jupyter/kernel.js', os.path.join(td, 'kernel.js'))
            shutil.copy('misc/logo-64x64.png', os.path.join(td, 'logo-64x64.png'))
            with open(os.path.join(td, 'kernel.json'), 'w') as f:
                json.dump(kernel_json, f, sort_keys=True)
            try:
                KS().install_kernel_spec(td, 'sos', user=self.user, replace=True, prefix=sys.exec_prefix)
                log.info('Use "jupyter notebook" to create or open SoS notebooks.')
            except:
                log.error("\nWARNING: Could not install SoS Kernel as %s user." % self.user)
        log.info('Run "python misc/patch_spyder.py" to patch spyder with sos support.')
        log.info('And "sos -h" to start using Script of Scripts.')

dest = '''\
Exploratory data analysis in computationally intensive disciplines such as computational
biology often requires one to exploit a variety of tools implemented in different programming
languages and analyzing large datasets on high performance computing systems (e.g. computer
clusters). On top of all the difficulties in exchanging data between languages and computing
systems and analyzing data on different platforms, it becomes challenging to keep track of
such fragmented workflows and reproduce prior analyses.

With strong emphases on readability, practicality, and reproducibility, we have developed
a workflow system called “Script of Scripts” (SoS) with a web front-end and notebook format
based on Jupyter. Major features of SoS for exploratory analysis include multi-language
support, explicit and automatic data exchange between running sessions (kernels) in
different languages, cell-specific kernel switch using frontend-UI or cell magics, 
a side-panel that allows scratch execution of statements, preview of files and expressions,
and line-by-line execution of statements in cells. In particular, variable and file preview
on the side panel makes it possible to trouble-shoot scripts in multiple languages without 
contaminating the main notebook or interrupting the logic flow of the analysis. For large-scale
data analysis, the SoS workflow engine provides a unified interface to executing and managing
tasks on a variety of computing platforms such as PBS/Torch/LSF/Slurm clusters and RQ and
Celery task queues. Specified files are automatically synchronized between file systems,
thus enabling a single workflow to utilize multiple remote computing environments.

Researchers will benefit from the SoS system the flexibility to use their preferred languages
and tools for tasks without having to worry about data flow, and can perform light interactive
analysis while executing heavy remote tasks simultaneous in the same notebook in a neat and
organized fashion. SoS is available at http://vatlab.github.io/SOS/ and is distributed freely
under a GPL3 license. A live Jupyter server and several docker containers are available for
testing and running SoS without a local installation.

Please refer to http://vatlab.github.io/SOS/ for more details on SoS.
'''

setup(name = "sos",
    version = __version__,
    description = 'Script of Scripts (SoS): an interactive, cross-platform, and cross-language workflow system for reproducible data analysis',
    long_description=dest,
    author = 'Bo Peng',
    url = 'https://github.com/vatlab/SOS',
    author_email = 'bpeng@mdanderson.org',
    maintainer = 'Bo Peng',
    maintainer_email = 'bpeng@mdanderson.org',
    license = 'GPL3',
    include_package_data = True,
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        ],
    packages = find_packages(),
    cmdclass={'install': InstallWithConfigurations},
    install_requires=[
          'psutil',
          # progress bar
          'tqdm',
          # for file lock
          'fasteners',
          'pyyaml',
          'pygments',
          # for jupyter notebook format conversion
          'nbformat',
          'nbconvert>=5.1.1',
          'ipython',
          'ipykernel',
          'notebook>=5.0.0',
          'ptpython',
          # for DAG
          'networkx',
          'pydotplus',
      ],
    entry_points= '''
[console_scripts]
sos = sos.__main__:main
sos-runner = sos.__main__:sosrunner

[sos_addons]
patch-spyder.parser = sos.addons.patch_spyder:get_patch_spyder_parser
patch-spyder.func   = sos.addons.patch_spyder:patch_spyder

[pygments.lexers]
sos = sos.converter:SoS_Lexer

[sos_targets]
dynamic = sos.target:dynamic
remote = sos.target:remote
local = sos.target:local
executable = sos.target:executable
sos_variable = sos.target:sos_variable
sos_step = sos.target:sos_step
env_variable = sos.target:env_variable
bundle = sos.target:bundle
R_library = sos.R.target:R_library
Py_Module = sos.Python3.target:Py_Module

[sos_actions]
execute_script = sos.actions:execute_script
sos_run = sos.actions:sos_run
fail_if = sos.actions:fail_if
warn_if = sos.actions:warn_if
stop_if = sos.actions:stop_if
download = sos.actions:download
run = sos.actions:run

bash = sos.Bash.actions:bash
csh = sos.Bash.actions:csh
tcsh = sos.Bash.actions:tcsh
zsh = sos.Bash.actions:zsh
sh = sos.Bash.actions:sh

perl = sos.actions:perl
ruby = sos.actions:ruby
node = sos.actions:node
JavaScript = sos.actions:JavaScript
report = sos.actions:report
pandoc = sos.actions:pandoc

python = sos.Python3.actions:python
python2 = sos.Python2.actions:python2
python3 = sos.Python3.actions:python3

R = sos.R.actions:R
Rmarkdown = sos.R.actions:Rmarkdown

docker_build = sos.docker.actions:docker_build [docker]
docker_commit = sos.docker.actions:docker_commit [docker]

[sos_taskengines]
process = sos.sos_task:BackgroundProcess_TaskEngine
rq = sos.rq.sos_task:RQ_TaskEngine [rq]
celery = sos.celery.sos_task:Celery_TaskEngine [celery]
pbs = sos.PBS.sos_task:PBS_TaskEngine

[sos_functions]
runfile = sos.jupyter.sos_executor:runfile


[sos_previewers]
*.pdf,1 = sos.jupyter.preview:preview_pdf
*.html,1 = sos.jupyter.preview:preview_html
*.csv,1 = sos.jupyter.preview:preview_csv
*.xls,1 = sos.jupyter.preview:preview_xls
*.xlsx,1 = sos.jupyter.preview:preview_xls
*.gz,1 = sos.jupyter.preview:preview_gz
*.txt,1 = sos.jupyter.preview:preview_txt
*.md,1 = sos.jupyter.preview:preview_md [md]
*.dot,1 = sos.jupyter.preview:preview_dot [dot]
imghdr:what,1 = sos.jupyter.preview:preview_img [image]
zipfile:is_zipfile,1 = sos.jupyter.preview:preview_zip
tarfile:is_tarfile,1 = sos.jupyter.preview:preview_tar
*,0 = sos.jupyter.preview:preview_txt

*.bam,1 = sos.bioinfo.preview:preview_bam [bam]

[sos_languages]
R = sos.R.kernel:sos_R [R]
Python2 = sos.Python2.kernel:sos_Python2
Python3 = sos.Python3.kernel:sos_Python3
Bash = sos.Bash.kernel:sos_Bash
JavaScript = sos.JavaScript.kernel:sos_JavaScript

[sos_converters]
sos-html.parser = sos.converter:get_script_to_html_parser
sos-html.func = sos.converter:script_to_html

sos-term.parser = sos.converter:get_script_to_term_parser
sos-term.func = sos.converter:script_to_term

sos-md.parser = sos.converter:get_script_to_markdown_parser
sos-md.func = sos.converter:script_to_markdown

sos-ipynb.parser = sos.jupyter.converter:get_script_to_notebook_parser
sos-ipynb.func = sos.jupyter.converter:script_to_notebook

ipynb-sos.parser = sos.jupyter.converter:get_notebook_to_script_parser
ipynb-sos.func = sos.jupyter.converter:notebook_to_script

ipynb-html.parser = sos.jupyter.converter:get_notebook_to_html_parser
ipynb-html.func = sos.jupyter.converter:notebook_to_html

ipynb-pdf.parser = sos.jupyter.converter:get_notebook_to_pdf_parser
ipynb-pdf.func = sos.jupyter.converter:notebook_to_pdf

ipynb-md.parser = sos.jupyter.converter:get_notebook_to_md_parser
ipynb-md.func = sos.jupyter.converter:notebook_to_md

''',
    extras_require = {
        ':sys_platform=="win32"': ['colorama'],
        'image':    ['wand'],
        'md':       ['markdown'],
        'R':        ['feather-format', 'pandas', 'numpy'],
        'rq':       ['rq', 'rq-dashboard'],
        'celery':   ['celery', 'flower'],
        'bam':      ['pysam'],
        'dot':      ['graphviz'],
        # docker-py is not working on windows 10 (as of Jan 2017)
        'docker':   ['docker'],
        'spyder':   ['spyder', 'jedi'],
    }
)
