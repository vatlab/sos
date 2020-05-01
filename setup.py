#!/usr/bin/env python
#
# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

import os
import sys

from setuptools import find_packages, setup
from setuptools.command.bdist_egg import bdist_egg

_py_ver = sys.version_info
if _py_ver.major == 2 or (_py_ver.major == 3 and
                          (_py_ver.minor, _py_ver.micro) < (6, 0)):
    raise SystemError(
        'sos requires Python 3.6 or higher. Please upgrade your Python {}.{}.{}.'
        .format(_py_ver.major, _py_ver.minor, _py_ver.micro))

# obtain version of SoS
with open('src/sos/_version.py') as version:
    for line in version:
        if line.startswith('__version__'):
            __version__ = eval(line.split('=')[1])
            break

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))


def get_long_description():
    with open(os.path.join(CURRENT_DIR, "README.md"), "r") as ld_file:
        return ld_file.read()


class bdist_egg_disabled(bdist_egg):
    """Disabled version of bdist_egg

    Prevents setup.py install performing setuptools' default easy_install,
    which it should never ever do.
    """

    def run(self):
        sys.exit(
            "Aborting implicit building of eggs. Use `pip install -U --upgrade-strategy only-if-needed .` to install from source."
        )


cmdclass = {
    'bdist_egg': bdist_egg if 'bdist_egg' in sys.argv else bdist_egg_disabled
}

setup(
    name="sos",
    version=__version__,
    description='Script of Scripts (SoS): an interactive, cross-platform, and cross-language workflow system for reproducible data analysis',
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    author='Bo Peng',
    url='https://github.com/vatlab/SoS',
    author_email='Bo.Peng@bcm.edu',
    maintainer='Bo Peng',
    maintainer_email='Bo.Peng@bcm.edu',
    license='3-clause BSD',
    include_package_data=True,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    packages=find_packages('src'),
    cmdclass=cmdclass,
    package_dir={'': 'src'},
    python_requires='>=3.6',
    install_requires=[
        'psutil',
        # progress bar
        'tqdm',
        # for file lock
        'fasteners',
        'pyyaml',
        'pygments',
        # for DAG, some version requires pydot, some requires pydotplus
        'networkx',
        'pydot',
        'pydotplus',
        'pexpect',
        # required by windows
        'ptyprocess',
        # for report regeneration
        'jinja2',
        # to execute workflow embedded in .ipynb files
        'nbformat',
        # zeromq for IPC
        'pyzmq',
    ],
    entry_points='''
[console_scripts]
sos = sos.__main__:main
sos-runner = sos.__main__:sosrunner


[pygments.lexers]
sos = sos.converter:SoS_Lexer

[sos_targets]
file_target = sos.targets:file_target
sos_tempfile = sos.targets:sos_tempfile
dynamic = sos.targets:dynamic
remote = sos.targets:remote
executable = sos.targets:executable
sos_variable = sos.targets:sos_variable
sos_step = sos.targets:sos_step
env_variable = sos.targets:env_variable
sos_targets = sos.targets:sos_targets
system_resource = sos.targets:system_resource

Py_Module = sos.targets_python:Py_Module
R_library = sos.targets_r:R_library

[sos_actions]
script = sos.actions:script
sos_run = sos.actions:sos_run
fail_if = sos.actions:fail_if
warn_if = sos.actions:warn_if
stop_if = sos.actions:stop_if
done_if = sos.actions:done_if
skip_if = sos.actions:skip_if
download = sos.actions:download
run = sos.actions:run

bash = sos.actions_bash:bash
csh = sos.actions_bash:csh
tcsh = sos.actions_bash:tcsh
zsh = sos.actions_bash:zsh
sh = sos.actions_bash:sh

node = sos.actions_javascript:node

julia = sos.actions_julia:julia

matlab = sos.actions_matlab:matlab
octave = sos.actions_matlab:octave

python = sos.actions_python:python
python2 = sos.actions_python:python2
python3 = sos.actions_python:python3

R = sos.actions_r:R
Rmarkdown = sos.actions_r:Rmarkdown

ruby = sos.actions_ruby:ruby

perl = sos.actions:perl
report = sos.actions:report
pandoc = sos.actions:pandoc

docker_build = sos.docker.actions:docker_build
singularity_build = sos.singularity.actions:singularity_build

[sos_taskengines]
process = sos.task_engines:BackgroundProcess_TaskEngine

[sos_workflowengines]
process = sos.workflow_engines:BackgroundProcess_WorkflowEngine

[sos_previewers]
*.pdf,1 = sos.preview:preview_pdf
*.html,1 = sos.preview:preview_html
*.csv,1 = sos.preview:preview_csv
*.xls,1 = sos.preview:preview_xls
*.xlsx,1 = sos.preview:preview_xls
*.gz,1 = sos.preview:preview_gz
*.txt,1 = sos.preview:preview_txt
*.md,1 = sos.preview:preview_md
*.dot,1 = sos.preview:preview_dot [dot]
*.svg,1 = sos.preview:preview_svg
imghdr:what,1 = sos.preview:preview_img
zipfile:is_zipfile,1 = sos.preview:preview_zip
tarfile:is_tarfile,1 = sos.preview:preview_tar
*,0 = sos.preview:preview_txt

[sos_converters]
sos-html = sos.converter:ScriptToHTMLConverter

''',
    # [sos_installers]
    # vim-syntax.parser = sos.install:get_install_vim_syntax_parser
    # vim-syntax.func = sos.install:install_vim_syntax
    extras_require={
        ':sys_platform=="win32"': ['colorama'],
        # faster hashlib
        ':sys_platform!="win32"': ['xxhash'],
        'dot': ['graphviz', 'pillow'],
    })
