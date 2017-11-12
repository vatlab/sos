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
with open('src/sos/_version.py') as version:
    for line in version:
        if line.startswith('__version__'):
            __version__ = eval(line.split('=')[1])
            break


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
        log.info('\nSoS is installed and configured to use with vim.')
        log.info('Use "set syntax=sos" to enable syntax highlighting.')

description = '''\
Computationally intensive disciplines such as computational biology often 
requires one to exploit a variety of tools implemented in different programming
languages, and to analyze large datasets on high performance computing systems.
Although scientific workflow systems are powerful in organizing and executing
large-scale data analysis processes, there are usually non-trivial learning
curve and engineering overhead in creating and maintaining such workflows,
making them unsuitable for data exploration and prototyping. To bridge the
gap between interactive analysis and workflow systems, we developed Script
of Scripts (SoS), a system with strong emphases on readability, practicality,
and reproducibility for daily computational research. For exploratory analysis
SoS provides a multi-language file format and scripting engine that centralizes
all computations, and creates dynamic report documents for publishing and
sharing. As a workflow engine, SoS provides an intuitive syntax to create
workflows in process-oriented, outcome-oriented and mixed styles, as well as
a unified interface to executing and managing tasks on a variety of computing
platforms with automatic synchronization of files between isolated systems.
In this paper we illustrate with real-world examples the use of SoS as both
interactive analysis tool and pipeline platform for all stages of methods
development and data analysis projects. In particular we demonstrate how SoS
can easily be adopted based on existing scripts and pipelines, yet resulting
in substantial improvement in terms of organization, readability and
cross-platform computation management.

Please refer to http://vatlab.github.io/SOS/ for more details on SoS.
'''

setup(name = "sos",
    version = __version__,
    description = 'Script of Scripts (SoS): an interactive, cross-platform, and cross-language workflow system for reproducible data analysis',
    long_description = description,
    author = 'Bo Peng',
    url = 'https://github.com/vatlab/SoS',
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
        'Operating System :: Microsoft :: Windows',
        'Intended Audience :: Information Technology',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
        ],
    packages = find_packages('src'),
    package_dir = {'': 'src'},
    # install vim syntax highlighting and other stuff
    cmdclass={'install': InstallWithConfigurations},
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
          'docker;platform_system!="Windows"',
      ],
    entry_points= '''
[console_scripts]
sos = sos.__main__:main
sos-runner = sos.__main__:sosrunner


[pygments.lexers]
sos = sos.converter:SoS_Lexer

[sos_targets]
file_target = sos.targets:file_target
dynamic = sos.targets:dynamic
remote = sos.targets:remote
executable = sos.targets:executable
sos_variable = sos.targets:sos_variable
sos_step = sos.targets:sos_step
env_variable = sos.targets:env_variable
sos_targets = sos.targets:sos_targets

[sos_actions]
script = sos.actions:script
sos_run = sos.actions:sos_run
fail_if = sos.actions:fail_if
warn_if = sos.actions:warn_if
stop_if = sos.actions:stop_if
download = sos.actions:download
run = sos.actions:run

perl = sos.actions:perl
ruby = sos.actions:ruby
report = sos.actions:report
pandoc = sos.actions:pandoc

docker_build = sos.docker.actions:docker_build

[sos_taskengines]
process = sos.sos_task:BackgroundProcess_TaskEngine

[sos_converters]
sos-html.parser = sos.converter:get_script_to_html_parser
sos-html.func = sos.converter:script_to_html

sos-term.parser = sos.converter:get_script_to_term_parser
sos-term.func = sos.converter:script_to_term

sos-md.parser = sos.converter:get_script_to_markdown_parser
sos-md.func = sos.converter:script_to_markdown
''',
    extras_require = {
        ':sys_platform=="win32"': ['colorama'],
        'dot':      ['graphviz'],
    }
)
