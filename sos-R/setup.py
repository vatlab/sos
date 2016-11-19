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

from setuptools import find_packages, setup

# obtain version of SoS-R
with open('sos_R/_version.py') as version:
    exec(version.read())

setup(name = "sos-R",
    version = __version__,
    description = 'A module that provides support for the R-language for Script of Scripts (SoS)',
    author = 'Bo Peng',
    url = 'https://github.com/BoPeng/SOS/sos-bioinfo',
    author_email = 'bpeng@mdanderson.org',
    maintainer = 'Bo Peng',
    maintainer_email = 'bpeng@mdanderson.org',
    license = 'GPL3',
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
    packages = find_packages(exclude=['test']),
    entry_points={
		'sos_languages': [
			'R = sos_R.kernel:sos_R [R]',
		],
        'sos_targets': [
            'R_library = sos_R.target:R_library'
        ],
        'sos_actions': [
            'R = sos_R.actions:R',
            'Rmarkdown = sos_R.actions:Rmarkdown',
        ],
	},
    extras_require = {
        # feather is used to convert between R and Python dataframes
        'R': ['feather', 'pandas', 'numpy']
    }
)
