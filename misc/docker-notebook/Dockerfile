# This file is part of Script of Scripts (sos), a workflow system
# for the execution of commands and scripts in different languages.
# Please visit https://github.com/bpeng2000/SOS for more information.
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

# SoS official docker image for latest version of SoS. Use command
# 
#     docker build -t mdabioinfo/sos-notebook:latest docker-notebook
#
# to build it.
#
FROM jupyter/datascience-notebook

MAINTAINER Bo Peng <bpeng@mdanderson.org>

RUN     mkdir /home/mdabioinfo
ENV     HOME  /home/mdabioinfo
WORKDIR ${HOME}

# let us also install some python and R/bioconductor packages

# because it is called a "FULL" container, let us install popular python
# and R packages
RUN     conda install -y numpy scipy pandas scikit-learn nltk
RUN     conda install -y -c https://conda.binstar.org/bpeng simuPOP

# R bioconductor packages that I would use from time to time
RUN     curl https://raw.githubusercontent.com/Bioconductor/bioc_docker/master/out/release_base/install.R -o /tmp/install.R
RUN     R -f /tmp/install.R

RUN     curl https://raw.githubusercontent.com/Bioconductor/bioc_docker/master/out/release_base/installFromBiocViews.R -o /tmp/installFromBiocViews.R
RUN     R -f /tmp/installFromBiocViews.R

RUN     curl https://raw.githubusercontent.com/Bioconductor/bioc_docker/master/out/release_core/installFromBiocViews.R -o /tmp/installCore.R
RUN     R -f /tmp/installCore.R

RUN     apt-get install -y libxml2 libxml2-dev
RUN     Rscript -e
"source('https://bioconductor.org/biocLite.R');biocLite('lumi')"
RUN     curl https://raw.githubusercontent.com/Bioconductor/bioc_docker/master/out/release_sequencing/installFromBiocViews.R -o /tmp/installSequencing.R
RUN     R -f /tmp/installSequencing.R

# these should be installed automatically by sos
RUN     conda install -y psutil
RUN     pip install sos


