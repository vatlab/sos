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

# SoS official docker image for latest version of SoS. Use command
# 
#     docker build -t mdabioinfo/sos-notebook:latest docker-notebook
#
# to build it.
#
FROM jupyter/datascience-notebook

MAINTAINER Bo Peng <bpeng@mdanderson.org>

USER    root

RUN     conda install -y feather-format -c conda-forge
# RUN     conda install -y spyder

RUN     pip install pip --upgrade
RUN     pip install xlrd docker
RUN     pip install markdown

RUN     pip install bash_kernel
RUN     python -m bash_kernel.install

RUN     pip install sos

COPY    doc /home/$NB_USER/work/doc
RUN     chown -R $NB_USER /home/jovyan/work

USER    $NB_USER
EXPOSE	8888


