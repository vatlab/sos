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

#       Tools
RUN     apt-get update
RUN     apt-get install -y graphviz
RUN     apt-get install -y texlive-xetex texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended

#       Bash
RUN     pip install bash_kernel
RUN     python -m bash_kernel.install

#       Octave
RUN     apt-get install -y octave
RUN     octave --eval 'pkg install -forge dataframe'
RUN     pip install octave_kernel
RUN     python -m octave_kernel.install

#       JavaScript
RUN     apt-get purge --auto-remove nodejs npm node
RUN     rm -rf ~/.nvm 
RUN     apt-get install -y nodejs-legacy npm
RUN     npm install -g ijavascript
RUN     npm rebuild
RUN     ijsinstall --spec-path=full

#        Julia
RUN     apt-get install -y libgmp3-dev
#RUN     apt-get install -y software-properties-common
#RUN     add-apt-repository -y ppa:staticfloat/juliareleases
#RUN     add-apt-repository -y ppa:staticfloat/julia-deps
#RUN     apt-get update
RUN     apt-get install -y julia
RUN     julia -e "ENV[\"JUPYTER\"]=\"$(which jupyter)\";Pkg.add(\"IJulia\")"
RUN     julia -e 'Pkg.add("Feather")'
RUN     julia -e 'Pkg.add("DataFrames")'
RUN     julia -e 'Pkg.add("NamedArrays")'

#       Python 2
RUN     conda create -n ipykernel_py2 python=2 ipykernel
RUN     /bin/bash -c "source activate ipykernel_py2"
RUN     python -m ipykernel install --user
RUN     /bin/bash -c "source deactivate ipykernel_py2"

# Bioinfo
RUN     pip install pysam

# SOS
RUN     pip install pip --upgrade
RUN     pip install xlrd docker
RUN     pip install markdown wand graphviz

RUN     conda install -y feather-format -c conda-forge
RUN     pip install nbformat --upgrade
RUN     pip install sos sos-notebook sos-r sos-julia sos-python sos-matlab sos-javascript sos-bash sos-bioinfo

#COPY    docs /home/$NB_USER/work/doc
RUN     chown -R $NB_USER /home/jovyan/work
RUN     chown -R $NB_USER /home/jovyan/.local

USER    $NB_USER
EXPOSE	8888


