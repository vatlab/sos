# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
# 
#     docker build -t mdabioinfo/sos-notebook:latest docker-notebook
#
# to build it.
#
FROM jupyter/datascience-notebook:1085ca054a5f

MAINTAINER Bo Peng <bpeng@mdanderson.org>

USER    root

#       Tools
RUN     apt-get update
RUN     apt-get install -y graphviz
RUN     apt-get install -y texlive-xetex texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended

RUN     apt-get install -y octave
#RUN     octave --eval 'pkg install -forge dataframe'

RUN     apt-get purge --auto-remove nodejs npm node
RUN     rm -rf ~/.nvm 
RUN     apt-get install -y nodejs-legacy npm

RUN     apt-get install -y libgmp3-dev
RUN     apt-get install -y software-properties-common
RUN     add-apt-repository -y ppa:staticfloat/juliareleases
#RUN     add-apt-repository -y ppa:staticfloat/julia-deps
RUN     apt-get update
RUN     apt-get install -y julia

# Install some packages for our examples
#RUN     R --slave -e 'install.packages(c("package1", "package2"), lib="/usr/local/lib/R/site-library")'
RUN     R --slave -e 'source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")'
RUN     R --slave -e 'install.packages("glmnet", repos="http://cran.us.r-project.org")'
RUN     pip install sklearn

USER    jovyan

#       Bash
RUN     pip install bash_kernel
RUN     python -m bash_kernel.install --user

#       Octave
RUN     pip install octave_kernel
RUN     python -m octave_kernel install --user

#       JavaScript
RUN     npm rebuild
RUN     npm install -g ijavascript
RUN     ijsinstall --spec-path=full

#       SPARQL kernel for testing
RUN     cd /tmp && git clone https://github.com/asanchez75/sparql-kernel.git && cd /tmp/sparql-kernel && python setup.py install
RUN     jupyter sparqlkernel install --user

#        Julia
RUN     julia -e "ENV[\"JUPYTER\"]=\"$(which jupyter)\";Pkg.add(\"IJulia\")"
RUN     julia -e 'Pkg.add("Feather")'
RUN     julia -e 'Pkg.add("DataFrames")'
RUN     julia -e 'Pkg.add("NamedArrays")'

#       Python 2
RUN     conda create -n ipykernel_py2 python=2 ipykernel
RUN     /bin/bash -c "source activate ipykernel_py2; python -m ipykernel install --user; source deactivate"

# Bioinfo
RUN     pip install pysam

# for testing
RUN     pip install xlsx2csv bs4

# SOS
RUN     pip install pip --upgrade
RUN     pip install xlrd docker
RUN     pip install markdown wand graphviz imageio pillow 


RUN     conda install -y feather-format -c conda-forge
RUN     pip install nbformat --upgrade
## trigger rerun for sos updates
ARG	DUMMY=unknown
RUN     DUMMY=${DUMMY} pip install sos sos-notebook sos-r sos-julia sos-python sos-matlab sos-javascript sos-bash sos-bioinfo --upgrade
RUN     python -m sos_notebook.install
RUN     pip install jupyterlab
RUN     jupyter labextension install jupyterlab-sos

COPY    examples /home/jovyan/examples

USER    root
RUN     chown -R jovyan /home/jovyan/examples
RUN     chown -R jovyan /home/jovyan/work
RUN     chown -R jovyan /home/jovyan/.local

USER    jovyan

#EXPOSE	8888
