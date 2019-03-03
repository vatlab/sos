# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
#
#     docker build -t mdabioinfo/sos-notebook:latest docker-notebook
#
# to build it.
#

# tag created in Fev 2019
FROM jupyter/datascience-notebook:83ed2c63671f

MAINTAINER Bo Peng <bpeng@mdanderson.org>

USER    root

#       Tools
RUN     apt-get update
RUN     apt-get install -y graphviz
RUN     apt-get install -y texlive-xetex texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended

RUN     apt-get install -y octave
RUN     octave --eval 'pkg install -forge dataframe'

RUN     apt-get install -y nodejs npm

RUN     apt-get install -y libgmp3-dev
RUN     apt-get install -y software-properties-common

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

#        Julia
RUN     julia -e "using Pkg;Pkg.add([\"Feather\", \"DataFrames\", \"NamedArrays\"])"

#       Python 2
RUN     conda create -n ipykernel_py2 python=2 ipykernel
RUN     /bin/bash -c "source activate ipykernel_py2; python -m ipykernel install --user; source deactivate"

#       Markdown kernel
RUN     pip install markdown-kernel
RUN     python -m markdown_kernel.install

# Bioinfo
RUN     pip install pysam

# SOS
RUN     pip install docker markdown wand graphviz imageio pillow nbformat

RUN     conda install -y feather-format -c conda-forge

## trigger rerun for sos updates
ARG	    DUMMY=unknown
RUN     DUMMY=${DUMMY} pip install sos sos-notebook sos-r sos-julia sos-python sos-matlab sos-javascript sos-bash sos-bioinfo --upgrade

RUN     python -m sos_notebook.install

# install the alpha version of jupyter lab to demonstrate jupyterlab-sos
RUN     pip install jupyterlab==1.0.0a1
RUN     jupyter labextension install transient-display-data@0.2.1
RUN     jupyter labextension install jupyterlab-sos@0.4.2

USER    jovyan
