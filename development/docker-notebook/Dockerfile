# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
#
#     docker build -t vatlab/sos-notebook:latest docker-notebook
#
# to build it.
#

# tag created in Fev 2019
FROM jupyter/datascience-notebook:r-3.6.3

MAINTAINER Bo Peng <Bo.Peng@bcm.edu>

USER    root

#       Tools
RUN     apt-get update
RUN     apt-get install -y graphviz zlib1g-dev libbz2-dev libcurl4-openssl-dev libssl-dev
RUN     apt-get install -y texlive-xetex texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended
RUN     apt-get install -y octave
RUN     octave --eval 'pkg install -forge dataframe'

RUN     apt-get install -y npm vim libgmp3-dev software-properties-common
RUN     apt-get install -y libtool libffi-dev ruby ruby-dev make  libzmq3-dev libczmq-dev

# Install some packages for our examples
RUN     conda install -c conda-forge -y r-arrow r-glmnet r-biocmanager
#  ruby
RUN     gem install ffi-rzmq
RUN     gem install iruby --pre
RUN     gem install daru nmatrix

RUN     iruby register --force

RUN     cd /home/jovyan; chown -R jovyan ../jovyan
USER    jovyan

RUN     pip install sklearn

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

RUN     conda install -y pyarrow -c conda-forge


# SOS
RUN     pip install docker markdown wand graphviz imageio pillow nbformat

## trigger rerun for sos updates
ARG	    DUMMY=unknown
RUN     DUMMY=${DUMMY} conda install -c conda-forge sos sos-notebook sos-r sos-julia sos-python sos-matlab sos-bash
RUN     pip install sos-ruby sos-javascript sos-bioinfo
RUN     conda install -c conda-forge jupyterlab-transient-display-data jupyterlab-sos
RUN     jupyter lab build --dev-build=False --minimize=False

#       Markdown kernel
RUN     pip install markdown-kernel
RUN     python -m markdown_kernel.install --prefix /opt/conda/

ENV     JUPYTER_ENABLE_LAB TRUE

USER    jovyan
