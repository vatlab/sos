# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
# 
#     docker build -t mdabioinfo/sos-notebook:latest docker-notebook
#
# to build it.
#

FROM mdabioinfo/sos-notebook:latest

MAINTAINER Bo Peng <bpeng@mdanderson.org>

USER    root

RUN     R --slave -e 'source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")'
RUN     R --slave -e 'install.packages("glmnet", repos="http://cran.us.r-project.org")'
RUN     pip install sklearn

USER    jovyan

#       SPARQL kernel for testing
RUN     cd /tmp && git clone https://github.com/asanchez75/sparql-kernel.git && cd /tmp/sparql-kernel && python setup.py install
RUN     jupyter sparqlkernel install --user

# for testing
RUN     pip install xlsx2csv bs4 xlrd

COPY    examples /home/jovyan/examples

USER    root
RUN     chown -R jovyan /home/jovyan/examples
RUN     chown -R jovyan /home/jovyan/work
RUN     chown -R jovyan /home/jovyan/.local

USER    jovyan
