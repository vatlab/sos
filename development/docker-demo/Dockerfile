# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
#
#     docker build -t vatlab/sos-notebook:latest docker-demo
#
# to build it.
#

FROM vatlab/sos-notebook:latest

MAINTAINER Bo Peng <Bo.Peng@bcm.edu>

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
