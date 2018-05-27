# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
# 
#     docker build -t mdabioinfo/sos:latest docker-base 
#
# to build it.
#
FROM continuumio/anaconda3

MAINTAINER Bo Peng <bpeng@mdanderson.org>

RUN     mkdir /home/mdabioinfo
ENV     HOME  /home/mdabioinfo
WORKDIR ${HOME}

RUN     apt-get update
RUN     apt-get install -y  gcc
# these should be installed automatically by sos
RUN     pip install nbformat --upgrade
RUN     pip install sos

