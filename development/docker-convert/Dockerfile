# Copyright (c) Bo Peng and the University of Texas MD Anderson Cancer Center
# Distributed under the terms of the 3-clause BSD License.

# SoS official docker image for latest version of SoS. Use command
# 
#     docker build -t vatlab/sos-convert:latest docker-convert
#
# to build it.
#
FROM vatlab/sos:latest

MAINTAINER Bo Peng <Bo.Peng@bcm.edu>

RUN     pip install sos-notebook

ENTRYPOINT ["sos", "convert"]

