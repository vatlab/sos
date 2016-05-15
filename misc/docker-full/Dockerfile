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
#     docker build -t mdabioinfo/sos:latest docker-base 
#
# to build it.
#
FROM mdabioinfo/sos:latest

MAINTAINER Bo Peng <bpeng@mdanderson.org>

WORKDIR ${HOME}

# these should be installed automatically by sos
RUN     apt-get update
RUN     apt-get install -y R-base
RUN     apt-get install -y libzmq3-dev libcurl4-openssl-dev 
RUN     Rscript -e "install.packages(c('rzmq','repr','IRkernel','IRdisplay'), repos = c('http://irkernel.github.io/', 'http://cran.us.r-project.org'))" 
RUN     Rscript -e "IRkernel::installspec()"
		
RUN     mkdir -p workspace/notebooks workspace/data .ipython/profile_default

RUN     (echo "c = get_config()" && \
         echo "c.NotebookApp.ip = '0.0.0.0'" && \
         echo "c.NotebookApp.open_browser = False") \
         > $HOME/.ipython/profile_default/ipython_notebook_config.py

EXPOSE  8888

RUN     jupyter notebook &
