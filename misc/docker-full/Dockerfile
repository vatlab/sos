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
#     docker build -t mdabioinfo/sos:latest docker-full 
#
# to build it.
#
FROM mdabioinfo/sos-base:latest

MAINTAINER Bo Peng <bpeng@mdanderson.org>

WORKDIR ${HOME}

# install the latest version of R-base because IRKernel supports only R 3.2
# or higher
RUN     echo "deb http://lib.stat.cmu.edu/R/CRAN/bin/linux/debian jessie-cran3/" >> /etc/apt/sources.list
RUN     apt-key adv --keyserver keys.gnupg.net --recv-key 6212B7B7931C4BB16280BA1306F90DE5381BA480
RUN     apt-get update
RUN     apt-get install -y  --force-yes libjpeg-dev r-base r-base-core r-recommended r-base-dev


RUN     apt-get install -y libzmq3-dev libcurl4-openssl-dev libssl-dev
RUN     Rscript -e "install.packages(c('rzmq', 'feather', 'devtools'), repos ='http://cran.us.r-project.org')" 
RUN     Rscript -e "devtools::install_github(paste0('IRkernel/', c('repr', 'IRdisplay', 'IRkernel')))"
RUN     Rscript -e "IRkernel::installspec()"

RUN     conda install -y feather-format -c conda-forge

RUN     mkdir -p workspace/notebooks workspace/data .ipython/profile_default

RUN     (echo "c = get_config()" && \
         echo "c.NotebookApp.ip = '0.0.0.0'" && \
         echo "c.NotebookApp.open_browser = False") \
         > $HOME/.ipython/profile_default/ipython_notebook_config.py

ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]

EXPOSE 80
CMD ["jupyter", "notebook", "--port=80","--no-browser"]

