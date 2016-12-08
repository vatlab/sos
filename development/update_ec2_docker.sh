#!/bin/bash

docker rm $(docker stop $(docker ps -a -q --filter="name=tmp"))
docker stop proxy
docker rm proxy

cd SOS
git pull
docker build -t mdabioinfo/sos-notebook:latest -f /home/ubuntu/SOS/development/docker-notebook/Dockerfile .

export TOKEN=$( head -c 30 /dev/urandom | xxd -p )
docker run --net=host -d -e CONFIGPROXY_AUTH_TOKEN=$TOKEN --name=proxy jupyter/configurable-http-proxy --default-target http://127.0.0.1:9999
docker run --net=host -d -e CONFIGPROXY_AUTH_TOKEN=$TOKEN --name tmp_sos -v /var/run/docker.sock:/docker.sock jupyter/tmpnb python orchestrate.py --image='mdabioinfo/sos-notebook' --command="jupyter notebook --NotebookApp.base_url={base_path} --ip=0.0.0.0 --port {port}"
