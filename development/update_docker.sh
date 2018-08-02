#/bin/bash
docker stop $(docker ps -a -q)
docker rm $(docker ps -a -q)
cd ~/SoS
git pull
cd development
docker build --build-arg DUMMY=`date +%s` docker-notebook -t mdabioinfo/sos-notebook:latest
docker rmi $(docker images | grep "^<none>" | awk '{print $3}')
export TOKEN=$( head -c 30 /dev/urandom | xxd -p )
docker run --net=host -d -e CONFIGPROXY_AUTH_TOKEN=$TOKEN --name=proxy    jupyterhub/configurable-http-proxy --default-target http://127.0.0.1:9999
docker run --net=host -d -e CONFIGPROXY_AUTH_TOKEN=$TOKEN --name=tmp_sos -v /var/run/docker.sock:/docker.sock jupyter/tmpnb python orchestrate.py --pool-size=25 --cull-timeout=900 --cull-max=3600 --image='mdabioinfo/sos-notebook' --container-user=jovyan --command='jupyter notebook --no-browser --port {port} --ip=0.0.0.0 --NotebookApp.base_url={base_path} --NotebookApp.port_retries=0 --NotebookApp.token="" --NotebookApp.disable_check_xsrf=True'
