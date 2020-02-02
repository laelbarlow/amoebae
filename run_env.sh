#!/usr/bin/env bash

# Generate a docker container from the docker image that was just built.
echo \nGenerating a container from the Docker image... &&\
#docker run -it --rm -p 8888:8888 -v ${PWD}/notebooks:/home/lael/notebooks amoebaeliteimage:latest &&\
#docker run -it --rm -p 8888:8888 amoebaeliteimage:latest &&\
#docker run --rm -p 8888:8888 -v $(PWD):/opt/notebooks continuumio/anaconda3 /bin/bash -c "/opt/conda/bin/jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir=/opt/notebooks --allow-root --no-browser"
docker run --rm -p 8888:8888 -v $(PWD):/opt/notebooks amoebaedockerimage:1.0 /bin/bash -c "/opt/conda/bin/jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir=/opt/notebooks --allow-root --no-browser"


