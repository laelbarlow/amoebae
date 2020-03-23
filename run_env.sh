#!/usr/bin/env bash

# Generate a docker container from the docker image that was just built.

## Works:
#printf "\nGenerating a container from the Docker image..." &&\
#docker run --rm -p 8888:8888 -v $(PWD):/opt/notebooks amoebaedockerimage:1.0 /bin/bash -c "/opt/conda/bin/jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir=/opt/notebooks/notebooks --allow-root --no-browser"

#printf "\nGenerating a container from the Docker image..." &&\
#docker run --rm -p 8888:8888 -v $(PWD):/amoebae amoebaedockerimage:1.0 /bin/bash -c "/opt/conda/bin/jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir=/amoebae/notebooks --allow-root --no-browser"

printf "\nGenerating a container from the Docker image..." &&\
docker run --rm -p 8888:8888 -v $(PWD):/opt/amoebae amoebaedockerimage:1.0 /bin/bash -c "/opt/conda/bin/jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir=/opt/amoebae/notebooks --allow-root --no-browser"
