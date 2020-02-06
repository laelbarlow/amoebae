#!/usr/bin/env bash

# Generate a docker container from the docker image that was just built and
# open a bash session.
printf "\nGenerating a container from the Docker image..." &&\
docker run -it -v $(PWD):/opt/notebooks amoebaedockerimage:1.0 /bin/bash


