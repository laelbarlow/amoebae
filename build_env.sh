#!/usr/bin/env bash
# This script is for generating a Docker image and container using a Docker
# file in the current directory which should be the amoebae git repository
# directory. 

SECONDS=0

# Clone the git repo.
#echo \nCloning git repo... &&\
#git clone https://laelbarlow@bitbucket.org/laelbarlow/amoebae-lite.git &&\
# Change directories into the cloned repo.
#echo \nChanging directories into git repo... &&\
#cd amoebae-lite &&\
# Prompt user to add the amoebae-lite path to the files for sharing via the
# Docker Desktop App GUI.
printf "Add this directory path to the Docker files that can be shared and
mounted to Docker containers via the Docker graphical user interface:\n" &&\
pwd &&\
printf "\nPress enter to continue" &&\
read EMPTYVAR &&\
# Change the owner of the files to be mounted to the container on your host machine.
#echo Changing permissions for current directory ${PWD}/notebooks... &&\
#sudo chown -R 1000:100 ${PWD}/notebooks &&\
# Build a docker image using the Dockerfile file in the amoebae-lite repo.
printf "\nBuilding Docker image from Dockerfile..." &&\
docker build -t amoebaedockerimage:1.0 .

ELAPSED="Generation of your Docker image took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"

echo $ELAPSED
