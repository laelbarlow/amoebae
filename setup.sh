#!/usr/bin/env bash
# Copyright 2018 Lael D. Barlow
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 

# Script for setting up AMOEBAE to run in a Docker container after cloning the
# git repository. This is a bash script so it does not depend on python3 being
# installed (for a setup.py script to work) and it does not depend on developer
# tools being installed on MacOS (for make to work).


# Start timing this script.
SECONDS=0


#######################################################
# Copy the example settings file.
if [ ! -f "settings.py" ]; then
    printf "\tWriting settings file.\n\n"
    cp settings.py.example settings.py
fi


#######################################################
# Make scripts executable.
chmod a+x amoebae
chmod a+x misc_scripts/*.py
chmod a+x misc_scripts/*.sh


#######################################################
# Check whether Docker is installed or not.
printf "\nDetermining whether Docker needs to be installed or not.\n\n"

if test "$(command -v docker)"; then
    printf "\tDocker is installed.\n\n"
    
else
    printf "\tDocker is not installed.\n\n"

# If Docker is not installed, detect what operating system is being used, and
# install Docker.

if [ "$(uname)" == "Darwin" ]; then
printf "\nDetected MacOS."
# Do something under Mac OS X platform
# Prompt user to install Docker Desktop for MacOS manually.
printf "\nInstall Docker Desktop for MacOS from the website:\n" &&\
printf "\n\thttps://docs.docker.com/docker-for-mac/install/\n\n" &&\
read -p "\nDone? (press enter to continue)" &&\
# Prompt user to add the amoebae repo path to the files for sharing via the
# Docker Desktop App GUI.
printf "If you have not already done so, add this directory path to the Docker
files that can be shared and mounted to Docker containers via the Docker
graphical user interface:\n\n" &&\
pwd &&\
read -p "\nDone? (press enter to continue)" &&\
# Prompt user to customize the resources provided to docker for running
# containers via the Docker Desktop App GUI.
printf "If you have not already done so, customize settings for defining number
of CPUs and amount of memory, swap space, and storage made available to Docker
containers that will run on your machine, via the Docker graphical user
interface. These can be accessed in Preferences > Resources.\n\n" &&\
read -p "\nDone? (press enter to continue)"

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux."
    # Do something under GNU/Linux platform

    # Install docker.
    # This should work on ubuntu:
    #sudo apt-get install -y \
    #    apt-transport-https \
    #    ca-certificates \
    #    curl \
    #    software-properties-common curl -fsSL https://yum.dockerproject.org/gpg | sudo apt-key add - sudo add-apt-repository \
    #    "deb https://apt.dockerproject.org/repo/ \
    #    ubuntu-$(lsb_release -cs) \
    #    main" sudo apt-get update
    #sudo apt-get -y install docker-engine # add current user to docker group so there is no need to use sudo when running docker
    #sudo usermod -aG docker $(whoami)

    printf "\nProcedures for installing Docker on linux vary. Consult the
    Docker documentation: https://docs.docker.com/install/"

    # ***Note: It may be necessary to customize memory and CPU usage options:
    # https://docs.docker.com/config/containers/resource_constraints/

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    printf "\nDetected Windows."
    # Do something under 32 bits Windows NT platform
    echo "AMOEBAE has not been tested on Windows, please contact the developer."
    exit 1

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    printf "\nDetected Windows."
    # Do something under 64 bits Windows NT platform
    echo "AMOEBAE has not been tested on Windows, please contact the developer."
    exit 1
fi

fi


#######################################################
# Build Docker image from Dockerfile.
# This generates a Docker image and container using a Docker
# file in the current directory which should be the amoebae git repository
# directory. 

# Change the owner of the files to be mounted to the container on your host machine.
#echo Changing permissions for current directory ${PWD}/notebooks... &&\
#sudo chown -R 1000:100 ${PWD}/notebooks &&\
# Build a docker image using the Dockerfile file in the amoebae-lite repo.
printf "\nBuilding Docker image from Dockerfile..." &&\
docker build -t amoebaedockerimage:1.0 .


#######################################################
# Run test code using pytest. This runs all the tests defined in files in the
# tests subdirectory.
#printf "\nTesting AMOEBAE source code using PyTest..." &&\
#pytest


#######################################################
# Report how much time it took for this script to run.
echo ""
ELAPSED="Setting up AMOEBAE took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED


