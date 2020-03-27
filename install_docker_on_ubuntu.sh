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

# This script is for installing docker on ubuntu, it may not work for all
# flavours of linux. Also, MacOS users will need to use the relevant desktop
# app for docker.

# Update the package manager.
sudo apt-get update &&\
# Install docker.
sudo apt install docker.io &&\
# Set docker to setup at startup.
sudo systemctl start docker &&\
sudo systemctl enable docker &&\
# Create a docker group.
#sudo groupadd docker &&\
# Add your user to the docker group.
sudo usermod -aG docker ${USER} &&\
# Prompt user to restart their computer.
printf "\n\nRestart your computer (or virtual machine) to begin using docker.\n\n
To test if docker works after restarting, enter the following command:\n\n
\tdocker run hello-world\n\n"
