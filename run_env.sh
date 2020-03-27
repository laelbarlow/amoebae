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

# Initiate optional command line options.
a_flag=''
b_flag=''
files=''
verbose='false'

# Define usage statement.
print_usage() {
  printf "Usage: 

This script is for running a docker container using a docker image generated
by using the setup.sh script. It loads the amoebae directory as a volume in
the container so that the amoebae directory can be accessed from both your
host machine and from within the docker container. Also, this script
optionally loads one additional directory as a volume in the container.

Always run from within the main directory of a clone of the amoebae
repository.

To load only the amoebae directory, and open the notebooks subdirectory with
jupyter:

    bash run_env.sh

To load the amoebae directory as well as a second directory with jupyter
notebooks using the -d option:

    bash run_env.sh -d /full/path/to/second/directory

To access the container via the command line instead of jupyter, use the -a
option (this can be combined with the -d option):

    bash run_env.sh -a

To exit the container, log out of the jupyter interface and type CTRL-C in the
terminal window, or, if logged in via the command line, then enter exit to exit
the container.

***When you exit the container, all files written to directories other than
the loaded directories will be erased.
  
"
}

# Handle optional command line options
while getopts 'abd:v' flag; do
  case "${flag}" in
    a) a_flag='true' ;;
    b) b_flag='true' ;;
    d) dirpath="${OPTARG}" ;;
    v) verbose='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

# Define name of docker image to use for running a container. This is the same
# name as the image generated using the setup.sh script.
DI=amoebaedockerimage:1.0
printf "\nDocker image used for running container: "$DI"\n"

# Define directories on host machine to load as volumes to docker container.
printf "\nDirectories loaded as volumes in container: " &&\
V1=${PWD}
V1N=/opt/amoebae
printf "\n\t$V1  loaded as  $V1N" &&\
V2=${dirpath:-''}
USEV2=' '
V2N=''
if test $V2; then 
USEV2=''
V2N=/opt/${basename -- "$V2"}
printf "\n\t$V2  loaded as  $V2N\n"
fi

# Generate a docker container from the docker image that was just built.
printf "\n\nGenerating a container from the Docker image...\n\n" &&\

if test $a_flag; then

docker run --rm -it \
    -v $(PWD):/opt/amoebae \
    ${USEV2:--v $V2:$V2N} \
    $DI /bin/bash

else

docker run --rm -p 8888:8888 \
    -v ${PWD}:/opt/amoebae \
    ${USEV2:--v $V2:$V2N} \
    $DI /bin/bash -c "/opt/conda/bin/jupyter notebook --ip=0.0.0.0 --port=8888 --notebook-dir=${V2N:-/opt/amoebae/notebooks} --allow-root --no-browser"

fi


