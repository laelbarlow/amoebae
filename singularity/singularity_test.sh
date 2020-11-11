#/bin/bash
#
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

# This script runs test code in the tests subdirectory using pytest (from
# within the singularity container. If singularity is not
# installed, then then it spins up a Vagrant virtual machine running Ubuntu
# (because singularity only runs on Linux). 

# Start timing this script.
SECONDS=0

# Determine name of current working directory.
repo_basename="$( basename $PWD)"


# Define a function to run a singularity exec command within a virtual
# machine.
singularity_test_in_vm () {
  # Spin up a Vagrant VM.
  vagrant up && \
  # Copy singularity container from the host to the VM.
  printf "\nCopying .sif file into virtual machine, if necessary.\n" && \
  vagrant ssh -c 'if [ ! -f "singularity.sif" ]; then cp '$repo_basename'/singularity.sif singularity.sif; fi' && \
  # Run the command to run PyTest.
  printf "\nRunning PyTest from within the Singularity container.\n" && \
  vagrant ssh -c 'cd '$repo_basename' && singularity exec -B /home/vagrant/'$repo_basename':/opt/'$repo_basename' ../singularity.sif pytest /opt/'$repo_basename'/tests' && \
  # Remove the copy of the .sif file.
  #printf "\nRemoving copy of .sif file." && \
  #vagrant ssh -c 'rm singularity.sif' && \
  # Shut down VM.
  vagrant halt
  # Erase VM (including copy of .sif file).
  #vagrant destroy
}


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS."
    # Do something under Mac OS X platform

    # Run singularity exec command from within a virtual machine.
    singularity_test_in_vm

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n"
    # Do something under GNU/Linux platform

    # Check whether singularity is installed, and if not then run singularity
    # in a virtual machine using vagrant.
    if test "$(command -v singularity)"; then
        printf "\nSingularity is installed.\n\n"

        singularity exec -B $repo_basename:/opt/$repo_basename singularity.sif pytest

    else
        printf "\nSingularity is not installed.\n"
        printf "\nRunning singularity in a virtual machine.\n\n"

        # Run singularity exec command from within a virtual machine.
        singularity_test_in_vm

    fi

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

#######################################################
# Report how much time it took for this script to run.
echo ""
ELAPSED="Testing AMOEBAE code from within a Singularity container took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

