#/bin/bash

# This script is for using an existing .sif file ("singularity.sif") to run
# AMOEBAE code from a linux machine. If run on MacOS, this script will run
# singularity within a virtual machine managed by Vagrant. 


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS.\n"
    # Do something under Mac OS X platform

    # Spin up a Vagrant VM.
    vagrant up && \
    # Copy singularity container from the host to the VM.
    vagrant ssh -c 'cp amoebae/singularity.sif singularity.sif' && \
    # Run the command to start Jupyter.
    vagrant ssh -c 'singularity exec singularity.sif \
                         jupyter notebook \
                         --ip=0.0.0.0 \
                         --port=8888 \
                         --notebook-dir=~/amoebae/notebooks \
                         --allow-root \
                         --no-browser' && \
    # Remove the copy of the .sif file.
    vagrant ssh -c 'rm singularity.sif' && \
    # Shut down VM.
    vagrant halt
    # Erase VM (including copy of .sif file).
    #vagrant destroy

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n"
    # Do something under GNU/Linux platform

    # Run jupyter using the singularity container.
    singularity exec singularity.sif \
        jupyter notebook --ip=0.0.0.0 \
                         --port=8888 \
                         --notebook-dir=notebooks \
                         --allow-root \
                         --no-browser


elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    printf "\nDetected Windows.\n"
    # Do something under 32 bits Windows NT platform
    echo "AMOEBAE has not been tested on Windows, please contact the developer."
    exit 1

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    printf "\nDetected Windows.\n"
    # Do something under 64 bits Windows NT platform
    echo "AMOEBAE has not been tested on Windows, please contact the developer."
    exit 1
fi
