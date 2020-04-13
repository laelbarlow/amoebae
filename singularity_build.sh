#/bin/bash

# This script builds a singularity container image using singularity. If run on
# MacOS, it will run singularity within a virtual machine managed using
# Vagrant. Only run this script if you have updated or customized the
# singularity.recipe file. Also, this requires the use of the sudo command, so
# you won't be able to run this on linux without being an administrator.


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS.\n"
    # Do something under Mac OS X platform

    # Spin up a Vagrant VM.
    vagrant up && \
    # Copy recipe file from the host to the VM.
    vagrant ssh -c 'cp amoebae/singularity.recipe singularity.recipe' && \
    # Run the singularity build script.
    vagrant ssh -c 'sudo singularity build singularity.sif singularity.recipe' && \
    # Move the .sif file into the directory synced with the host. This will
    # overwrite any existing file with the same name in that directory.
    vagrant ssh -c 'mv singularity.sif amoebae/singularity.sif' && \
    # Shut down VM.
    vagrant halt
    # Erase VM.
    #vagrant destroy

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n"
    # Do something under GNU/Linux platform

    # Use singularity directly.
    sudo singularity build singularity.sif singularity.recipe


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
