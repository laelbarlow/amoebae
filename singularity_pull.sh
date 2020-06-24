#/bin/bash

# This script pulls the latest pre-built singularity container image from the
# singularity library
# (https://cloud.sylabs.io/library/_container/5e8ca8fff0f8eb90a8a7b60d).
# If this is run on MacOS, then it spins up a Vagrant virtual machine running
# Ubuntu (because singularity only runs on Linux).

# Start timing this script.
SECONDS=0

singularity_pull_in_vm () {
  # Spin up a Vagrant VM.
  vagrant up && \
  # Pull the latest container from the singularity library into the
  # singularity container.
  vagrant ssh -c 'singularity pull singularity.sif library://laelbarlow/default/amoebae:latest' && \
  # Move the .sif file into the directory synced with the host. This will
  # overwrite any existing file with the same name in that directory.
  vagrant ssh -c 'mv singularity.sif amoebae/singularity.sif' && \
  # Shut down VM.
  vagrant halt
  # Erase VM.
  #vagrant destroy
}


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS."
    # Do something under Mac OS X platform

    # Pull pre-built singularity container from within a virtual machine.
    singularity_pull_in_vm

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n"
    # Do something under GNU/Linux platform

    # Check whether singularity is installed, and if not then run singularity
    # in a virtual machine using vagrant.
    if test "$(command -v singularity)"; then
        printf "\nSingularity is installed.\n\n"

        # Pull the latest container from the singularity library.
        singularity pull singularity.sif library://laelbarlow/default/amoebae:latest
        
    else
        printf "\nSingularity is not installed.\n"
        printf "\nRunning singularity in a virtual machine.\n\n"

        # Run singularity pull command from within a virtual machine.
        singularity_pull_in_vm

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
ELAPSED="Downloading a Singularity container for AMOEBAE took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

