#/bin/bash

# This script is for using an existing .sif file ("singularity.sif") to run
# AMOEBAE code from a linux machine. If run on MacOS, this script will run
# singularity within a virtual machine managed by Vagrant. 

singularity_jupyter_in_vm () {
  # Spin up a Vagrant VM.
  vagrant up && \
  # Copy singularity container from the host to the VM.
  vagrant ssh -c 'cp amoebae/singularity.sif singularity.sif' && \
  # Run the command to start Jupyter.
  vagrant ssh -c 'singularity exec \
                       -B /home/vagrant/amoebae:/opt/amoebae \
                       singularity.sif \
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
}


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS.\n"
    # Do something under Mac OS X platform

    # Run jupyter using singularity from within a virtual machine.
    singularity_jupyter_in_vm

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n"
    # Do something under GNU/Linux platform

    # Check whether singularity is installed, and if not then run singularity
    # in a virtual machine using vagrant.
    if test "$(command -v singularity)"; then
        printf "\nSingularity is installed.\n\n"

        # Run jupyter using the singularity container.
        singularity exec \
            -B $PWD:/opt/amoebae \
            singularity.sif \
            jupyter notebook \
            --ip=0.0.0.0 \
            --port=8888 \
            --notebook-dir=notebooks \
            --allow-root \
            --no-browser
        
    else
        printf "\nSingularity is not installed.\n"
        printf "\nRunning singularity in a virtual machine.\n\n"

        # Run jupyter with singularity from within a virtual machine.
        singularity_jupyter_in_vm
    fi

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

