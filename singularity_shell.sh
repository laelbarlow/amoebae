#/bin/bash

# This script is for using an existing .sif file ("singularity.sif") to run
# code from a linux machine. If run on MacOS, this script will run singularity
# within a virtual machine managed by Vagrant. 


# Determine name of current working directory.
repo_basename="$( basename $PWD)"

# Define a function to run a singularity shell command within a virtual
# machine.
singularity_shell_in_vm () {
  # Spin up a Vagrant VM.
  vagrant up && \
  # Copy singularity container from the host to the VM.
  printf "\nCopying .sif file into virtual machine.\n" && \
  vagrant ssh -c 'cp '$repo_basename'/singularity.sif singularity.sif' && \
  # Run the command to open a shell session.
  printf "\nStarting shell session in Singularity container.\n\n\t***Type exit to exit the session when you are done.\n\n" 
  vagrant ssh -c 'cd '$repo_basename' && singularity shell -B /home/vagrant/'$repo_basename':/opt/'$repo_basename' ../singularity.sif' && \
  # Remove the copy of the .sif file.
  printf "\nRemoving copy of .sif file." && \
  vagrant ssh -c 'rm singularity.sif' && \
  # Shut down VM.
  vagrant halt
  # Erase VM (including copy of .sif file).
  #vagrant destroy
}


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS.\n"
    # Do something under Mac OS X platform

    # Run shell session using singularity from within a virtual machine.
    singularity_shell_in_vm

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n"
    # Do something under GNU/Linux platform

    # Check whether singularity is installed, and if not then run singularity
    # in a virtual machine using vagrant.
    if test "$(command -v singularity)"; then
        printf "\nSingularity is installed.\n\n"

        # Open a shell session in the singularity container.
        singularity shell \
            -B $PWD:/opt/$repo_basename \
            singularity.sif         

    else
        printf "\nSingularity is not installed.\n"
        printf "\nRunning singularity in a virtual machine.\n\n"

        # Run shell session with singularity from within a virtual machine.
        singularity_shell_in_vm
    fi

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    printf "\nDetected Windows.\n"
    # Do something under 32 bits Windows NT platform
    echo "This code has not been tested on Windows, please contact the developer."
    exit 1

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    printf "\nDetected Windows.\n"
    # Do something under 64 bits Windows NT platform
    echo "This code has not been tested on Windows, please contact the developer."
    exit 1
fi

