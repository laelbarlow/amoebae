#/bin/bash

# This script is for using an existing .sif file ("singularity.sif") to run
# AMOEBAE code from a MacOS machine. This involves spinning up a Vagrant
# virtual machine running Ubuntu (because singularity only runs on Linux).

# Spin up a Vagrant VM.
vagrant up && \
# Copy recipe file from the host to the VM.
vagrant ssh -c 'cp amoebae/singularity.recipe singularity.recipe' && \
# Copy singularity build script from the host to the VM.
vagrant ssh -c 'cp amoebae/singularity_build.sh singularity_build.sh' && \
# Run the singularity build script.
vagrant ssh -c '/bin/bash singularity_build.sh' && \
# Move the .sif file into the directory synced with the host. This will
# overwrite any existing file with the same name in that directory.
vagrant ssh -c 'mv singularity.sif amoebae/singularity.sif' && \
# Shut down VM.
vagrant halt
# Erase VM.
#vagrant destroy
