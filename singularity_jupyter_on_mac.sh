#/bin/bash

# This script is for using an existing .sif file ("singularity.sif") to run
# AMOEBAE code from a MacOS machine. This involves spinning up a Vagrant
# virtual machine running Ubuntu (because singularity only runs on Linux).

# Spin up a Vagrant VM.
vagrant up && \
# Copy files from the host to the VM.
vagrant ssh -c 'cp amoebae/singularity.sif singularity.sif' && \
vagrant ssh -c 'cp amoebae/singularity_jupyter.sh singularity_jupyter.sh' && \
# Run script to start Jupyter.
vagrant ssh -c '/bin/bash singularity_jupyter.sh' && \
# Remove the copy of the .sif file.
vagrant ssh -c 'rm singularity.sif' && \
# Shut down VM.
vagrant halt
# Erase VM (including copy of .sif file).
#vagrant destroy


