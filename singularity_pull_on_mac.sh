#/bin/bash

# This script is for downloading a pre-built singularity container on a MacOS
# machine. This involves spinning up a Vagrant virtual machine running Ubuntu
# (because singularity only runs on Linux).

# Spin up a Vagrant VM.
vagrant up && \
# Run the singularity pull script.
vagrant ssh -c '/bin/bash singularity_pull.sh' && \
# Move the .sif file into the directory synced with the host. This will
# overwrite any existing file with the same name in that directory.
vagrant ssh -c 'mv singularity.sif amoebae/singularity.sif' && \
# Shut down VM.
vagrant halt
# Erase VM.
#vagrant destroy
