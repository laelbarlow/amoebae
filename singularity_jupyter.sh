#/bin/bash

# This script is for using an existing .sif file ("singularity.sif") to run
# AMOEBAE code from a linux machine. The path to the notebook directory will
# change depending on whether this script is being run from inside a virtual
# machine or not. 

# Set directory in which to run jupyter.
NOTEBOOKDIR=notebooks
if [ ! -f "amoebae" ]; then
    printf "\tRunning singularity from within a virtual machine.\n\n"
    NOTEBOOKDIR=~/amoebae/notebooks
fi

# Run jupyter using the singularity container.
singularity exec singularity.sif \
    jupyter notebook --ip=0.0.0.0 \
                     --port=8888 \
                     --notebook-dir=${NOTEBOOKDIR} \
                     --allow-root \
                     --no-browser
