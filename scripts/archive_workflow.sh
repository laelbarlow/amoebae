#!/usr/bin/env bash

# Activate python virtual environment.
source scripts/workflow_python_env_definition.sh

# Define a timestamp for the archive file.
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`

# Run snakemake in python virtual environment.
snakemake --archive amoebae_example_workflow_$TIMESTAMP.tar.gz --cores 1

# Deactivate python virtual environment.
source scripts/deactivate_python_env.sh

