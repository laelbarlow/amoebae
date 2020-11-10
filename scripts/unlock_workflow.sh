#!/usr/bin/env bash

# Activate python virtual environment.
source scripts/workflow_python_env_definition.sh

# Run snakemake in python virtual environment.
snakemake --unlock --cores 1

# Deactivate python virtual environment.
source scripts/deactivate_python_env.sh
