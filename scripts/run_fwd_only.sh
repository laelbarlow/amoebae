#!/usr/bin/env bash

# Determine what Snakemake profile to use.
source scripts/determine_snakemake_profile.sh

# Activate python virtual environment.
source scripts/workflow_python_env_definition.sh

# Run snakemake in python virtual environment.
snakemake get_fwd_hit_seqs -j 100 --use-singularity --profile $snakemake_profile

# Deactivate python virtual environment.
source scripts/deactivate_python_env.sh
