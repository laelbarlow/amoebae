#!/usr/bin/env bash

# Determine what Snakemake profile to use.
source scripts/determine_snakemake_profile.sh

# Activate python virtual environment.
source scripts/workflow_python_env_definition.sh

# Run snakemake in python virtual environment.
#snakemake plot_workflow -j 100 --use-singularity --profile $snakemake_profile \
python3 --version && \
snakemake --version && \
snakemake -n  --snakefile workflow/Snakefile && \
snakemake plot_workflow -j 100 --use-conda --profile $snakemake_profile \
--snakefile workflow/Snakefile --verbose

# Deactivate python virtual environment.
source scripts/deactivate_python_env.sh
