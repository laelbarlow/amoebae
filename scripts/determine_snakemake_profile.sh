#!/usr/bin/env bash

# Determine what string to use to define Snakemake profile in current
# environment.

snakemake_profile=""

if test "$(command -v sbatch)"; then
    echo SLURM job scheduler detected.
    snakemake_profile="slurm"

elif qstat -help; then
    printf "(SUN/Univa) Grid Engine job scheduler detected."
    snakemake_profile="sge"

elif pbs-config --version; then
    printf "PBS-TORQUE job scheduler detected."
    snakemake_profile="pbs-torque"

else
    echo Unable to detect a job scheduler.
    #exit 1

fi
