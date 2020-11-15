#!/bin/bash
# properties = {properties}

# Load modules.
module load  \
nixpkgs/16.09 \
gcc/7.3.0 \
mesa/18.3.6 \
blast+/2.9.0 \
hmmer/3.2.1 \
muscle/3.8.31 \
exonerate/2.4.0

## Set up a python virtual environment with necessary python modules.
#virtualenv --no-download $SLURM_TMPDIR/env
#source $SLURM_TMPDIR/env/bin/activate
#pip install --no-index --upgrade pip
#pip install --no-index \
#pandas==1.1.3 \
#matplotlib==3.3.2 \
#biopython==1.77 \
#gffutils==0.9 \
#graphviz==0.13.2 \
#snakemake==5.2.4

# Execute job.
{exec_job}
