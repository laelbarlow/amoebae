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

# Execute job.
{exec_job}
