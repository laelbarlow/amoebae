#/bin/bash

# This script builds a singularity container image. This simply runs a
# singularity command, and it must be run on a Linux machine.

sudo singularity build singularity.sif singularity.recipe

