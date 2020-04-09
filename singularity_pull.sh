#/bin/bash

# This script pulls a pre-built singularity container image from the library.
# This simply runs a singularity command, and it must be run on a Linux
# machine.

singularity pull singularity.sif library://laelbarlow/default/amoebae:0.0.0

