#!/usr/bin/env bash

# Determine which environment management options to use.
env_options=""

if test "$(command -v conda)"; then
  env_options="--use-conda"

else
  env_options="--use-conda --use-singularity"

fi

