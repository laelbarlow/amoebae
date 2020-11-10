#!/usr/bin/env bash

# Deactivate python virtual environment.
if test "$(command -v conda)"; then
  # Deactivate Conda environment.
  conda deactivate
else
  # Assume a Python venv environment was used instead of Conda.
  deactivate
fi

