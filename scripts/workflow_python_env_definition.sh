#!/usr/bin/env bash

# Script for building and accessing virtual environments for running snakemake.
# If you are running this on a cluster that uses the Conda package manager,
# make sure that you have activated a Conda environment using an appropriate
# version of Python.

# Determine snakemake profile name.
source scripts/determine_snakemake_profile.sh

# Determine whether packages can be loaded using the 'module' command.
if test "$(command -v module)"; then
    # Use the 'module' command to load Python and Singularity modules.

    # Load appropriate version of Python.
    # (tested using 3.7.4).
    module load python/3.7

    # Load an appropriate version of Singularity.
    module load singularity/3.6

else
    # Print the versions.
    echo ""
    echo Using this version of Python:
    python3 --version
    echo ""
    echo Using this version of Singularity:
    singularity --version

fi

# Check that Python can be used.
if test "$(command -v python3)"; then
    echo ""
else
    echo Unable to access Python3.
    exit 1
fi

# Check that Singularity can be used.
if test "$(command -v singularity)"; then
    echo ""
else
    echo Unable to access Singularity.
    #exit 1
fi

# Define path to directory to contain python virtual environment files
# (installation of dependencies)
DIR=~/env_amoebae_snakemake_workflow

# Build the virtual environment, if necessary.

if test "$(command -v conda)"; then
  echo Conda is available.

  # Install dependencies in virtual environment using conda.
  if ! test "$(conda info --envs | grep "conda_env_amoebae_snakemake_workflow")"; then
    #conda create --name conda_env_amoebae_snakemake_workflow python=3.7.4
    #source activate conda_env_amoebae_snakemake_workflow

    # Check that mamba is installed.
    if ! test "$(command -v mamba)"; then
        echo Installing mamba.
        # Install Mamba for making snakemake install easier.
        conda install -c conda-forge mamba
    fi

    echo ""
    echo Setting up environment for snakemake with mamba.
    echo ""


    # Create conda environment for snakemake with mamba.
    mamba create -c conda-forge -c bioconda \
        -n conda_env_amoebae_snakemake_workflow \
            snakemake \
            matplotlib \
            numpy \
            graphviz \
            requests \
            pulp \

    echo Environment set up. Re-run to activate.
    exit 1

  else
    # Otherwise, just activate the environment.
    echo Attempting to activate conda environment from script.

    # Get the path to the conda base environment.
    conda_base_env=$(conda info | grep -i 'base environment' | tr -s ' ' | cut -d' ' -f5)
    # Source conda from the base environment so that it can be used in a
    # subshell.
    source $conda_base_env/etc/profile.d/conda.sh
    # Activate the conda environment.
    conda activate conda_env_amoebae_snakemake_workflow

    # Add profile path to $PATH.
    export PATH=$PATH:~/.config/snakemake/$snakemake_profile

  fi

else
  # Install dependencies in virtual environment using pip.
  if [ ! -d "$DIR" ]; then
    python3 -m venv $DIR 
    source $DIR/bin/activate
    pip install --upgrade pip
    pip install setuptools wheel
    pip install requests
    pip install \
    numpy==1.18.1 \
    biopython==1.77 \
    gffutils==0.10.1 \
    pypdf2==1.26.0 \
    reportlab==3.5.51 \
    pandoc==2.0a4 \
    pandas==1.1.3 \
    matplotlib==3.3.2 \
    graphviz==0.14.2 \
    snakemake==5.26.1

  else
    # Otherwise, just activate the environment.
    source $DIR/bin/activate
  fi
fi


