#!/usr/bin/env bash

# Script for setting up amoebae_example_workflow.

# Start timing this script.
SECONDS=0

# Set to exit immediately if any commands return an error.
set -e

# Detect job scheduler.
source scripts/determine_snakemake_profile.sh

# Determine whether packages can be loaded using the 'module' command.
if test "$(command -v module)"; then
    # Use the 'module' command to load Python module.

    # Load appropriate version of Python.
    module load python/3.7

else
    # Print the version.
    echo ""
    echo Using this version of Python3:
    python3 --version

fi

# Check that Python can be used.
if test "$(command -v python3)"; then
    echo ""
else
    echo Unable to access Python3.
    exit 1
fi

# Define path to python environment directory.
PYENVDIR=~/env_amoebae_workflow_setup_py

# Set up a python environment with cookiecutter.
if [ ! -d "$PYENVDIR" ]; then
  #virtualenv --no-download $PYENVDIR 
  python3 -m venv $PYENVDIR 
  source $PYENVDIR/bin/activate
  pip install --upgrade pip
  pip install cookiecutter
  pip install -U pytest
else
  source $PYENVDIR/bin/activate
fi


if [ "$(uname)" == "Darwin" ]; then
    printf "\nDetected MacOS."
    # Do something under Mac OS X platform
    echo "This software has not been tested on MacOS, as it has been designed to run on linux clusters."
    exit 1

elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    printf "\nDetected Linux.\n\n"
    # Do something under GNU/Linux platform
 
    # Determine what job scheduler is in use, and set up the appropriate
    # SnakeMake profile.
    if [ "$snakemake_profile" == "slurm" ]; then
        printf "\nSetting up config files for running SnakeMake with the SLURM job scheduler.\n\n"
        # Install configuration files for running on a SLURM cluster
        # (https://github.com/Snakemake-Profiles/slurm).
        if [ ! -d ~/.config/snakemake/slurm ]; then
            mkdir -p ~/.config/snakemake
            cd ~/.config/snakemake
            cookiecutter https://github.com/Snakemake-Profiles/slurm.git
            cd -
            # Now you can run snakemake as: snakemake --profile slurm ...

            # Copy amoebae_cluster_config.yaml example file from resources directory to
            # profile files.
            cp resources/example_slurm_amoebae_cluster_config.yaml \
            ~/.config/snakemake/slurm/amoebae_cluster_config.yaml

            # Update value of the amoebae_cluster_config variable in the submission
            # python file with the name of the amoebae_cluster_config.yaml file.
            sed -i 's/amoebae_cluster_config = \"\"/amoebae_cluster_config = \"amoebae_cluster_config.yaml\"/g' \
                ~/.config/snakemake/slurm/slurm-submit.py 
        fi

    elif [ "$snakemake_profile" == "sge" ]; then
        printf "\nSetting up config files for running SnakeMake with the (SUN/Univa) Grid Engine job scheduler.\n\n"

        # Install configuration files for running on a cluster with the Grid
        # Engine scheduler.
        if [ ! -d ~/.config/snakemake/sge ]; then
            mkdir -p ~/.config/snakemake
            cd ~/.config/snakemake
            # Make profile directory using cookiecutter to install from GitHub.
            cookiecutter https://github.com/Snakemake-Profiles/sge.git
            # Make Python scripts in profile executable.
            chmod a+x ~/.config/snakemake/sge/*.py
            # Change back to previous directory.
            cd -
            # Now you can run snakemake as: snakemake --profile sge ...

        fi

    elif [ "$snakemake_profile" == "pbs-torque" ]; then
        # Install configuration files for running on a cluster with the
        # PBS-TORQUE scheduler.
        if [ ! -d ~/.config/snakemake/pbs-torque ]; then
            printf "\nSetting up config files for running SnakeMake with the PBS-TORQUE job scheduler.\n\n"

            mkdir -p ~/.config/snakemake
            cd ~/.config/snakemake
            # Make profile directory using cookiecutter to install from GitHub.
            cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git
            # Make Python scripts in profile executable.
            chmod a+x ~/.config/snakemake/pbs-torque/*.py
            # Change back to previous directory.
            cd -
            # Now you can run snakemake as: snakemake --profile sge ...

        else
            printf "\nA profile already exists for running SnakeMake with the PBS-TORQUE job scheduler.\n\n"

        fi

    else
        printf "\nUnable to automatically set up config files for submitting
        jobs through the job scheduler on this system. You must do this
        manually.\n\n"
    fi


elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    printf "\nDetected Windows."
    # Do something under 32 bits Windows NT platform
    echo "This software has not been tested on Windows, as it has been designed to run on linux clusters."
    exit 1

elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    printf "\nDetected Windows."
    # Do something under 64 bits Windows NT platform
    echo "This software has not been tested on Windows, as it has been designed to run on linux clusters."
    exit 1
fi



#######################################################

# Deactivate and remove Python environment.
deactivate
rm -rf $PYENVDIR


#######################################################
# Report how much time it took for this script to run.
echo ""
ELAPSED="Setting up this workflow took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED


