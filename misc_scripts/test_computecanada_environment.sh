#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5000M
#SBATCH --time=00:05:00
#SBATCH --account=def-dacks
#SBATCH --mail-user=lael@ualberta.ca

# This script tests whether the environment (dependencies, modules, etc.) can
# be set up properly on a computecanada account. To use, just navigate to the
# amoebae/misc_scripts.py directory on the server (after cloning the amoebae
# repository) and run the script ( ./test_computecanada_environment.sh ).

# Load python.
module load python/3.7.0
# Load this because it has pandas and numpy.
module load scipy-stack 

# Load software to use with AMOEBAE.
module load nixpkgs/16.09
module load gcc/7.3.0
module load blast+/2.7.1
module load hmmer/3.1b2
module load iq-tree/1.5.5
#module load mrbayes/3.2.6
module load muscle/3.8.31

# Make a virtual environment for python and install python modules.
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
pip install biopython --no-index
pip install ete3 --no-index
pip install gffutils
pip install PyPDF2 
pip install reportlab


# Run amoebae commands to check whether all the subprocess calls and import statements will work.
amoebae check_depend
amoebae check_imports


