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

# Load modules that are dependencies for ete3.
module load qt/5.11.3

# Make a virtual environment for python and install python modules.
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
pip install biopython --no-index
pip install gffutils
pip install PyPDF2 
pip install reportlab

# Install dependencies of ete3
pip install six
pip install PyQt5

# Install ete3.
pip install ete3 --no-index

# Write a python script with an import statement. 
PYTHONSCRIPT="test_imports.py"
echo #/usr/bin/env python3 > $PYTHONSCRIPT
echo from PyQt5 import QtGui, QtCore >> $PYTHONSCRIPT
PYTHONOUTPUT="test_imports_output.txt"
python3 $PYTHONSCRIPT > $PYTHONOUTPUT

# Run amoebae commands to check whether all the subprocess calls and import statements will work.
#amoebae check_depend
#amoebae check_imports

# Shut down the virtual environment.
deactivate


