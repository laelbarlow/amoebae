#!/bin/bash
# This script writes a script with which to run AMOEBAE commands with. The
# output script must then be submitted via the sbatch command.

# Define input amoebae command to be written to output script file.
AMOEBAECOMMAND="$@"

# Define a timestamp function
timestamp() {
    date +"%Y%m%d%H%M%S"
}

TEXT=$"#!/bin/bash
#SBATCH --ntasks=1              # number of MPI processes (mb should be faster with 16, but for some reason when 16 is specified the job times out without starting on cedar).
#SBATCH --mem-per-cpu=125000M   # maximum on cedar?
#SBATCH --time=3:00:00
#SBATCH --account=def-dacks
#SBATCH --mail-user=lael@ualberta.ca
#SBATCH --mail-type=END

# Load python.
#module load python/3.7.0
# Load older version so that it will work with PyQt
module load python/3.6.3

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
#module load qt/5.11.3
# Load an older version of qt which actually works (with python 3.6.3), because
# version 5.11.3 has a bug.
module load qt/5.10.1

# Make a virtual environment for python and install python modules.
virtualenv --no-download \$SLURM_TMPDIR/env
source \$SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
pip install biopython --no-index
#pip install gffutils
pip install PyPDF2 
pip install reportlab

# Install dependencies of ete3
pip install six

# Install ete3.
pip install ete3 --no-index

### Run amoebae commands to check whether all the subprocess calls and import statements will work.
##printf \"\n\n\nChecking dependencies of amoebae:\n\"
##amoebae check_depend
##printf \"\n\n\nChecking all import statements for amoebae:\n\"
##amoebae check_imports

#*** AMOEBAE command here:

amoebae $AMOEBAECOMMAND

#***

# Shut down the virtual environment.
deactivate

"

# Define a path to an output script.
TIMESTAMP=$(date +%Y%m%d%H%M%S)
OUTPATH=cc_amoebae_$TIMESTAMP.sh

# Write script text to output path.
printf "$TEXT" > $OUTPATH 

printf "
Now submit the job with the sbatch command.

For example:
    sbatch $OUTPATH
    
"

