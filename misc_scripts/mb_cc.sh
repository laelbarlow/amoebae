#!/bin/bash
# This script writes script with which to run phylogenetic analyses with
# MrBayes for all .mb.nex (nexus alignment) files in the current working
# directory. These scripts must then each be submitted via the sbatch command.
#
# For example:
#     module load mrbayes/3.2.6
#     for VAR in *_mb.sh; do sbatch $VAR; done

# Load MrBayes (for some reason it doesn't load properly if you load it for each output script).
#module load iq-tree/1.5.5

for FILE in *mb.nex
do 
NEWDIR=$FILE'_cedar'
mkdir $NEWDIR
NEWFILE=$NEWDIR'/'$FILE
cp $FILE $NEWFILE
TEXT=$"#!/bin/bash
#SBATCH --ntasks=8              # number of MPI processes (mb should be faster with 16, but for some reason when 16 is specified the job times out without starting on cedar).
#SBATCH --mem-per-cpu=125000M   # maximum on cedar?
#SBATCH --time=24:00:00
#SBATCH --account=def-dacks
#SBATCH --mail-user=lael@ualberta.ca
#SBATCH --mail-type=END

cd $NEWDIR
mpirun -np \$SLURM_NTASKS mb $FILE > $FILE'_log.txt' 
"
printf "$TEXT" > $FILE'_mb.sh' 
done

printf "
Now make sure that the IQ-tree module is loaded, and submit each output script
via sbatch.

For example:
    module load mrbayes/3.2.6
    for VAR in *_mb.sh; do sbatch \$VAR; done
    
"

