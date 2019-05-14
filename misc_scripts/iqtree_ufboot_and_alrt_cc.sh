#!/bin/bash
# This script writes script with which to run phylogenetic analyses with
# IQ-tree for all .phy (phylip alignment) files in the current working
# directory. These scripts must then each be submitted via the sbatch command.
#
# For example:
#     for VAR in *_iqtree_ufboot_and_alrt.sh; do sbatch $VAR; done

# Load IQ-tree (for some reason it doesn't load properly if you load it for each output script).
#module load iq-tree/1.5.5

for FILE in *.phy
do 
TEXT=$"#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100000M
#SBATCH --time=01:00:00
#SBATCH --account=def-dacks
#SBATCH --mail-user=lael@ualberta.ca
#SBATCH --mail-type=END

mkdir $FILE'_cedar'
iqtree -s $FILE -bb 1000 -wbt -alrt 1000 -m MFP -nt AUTO -pre $FILE'_cedar/output'
"
printf "$TEXT" > $FILE'_iqtree_ufboot_and_alrt.sh' 
done

printf "
Now make sure that the IQ-tree module is loaded, and submit each output script
via sbatch.

For example:
    module load iq-tree/1.5.5
    for VAR in *_iqtree_ufboot_and_alrt.sh; do sbatch \$VAR; done
    
"

