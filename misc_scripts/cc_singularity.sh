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
#SBATCH --mem-per-cpu=16000M   
#SBATCH --time=0:10:00
#SBATCH --account=def-dacks
#SBATCH --mail-user=lael@ualberta.ca
#SBATCH --mail-type=END

# Import singularity.
module load singularity/3.5

# Run amoebae commands to check whether all the subprocess calls and import statements will work.
printf \"\n\n\nChecking dependencies of amoebae:\n\"
singularity exec -B /home -B /project -B /scratch -B /localscratch singularity.sif echo \$PATH

singularity exec -B /home -B /project -B /scratch -B /localscratch singularity.sif /home/lael/projects/def-dacks/lael/amoebae/amoebae check_depend

printf \"\n\n\nChecking all import statements for amoebae:\n\"
singularity exec -B /home -B /project -B /scratch -B /localscratch singularity.sif /home/lael/projects/def-dacks/lael/amoebae/amoebae check_imports

printf \"\n\n\nRunning amoebae command:\n\"
#*** AMOEBAE command here:

singularity exec -B /home -B /project -B /scratch -B /localscratch singularity.sif /home/lael/projects/def-dacks/lael/amoebae/amoebae $AMOEBAECOMMAND

#***
"

# Define a path to an output script.
TIMESTAMP=$(date +%Y%m%d%H%M%S)
OUTPATH=cc_singularity_$TIMESTAMP.sh

# Write script text to output path.
printf "$TEXT" > $OUTPATH 

printf "
Note: If you are running something that uses iqtree for tree searching, then
you may need to specify 2 threads in the output .sh script before submitting.

Submit the job with the sbatch command.

For example:
    sbatch $OUTPATH
    
"

