#!/bin/bash
# This script writes a script with which to run notebooks on computecanada
# servers. This may be similar for other clusters that use the SLURM job
# scheduler. The output script must then be submitted via the sbatch command to
# run the notebook.

# Usage:
#     write_notebook_slurm_script.sh <notebook file name (.ipynb)>

# Define the filename (including extension) of the jupyter notebook to be run.
NBFILE="$1"

if [ ! ${NBFILE: -6} == ".ipynb" ];then
    echo Error: Input file name must have the extension .ipynb
    exit 1
fi

# Define a timestamp function
timestamp() {
    date +"%Y%m%d%H%M%S"
}

TEXT=$"#!/bin/bash
#SBATCH --ntasks=1              
#SBATCH --mem-per-cpu=8000M   
#SBATCH --time=0:05:00 # More time will be required!
#SBATCH --account=def-dacks # Update this!
#SBATCH --mail-user=lael@ualberta.ca # Update this too!
#SBATCH --mail-type=END

## Change directories into the amoebae directory (parent of current notebooks
## directory).
#cd ..

# Import singularity.
module load singularity/3.5

# Add current directory (the main amoebae directory) to the \$PATH so that the
# amoebae script can be accessed.
export PATH=\"\$PATH:\$\( dirname \$PWD\)\"

# Run the notebook.
singularity exec -B /home -B /project -B /scratch -B /localscratch -B \$\( dirname \$PWD\) \
    ../singularity.sif jupyter nbconvert \
        --to notebook \
        --allow-errors \
        --inplace \
        --execute \
        --ExecutePreprocessor.timeout=None \
        $NBFILE


#***
"

# Define a path to an output script.
TIMESTAMP=$(date +%Y%m%d%H%M%S)
OUTPATH=${NBFILE}_$TIMESTAMP.sh

# Write script text to output path.
printf "$TEXT" > $OUTPATH 

printf "
Note: Depending on what software you are using, you may need to specify 2
threads in the output .sh script before submitting (set \"--ntasks=1\" to
\"--ntasks=2\"). Also, check that the email address is correct.

Submit the job with the sbatch command.

For example:
    sbatch $OUTPATH
    
"

