#/bin/bash

# Download FASTA files (and GFF3 annotation files) for genomes, transcriptomes,
# etc.

# Start timing this script.
SECONDS=0

# Parse command line arguments.
input_csv_path=$1
temp_data_dir=$2

# Make directory to contain genome files.
mkdir $temp_data_dir

# Parse CSV file in resources directory, which lists files to download.
OLDIFS=$IFS
IFS=','
[ ! -f $input_csv_path ] && { echo "$input_csv_path file not found"; exit 99; }
while read filename delimiter position compression location blank
do
    if [ ! "$filename" == "Filename" ]; then
    if [ ! "$filename" == "" ]; then

        # Report identified CSV contents.
        echo "Name : $filename"
        echo "Header delimiter : $delimiter"
        echo "Sequence ID position : $position"
        echo "Compression type : $compression"
        echo "Location : $location"

        # If no location is specified, then assume that the file is already
        # present in a local directory.
        if [ "$location" == "" ]; then
            :
        else
            # Define destination path for downloaded file.
            destination=$temp_data_dir/$filename

            # Download file and decompress, if necessary.
            if [ "$compression" == "None" ]; then
                curl $location --output $destination
            elif [ "$compression" == "gzip" ]; then
                #curl "${location%?}" --output $destination'.gz' && \
                curl $location --output $destination'.gz' && \
                gunzip $destination'.gz'
            else
                printf "\nCould not decompress file.\n"
                exit 3
            fi
        fi
    fi
    fi

done < $1
IFS=$OLDIFS


#######################################################
# Report how much time it took for this script to run.
echo ""
ELAPSED="Downloading FASTA files took the following amount of time: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

