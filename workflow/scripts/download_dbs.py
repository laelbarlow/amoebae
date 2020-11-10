#!/bin/env python3
"""Script for downloading FASTA files (genomes, proteomes, transcriptomes)
and/or GFF3 annotation files listed in a CSV file.

Usage:
    download_db.py \
        /path/to/csv_file/listing/databases \
        /path/to/destination/directory \

"""
import sys
import os
import requests
import time
import subprocess
from get_filenames_from_csv import get_database_filenames_from_csv

# Parse command line arguments.
cmdln = sys.argv
input_csv_path = cmdln[1]
temp_data_dir = cmdln[2]

# Get a list of filenames for genomes/proteomes/transcriptome files to be
# downloaded.
database_names_ncbi = get_database_filenames_from_csv(input_csv_path, 'ncbi')

# Make directory to contain genome files.
if not os.path.isdir(temp_data_dir):
    os.mkdir(temp_data_dir)

# Parse CSV file in resources directory, which lists files to download.
with open(input_csv_path) as infh:
    for i in infh:
        if not i.startswith('Filename') and not i.startswith('\n') and not\
        i.startswith(','):

            s = i.strip().split(',')
            filename = s[0]
            delimiter = s[1]
            position = s[2]
            compression = s[3]
            location = s[4]

            # Report identified CSV contents.
            #print('filename')
            #print(filename)
            #print('delimiter')
            #print(delimiter)
            #print('position')
            #print(position)
            #print('compression')
            #print(compression)
            #print('location')
            #print(location)

            # If no location is specified, then assume that the file is already
            # present in a local directory.
            #if location.strip() != "":
            if filename in database_names_ncbi:
                
                # Check that a location is provided for where to download from.
                assert location.strip() != '', """Error: No remote location
                provided for file %s in CSV sheet %s.""" % (filename, input_csv_path)

                # Define destination path for downloaded file.
                destination = os.path.join(temp_data_dir, filename)

                # Download file and decompress, if necessary.
                if compression is None:
                    subprocess.call(['curl', location, '--output', destination])
                elif compression == 'gzip':
                    subprocess.call(['curl', location, '--output', destination + '.gz']) 
                    subprocess.call(['gunzip', destination + '.gz'])
                else:
                    print("\nCould not decompress file.")






