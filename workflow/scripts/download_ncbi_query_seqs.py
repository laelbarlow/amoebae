#!/bin/env python3
"""Script for downloading FASTA sequences defined in a CSV file.

Usage:
    download_ncbi_query_seqs.py \
        /path/to/csv_file/listing/queries \
        /path/to/single_fasta/destination/directory \
        /path/to/multi_fasta/destination/directory

"""
import sys
import os
import requests
import time
from get_filenames_from_csv import get_query_file_dict_from_csv

# Parse command line arguments.
cmdln = sys.argv
input_csv_path = cmdln[1]
seq_destination_dir_path = cmdln[2]
ali_destination_dir_path = cmdln[3]


############################################################
############################################################
# Parse input CSV file to define a dictionary with desination filenames as
# values and NCBI sequence accessions keys.
query_dict = get_query_file_dict_from_csv(input_csv_path)

# Filter out query filenames for which there is no NCBI accession, because
# these query files do not need to be downloaded as they are already present in
# a local directory.
query_dict = {k: v for k, v in query_dict.items() if v != []}


############################################################
############################################################
# Download single-sequence queries.

# Make a new temporary directory to store sequence files.
temp_query_dir_name = seq_destination_dir_path
if not os.path.isdir(temp_query_dir_name):
    os.mkdir(temp_query_dir_name)

# Loop over keys in the query_dict dictionary.
for filename in query_dict.keys():
    accessions = query_dict[filename]
    # Only download single-sequence queries.
    if len(accessions) == 1:
        filepath = os.path.join(temp_query_dir_name, filename)
        # Only download sequences that have not already been downloaded.
        if not os.path.isfile(filepath):
            try:
                # Download the sequence from the NCBI Protein database.
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=' + ','.join(accessions)
                r = requests.get(url)
                with open(filepath, 'w') as o:
                    o.write(r.text)
            except:
                print("\nError: Could not access NCBI server.\n")
        # Check that the sequence was actually downloaded.
        assert os.path.isfile(filepath), """The sequence with the following accession could not be downloaded from NCBI: %s\n
        Try re-running this script.""" % accession
    # Wait one second between requests so that the API does not throw an error
    # and just write an error message to the output FASTA file paths.
    time.sleep(1)


############################################################
############################################################
# Download multiple-sequence queries.

# Make a new temporary directory to store sequence files.
temp_alignment_dir_name = ali_destination_dir_path
if not os.path.isdir(temp_alignment_dir_name):
    os.mkdir(temp_alignment_dir_name)

# Loop over keys in the query_dict dictionary.
for filename in query_dict.keys():
    accessions = query_dict[filename]
    # Only download multi-sequence queries.
    if len(accessions) > 1:
        filepath = os.path.join(temp_alignment_dir_name, filename)
        # Only download sequences that have not already been downloaded.
        if not os.path.isfile(filepath):
            try:
                # Download the sequence from the NCBI Protein database.
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=' + ','.join(accessions)
                r = requests.get(url)
                with open(filepath, 'w') as o:
                    o.write(r.text)
            except:
                print("\nError: Could not access NCBI server.\n")
        # Check that the sequence was actually downloaded.
        assert os.path.isfile(filepath), """The sequence with the following accession could not be downloaded from NCBI: %s\n
        Try re-running this script.""" % accession
        # Check that the file written is in FASTA format.
        first_char = None
        with open(filepath, 'r') as ofh:
            first_line = ofh.readline()
            first_char = first_line[0]
        if not first_char == '>':
            print("\nError: The first character of the resulting file is not >.\n")
            # Remove output file (Snakemake may try to re-run the script, and
            # it might work then).
            os.remove(filepath)

    # Wait one second between requests so that the API does not throw an error
    # and just write an error message to the output FASTA file paths.
    time.sleep(1)

