#!/usr/bin/env python3
"""
Script for downloading FASTA sequences defined in a CSV file using the NCBI E-utilities.

Usage:
    download_ncbi_query_seqs.py \
        /path/to/csv_file/listing/queries \
        /path/to/single_fasta/destination/directory \
        /path/to/multi_fasta/destination/directory

The CSV file is parsed to obtain a mapping of filenames to lists of NCBI protein
accessions. Sequences for single-accession queries are stored in the directory
specified for single_fasta; multi-accession queries are stored in the directory
specified for multi_fasta.

Requires:
    - Python 3
    - requests (pip install requests)
"""

import sys
import os
import certifi
import requests
import time

# Custom function for parsing the CSV file. Expects a dictionary where
# keys = filename, values = list of accession IDs.
from get_filenames_from_csv import get_query_file_dict_from_csv

# Parse command line arguments
cmdln = sys.argv
if len(cmdln) < 4:
    sys.exit(
        "\nError: Missing arguments.\n\n"
        "Usage:\n"
        "    download_ncbi_query_seqs.py "
        "/path/to/csv_file/listing/queries "
        "/path/to/single_fasta/destination/directory "
        "/path/to/multi_fasta/destination/directory\n"
    )

input_csv_path = cmdln[1]
seq_destination_dir_path = cmdln[2]
ali_destination_dir_path = cmdln[3]

# Fetch dictionary: { filename : [accession_1, accession_2, ...] }
query_dict = get_query_file_dict_from_csv(input_csv_path)

# Filter out any entries that do not contain NCBI accessions
query_dict = {fname: accs for fname, accs in query_dict.items() if accs}

# Define NCBI E-utilities URL and common parameters
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
COMMON_PARAMS = {
    "db": "protein",
    "rettype": "fasta",
    "retmode": "text"
}

########################################################################
# Download single-sequence queries
########################################################################

# Ensure destination directory for single sequences exists
temp_query_dir_name = seq_destination_dir_path
if not os.path.isdir(temp_query_dir_name):
    os.mkdir(temp_query_dir_name)

# Loop over queries. If there is exactly one accession, treat it as a single-sequence query
for filename, accessions in query_dict.items():
    if len(accessions) == 1:
        filepath = os.path.join(temp_query_dir_name, filename)
        if not os.path.isfile(filepath):
            try:
                # Construct parameters for request
                params = COMMON_PARAMS.copy()
                params["id"] = ",".join(accessions)

                response = requests.get(EFETCH_URL, params=params, verify=certifi.where())
                # Raise an error if the request was not successful (e.g., status code != 200)
                response.raise_for_status()

                # Write the response text to the output file
                with open(filepath, "w") as outfile:
                    outfile.write(response.text)

            except requests.exceptions.RequestException as e:
                print(f"\nError: Could not download {accessions[0]} from NCBI. Details: {e}\n")

        # Validate that the file was successfully created and appears to contain FASTA data
        assert os.path.isfile(filepath), (
            f"The sequence with the following accession could not be "
            f"downloaded from NCBI: {accessions[0]}\nTry re-running this script."
        )
        with open(filepath, "r") as infh:
            contents = infh.read(1024)  # read just a bit
            assert contents.startswith(">"), (
                f"The output file for {accessions[0]} does not appear to "
                f"contain a FASTA sequence.\nYou may need to try re-running this script."
            )

        # Delay between requests to respect NCBI usage guidelines
        time.sleep(1)

########################################################################
# Download multiple-sequence queries
########################################################################

# Ensure destination directory for multiple sequences (alignments) exists
temp_alignment_dir_name = ali_destination_dir_path
if not os.path.isdir(temp_alignment_dir_name):
    os.mkdir(temp_alignment_dir_name)

# Loop over queries. If there are multiple accessions, treat them as multi-sequence queries
for filename, accessions in query_dict.items():
    if len(accessions) > 1:
        filepath = os.path.join(temp_alignment_dir_name, filename)
        if not os.path.isfile(filepath):
            try:
                # Construct parameters for request
                params = COMMON_PARAMS.copy()
                params["id"] = ",".join(accessions)

                response = requests.get(EFETCH_URL, params=params, verify=certifi.where())
                response.raise_for_status()

                # Write the response text to the output file
                with open(filepath, "w") as outfile:
                    outfile.write(response.text)

            except requests.exceptions.RequestException as e:
                print(f"\nError: Could not download {','.join(accessions)} from NCBI. Details: {e}\n")

        # Validate that the file was successfully created and starts with '>'
        assert os.path.isfile(filepath), (
            f"The sequences with the following accessions could not be "
            f"downloaded from NCBI: {accessions}\nTry re-running this script."
        )
        with open(filepath, "r") as ofh:
            first_line = ofh.readline()
            if not first_line.startswith(">"):
                print(
                    f"\nError: The resulting file for accessions "
                    f"{','.join(accessions)} does not start with '>'.\n"
                )
                # Remove the file so the script can attempt re-download
                os.remove(filepath)

        # Delay between requests
        time.sleep(1)

