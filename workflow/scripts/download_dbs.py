#!/bin/env python3
"""
Script for downloading FASTA files (genomes, proteomes, transcriptomes)
and/or GFF3 annotation files listed in a CSV file, using Python's requests
library (rather than curl).

Usage:
    download_db.py \
        /path/to/csv_file/listing/databases \
        /path/to/destination/directory \
"""

import sys
import os
import certifi
import requests
import time
import gzip
import shutil
from get_filenames_from_csv import get_database_filenames_from_csv

# Parse command line arguments.
cmdln = sys.argv
if len(cmdln) < 3:
    sys.exit(
        "\nError: Missing arguments.\n\n"
        "Usage:\n"
        "    download_db.py "
        "/path/to/csv_file/listing/databases "
        "/path/to/destination/directory\n"
    )

input_csv_path = cmdln[1]
temp_data_dir = cmdln[2]

# Get a list of filenames for genomes/proteomes/transcriptome files
# that need to be downloaded from NCBI or another remote source.
database_names_ncbi = get_database_filenames_from_csv(input_csv_path, 'ncbi')

# Create destination directory if it does not exist.
if not os.path.isdir(temp_data_dir):
    os.mkdir(temp_data_dir)

# Open the CSV file and iterate over each line to identify files that require downloading.
with open(input_csv_path, "r") as infh:
    for line in infh:
        line = line.strip()
        # Skip header lines, empty lines, or lines starting with ','.
        if (not line) or line.startswith('Filename') or line.startswith(','):
            continue

        # Each line is assumed to have the following format:
        # filename,delimiter,position,compression,location
        s = line.split(',')
        filename = s[0]
        delimiter = s[1]
        position = s[2]
        compression = s[3]
        location = s[4]

        # If the file is listed among those needing download from NCBI (or remote source),
        # proceed with downloading.
        if filename in database_names_ncbi:
            # Verify that a remote location is specified.
            assert location.strip() != '', (
                f"\nError: No remote URL provided for file '{filename}' in CSV {input_csv_path}.\n"
            )

            destination = os.path.join(temp_data_dir, filename)

            # If there is no compression specified, download the file directly.
            if not compression:
                try:
                    print(f"Downloading {filename} to {destination} ...")
                    r = requests.get(location, stream=True, verify=certifi.where())
                    r.raise_for_status()

                    # Write file in chunks to handle large downloads
                    with open(destination, 'wb') as outfh:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                outfh.write(chunk)

                except requests.exceptions.RequestException as e:
                    print(f"\nError: Could not download {location}. Details: {e}\n")
                    continue

            # If the file is gzip-compressed, download to temporary .gz file, then decompress.
            elif compression.lower() == 'gzip':
                gz_destination = destination + '.gz'
                try:
                    print(f"Downloading and decompressing {filename} (gzip) to {destination} ...")
                    r = requests.get(location, stream=True, verify=certifi.where())
                    r.raise_for_status()

                    # Write gzipped file in chunks
                    with open(gz_destination, 'wb') as gz_outfh:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                gz_outfh.write(chunk)

                    # Decompress the .gz file
                    with gzip.open(gz_destination, 'rb') as f_in, open(destination, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                    # Remove the temporary .gz file
                    os.remove(gz_destination)

                except requests.exceptions.RequestException as e:
                    print(f"\nError: Could not download {location}. Details: {e}\n")
                    continue
                except OSError as oe:
                    print(f"\nError decompressing file {gz_destination}. Details: {oe}\n")
                    continue

            else:
                print("\nCould not decompress file because compression format is unsupported.\n")
                continue

            # Validate that the file was indeed downloaded
            assert os.path.isfile(destination), (
                f"Failed to download file {filename} from {location}. "
                "Check the URL or try again if there was a temporary server issue."
            )

            # Wait briefly between downloads to reduce the risk of server issues or rate limits
            time.sleep(1)

