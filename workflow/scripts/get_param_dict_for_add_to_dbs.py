#!/bin/env python3
"""Module defining functions for parsing information about how to extract
sequence IDs from FASTA headers in each FASTA file listed in the
resources/genomes.csv file.
"""
import sys
import os

def get_param_dict_for_add_to_dbs(input_csv_path):
    """Parse input CSV file to define a dictionary with desination filenames as
    keys and tuples of parameter strings as values.
    """
    param_dict = {}
    with open(input_csv_path) as infh:
        # Iterate over lines in CSV file.
        lnum = 0
        for i in infh:
            lnum += 1
            # Ignore irrelevant lines.
            if not i.startswith('Filename') \
                    and not i.startswith('\n') \
                    and not i.startswith(','):
                # Parse info in line.
                spliti = i.strip().split(',')
                filename = str(spliti[0])
                delimiter = str(spliti[1])
                seq_id_pos = str(spliti[2])

                # Check that the filename extension is acceptable.
                assert filename.rsplit('.', 1)[1] in ['faa', 'fna', 'gff3'],\
                """Error: Input file %s line %s: Data filename must have extension
                faa, fna, or gff3 not %s.""" % (input_csv_path, str(lnum), filename)

                # Check that the filename was not listed previously.
                assert filename not in param_dict.keys(), """Error: Data
                filename %s is listed more than once in CSV file %s.""" % \
                (filename, input_csv_path)

                # If delimiter is a space character, then put a backslash in
                # front of it.
                if delimiter == ' ':
                    delimiter = '\ '

                # Add info to parameter dictionary.
                param_dict[filename] = (delimiter, seq_id_pos)

    # Return the dict.
    return param_dict
