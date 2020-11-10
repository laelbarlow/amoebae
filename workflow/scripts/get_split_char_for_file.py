#!/bin/env python3
"""Script for retrieving the text character to split FASTA headers on for
identifying sequence IDs.
"""
import sys
import os
from get_param_dict_for_add_to_dbs import get_param_dict_for_add_to_dbs

if __name__ == '__main__':
    # Parse command line arguments.
    cmdln = sys.argv
    database_csv = cmdln[1]
    sequence_filename = cmdln[2]

    # Get a dictionary of parameters for the add_to_dbs rule with filenames as keys
    # and parameter tuples as values.
    add_to_dbs_params = get_param_dict_for_add_to_dbs(database_csv)

    # Get split character for filename.
    split_char = add_to_dbs_params[os.path.basename(sequence_filename)][0]

    # Print split character.
    print(split_char)

