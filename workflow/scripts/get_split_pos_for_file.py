#!/bin/env python3
"""Script for retrieving the position/field number for identifying sequence IDs
in a FASTA sequence header after splitting on a particular character.
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

    # Get split position for filename.
    split_pos = add_to_dbs_params[os.path.basename(sequence_filename)][1]

    # Print split character.
    print(split_pos)

