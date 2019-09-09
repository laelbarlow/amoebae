#!/usr/bin/env python3
# Copyright 2018 Lael D. Barlow
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
"""
For generating text to use as taxon names for figures generated using boxshade.
Export output as pdf, then copy and paste into the boxshade output using Adobe 
illustrator.

Usage:
    get_residues_numbers_for_position.py <nexus alignment file path> <number of
    the position that alignment figure starts at>

"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from module_afa_to_nex import delete_extra_mesquite_lines
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
import module_dacksify_pos_hmmer_hits
import settings

command_line_list = sys.argv
infp = str(command_line_list[1])
position_num = str(command_line_list[2])


def get_residue_numbers_for_position(infp, position_num):
    """Takes a filepath and number for position of interest, and writes a file
    with residues numbers for all the taxa at that position.
    """
    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(infp)

    # Define the name of the output file.
    outfilename = infp.replace('.nex', '_taxon_names.txt')

    with open(infp) as infh, open(outfilename, 'w') as o:
        # Check that the input file has the filename extension ".nex".
        assert infp.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')

        # Get length of sequences and number of sequences in alignment.
        seq_len = alignment.get_alignment_length()
        num_seqs = len(alignment)

        # get a list of columns as strings in the original alignment.
        columns = [alignment[:, col] for col in range(seq_len)] 

        # get a list of records from the alignment.
        records = []
        for record in alignment:
            records.append(record)

        # Iterate over residues in column, and get the residue number.
        res_nums = ''
        tax_num = 0
        pos_index = int(position_num) -1
        print('\nPosition ' + position_num + ' of alignment\n' + '*')
        for r in columns[pos_index]:
            tax_num += 1

            # Get formatted taxon name.
            rec_id = records[tax_num -1].id
            #split = rec_id.strip().split(' ')
            #name = split[1].split('_')[0]
            #unabbrev_name = module_dacksify_pos_hmmer_hits.find_unabbrev_name(name, tablefilepath)
            #taxon_name = unabbrev_name + ' ' + split[1].split('_')[1] + ' ' + split[0]
            taxon_name = rec_id

            # Get the residue number for r.
            full_seq = str(records[tax_num -1].seq)
            seq_up_to_and_incl_r = full_seq[0:pos_index + 1]
            seq_up_to_and_incl_r_no_gaps = seq_up_to_and_incl_r.replace('-', '')
            res_num = len(seq_up_to_and_incl_r_no_gaps)

            # Print relevant info for checking that it's working.
            print(full_seq[pos_index:pos_index + 11])

            # Write taxon name and residue number to output file (justified
            # text).
            string = "{:<43s}{:>0d}".format(taxon_name, res_num) # Adjust if alignment is off.
            o.write(string + '\n')


if __name__ == '__main__':
    get_residue_numbers_for_position(infp, position_num)
