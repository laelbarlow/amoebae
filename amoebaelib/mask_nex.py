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
"""Mask a nex alignment by adding a new taxon, "MASK" sequence.

This is only a rough mask that is applied based on very simple criteria.

It would be good to eventually add more sophistocated criteria, especially
criteria based on similarity (according to a scoring matrix) rather than just
identity, or number of sequences without gaps.

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
from afa_to_nex import delete_extra_mesquite_lines
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections


def apply_mask_criteria(column):
    """Apply simple masking criteria to a single column, and return '-' if the
    column does not meet the criteria, and 'I' if it does.
    """
    # Return '-' by default.
    mask_char = '-'

    # Get column features.
    num_seqs = len(column)
    half_num_seqs = num_seqs / 2
    num_gaps_in_col = column.count('-')
    column_no_gaps = column.replace('-', '')

    # Check that the column is not entirely composed of gaps.
    #assert not column_no_gaps == '', "Error: Empty positions in input alignment."

    if column_no_gaps == '':
        return mask_char

    elif not column_no_gaps == '':

        most_common_residue = collections.Counter(column_no_gaps).most_common(1)[0]
        most_common_residue_count = most_common_residue[1]
        percent_identity = most_common_residue_count * 100 / num_seqs 

        # If less than half the sequences have a gap at this position of the
        # alignment, then include the position. 
        #if num_gaps_in_col < half_num_seqs:
        #    mask_char = 'I'

        if num_gaps_in_col < (num_seqs * 0.30):
            mask_char = 'I'

        # If percent identity is at least 50, then include position.
        if percent_identity >= 50:
            mask_char = 'I'

        return mask_char


def mask_alignment(alignment):
    """Takes an alignment object and adds a 'MASK' sequence using certain
    criteria for inclusion of positions.
    """
    # Get length of sequences and number of sequences in alignment.
    seq_len = alignment.get_alignment_length()
    num_seqs = len(alignment)

    # get a list of columns as strings in the original alignment.
    columns = [alignment[:, col] for col in range(seq_len)] 

    # Iterate over columns and make a mask sequence to append.
    mask_seq = ''
    for col in columns:
        mask_char = apply_mask_criteria(col) 
        mask_seq = mask_seq + mask_char

    # Generate a SeqRecord object with the mask_seq.
    #empty_mask_seq = Seq(mask_seq, IUPAC.protein)
    empty_mask_seq = Seq(mask_seq)
    empty_mask_rec = SeqRecord(empty_mask_seq, id='MASK', name='MASK')

    # Add the mask_seq sequence to the alignment.
    masked_alignment = alignment
    masked_alignment.append(empty_mask_rec)

    return masked_alignment


def mask_nex(infilepath, outfilepath=None):
    """Takes a filepath and adds a MASK sequence.
    """
    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(infilepath)

    # Define the name of the output file.
    if outfilepath is None:
        outfilepath = infilepath.replace('.nex', '.mask.nex')

    with open(infilepath) as infh, open(outfilepath, 'w') as o:
        # Check that the input file has the filename extension ".nex".
        assert infilepath.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')

        masked_alignment = mask_alignment(alignment)
        AlignIO.write(masked_alignment, o, 'nexus')

