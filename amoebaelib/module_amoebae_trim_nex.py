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
"""Trim an nex alignment using the last taxon "MASK" sequence as a guide.

Usage:
    trim_nex.py <path to nex alignment with MASK>
"""
import sys
import os
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from module_afa_to_nex import delete_extra_mesquite_lines

## Increase recursion limit for python. Potentially dangerous, but necessary for
## large alignments, like for Giantin and lava lamp proteins!
#sys.setrecursionlimit(10000)

# Recursive function to trim alignment.
# Base case is when there are no '-' characters in the MASK (last sequence in
# the alignment).
def trim_alignment(alignment):
    """Recursively trims an alignment by removing positions with 'I' in the
    MASK sequence.
    """
    # Increase recursion limit for python. Potentially dangerous, but necessary for
    # large alignments, like for Giantin and lava lamp proteins!
    sys.setrecursionlimit(10000)

    # Define mask row.
    mask_seq = alignment[-1].seq

    # If all the positions to trim have been removed, then return the input
    # alignment.
    if '-' not in mask_seq:
        return alignment[:-1]

    # If not, then remove the first position that was not masked for inclusion. 
    else:
        alignment2 = None
        num = -1
        for position in mask_seq:
            num += 1
            if position == '-':
                alignment2 = alignment[:, :num] + alignment[:, num + 1:]
                break

        return trim_alignment(alignment2)


def trim_nex(infilepath, outfilename):
    """Takes a filepath for a nexus alignment, and writes a new file containing
    a trimmed alignment.
    """
    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(infilepath)

    # Define the name of the output file.
    #outfilename = infilepath.replace('.nex', '.trim.nex')

    with open(infilepath) as infh, open(outfilename, 'w') as o:
        # Check that the input file has the filename extension ".nex".
        assert infilepath.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')
        mask_seq = alignment[-1]
        
        # Check that the MASK is present.
        assert mask_seq.id == 'MASK', "Error: Could not identify MASK sequence."

        trimmed_alignment = trim_alignment(alignment)
        AlignIO.write(trimmed_alignment, o, 'nexus')


