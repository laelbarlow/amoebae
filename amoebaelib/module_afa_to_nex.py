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
"""This module is for afa_to_nex.py script and others, which use the same
functions.
"""
import sys
import os
import shutil
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from module_amoebae_get_datatype import get_dbtype
from Bio.Seq import Seq


def determine_alphabet(filepath):
    """Automatically determine what alphabet to use for a given file."""
    dtype = get_dbtype(filepath)
    if dtype == 'prot':
        return Gapped(IUPAC.protein)
    if dtype == 'nucl':
        return Gapped(IUPAC.ambiguous_dna)


def afa_to_nex(infilepath, outfilepath=None):
    """Takes an aligned fasta file (protein) and writes a nexus format file.
    """
    # Delete extra mesquite lines if present (otherwise they cause errors).
    # These may have been added if the nex file was opened and saved in
    # mesquite. IS THIS NECESSARY???
    delete_extra_mesquite_lines(infilepath)

    # Convert afa to nex.
    if outfilepath is None:
        outfilepath = os.path.dirname(infilepath)
    inhandle = open(infilepath)
    alph = determine_alphabet(infilepath) # Get dbtype
    alignment = AlignIO.read(inhandle, "fasta",
        alphabet=alph)
    for record in alignment:
        x = record.description
        record.id = x
    outhandle = open(outfilepath, "w")
    AlignIO.write(alignment, outhandle, "nexus")
    inhandle.close()
    outhandle.close()

def nex_to_afa(infilepath, outfilepath):
    """Takes a nexus file (protein) and writes an aligned fasta file.
    """
    # Delete extra mesquite lines if present (otherwise they cause errors).
    # These may have been added if the nex file was opened and saved in
    # mesquite.
    delete_extra_mesquite_lines(infilepath)

    # Convert nex to afa.
    inhandle = open(infilepath)
    alph = determine_alphabet(infilepath) # Get dbtype
    alignment = AlignIO.read(inhandle, "nexus",
        alphabet=alph)
    outhandle = open(outfilepath, "w")
    AlignIO.write(alignment, outhandle, "fasta")
    inhandle.close()
    outhandle.close()

def delete_extra_mesquite_lines(infilepath):
    """Remove extra blocks added by Mesquite, if present.  This prevents
    Biopython from raising an error about an unmatched 'end' in one of the
    nexus blocks.  The extra blocks are not necessary anyway.
    """
    inhandle = open(infilepath, 'r')
    lines = inhandle.readlines()
    inhandle.close()
    inhandle = open(infilepath, 'w')
    x = False
    for line in lines:
        if line == 'BEGIN ASSUMPTIONS;\n':
            x = True
        if not x:
            inhandle.write(line)
    inhandle.close()

def nex_to_phylip(infilepath, outfilepath):
    """Takes nexus format file and makes a phylip format file.
    Warning: Taxon names must be no more than 10 alphanumeric characters.
    """
    # Delete extra mesquite lines if present (otherwise they cause errors).
    # These may have been added if the nex file was opened and saved in
    # mesquite.
    delete_extra_mesquite_lines(infilepath)

    # Convert nex to phylip.
    inhandle = open(infilepath)
    alph = determine_alphabet(infilepath) # Get dbtype
    alignment = AlignIO.read(inhandle, "nexus",
        alphabet=alph)
    outhandle = open(outfilepath, "w")
    AlignIO.write(alignment, outhandle, "phylip-sequential")
    inhandle.close()
    outhandle.close()

def nex_to_mbnex(infilepath, outfilepath, mrbayescodeblocks):
    """Takes a nexus format alignment file and converts to a format for use in
    a phylogenetic analysis using MrBayes.
    """
    # Delete extra mesquite lines if present (otherwise they cause errors).
    # These may have been added if the nex file was opened and saved in
    # mesquite.
    delete_extra_mesquite_lines(infilepath)

    # Parse alignment sequences and replace '?' with 'X'.
    infilepath2 = infilepath + '_TEMP.nex'
    with open(infilepath) as infh, open(infilepath2, 'w') as o:
        alignment = AlignIO.read(infh, "nexus")
        for i in alignment:
            #i.seq.seq = i.seq.seq.replace('?', 'X')
            i.seq = Seq(str(i.seq).replace('?', 'X'))

        AlignIO.write(alignment, o, "nexus")

    # Modify input file to make output file.
    with open(infilepath2) as infh, open(outfilepath, 'w') as o:
        for i in infh:
            # Remove "gap" and change "missing" to "-".
            if i.startswith('\tformat'):
                newi = i.replace('missing=?', 'missing=-').replace('gap=-', '')
                o.write(newi)

            # Append MrBayes code blocks.
            elif i.startswith('end;'):
                o.write(i + '\n')
                o.write(mrbayescodeblocks)

            else:
                o.write(i)

    # Remove temporary file.
    os.remove(infilepath2)



