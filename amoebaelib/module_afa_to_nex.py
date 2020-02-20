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
import glob
import subprocess
import settings


def determine_alphabet(filepath):
    """Automatically determine what alphabet to use for a given file."""
    dtype = get_dbtype(filepath)
    if dtype == 'prot':
        return Gapped(IUPAC.protein)
    if dtype == 'nucl':
        return Gapped(IUPAC.ambiguous_dna)


def afa_to_nex(infilepath, outfilepath=None):
    """Takes an aligned fasta file (protein) and writes a nexus format file.
    This function replaces the record.id with the record.description in the
    nexus files, otherwise only the accession (before the first space character) is
    written to the nexus file records.
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

    # Move the original file to a new path temporarily.
    temp_copy_path = infilepath + '_temp_copy'
    os.rename(infilepath, temp_copy_path)

    # Write a new file.
    inhandle = open(infilepath, 'w')
    x = False
    for line in lines:
        if line == 'BEGIN ASSUMPTIONS;\n':
            x = True
        if not x:
            inhandle.write(line)
    inhandle.close()

    assert os.path.isfile(infilepath)
    # Remove the temporary copy of the input file.
    os.remove(temp_copy_path)


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


def align_fa(infilepath, outfilepath, aamatrix_path, program='muscle'):
    """Calls MUSCLE or MAFFT."""
    # Run MUSCLE by default.
    if program is None:
        program='muscle'

    # Check that an alignment program was specified.
    #assert program == 'muscle' or program == 'mafft', 'Error: Must specify\
    #alignment program as either muscle or mafft.'

    if program == 'muscle':
        # Call MUSCLE with default options.
        subprocess.call(["muscle", "-in", infilepath, "-out", outfilepath])

    # Run MAFFT if specified using somewhat default parameters.
    elif program == 'mafft':
        # Get the matrix file path with user's home directory. 
        aamatrix_path = os.path.expanduser(aamatrix_path)
        with open(outfilepath, 'w') as o:
            # The --localpair option may not be optimal. There are two other
            # options suggested for accurate alignments in the mafft
            # documentation. The --genafpair option is what Surkont et al.
            # used.
            stdoutdata = subprocess.call(['mafft', '--maxiterate', '1000',
                '--ep', '0', '--genafpair', '--aamatrix', aamatrix_path, infilepath], stdout=o)

    # Run MAFFT with variable parameters (modify this code).
    elif program == 'altmafft':
        # Get the matrix file path with user's home directory. 
        aamatrix_path = os.path.expanduser(aamatrix_path)
        with open(outfilepath, 'w') as o:
            # The --localpair option may not be optimal. There are two other
            # options suggested for accurate alignments in the mafft
            # documentation.
            stdoutdata = subprocess.call(['mafft', '--maxiterate', '1000',
                '--localpair', '--aamatrix', aamatrix_path, infilepath], stdout=o)


def align_all_fa(indirpath=None, outdirpath=None):
    """Calls the align_fa_muscle function on all .fa files in a given
    directory, which is the cwd by default.  And, outputs to a given output
    directory, which is also cwd by default.
    """
    cwd = os.getcwd()
    if indirpath == None:
        indirpath = cwd
    if outdirpath == None:
        outdirpath = cwd
    for f in [x for x in glob.glob(os.path.join(indirpath, '*.fa'))]:
        outfilename = os.path.basename(f).replace('.fa', '.afaa')
        align_fa(f, os.path.join(outdirpath, outfilename))


def align_one_fa(infilepath, outdirpath=None, program=None, aamatrix_path=None,
        conv_to_nex=None):
    """Calls the align_fa_muscle function on a .fa file.  And, outputs to a
    given output directory, which is the cwd by default.
    """
    #if not infilepath.rsplit('.', 1)[1] == 'fa':
    #    print('\n*** Warning: The file specified does not have the extension fa.')
    if outdirpath == None:
        outdirpath = os.path.dirname(infilepath)
    outfilename = os.path.basename(infilepath).rsplit('.', 1)[0] + '.afaa'
    outfilepath = os.path.join(outdirpath, outfilename)
    # Align with muscle with default settings and aamatrix.
    align_fa(infilepath, outfilepath, aamatrix_path,\
            program)

    # Optionally convert output file to nex and delete afa.
    if conv_to_nex:
        o = os.path.join(outdirpath, outfilename)
        outfilepath2 = outfilepath.rsplit('.', 1)[0] + '.nex'
        afa_to_nex(outfilepath, outfilepath2)
        os.remove(outfilepath)
        outfilepath = outfilepath2
    
    # Return path to output file.
    return outfilepath


