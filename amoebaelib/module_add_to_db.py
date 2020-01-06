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
"""Module for script: amoebae.
"""
# Import built-in modules.
#import argparse
import sys
import os
import subprocess
import re
#import settings
import shutil
#import glob
import time
#import pandas as pd

# Import modules from installed libraries/packages.
from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped

from module_amoebae_get_datatype import get_dbtype


# Define functions to be used in amoebae.

def make_easel_index(outfp):
    """Use the esl-sfetch miniapp from the Easel package associated with the
    HMMER3 software package to index a fasta file.
    """
    subprocess.call(['esl-sfetch', '--index', outfp])


def convert_headers(infp, outfp, split_char=' ', split_pos='0',
        skip_header_reformat=False, auto_extract_seqs=False):
    """For converting headers in fasta files, and writing a new fasta file with
    appropriate filename extension.

    ***If told to split on "empty" (empty string), then won't split at all (for
    headers that are composed only of accession).
    """
    # Check that the outfilepath doesn't already exist.
    assert not os.path.isfile(outfp), """Error: Specified output path already
    exists: %s""" % outfp
     
    if not skip_header_reformat: 
        # If automatically identifying accessions, then figure out what position in
        # the split headers the accessions are in.
        header_split_position_for_auto_extract = None
        split_chars = re.compile(r'(\ |\||:|;|\t)')
        if auto_extract_seqs:
            for position in range(0, 9):
                # Get list of accessions at this position of the split headers.
                accs = []
                with open(infp) as infh:
                    for seq in SeqIO.parse(infp, 'fasta'):
                        header = seq.description
                        acc = split_chars.split(header)[position]
                        accs.append(acc)
                # See to what extent the accessions are unique.
                num_identical_accs = len(accs) - len(set(accs))
                percent_identical_accs = (num_identical_accs / len(accs)) * 100
                if percent_identical_accs < 5:
                    header_split_position_for_auto_extract = position
                    break
            # Check that a position with unique (enough) accessions was found.
            assert header_split_position_for_auto_extract is not None, """Could
            not find a position with unique accessions after splitting on
            spaces and pipe characters."""

        # Convert all the headers by manipulating sequence objects.
        accs = []
        with open(infp) as infh, open(outfp, 'w') as o:
            for seq in SeqIO.parse(infp, 'fasta'):
                header = seq.description # Should contain seq.id as well.
                acc = None
                new_header = None
                # Extract accession from header.
                if auto_extract_seqs:
                    # Extract acc from header as determined above.
                    acc =\
                    split_chars.split(header)[header_split_position_for_auto_extract]
                    print(acc)
                else:
                    if split_char == 'empty':
                        acc = header
                    else:
                        acc = header.split(split_char)[int(split_pos)]

                # Modify the accession if a previous sequence has the same
                # accession.
                alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                if acc in accs:
                    modified = False
                    old_acc = acc
                    for letter in alphabet:
                        modified_acc = acc + '_' + letter
                        if modified_acc not in accs:
                            acc = modified_acc
                            modified = True
                            break
                    if not modified:
                        for letter1 in alphabet:
                            for letter2 in alphabet:
                                modified_acc = acc + '_' + letter1 + letter2
                                if modified_acc not in accs:
                                    acc = modified_acc
                                    modified = True
                                    break
                            if modified:
                                break
                    assert modified, """Could not modify the accession number
                    for sequence with duplicate accession number %s.""" % acc

                    # Report change to accession number.
                    print("""Accession %s already in set for writing to output,
                    so modified to %s.""" % (old_acc, acc))

                # Add accession to list.
                accs.append(acc)
                # Define new header, and change string object id attribute.
                new_header = acc + ' ' + '\"' + header + '\"'
                seq.id = new_header
                # Remove additional sequence description attribute if present.
                seq.description = ''
                # Write sequence with modified header to output file.
                SeqIO.write(seq, o, 'fasta')

        ## Check that all the accessions are unique within the input file.
        #assert len(accs) == len(set(accs)), """Error: Some accessions in input
        #fasta file (%s) are identical.""" % os.path.basename(infp)

        # Check that all the accessions are unique within the output file.
        assert len(accs) == len(set(accs)), """Error: Some accessions for
        output fasta file are identical."""

    else:
        # Just copy the file as-is.
        shutil.copyfile(infp, outfp)


def make_blast_db(infp):
    """Makes blastable database for given fasta file (in same directory).
    """
    # Determine data type from filename extension.
    dbtype = get_dbtype_from_file_exten(infp)

    # Make the file blast-able with makeblastdb.
    subprocess.call(["makeblastdb", "-in", infp, "-dbtype",
        dbtype])


def get_corr_fasta_exten(infp):
    """Determines the correct filename extension that a given fasta file should
    have.
    """
    # Automatically determine what the dbtype is.
    dbtype = get_dbtype(infp)
    assert dbtype == 'prot' or dbtype == 'nucl', """Error: Could not determine
    data type in file " + f + " using get_datatype module."""

    # Determine appropriate filename extension to use.
    exten = None
    if dbtype == 'prot':
        exten = 'faa'
    elif dbtype == 'nucl':
        exten = 'fna'

    return exten


def get_dbtype_from_file_exten(infp):
    """Determine data type from fasta filename extension.
    """
    dbtype = None
    exten = infp.rsplit('.', 1)[1]
    if exten == 'faa':
        dbtype = 'prot'
    elif exten == 'fna':
        dbtype = 'nucl'
    elif exten == 'hmmdb':
        dbtype = 'prot'
    elif exten == 'sql':
        dbtype = 'annotations'
    assert dbtype is not None, """Could not determine datase type based on
    filename extension: %s""" % exten
    return dbtype


def update_csv(outfp, csv_file):
    """Appends a line to the given spreadsheet with info about the given fasta
    file added to a directory.
    """
    # Determine data type from filename extension.
    dbtype = get_dbtype_from_file_exten(outfp)

    # Get current date.
    cur_date = time.strftime("%Y/%m/%d")

    # Get abbreviation to add.
    fp_without_ext = os.path.basename(outfp.rsplit('.',1)[0])

    # Check whether file has a newline character at the end so a new line can
    # start correctly.
    newline = False
    if open(csv_file, 'r').read()[-1] == '\n':
        newline = True

    # Write to file.
    with open(csv_file, 'a') as csvh:
        x = '\n'
        if newline:
            x = ''
        csvh.write(x + os.path.basename(outfp) + ',,,,' + fp_without_ext + ',,,' +\
                        dbtype + ',,,' + cur_date +\
        '\n')

