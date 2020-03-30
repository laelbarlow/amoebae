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
#import re
#import settings
#import shutil
#import glob
import time
import pandas as pd

# Import modules from installed libraries/packages.
from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped

import amoebae_m
from module_amoebae_get_datatype import get_dbtype


# Define functions to be used in amoebae.


def get_query_file_type(query_file):
    """Takes an input query file and determines whether it is a FASTA file or
    an HMM file. Returns either "FASTA" or "HMM".
    """
    file_type = None
    with open(query_file) as infh:
        first_line = infh.readline().strip()
        if first_line.startswith('>'):
            file_type = 'FASTA'
        elif first_line.startswith('HMM'):
            file_type = 'HMM'

    # Check that it worked.
    assert file_type is not None, """Error: Input file type could not be
    determined: %s""" % query_file

    # Return file type.
    return file_type


def get_hmm_datatype(query_file):
    """Takes an HMM file (HMMer3 software package) and determines what data
    type it has (i.e., generated from an amino acid or nucleic acid alignment).

    Returns either "prot" or "nucl".
    """
    datatype = None
    with open(query_file) as infh:
        for i in infh:
            if i.startswith('ALPH'):
                dname = i.strip().split('  ')[1]
                if dname == 'amino':
                    datatype = 'prot'
                elif dname == 'DNA':
                    datatype = 'nucl'
                break

    # Check that it worked.
    assert datatype is not None, """Error: Data type could not be
    determined for input file: %s""" % query_file

    # Return the data type.
    return datatype


def get_mod_query_path(query_file, filetype, datatype, query_dir):
    """Define a new query path.
    """
    exten = None
    if filetype == 'fa': 
        if datatype == 'prot':
            exten = 'faa'
        elif datatype == 'nucl':
            exten = 'fna'
    elif filetype == 'afa':
        if datatype == 'prot':
            exten = 'afaa'
        elif datatype == 'nucl':
            extenn = 'afna'

    # Check that it worked.
    assert exten is not None, """Error: New extension could not be
    determined for input file: %s""" % query_file

    # Define basename with new extension.
    query_file_new_exten = os.path.basename(query_file).rsplit('.', 1)[0] +\
    '.' + exten

    # Get new path.
    new_query_file_path = os.path.join(query_dir, query_file_new_exten)

    # Return the new path.
    return new_query_file_path


def is_single_fasta(query_file):
    """Return True if the given fasta file contains only one sequence.
    """
    single = False
    num_seqs = 0
    with open(query_file) as infh:
        for seq in SeqIO.parse(query_file, 'fasta'):
            num_seqs += 1
    if num_seqs == 1:
        single = True

    return single


def fasta_seqs_all_same_len(query_file):
    """Returns True if the sequences in a given fasta file are all the same
    length.
    """
    seq_lengths = []
    with open(query_file) as infh:
        for seq in SeqIO.parse(infh, 'fasta'):
            seq_lengths.append(len(seq))
    # If the sequences are the same length, then a nonredundant list of
    # sequence lengths will only contain one element.
    if len(set(seq_lengths)) == 1:
        return True
    else:
        return False

def update_query_csv(csv_file, mod_query_path, datatype):
    """Appends a line to the given spreadsheet with info about the given fasta
    file added to a directory.
    """
    # Define column headers.
    headers = ['Filename',
               'Query title',
               'Query source description',
               'Query taxon (species if applicable)',
               'Query database filename (if applicable)',
               'File type',
               'Data type',
               'Date added',
               'Citation'
               ]

    # Make new csv file if necessary.
    if not os.path.isfile(csv_file):
        #df = pd.DataFrame(columns=headers)
        #df.to_csv(csv_file) #, index_label='Filename')
        with open(csv_file, 'w') as o:
            o.write(','.join(headers)) # + '\n')

    # Load dataframe from csv file.
    df = pd.read_csv(csv_file) #, index_col='Filename')

    # Get current date.
    cur_date = time.strftime("%Y/%m/%d")

    # Get query filename and extension.
    full = os.path.basename(mod_query_path).rsplit('.', 1)
    filename = full[0]
    exten = full[1]

    # Get query basename.
    query_basename = os.path.basename(mod_query_path)

    # Extract info from query filename.
    query_title = '?'
    taxon = '?'
    species = '?'
    if len(query_basename.split('_')) > 2:
        # Get query title name.
        query_title = query_basename.split('_')[0]

        # Get query taxon name.
        # Assumes that there is an accession or something else after the "taxon" in
        # the query file name.
        taxon = amoebae_m.get_query_taxon_from_filename(query_basename)

        # Get species based on taxon.
        species = amoebae_m.get_species_from_db_csv(taxon)
    else:
        ## Print warning.
        #print("""Warning: Could not identify query title or database/taxon name
        #in input filename.""")
        ## Just use the whole filename minus the filename extension.
        #query_title = query_basename.rsplit('.')[0]

        query_title = query_basename.split('_')[0]

    # Initiate dataframe for line to append.
    new_row = pd.DataFrame(columns=headers)

    # Get database filename.
    db_filename = amoebae_m.get_db_filename_for_query_from_db_csv(taxon)

    # Add info to new row.
    new_row.loc[0] = ['???'] * len(headers)

    new_row.loc[0]['Filename'] = query_basename
    new_row.loc[0]['Query title'] = query_title
    new_row.loc[0]['Query source description'] = taxon
    new_row.loc[0]['Query taxon (species if applicable)'] = species
    new_row.loc[0]['Data type'] = datatype
    new_row.loc[0]['File type'] = exten
    new_row.loc[0]['Date added'] = cur_date
    new_row.loc[0]['Citation'] = '?'
    new_row.loc[0]['Query database filename (if applicable)'] =\
    db_filename
    # Check that it worked.
    assert not '???' in new_row.loc[0], """Could not add all the necessary info
    to the query info spreadsheet for query file:\n\t%s""" % mod_query_path

    # Append new row to dataframe.
    df = df.append(new_row, ignore_index=True)

    # Re-order columns in output dataframe.
    df = df[new_row.columns]

    # Write updated dataframe to csv file.
    #df[1:].to_csv(csv_file) 
    df.to_csv(csv_file, index=False) 

    # Report activity:
    #print('Information added to spreadsheet %s:' % os.path.basename(csv_file))
    #print('\tFilename: ' +  query_basename)
    #print('\tQuery title: ' + query_title)
    #print('\tQuery source description: ' + taxon)
    #print('\tQuery taxon (species if applicable): ' + species)
    #print('\tData type: ' + datatype)
    #print('\tFile type: ' + exten)
    #print('\tDate added: ' + cur_date)
    #print('\tCitation: ' + '?')
    #print('\tQuery database filename (if applicable): ' + db_filename)


