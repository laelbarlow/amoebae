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
"""Module for amoebae script with functions for running hmmscan (from the
HMMer3 software package).
"""
# Import built-in modules.
#import argparse
import sys
import os
import subprocess
#import re
#import shutil
import glob
import time
import pandas as pd

# Import modules from installed libraries/packages.
#from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped

#from get_datatype import get_dbtype
#from paralogue_counter import get_seq_obj_from_srch_res_csv_info
from nex_to_hmm import nex_to_hmm

#import settings

# Define functions to be used in amoebae.

def all_nex_to_hmm(indirpath):
    """Converts all nex files in a given directory to .hmm (writes .hmm files).
    """
    # Loop over files.
    for nex in glob.glob(os.path.join(indirpath, '*.nex')):
        # Convert nex to hmm.
        nex_to_hmm(nex)


def get_hmm_paths(indirpath):
    """Returns a list of all .hmm files in a given directory.
    """
    return glob.glob(os.path.join(indirpath, '*.hmm'))

def cat_hmm(hmm_paths, outfilepath):
    """Takes a list of .hmm file paths, and outputs a file that is a
    concatenation of these, and returns its path.
    """
    outfilehandle = open(outfilepath, 'w')
    for f in hmm_paths:
        with open(f) as h:
            for l in h:
                outfilehandle.write(l)
    outfilehandle.close()

def call_hmmpress(infilepath):
    """Calls hmmpress on a file to make a concatenation of .hmm files into an
    HMM database for searching with hmmscan.
    """
    os.chdir(os.path.dirname(infilepath))
    subprocess.call(['hmmpress', infilepath])

def make_hmm_db(indirpath, hmmdbname=None, outdirpath=None):
    """Takes .hmm files in a given input directory, and outputs a searchable
    HMM database (files) in a given output directory.
    """
    # If no hmm db name provided, then set to default name.
    if hmmdbname is None:
        hmmdbname = 'hmmdb'

    # Append a timestamp to the database filename.
    cur_time = time.strftime("%Y_%m_%d_%H_%M_%S")
    hmmdbname = hmmdbname + '_' + cur_time + '.hmmdb'

    # If no out dir path specified, then just use in dir path.
    if outdirpath is None:
        outdirpath = indirpath

    # Convert .nex alignment files to HMMs (.hmm).
    all_nex_to_hmm(indirpath)

    # Get the .hmm file paths.
    hmms = get_hmm_paths(indirpath)

    # Concatenate .hmm files.
    hmmdbpath = os.path.join(outdirpath, hmmdbname)
    chmm = cat_hmm(hmms, hmmdbpath) 

    # Compress the concatenated file.
    call_hmmpress(hmmdbpath)


def determine_if_is_hmmdb(infp):
    """Return True if the given file is an HMM database (generated using
    hmmpress from the HMMer3 software package), and return False otherwise.
    """
    #if open(infp, 'r').read().startswith('HMMER3/f'):
    if open(infp, 'r').readline().startswith('HMMER3/f'):
        return True
    else:
        return False

