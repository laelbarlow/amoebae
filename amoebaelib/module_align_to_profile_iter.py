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
"""Script for iteratively aligning sequences from a fasta file one-by-one with
an existing alignment using MUSCLE's -profile option.

********** Ideally, the fastas should be added to the alignment in decreasing order of
how well they align to the HMM (based on score).

Usage:
    align_to_profile_iter.py <nex starting alignment> <fasta file with seqs to
    align to nex alignment>
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import subprocess
import argparse
from Bio import SeqIO
from module_afa_to_nex import nex_to_afa, afa_to_nex
import shutil



def split_fasta(infilepath, temp_subdirname):
    """Takes a fasta files and writes new fasta files in a subdirectory for
    each sequence in the fasta file.
    """
    # Loop through the sequences in the input fasta file
    with open(infilepath) as i:
        fasta_sequences = SeqIO.parse(i,'fasta')

        num = 0
        for seq in fasta_sequences:
            num += 1
            seqfilename = os.path.join(temp_subdirname, str(num) + '.fa')
            with open(seqfilename, 'w') as seqf:
                SeqIO.write([seq],seqf,'fasta')

def align_iteratively(inalignfilepath, temp_subdirname, outfilepath):
    """Takes an afa file and aligns it one-by-one to sequences in files in the
    temporary subdirectory, and outputs the final alignment in fasta format.
    """
    # Make a list of the input fa files, and sort it numerically based on
    # filename.
    file_list = glob.glob(temp_subdirname + '/*.fa')
    l = lambda x: int(os.path.basename(x).split('.')[0])
    file_list_sorted = sorted(file_list, key=l)

    # Loop through the sorted list of files and align each sequence to the
    # input alignment.
    for f in file_list_sorted:
        print('\n\n' + f + '\n\n')
        # Call MUSCLE with default options.
        subprocess.call(["muscle", "-profile", "-in1", inalignfilepath, "-in2", f,
            "-out", outfilepath])
    
def do_align_iteratively(innexpath, infapath, outnexpath=None):
    """Do the whole process
    """
    # Make a temporary subdirectory for storing individual sequence files.
    temp_subdirname = infapath.rsplit('.', 1)[0] + '_temp'
    os.mkdir(temp_subdirname)

    # Split the input fasta file.
    split_fasta(infapath, temp_subdirname)

    # Name output file.
    infapath_basename = os.path.basename(infapath)
    outafapath = None
    if outnexpath is None:
        outafapath = innexpath.rsplit('.', 1)[0] + '+' + input('Name: ') + '.afaa'
    else:
        outafapath = outnexpath.rsplit('.', 1)[0] + '.afaa'

    # Convert nex to afa.
    nex_to_afa(innexpath, outafapath) 

    # Align seqs in tempdir to afa iteratively.
    align_iteratively(outafapath, temp_subdirname, outafapath)

    # Delete tempdir.
    shutil.rmtree(temp_subdirname)

    # Define outnex path.
    if outnexpath is None:
        outnexpath = outafapath.rsplit('.', 1)[0] + '.nex'

    # Convert outafa to nex format.
    afa_to_nex(outafapath, outnexpath)

    # Delete afa file.
    os.remove(outafapath)

    #print(outnexpath)
