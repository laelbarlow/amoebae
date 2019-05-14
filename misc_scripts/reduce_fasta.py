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
"""Removes all but the first instance of a sequence that came from each genome
(for working with results from get_top_n_... script).

Or, writes a reduced fasta file from the first input fasta file, which only
includes sequences with assessions that are present in the second input fasta
file.
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


# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Outputs positive hits and sequences.""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('infp1', help='fasta file to be reduced.')
parser.add_argument('infp2', help='''fasta file that has the more limited set
of sequences.''')


args = parser.parse_args()

def reduce_fasta(infilepath):
    # Loop through the sequences in the input fasta file
    seqfilename = infilepath.replace('.fa', '_reduced.fa')
    with open(infilepath) as i, open(seqfilename, 'w') as seqf:
        fasta_sequences = SeqIO.parse(i,'fasta')
        num = 0
        db_list = []
        for seq in fasta_sequences:
            num += 1
            split_header = seq.description.split(' ', 1)
            acc = split_header[0]
            dbname = split_header[1]
            if dbname not in db_list:
                SeqIO.write([seq],seqf,'fasta')
                db_list.append(dbname)
            else:
                pass


def reduce_fasta2(infilepath1, infilepath2):
    # Loop through the sequences in the input fasta file
    seqfilename = infilepath1.rsplit('.', 1)[0] + '_reduced.fa'
    with open(infilepath1) as i, open(infilepath2) as i2, open(seqfilename, 'w') as seqf:
        # Compile a reduced list of accessions (for sequences to keep in
        # output).
        reduced_acc_list = []
        for x in SeqIO.parse(i2, 'fasta'):
            reduced_acc_list.append(x.id)

        # Iterate through the first input fasta, and only write appropriate
        # sequences to output.
        db_list = []
        fasta_sequences = SeqIO.parse(i,'fasta')
        for seq in fasta_sequences:
            acc = seq.id
            if acc in reduced_acc_list:
                SeqIO.write([seq],seqf,'fasta')
            else:
                pass

    

if __name__ == '__main__':
    #reduce_fasta(args.infp1)
    reduce_fasta2(args.infp1, args.infp2)
