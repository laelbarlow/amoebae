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
"""This script is for writing fasta file that does not have the header
reformatting that is in the input fasta file written with the amoebae add_to_dbs
command.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#import glob
#import subprocess
import argparse
from Bio import SeqIO
#from afa_to_nex import nex_to_afa, afa_to_nex
#import shutil



def undo_formatting(infilepath):
    """Takes a fasta file and writes a new fasta file with unformatted
    (original) headers.
    """
    # Define output file path.
    outfilepath = infilepath.rsplit('.', 1)[0] + '_unformatted.'\
            + infilepath.rsplit('.', 1)[1]

    # Loop through the sequences in the input fasta file
    with open(infilepath) as i, open(outfilepath, 'w') as seqf:
        fasta_sequences = SeqIO.parse(i,'fasta')
        num = 0
        for seq in fasta_sequences:
            num += 1
            # Identify original header.
            header = seq.description
            original_header = header.split('\"', 1)[1][:-1]

            # Change header back to original.
            seq.description = original_header

            # Write output. 
            SeqIO.write([seq],seqf,'fasta')


if __name__ == '__main__':
    # Set up argument parsing.
    parser = argparse.ArgumentParser(
        description = """Writes a fasta file with unformatted headers.""",
        epilog = r"""epilog""")
    
    # Positional arguments.
    parser.add_argument('infp1', help='fasta file to be processed.')
    
    args = parser.parse_args()

    #do_split_fasta(args.infp1)
    #do_split_fasta2(args.infp1)
    undo_formatting(args.infp1)
