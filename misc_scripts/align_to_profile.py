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
from afa_to_nex import nex_to_afa, afa_to_nex
import shutil

from align_to_profile_iter import do_align_iteratively




if __name__ == '__main__':
    # Set up argument parsing.
    parser = argparse.ArgumentParser(
        description = """Outputs positive hits and sequences.""",
        epilog = r"""epilog""")
    
    # Positional arguments.
    parser.add_argument('infp1', help='nex alignment used to make an HMM.')
    parser.add_argument('infp2', help='fasta file with sequences to be aligned\
            iteratively.')
    
    # Optional arguments.
    #parser.add_argument('-p', '--program', help='Program (either muscle or mafft)',
    #        default='muscle')
    #parser.add_argument('-m', '--aamatrix', help='Amino acid matrix (NCBI format)',
    #        default=None)
    
    args = parser.parse_args()

    do_align_iteratively(args.infp1, args.infp2)
