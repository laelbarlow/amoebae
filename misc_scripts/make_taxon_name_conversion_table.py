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
"""Takes an alignment before simplifying for input to mrbayes format (or other
formats), and after. Then makes a table for converting between the simplified
and unsimplified names.

This script assumes that the sequences are listed in the same order in both
input files.

"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
from afa_to_nex import nex_to_afa, afa_to_nex

command_line_list = sys.argv

infp1 = str(command_line_list[1]) # pre-simplification nex
infp2 = str(command_line_list[2]) # post-simplification nex
#infp3 = str(command_line_list[3]) 

# Convert nex files to afa so they can be easily parsed.
infp1a = infp1.replace('.nex', '.afa')
nex_to_afa(infp1, infp1a)

infp2a = infp2.replace('.nex', '.afa')
nex_to_afa(infp2, infp2a)

# Make table file for input to Name_replace.pl.
outfile_name = infp2.replace('.nex', '.table')

with open(infp1a) as infh1, open(infp2a) as infh2, open(outfile_name, 'w') as o:
    # Generate a list of headers from infp1.
    infp1_headers = []
    for r1 in SeqIO.parse(infh1, 'fasta'):
        infp1_headers.append(r1.description)
    
    # Generate a list of headers from infp2.
    infp2_headers = []
    for r2 in SeqIO.parse(infh2, 'fasta'):
        infp2_headers.append(r2.description)

    # Write a table file with alternating new and old headers, starting with
    # new.
    n = 0
    for i in infp1_headers:
        o.write(infp2_headers[n] + '\n')
        o.write(i + '\n')
        n += 1

# Remove afa files.
os.remove(infp1a)
os.remove(infp2a)


