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
"""For appending afa files in a given directory, so they can be viewed as a
single alignment file (if there is a limitation on supplementary files uploaded
with a manuscript, for example).
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from afa_to_nex import afa_to_nex

command_line_list = sys.argv
indp = str(command_line_list[1])

# Find the longest alignment.
files = glob.glob(indp + '/*.afa')
lengths = []
for f in files:
    #print(os.path.basename(f))
    l = None
    num = 0
    for record in SeqIO.parse(f, "fasta"):
        num += 1
        if num == 1:
            l = len(record.seq)
            lengths.append(l)
            #print(l)
        else:
            break

max_len = max(lengths)
#print(max_len)


# Write a new alignment with sequences the same length.
files = glob.glob(indp + '/*.afa')
with open(indp + '/output.afa', 'w') as o:
    for f in files:
        alignment_name = os.path.basename(f).rstrip('.afa')
        #print(alignment_name)
        o.write('>***Begin alignment ' + alignment_name + '\n' + '-' * max_len + '\n')
        for record in SeqIO.parse(f, "fasta"):
            to_append = max_len - len(record.seq)
            record.seq = record.seq + '-' * to_append 
            SeqIO.write(record, o, "fasta")
        o.write('>***End alignment' + alignment_name + '\n' + '-' * max_len + '\n')

