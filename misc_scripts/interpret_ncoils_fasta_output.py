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
"""Interpret the output of NCOILS.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import settings
from Bio import SeqIO
from misc_functions import get_abbrev_fa_record_text_from_obj
import re

command_line_list = sys.argv
ncoils_out = str(command_line_list[1])

def interpret_ncoils_fasta(infp):
    outfp = infp.replace('.fa', '_interpreted.csv')
    with open(infp) as infh, open(outfp, 'w') as o:
        o.write('Record description,Sequence length,Number of coiled-coil\
residues,Percent coiled-coil residues,Coiled-coil segment lengths (aa)\n')
        for record in SeqIO.parse(infh, 'fasta'):
            # Calculate relevant values.
            seq = str(record.seq)
            seq_len = len(seq)
            num_x = seq.count('x')
            percent_x = (num_x * 100)/seq_len

            # compile a regular expression for identifying segments composed of
            # consecutive 'x' characters.
            xseg = re.compile(r'x+')

            # Get a list of lengths for strings of 'x's.
            xsegs = [len(x) for x in xseg.findall(seq)]
            xsegs2 = sorted(xsegs, reverse=True)
            xsegs3 = [str(x) for x in xsegs2]

            # Write values to file.
            o.write(record.description + ',' + str(seq_len) + ',' + str(num_x)\
                    + ',' + str(percent_x) + ',' + ','.join(xsegs3) + '\n')


if __name__ == '__main__':
    interpret_ncoils_fasta(ncoils_out)

