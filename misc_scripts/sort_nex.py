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
"""Takes an nex alignment and sorts the taxa according to their names. 
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import subprocess
from Bio import SeqIO
from module_afa_to_nex import nex_to_afa, afa_to_nex, delete_extra_mesquite_lines
import re

command_line_list = sys.argv
infp1 = str(command_line_list[1])

def sort_headers(headers):
    """Sort a list of headers.
    """
    # rev = lambda x: x[::-1] # Returns reverse of header string.
    prot = lambda x: x.split(' ')[-1].split('_')[-1] # Returns the protein name.

    # for i in headers:
    #     print(prot(i))
    sorted_headers = headers
    sorted_headers.sort(key=prot)

    assert len(headers) == len(sorted_headers), 'Error: Different number of\
 headers in original and sorted list.'

    return sorted_headers
   

# Convert input alignments to temporary afa files.
infp1afa = infp1.replace('.nex', '.afa_temp')
nex_to_afa(infp1, infp1afa)

# Get the order of IDs in the afa file.
ids_in_order = []
for record in SeqIO.parse(infp1afa, 'fasta'):
    ids_in_order.append(record.description)

# Get new nex name.
new_afa = infp1afa.replace('.afa_temp', '.sort.afa')
new_nex = new_afa.replace('.afa', '.nex')

# Write new afa file in new order.
with open(infp1afa) as infp1afah, open(new_afa, 'w') as new_afah:
    # Make a list of sequence objects in the first alignment.
    seq_obj_list = []
    for record in SeqIO.parse(infp1afah, 'fasta'):
        seq_obj_list.append(record)

    # Get a list of headers.
    headers = []
    for record in seq_obj_list:
        headers.append(record.description)

    # Sort header list.    
    sorted_headers = sort_headers(headers)

    # Write shared sequences to new fasta file.
    for i in sorted_headers:
        for r in seq_obj_list:
            if i == r.description:
                SeqIO.write(r, new_afah, 'fasta')


# Convert afa to nex.
afa_to_nex(new_afa, new_nex)

# Delete temporary files.
os.remove(infp1afa)
os.remove(new_afa)

