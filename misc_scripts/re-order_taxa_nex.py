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
"""Takes two nex alignments and makes the taxa in the second one in the same
order as in the first one.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import subprocess
from Bio import SeqIO
from afa_to_nex import nex_to_afa, afa_to_nex, delete_extra_mesquite_lines

command_line_list = sys.argv
infp1 = str(command_line_list[1])
infp2 = str(command_line_list[2])

# Convert input alignments to temporary afa files.
infp1afa = infp1.replace('.nex', '.afa_temp')
nex_to_afa(infp1, infp1afa)
infp2afa = infp2.replace('.nex', '.afa_temp')
nex_to_afa(infp2, infp2afa)

# Get the order of IDs in the first afa file.
ids_in_order = []
for record in SeqIO.parse(infp1afa, 'fasta'):
    ids_in_order.append(record.description)

# Get new nex name.
new_afa = infp2.replace('.nex', '_re-ordered.afa')
new_nex = new_afa.replace('.afa', '.nex')

# Write new afa file in new order.
with open(new_afa, 'w') as new_afah, open(infp2afa) as infp2afah:
    # Make a list of sequence objects in the second alignment.
    seq_obj_list = []
    for record in SeqIO.parse(infp2afah, 'fasta'):
        seq_obj_list.append(record)

    # Write shared sequences to new fasta file.
    written_ids = []
    for i in ids_in_order:
        for r in seq_obj_list:
            if i == r.description:
                SeqIO.write(r, new_afah, 'fasta')
                written_ids.append(r.description)

    # Make sure any remaining sequences are still written.
    for i in seq_obj_list:
        if i.description not in written_ids:
            SeqIO.write(i, new_afah, 'fasta')

# Convert afa to nex.
afa_to_nex(new_afa, new_nex)

# Delete temporary files.
os.remove(infp1afa)
os.remove(infp2afa)
os.remove(new_afa)

