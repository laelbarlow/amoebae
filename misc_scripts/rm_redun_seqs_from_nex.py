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
"""For removing sequences with redundant accessions from nex alignments,
without realigning.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
from afa_to_nex import nex_to_afa, afa_to_nex
#from afa_to_fa import afa_to_fa

command_line_list = sys.argv

infp = str(command_line_list[1])

afafp = infp.replace('.nex', '.afa')

nex_to_afa(infp, afafp)

afafp2 = afafp.replace('.afa', '_nonredun.afa')

with open(afafp) as i, open(afafp2, 'w') as o:
    acc_list = []
    for seq in SeqIO.parse(i, 'fasta'):
        acc = seq.description.split(' ', 1)[0]
        if acc not in acc_list:
            SeqIO.write(seq, o, 'fasta') 
        else:
            print('Removing redundant sequence: ' + seq.description)
        acc_list.append(acc)

os.remove(afafp)

nexfp2 = afafp2.replace('.afa', '.nex')

afa_to_nex(afafp2, nexfp2) 

os.remove(afafp2)



