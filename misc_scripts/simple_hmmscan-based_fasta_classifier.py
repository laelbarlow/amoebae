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

"""Takes a fasta file and an hmmscan output file (tabular output), and outputs
fasta files with sequences that score best against each of the HMMs.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
import re


command_line_list = sys.argv
infasta = str(command_line_list[1])
inhmmscan = str(command_line_list[2])

gaps = re.compile(r' +')

# Get classifications for each sequence in the fasta file from the hmmscan
# output.
classification_dict = {}
with open(infasta) as infh1:
    for record in SeqIO.parse(infh1, 'fasta'):
        acc = record.id
        with open(inhmmscan) as infh2:
            for i in infh2:
                if not i.startswith('#'):
                    acc2 = gaps.split(i)[2]
                    hmm = gaps.split(i)[0]
                    if acc2 == acc:
                        if not hmm in classification_dict.keys():
                            classification_dict[hmm] = []
                        # Add accession to list for hmm in dict.
                        classification_dict[hmm] = classification_dict[hmm] + [acc]
                        # Take the first instance.
                        break


# Check that the lists are mutually exclusive.
#print('Checking that the classification groups are mutually exclusive.')
#for hmm1 in classification_dict.keys():
#    set1 = set(classification_dict[hmm])

    #for hmm2 in classification_dict.keys():
    #    if hmm2 is not hmm1:
    #        print('\tComparing ' + hmm1 + ' to ' + hmm2)
    #        set2 = set(classification_dict[hmm2])
    #        print('\tAccessions found in both lists:')
    #        print('\t\t' + str(set1.intersection(set2)))
    #        #print(hmm1)
    #        #print(hmm2)
    #        #print(set1.intersection(set2))
    #        #assert len(set1.intersection(set2)) == 0, "Error"


# loop through keys in classification dict, and write corresponding files with
# fasta sequences.
for hmm in classification_dict.keys():
    print('\n')
    print(hmm)
    outfp = infasta + '_' + hmm + '_matches.fa'
    with open(infasta) as infh1, open(outfp, 'w') as o:
        for record in SeqIO.parse(infh1, 'fasta'):
            if record.id in classification_dict[hmm]:
                print('\t' + record.id)
                SeqIO.write(record, o, 'fasta')


