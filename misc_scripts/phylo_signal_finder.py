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

"""Ranks sequences in a fasta file according to how likely they are to contain
stronger phylogenetic signal for a specific clade represented by an HMM. This
is for comparing two HMMs/clades, and assumes that if the sequences can be
easily classified by their relative scores against the two HMMs, then this will
correspond to similarly favourable phylogenetic resolution.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
import re
import pylab
import subprocess


command_line_list = sys.argv
inhmmdb = str(command_line_list[1])
infasta = str(command_line_list[2])

gaps = re.compile(r' +')

def get_score_diff_from_dict(hmm1, hmm2, acc_sub_dict):
    """Determine the difference in scores to be used for ranking sequences.
    """
    score_diff = None

    if len(acc_sub_dict.keys()) == 2:
        score_diff = float(acc_sub_dict[hmm1]) - float(acc_sub_dict[hmm2])

    elif len(acc_sub_dict.keys()) == 1:
        if list(acc_sub_dict.keys())[0] == hmm1:
            score_diff = float(acc_sub_dict[hmm1])
        elif list(acc_sub_dict.keys())[0] == hmm2:
            score_diff = 0 - float(acc_sub_dict[hmm2])

    assert score_diff != None, "Error"

    return score_diff
    

# Check that there are not identical accessions in the input fasta file.


# Reduce redundancy in input fasta file using cd-hit.
#nonredunfasta = infasta + '_cd-hit_90.fa'
#subprocess.call(['cd-hit', '-i', infasta, '-o', nonredunfasta, '-c', '0.90'])
nonredunfasta = infasta # TEMPORARY.

# Run hmmscan.
hmmscanout = nonredunfasta + '_hmmscan_out.txt'
subprocess.call(['hmmscan', '--tblout', hmmscanout, inhmmdb, nonredunfasta])

# Get HMM names.
hmmlist = []
with open(hmmscanout) as infh1:
    for i in infh1:
        if not i.startswith('#'):
            hmmname = gaps.split(i)[0] 
            if hmmname not in hmmlist:
                hmmlist.append(hmmname)
            if len(hmmlist) == 2:
                break
hmm1 = hmmlist[0]
hmm2 = hmmlist[1]


# Construct dict of scores.
acc_dict = {}
with open(hmmscanout) as infh1:
    for i in infh1:
        if not i.startswith('#'):
            acc = gaps.split(i)[2]
            hmmname = gaps.split(i)[0] 
            score = gaps.split(i)[5]

            if acc not in acc_dict.keys():
                acc_dict[acc] = {}

            acc_dict[acc][hmmname] = score


# Get score differences for accessions.
acc_score_dict = {}
for acc in acc_dict.keys():
    acc_sub_dict = acc_dict[acc]
    score_diff = get_score_diff_from_dict(hmm1, hmm2, acc_sub_dict)
    acc_score_dict[acc] = score_diff


# Rank list of accessions based on score differences. 
ranked_accs = sorted(acc_score_dict.keys(), key=lambda x: acc_score_dict[x], reverse=True)

for acc in sorted(acc_score_dict.keys(), key=lambda x: acc_score_dict[x], reverse=True):
    print(str(acc_score_dict[acc]) + ' ' + acc)

# Write fasta files for ranked sequences for hmm1 and hmm2 (don't want reduced
# redundancy).
for hmmname in hmmlist:
    outfasta = infasta + '_' + hmmname + '_ranked.fa' 
    with open(outfasta, 'w') as o:
        if hmmname == hmm1:
            for acc in ranked_accs:
                if acc_score_dict[acc] >= 0:
                    with open(infasta) as infh1:
                        for record in SeqIO.parse(infh1, 'fasta'):
                            if record.id == acc:
                                SeqIO.write(record, o, 'fasta')
                    
        elif hmmname == hmm2:
            for acc in reversed(ranked_accs):
                if acc_score_dict[acc] <= 0:
                    with open(infasta) as infh1:
                        for record in SeqIO.parse(infh1, 'fasta'):
                            if record.id == acc:
                                SeqIO.write(record, o, 'fasta')



# Make histogram of scores.
score_diffs = [acc_score_dict[x] for x in acc_score_dict.keys()]
pylab.hist(score_diffs, bins=50) #specify the number of bins for the histogram
pylab.title("Score differences\n" + hmm1 + '\nvs\n' + hmm2)
pylab.xlabel("Difference in HMM sequence scores")
pylab.ylabel("Count")
pylab.show() #can do this instead of the savefig method if just want to view
#pylab.savefig(infilepath.rstrip('.fa') + '_seq_len_histo.output.pdf') #works for pdf or png
            



