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
"""This script is for adding fasta sequences that are not already contained in
a fasta file.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
import argparse
from misc_functions import get_fa_record_text_from_obj

# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Adds only new sequences.""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('old_fa', help='old fasta file path')
parser.add_argument('new_fa', help='new fasta file path')

args = parser.parse_args()

def add_new_seqs(old_fa, new_fa, old_fa2=None):
    """Adds sequences from new_fa to old_fa, if they are not already contained
    in old_fa.
    """
    mode = 'a'
    if old_fa2 is None:
        old_fa2 = old_fa
    else:
        mode = 'w'

    seq_dict = {}
    for record in SeqIO.parse(old_fa, 'fasta'):
        #print(record.description)
        acc = record.description.split(' ', 1)[0]
        # this may have to change if there is no protein name after the db_name.
        print('record.description = ' + record.description)
        db_name = record.description.split(' ', 1)[1].split('_', 1)[0] 
        
        if db_name not in seq_dict.keys():
            seq_dict[db_name] = [acc]
        else:
            seq_dict[db_name].append(acc)
    
    old_fa_handle = open(old_fa2, mode)
    old_fa_handle.write('\n')
    added_seq_count = 0
    for record in SeqIO.parse(new_fa, 'fasta'):
        acc = record.description.split(' ', 1)[0]
        # this may have to change if there is no protein name after the db_name.
        db_name = record.description.split(' ', 1)[1].split('_', 1)[0] 
        p = False
        if db_name in seq_dict.keys():
            if acc not in seq_dict[db_name]:
                p = True
        else:
            p = True

        if p:
            added_seq_count += 1
            old_fa_handle.write(get_fa_record_text_from_obj(record))

    print('\nAdded ' + str(added_seq_count) + ' new sequences.')

if __name__ == '__main__':
    add_new_seqs(args.old_fa, args.new_fa)
    
