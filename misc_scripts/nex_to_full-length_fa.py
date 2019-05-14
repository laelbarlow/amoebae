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
"""This script is for reading a nexus alignment file with trimmed sequences and
writing a fasta file with corresponding full-length sequences.

Usage:
    nex_to_full-length_fa.py <nexus file>

Note: The sequence names must be formatted as: ACCESSION DBNAME_OTHERINFO
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import argparse
import settings
from module_get_fas_from_db_dir import get_fas_from_db_dir
from Bio import AlignIO
import module_afa_to_nex

# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Gets full-length fasta sequences for sequences in nexus alignment.""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('infp1', help='infilepath')

# Optional arguments.
#parser.add_argument('-f', '--fwdeval', help='evalue cutoff for fHMMer hits', default=0.05)
#parser.add_argument('-r', '--reveval', help='evalue cutoff for rBLAST hits', default=0.05)
#parser.add_argument('-s', '--getseqs', help='write sequences to a file.',
#        action='store_true')
#parser.add_argument('-p', '--prot_name', help='name of protein that was \
#searched for', default=None)

args = parser.parse_args()
    

def get_out_name(infilepath):
    """Modifies the input file name to make the output file name.
    """
    return infilepath.replace('.nex', '_full_length.fa')

def write_pos_seqs(infilepath, dbdirpath, outfilepath, prot_name=None):
    """Uses the get_fas_from_db_dir module to get all the sequences
    corresponding to hits in a nexus alignment file.
    """
    with open(infilepath) as infh, open(outfilepath, 'w') as o:
        id_list = []
        alignment = AlignIO.read(infh, 'nexus')
        for i in alignment:
            print(i.id)
            acc = None
            dbname = None
            if ' ' in i.id:
                splitid = i.id.split(' ')
                acc = splitid[0]
                if '_' in splitid[1]:
                    dbname = splitid[1].split('_', 1)[0]
                else:
                    dbname = splitid[1]
            else:
                splitid = i.id.split('_', 1)
                acc = splitid[0]
                if '_' in splitid[1]:
                    dbname = splitid[1].split('_', 1)[0]
                else:
                    dbname = splitid[1]
            acc_list = [acc]
            id_list.append([dbname, acc_list])
         
        for l in id_list:
             o.write(get_fas_from_db_dir(l[0], l[1], dbdirpath,
                 prot_name))


if __name__ == '__main__':
    dbdirpath = settings.dbdirpath
    outfp = get_out_name(args.infp1)
    prot_name = input('\nProtein name: ')

    # Remove any additional mesquite blocks.
    module_afa_to_nex.delete_extra_mesquite_lines(args.infp1)

    # Get necessary sequences and write to file.
    write_pos_seqs(args.infp1, dbdirpath, outfp, prot_name)


