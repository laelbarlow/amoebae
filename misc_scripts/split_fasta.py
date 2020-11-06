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
"""This script is for splitting fasta files into individual files for each
sequence.  It is assumed that the sequence headers will be in the format
>ACCESSION DBNAME_PROTNAME

Advantage over split_fast_into_query_files.py is that it reads the database of
origin from each sequence header, instead of relying on input.

Usage:
    split_fasta.py <fasta file> 
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import subprocess
import argparse
from Bio import SeqIO
from afa_to_nex import nex_to_afa, afa_to_nex
import shutil


# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Outputs positive hits and sequences.""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('infp1', help='fasta file to be split.')

# Optional arguments.
#parser.add_argument('-p', '--program', help='Program (either muscle or mafft)',
#        default='muscle')
#parser.add_argument('-m', '--aamatrix', help='Amino acid matrix (NCBI format)',
#        default=None)

args = parser.parse_args()



def split_fasta(infilepath, subdirname):
    """Takes a fasta files and writes new fasta files in a subdirectory for
    each sequence in the fasta file.
    """
    # Loop through the sequences in the input fasta file
    with open(infilepath) as i:
        fasta_sequences = SeqIO.parse(i,'fasta')

        print('\nSplitting fasta file: ' + infilepath + '\n')

        num = 0
        for seq in fasta_sequences:
            num += 1
            split_header = seq.description.split(' ', 1)
            acc = split_header[0]
            acc = seq.description
            #dbname = split_header[1].split('_', 1)[0]
            #protname = split_header[1].split('_', 1)[1]

            #filename = dbname + '_' + acc[:25] + '_' + protname
            #filename = dbname + '_' + protname + '_' + acc[:25]
            #filename = 'QcSNARE2' + '_' + acc[:25].replace('|', '_').replace(' ', '_').replace('/', '_') + '_2'
            filename = acc[:25].replace('|', '_').replace(' ', '_').replace('/', '_') + '_2'

            print('Writing file: ' + filename + '.fa')
            
            seqfilename = os.path.join(subdirname, filename + '.fa')
            with open(seqfilename, 'w') as seqf:
                SeqIO.write([seq],seqf,'fasta')

    
def do_split_fasta(infapath):
    """Do the whole process
    """
    # Make a subdirectory for storing individual sequence files.
    subdirname = infapath.rsplit('.', 1)[0] + '_queries'
    os.mkdir(subdirname)

    # Split the input fasta file.
    split_fasta(infapath, subdirname)


def split_fasta2(infilepath, subdirname):
    """Takes a fasta files and writes new fasta files in a subdirectory for
    each sequence in the fasta file.
    """
    # Loop through the sequences in the input fasta file
    with open(infilepath) as i:
        fasta_sequences = SeqIO.parse(i,'fasta')

        print('\nSplitting fasta file: ' + infilepath + '\n')

        num = 0
        for seq in fasta_sequences:
            num += 1
            split_header = seq.description.split(' ', 1)
            acc = split_header[0]
            description = split_header[1].split('_', 1)[0]
            #protname = split_header[1].split('_', 1)[1]

            filename = 'Mbalamuthi' + '_' + acc[:25] # + '_' + protname

            print('Writing file: ' + filename + '.fa')
            
            seqfilename = os.path.join(subdirname, filename + '.fa')
            with open(seqfilename, 'w') as seqf:
                SeqIO.write([seq],seqf,'fasta')

    
def do_split_fasta2(infapath):
    """Do the whole process
    """
    # Make a subdirectory for storing individual sequence files.
    subdirname = infapath.split('.', -1)[0] + '_queries'
    os.mkdir(subdirname)

    # Split the input fasta file.
    split_fasta2(infapath, subdirname)


def split_fasta3(infilepath, subdirname):
    """Takes a fasta files and writes new fasta files in a subdirectory for
    each sequence in the fasta file.
    """

    # Define dictionary with query titles as keys and uniprot IDs as values.
    ne_query_titles = {}
    #with open('query_title_id.csv', 'r', encoding='utf8') as infh:
    with open('query_title_id.csv', 'r') as infh:
        for i in infh:
            s = i.split(',')
            query_title = s[0].strip()
            seq_id = s[1].strip()
            ne_query_titles[query_title] = seq_id


    print('\nSplitting fasta file: ' + infilepath + '\n')

    for title in ne_query_titles.keys():
        print('\n' + title)
        query_title = None

        # Loop through the sequences in the input fasta file
        with open(infilepath) as i:
            fasta_sequences = SeqIO.parse(i,'fasta')
            num = 0
            for seq in fasta_sequences:
                num += 1
                split_header = seq.description.split(' ', 1)
                #acc = split_header[0]
                #acc = seq.description
                acc = seq.description.split('|')[1]
                print('\t' + acc)
                organism = seq.description.split(' ')[0].rsplit('_', 1)[1]
                #dbname = split_header[1].split('_', 1)[0]
                #protname = split_header[1].split('_', 1)[1]

                if ne_query_titles[title] == acc:
                    print('\t' + 'match')
                    query_title = title

                    #filename = dbname + '_' + acc[:25] + '_' + protname
                    #filename = dbname + '_' + protname + '_' + acc[:25]
                    #filename = 'QcSNARE2' + '_' + acc[:25].replace('|', '_').replace(' ', '_').replace('/', '_') + '_2'
                    #filename = acc[:25].replace('|', '_').replace(' ', '_').replace('/', '_') + '_2'
                    filename = query_title + '_' + acc + '_' + organism

                    print('Writing file: ' + filename + '.fa')
                    
                    seqfilename = os.path.join(subdirname, filename + '.faa')
                    with open(seqfilename, 'w') as seqf:
                        SeqIO.write([seq],seqf,'fasta')

                if query_title is not None:
                    break

        assert query_title is not None, """Error: No sequence found for
        query title %s and ID %s.""" % (title, ne_query_titles[title])


    
def do_split_fasta3(infapath):
    """Do the whole process
    """
    # Make a subdirectory for storing individual sequence files.
    subdirname = infapath.rsplit('.', 1)[0] + '_queries'
    os.mkdir(subdirname)

    # Split the input fasta file.
    split_fasta3(infapath, subdirname)


if __name__ == '__main__':

    do_split_fasta(args.infp1)
    #do_split_fasta2(args.infp1)
    #do_split_fasta3(args.infp1)
