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
"""Takes a spreadsheet listing gene/protein names in one column and
accessions/headers in another column, and a file with a list of accessions (for
sequences in a specific clade of a tree), and a name for all the sequences
represented by the list of accessions. Outputs a new version of the input
spreadsheet with updated names. 

Usage:
    update_gene_classifications.py <input spreadsheet path> <accesions list
    filepath> <protein/gene name>

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import argparse


# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Updates gene/protein classifications (names)""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('infp1', help='input spreadsheet path')
parser.add_argument('infp2', help='accessions list filepath')
parser.add_argument('name', help='Protein/Gene name')

args = parser.parse_args()


def update_gene_names(infp1, infp2, name, outfp):
    # Settings.
    column_with_names = 1
    column_with_accessions = 5

    # Loop through input files and write an updated file.
    with open(infp1) as infh1, open(infp2) as infh2, open(outfp, 'w') as o:
        accessions = []
        for j in infh2:
            accessions.append(j.strip())
        for i in infh1:
            if not i.startswith('Protein'):
                spliti = i.split(',')
                original_name = spliti[column_with_names -1]
                accession = spliti[column_with_accessions -1]
                for a in accessions:
                    if a.startswith(accession):
                        spliti[column_with_names -1] = name
                o.write(','.join(spliti))
            else:
                o.write(i)

                
if __name__ == '__main__':

    # Define output filepath.
    outfilepath = args.infp1 + '_updated.csv'

    # Write updated file.
    update_gene_names(args.infp1, args.infp2, args.name, outfilepath)

