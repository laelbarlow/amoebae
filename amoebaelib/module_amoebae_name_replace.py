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
"""Module defining functions for coding and decoding taxon/sequence names in
alignments and tree files for input to phylogenetics software. Based on the
name_replace.pl script written by Jed Barlow.
"""

import sys
from ete3 import Tree
from module_afa_to_nex import nex_to_afa, afa_to_nex


def generate_code_name(number):
    """Take the number of the sequence in the input file, and output an
    alphabetical code that uniquely corresponds to that number.
    """
    # Check that the number is not greater than the maximum possible to code
    # using this function.
    assert number <= 99999999, """Given sequence number exceeds maximum limit
    for coding function."""

    # Define dictionary for converting digits to letters.
    conversion_dict = {'0': 'Z',
                       '1': 'W',
                       '2': 'R',
                       '3': 'P',
                       '4': 'D',
                       '5': 'J',
                       '6': 'K',
                       '7': 'L',
                       '8': 'V',
                       '9': 'Q'
                       }

    # Convert digits to letters.
    codename = ''
    for i in str(number):
        codename = codename + conversion_dict[i]

    # Add prefix.
    prefix = '' + ('Z' * (8-len(codename)))
    codename = prefix + codename

    # Return codename
    return codename


def write_afa_with_code_names(infile, outfile, tablefile):
    """Take an aligned (or unaligned) fasta file, and write a fasta file with
    coded names, and a table for converting names back.
    """
    original_names = []
    with open(infile) as infh, open(outfile, 'w') as o, open(tablefile, 'w') as tableh:
        inum = -1 
        for i in infh:
            if i.startswith('>'):
                inum += 1
                original_name = i[1:].strip()
                original_names.append(original_name)

                # Define code name.
                code_name = generate_code_name(inum)

                # Write coded name to output sequence file.
                o.write('>' + code_name + '\n')

                # Write values to table file.
                tableh.write(code_name + '\n')
                tableh.write(original_name + '\n')

            else:
                # Write line with sequence to output sequence file.
                o.write(i)

    # Check for redundant sequence names.
    assert len(original_names) == len(list(set(original_names))), """Redundant
    sequence names in input fasta file."""


def write_nex_with_code_names(infile, outfile, tablefile):
    """Take a nexus file, and write a nexus file with
    coded names, and a table for converting names back.
    """
    # Convert nexus to afa.
    afa_path = infile + '_TEMP.afa'
    nex_to_afa(infile, afa_path)

    # Code names.
    afa_coded_path = afa_path.rsplit('.', 1)[0] + '.C.afa'
    write_afa_with_code_names(afa_path, afa_coded_path, tablefile) 

    # Convert afa to nexus.
    afa_to_nex(afa_coded_path, outfile)

    # Remove intermediate afa files.
    os.remove(afa_path)
    os.remove(afa_coded_path)


def get_conversion_dict_from_table(tablefile):
    """Take a conversion table and return a dictionary for uncoding names.
    """
    conv_dict = {}
    with open(tablefile) as tablefileh:
        oddlist = []
        evenlist = []
        inum = 0
        for i in tablefileh:
            if not i.startswith('\n'):
                inum += 1
                if inum % 2 == 0:
                    evenlist.append(i.strip())
                else:
                    oddlist.append(i.strip())
    for code, original in zip(oddlist, evenlist):
        conv_dict[code] = original

    return conv_dict


def write_newick_tree_with_uncoded_names(infile, outfile, tablefile,
        quoted_names=False):
    """Take a text file with a newick tree, and write a new one with names
    converted back to original names (given values in a conversion table), and
    optionally put quotation marks around the names in the output tree.
    """
    # Generate a dictionary for converting names.
    conv_dict = get_conversion_dict_from_table(tablefile)

    ## Look for each code in the input tree, and replace with original name.
    #tree_string = None
    #with open(infile) as infh:
    #    tree_string = infh.readline()

    #for code in conv_dict:
    #    if quoted_names:
    #        tree_string = tree_string.replace(code, '"' + conv_dict[code] + '"')
    #    else:
    #        tree_string = tree_string.replace(code, conv_dict[code])

    #t1 = Tree(infile)
    t1 = Tree(infile, format=1)
    t2 = t1.copy()
    for node in t2.traverse():
        if node.is_leaf():
            for x in conv_dict.keys():
                if node.name.strip('\'').replace(' ', '_') == x.strip('\'').replace(' ', '_'):
                    node.name = conv_dict[x]

    # Write uncoded tree to output file.
    #with open(outfile, 'w') as o:
        #o.write(tree_string)
        #o.write(t2)
    t2.write(outfile=outfile, format=1)


def write_newick_tree_with_coded_names(infile, outfile, tablefile,
        quoted_names=False):
    """Take a newick file and replace full taxon names with coded names from a
    table.
    """
    # Generate a dictionary for converting names.
    conv_dict = get_conversion_dict_from_table(tablefile)

    # Look for each original name in the input tree and replace with the coded
    # name for the output tree.
    tree_string = None
    with open(infile) as infh:
        tree_string = infh.readline()

    for code in conv_dict.keys():
        tree_string = tree_string.replace(conv_dict[code], code)

    # Write uncoded tree to output file.
    with open(outfile, 'w') as o:
        o.write(tree_string)


