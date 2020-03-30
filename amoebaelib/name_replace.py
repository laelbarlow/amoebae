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
import os
import re
import subprocess
from ete3 import Tree
from afa_to_nex import nex_to_afa, afa_to_nex, delete_extra_mesquite_lines, nex_to_phylip, nex_to_mbnex
from string import Template

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
    print('\n\n\n')
    print(open(infile, 'r').read())
    print('\n\n\n')
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
        #tree_string = infh.readline()
        tree_string = infh.read()

    for code in conv_dict.keys():
        tree_string = tree_string.replace(conv_dict[code], code)

    # Write uncoded tree to output file.
    with open(outfile, 'w') as o:
        o.write(tree_string)


def codenames_nex(infp, aamodelpr=None):
    """Takes a path to a multiple sequence alignment file in nexus format and
    writes a new file with coded names.
    """
    # Set default amino acid substitution model to list in MrBayes code block.
    if aamodelpr is None:
        aamodelpr = 'LG4X'

    # Find the number of taxa in the original file.
    ntaxre = re.compile(r"(NTAX|ntax)=\d+")
    ntaxa = None
    with open(infp) as f:
        for l in f:
            if ntaxre.search(l):
                s = ntaxre.search(l)
                ntaxa = int(s.group().split('=')[1])
                break
    assert ntaxa is not None, 'Error: Could not identify number of taxa in input alignment.'
    
    # Make a name for afa file.
    tempafa1 = infp.replace('.nex', '.afa_temp')
    
    # Write an afa file from the nex file.
    nex_to_afa(infp, tempafa1)
    
    # Make name for temporary afa file with new taxon names.
    tempafa2 = tempafa1.replace('.afa_temp', '_C.afa_temp')
    
    # Make a table output file name.
    tablefilename = tempafa2.replace('.afa_temp', '.table')
     
    # Call name_replace.pl (needs to be in $PATH).
    #subprocess.call(['name_replace.pl', '-f', tempafa1, tempafa2, tablefilename]) 
    
    # Code names with module.
    write_afa_with_code_names(tempafa1, tempafa2, tablefilename)
    
    # Check that the number of taxa in the input alignment is the same as the
    # number listed in the output conversion table file.
    num_lines = sum(1 for line in open(tablefilename))
    ntaxa_out = int(num_lines / 2)
    assert ntaxa_out == ntaxa, 'Error: Different number of taxa in input alignment\
     compared to output conversion table.'
    
    # Check for redundant sequence names in output table file.
    lines = [line for line in open(tablefilename)]
    old_names = [line.strip() for line in lines[1::2]]
    old_name_nonredun = list(set(old_names))
    assert len(old_names) <= len(old_name_nonredun), 'Error: Redundant sequences\
     identified in input alignment.'
    
    for name in old_names:
        assert not name.endswith('.copy'), 'Error: duplicate sequences detected.'
    
    # Check that number of sequences is the same in input and in table.
    assert len(old_names) == ntaxa, 'Error: Number of taxa in table is different\
     from number of taxa in input alignment.'
    
    # Make name for new nex file with new taxon names.
    newnex = tempafa2.replace('.afa_temp', '.nex')
     
    # Convert afa file to nex.
    afa_to_nex(tempafa2, newnex)
    
    # Convert nex output file to phylip.
    newphylip = newnex.rsplit('.', 1)[0] + '.phy'
    nex_to_phylip(newnex, newphylip) 
    
    
    # Convert nex output file to nexus for MrBayes:
    
    # Define name for output file.
    mbnex = newnex.rsplit('.', 1)[0] + '.mb.nex'
    
    # Define code blocks for mrbayes input nexus alignment files, using output name.
    # Note: There is probably a better way to handle the variables for MrBayes, but
    # the code blocks can be manually edited after this file is produced.
    mrbayescodeblocks_template =Template(\
    """
    BEGIN mrbayes;
    
       set autoclose=yes;
       log start filename = mbLogFile.log.txt;
       prset aamodelpr=fixed($aamodelpr_var);   lset rates=gamma Ngammacat=4;
       mcmc ngen=10000000 stoprule=Yes stopval=0.01 printfreq=10000 samplefreq=1000 nchains=4 savebrlens=yes filename=$filename_var1;
       log stop;
        quit;
    
    END;
    
    
    
    [
    Begin mrbayes;
    log start filename=$filename_var2.con.log replace;
    sumt filename=$filename_var1 relburnin=yes burninfrac=0.25 contype=allcompat;
    end;
    ]
    """)
    
    mrbayescodeblocks = mrbayescodeblocks_template.substitute(
                                                    aamodelpr_var = aamodelpr,
                                                    filename_var1 = os.path.basename(mbnex),
                                                    filename_var2 = os.path.basename(mbnex.rsplit('.', 1)[0])
                                                    )
    
    # Use code blocks and output name to write mrbayes nexus file.
    nex_to_mbnex(newnex, mbnex, mrbayescodeblocks)
    
     
    # # Delete unnecessary files.
    os.remove(tempafa1)
    os.remove(tempafa2)


