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
"""Module that defines functions for mapping bootstrap values onto a MrBayes
tree topology.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
import subprocess
import glob
import re
from ete3 import Tree
from settings import raxmlname
from name_replace import write_newick_tree_with_uncoded_names


def reformat_combined_supports(tree_string):
    cs = re.compile(r'\)\d+:')
    for instance in cs.findall(tree_string):
        #print(instance)
        # Define a reformatted support string.
        supcomb = instance[1:-1]
        boot = supcomb[-3:].lstrip('0')
        prob = str(int(supcomb[:-3])/100)
        supcomb2 = prob + '/' + boot

        # Replace the instance with a reformatted support string in the
        # tree_string.
        tree_string = tree_string.replace(instance, instance[0] +\
                supcomb2 + instance[-1])

    # Return modified tree string.
    return tree_string


def combine_supports(boot_newick, prob_newick, combined_figtree_newick):
    """Takes two newick files and writes another newick file that when opened
    in figtree will show posterior probabilities and bootstrap values together.

    ***Should add a function to check whether the two trees being compared are
    rooted on the same node, because that might affect the accuracy of
    comparisons.
    """
    # Define function for adding zeros as necessary.
    get_3_digit = lambda x: '0'*(3 - len(x)) + x
    
    # Parse the input newick tree files.
    boot_newick_tree = Tree(boot_newick)
    prob_newick_tree = Tree(prob_newick)

    # Root trees on the same node arbitrarily.
    arbitrary_leaf_node = prob_newick_tree.get_leaves()[0]
    prob_newick_tree.set_outgroup(arbitrary_leaf_node.name)
    boot_newick_tree.set_outgroup(arbitrary_leaf_node.name)

    # Iterate through the rooted trees matching nodes and combining the support
    # values onto a single tree.
    for n1 in prob_newick_tree.traverse():
        found_boot = False
        if len(n1.get_leaf_names()) > 1:
            for n2 in boot_newick_tree.traverse():
                if set(n1.get_leaf_names()) == set(n2.get_leaf_names()):
                    found_boot = True

                    mb_support = str(n1.support)[:-2]
                    raxml_support = get_3_digit(str(n2.support)[:-2])

                    combined_support = mb_support + raxml_support
                    n1.support = int(combined_support)

        # Make sure that all the right nodes were identified.
        if not found_boot:
            if len(n1.get_leaves()) == 1:
                found_boot = True
        assert found_boot, "Error: could not identify one of the nodes in the"

    # Write the newick tree file to a temporary output file.
    temp_file = combined_figtree_newick + '_temp.newick'
    prob_newick_tree.write(outfile=temp_file)

    # Reformat supports in temporary file with tree so that they will be
    # displayed properly in FigTree.
    with open(temp_file) as infh, open(combined_figtree_newick, 'w') as o:
        for i in infh:
            # Call function to reformat tree string.
            new_tree_string = reformat_combined_supports(i)
            # Write reformatted tree string.
            o.write(new_tree_string)

    # Remove temporary file.
    os.remove(temp_file)


def mbcontre_to_newick_w_probs(intreepath, outtreepath):
    # Construct a taxon number: name dict.
    tax_dict = {}
    x = re.compile(r'\t\t\d')
    with open(intreepath) as infh:
        for i in infh:
            if x.search(i):
                split = i.split('\t')
                num = split[2]
                name = split[3].replace(',', '').strip()
                tax_dict[num] = name

    # Get the line with the tree as a string.
    intreestring = ''
    with open(intreepath) as infh:
        for i in infh:
            if i.startswith('   tree con_all_compat'):
                intreestring = i

    # Cleanup string.
    intreestring2 = '(' + intreestring.split('(', 1)[1]

    # Define REs to identify unnecessary information in tree string.
    extra1_external = re.compile(r'\d\[&prob=[\w\s,.+-={}_()"]+\]')
    extra1 = re.compile(r'\[&prob=[\w\s,.+-={}_()"]+\]')
    extra2 = re.compile(r'\[&length_mean=[\w\s,.+-={}_()"%]+\]')

    intreestring3 = intreestring2

    # Remove terminal node probability labels.
    for i in extra1_external.findall(intreestring3):
        intreestring3 = intreestring3.replace(i, i[0])

    # Reformat internal node probability labels.
    for i in extra1.findall(intreestring3):
        # Calculate a replacement percent probability value.
        replacement = str(round(float(i.split('=', 1)[1].split(',', 1)[0]) * 100))
        #if len(replacement) < 3:
        #    replacement = '0' + replacement

        # Replace long annotation with replacement percent probability.
        intreestring3 = intreestring3.replace(i, replacement)

    # Remove unnecessary branch length annotation.
    intreestring3 = extra2.sub('', intreestring3)
    
    # Construct a list of taxon numbers in string.
    z = re.compile(r'[(|,]\d+:')
    tax_nums = z.findall(intreestring3)

    # Replace taxon numbers with names.
    intreestring4 = intreestring3
    for key in tax_dict.keys():
        for tax_num in tax_nums:
            if key == tax_num[1:-1]:
                intreestring4 = intreestring4.replace(tax_num, tax_num[0] + tax_dict[key] + tax_num[-1])

    # Replace scientific notation numbers to decimal form so RAxML can parse.
    # (May not be necessary)
    #y = re.compile(r':[\w.-]+')
    y = re.compile(r':[\w.-]+\+?\d+')
    numbers = y.findall(intreestring4)
    intreestring5 = intreestring4
    for num in numbers:
        numx = num.strip(':')
        repl = '{:.7f}'.format(float(numx))
        intreestring5 = intreestring5.replace(numx, repl)

    # Write final tree string to output file.
    with open(outtreepath, 'w') as o:
        o.write(intreestring5)

def contre_to_newick(intreepath, outtreepath):
    """Takes a .con.tre file output from MrBayes and writes the same tree in
    newick format.
    """
    # Construct a taxon number: name dict.
    tax_dict = {}
    x = re.compile(r'\t\t\d')
    with open(intreepath) as infh:
        for i in infh:
            if x.search(i):
                split = i.split('\t')
                num = split[2]
                name = split[3].replace(',', '').strip()
                tax_dict[num] = name

    # Get the line with the tree as a string.
    intreestring = ''
    with open(intreepath) as infh:
        for i in infh:
            if i.startswith('   tree con_all_compat'):
                intreestring = i

    # Cleanup string.
    intreestring2 = '(' + intreestring.split('(', 1)[1]

    # Remove unnecessary information from tree string.
    extra1 = re.compile(r'\[&prob=[\w\s,.+-={}_()"]+\]')
    extra2 = re.compile(r'\[&length_mean=[\w\s,.+-={}_()"%]+\]')
    intreestring3 = extra1.sub('', extra2.sub('', intreestring2))

    # Construct a list of taxon numbers in string.
    z = re.compile(r'[(|,]\d+:')
    tax_nums = z.findall(intreestring3)

    # Replace taxon numbers with names.
    intreestring4 = intreestring3
    for key in tax_dict.keys():
        for tax_num in tax_nums:
            if key == tax_num[1:-1]:
                intreestring4 = intreestring4.replace(tax_num, tax_num[0] + tax_dict[key] + tax_num[-1])

    # Replace scientific notation numbers to decimal form so RAxML can parse.
    #y = re.compile(r':[\w.-]+')
    y = re.compile(r':[\w.-]+\+?\d+')
    numbers = y.findall(intreestring4)
    intreestring5 = intreestring4
    for num in numbers:
        numx = num.strip(':')
        repl = '{:.7f}'.format(float(numx))
        intreestring5 = intreestring5.replace(numx, repl)
    
    # Write final tree string to output file.
    with open(outtreepath, 'w') as o:
        o.write(intreestring5)


