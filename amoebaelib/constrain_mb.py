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

"""Contains functions for adding contstraint commands to MrBayes input files.
"""
from datapaths import DataPaths
import sys
import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
from afa_to_nex import delete_extra_mesquite_lines, afa_to_nex, nex_to_afa, nex_to_phylip
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections
import argparse
import subprocess
import time
import datetime
import glob
from ete3 import Tree
import re
from trim_nex import trim_nex
from paralogue_counter import get_seq_obj_from_srch_res_csv_info
import pandas as pd
import column_header_lists

def get_taxon_number_dict(alignment):
    """Take a nexus alignment file and return a dictionary with taxon names as
    keys and their numbers as values.
    """
    taxon_number_dict = {}
    with open(alignment) as infh:
        started = False
        taxon_num = 0
        for i in infh:
            if i.startswith('matrix') or i.startswith('\tMATRIX'):
                started = True
            if i.startswith(';'):
                break
            
            if started and not i.startswith('matrix') and not i.startswith('\tMATRIX'):
                taxon_num += 1
                taxon_name = i.rsplit(' ', 1)[0].strip()
                taxon_number_dict[taxon_name] = taxon_num

    assert taxon_number_dict != {}, """Could not generate a dictionary of taxon
    numbers from nexus alignment file."""

    return taxon_number_dict


def constrain_mb_with_tree(alignment, tree, out_alignment=None):
    """Take a nexus alignment and newick tree and write a new alignment file
    with MrBayes constraint commands added to constrain the tree search to only
    those topologies that match input tree topology (at all nodes).
    """
    # Set output file path.
    if out_alignment is None:
        out_alignment = alignment.rsplit('.', 1)[0] + '_constrained_mb.nex'

    # Check that output file does not exist already.
    assert not os.path.isfile(out_alignment), """Specified output file already
    exists: %s""" % out_alignment

    # Check that output alignment path is not over 99 characters long.
    assert len(out_alignment) > 99, """Alignment file name too long."""
    
    # Get taxon number dict for converting names to numbers.
    taxon_number_dict = get_taxon_number_dict(alignment)

    # Initiate construction of a string of all constraint commands.
    command_prefix = '   constraint '
    constraint_commands = ''

    # Parse input tree using ete3.
    print('\n\n\n')
    print('Topology constraint for mb:')
    print(tree)
    print('\n\n\n')
    t1 = Tree(tree)

    # Count number of leaf nodes in tree and get list of all the leaf names.
    total_leafs = []
    for node in t1.traverse():
        if node.is_leaf():
            total_leafs.append(node.name)
    total_leaf_count = len(total_leafs)

    # Check that the number of keys in the taxon number dict is the same as the
    # number of leaves in the tree (not always necessary).
    #assert total_leaf_count == len(taxon_number_dict.keys()), """Apparently
    #different numbers of sequence names in tree compared to alignment.""" 

    # Check that every leaf name can be coded using the dictionary.
    for node in t1.traverse():
        if node.is_leaf():
            found_match = False
            for key in taxon_number_dict.keys():
                if node.name.strip('\'') == key:
                    found_match = True
            assert found_match, """Could not identify key in dict corresponding
            to the leaf name: %s""" % node.name

    # Iterate over all internal nodes in tree topology.
    node_num = 0
    internal_node_list = []
    for node in t1.traverse():
        if not node.is_leaf():
            node_num += 1
            # Get list of all daughter leaf node names.
            leaf_names = []
            for dnode in node.iter_descendants():
                if dnode.is_leaf():
                    leaf_names.append(dnode.name.strip('\''))

            # Get list as corresponding numbers.
            taxon_numbers = [str(taxon_number_dict[x]) for x in leaf_names]

            # Construct constraint command.
            internal_node_name = 'Node' + str(node_num)
            internal_node_list.append(internal_node_name)
            command = command_prefix + internal_node_name + ' -1 = ' + ' '.join(taxon_numbers) + ';\n'

            # Append new command to constraint command string.
            constraint_commands = constraint_commands + command

    # Add prset command to commands.
    constraint_commands = constraint_commands + '   prset topologypr=constraints(' +\
            ', '.join(internal_node_list) + ');\n'

    # Write output file with constraint commands inserted into the appropriate
    # MrBayes code block.
    with open(alignment) as infh, open(out_alignment, 'w') as o:
        prset_line = re.compile(r'^ +prset')
        inserted = False
        insert = False
        insertnum = 0
        for i in infh:
            #if i.startswith('   prset'):
            if prset_line.match(i):
                insert = True
                o.write(i)
            else:
                if 'filename=' in i:
                    o.write(i.replace('.mb.', '.mb_constrained_mb.'))
                else:
                    o.write(i)

            if insert and insertnum == 0:
                o.write(constraint_commands)
                inserted = True
                insert = False
                insertnum = 1

        # Check that the constraint commands were inserted.
        assert inserted, """Constraint commands were not inserted."""

    # Return main output path.
    return out_alignment

            

