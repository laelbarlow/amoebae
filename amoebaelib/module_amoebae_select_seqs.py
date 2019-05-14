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

"""Contains functions for modifying sequence selection for a given phylogenetic
tree.
"""
import settings
import sys
import os
import numpy as np
import collections
import argparse
import subprocess
import time
import datetime
import glob
import shutil
import re
import statistics
import pandas as pd
from random import shuffle

from module_amoebae_trim_nex import trim_nex
from module_paralogue_counter import get_seq_obj_from_srch_res_csv_info
import module_amoebae_column_header_lists
from module_amoebae_phylo_clas import ModelInfoFromCSV,\
get_clade_name_from_model, code_names_in_ali, quote_tree, code_tree,\
uncode_tree, uncode_tree_obj
from module_paralogue_counter import add_seq_to_alignment3

# Import functions for working with sequences, alignments, and trees.
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
from module_afa_to_nex import delete_extra_mesquite_lines, afa_to_nex, nex_to_afa, nex_to_phylip
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace


def get_clade_name_from_model2(type_seq_name, type_seqs_dict):
    """Return clade name for given type sequence name and type seqs info file
    (from models directory).
    """
    # Get clade name.
    clade_name = '?'
    for clade in type_seqs_dict.keys():
        if type_seq_name == type_seqs_dict[clade]:
                clade_name = clade
                break
        elif type_seq_name == type_seqs_dict[clade].replace(' ', '_'):
            clade_name = clade
            break
    assert clade_name != '?', """Could not identify clade name for type
    sequence name %s.""" % type_seq_name

    return clade_name


def get_nodes_from_actual_tree_obj(t3, nodes_of_interest):
    """Get nodes in the given tree object that correspond to the nodes in a
    given list of nodes from another tree object, and return as a list.
    """
    # Get original root node for tree object.
    original_outgroup = list(t3.traverse())[1]
    # Compile a list of corresponding nodes.
    new_node_list = []
    for i in nodes_of_interest:
        i_leaf_names = set([x.name for x in i.get_leaves()])
        corresponding_node = None
        # Traverse the tree to find node.
        for node in t3.traverse():
            node_leaf_names = [x.name for x in node.get_leaves()]
            # If the node has the same set of leaf names, then it is the
            # corresponding node.
            if node_leaf_names == i_leaf_names:
                corresponding_node = node
                break
        if corresponding_node is None:
            # Generate a list of three different rootings for the tree, at
            # least one of which should allow identification of the clade of
            # interest.
            roots = []
            roots.append(t3.get_midpoint_outgroup())
            t3.set_outgroup(t3.get_midpoint_outgroup())
            side1 = list(t3.traverse())[1]
            side2 = list(t3.traverse())[2]
            for side in [side1, side2]:
                side_longest_length = 0
                side_longest_leaf = None
                for leaf in side.get_leaves():
                    length = side.get_distance(leaf)
                    if length > side_longest_length:
                        side_longest_length = length
                        side_longest_leaf = leaf
                # Check that a leaf was identified.
                assert side_longest_leaf is not None
                # Add the longest leaf on this side of the midpoint to the list
                # of nodes on which to try rooting the tree.
                roots.append(side_longest_leaf)

            # Try finding the corresponding node again.
            for root_node in roots:
                t3.set_outgroup(root_node)
                for node in t3.traverse():
                    node_leaf_names = set([x.name for x in node.get_leaves()])
                    # If the node has the same set of leaf names, then it is the
                    # corresponding node.
                    if node_leaf_names == i_leaf_names:
                        corresponding_node = node
                        break

            # Check that the corresponding node was found.
            assert corresponding_node is not None

            # Root on the corresponding node identified.
            t3.set_outgroup(corresponding_node)

        # Check that the corresponding node was found.
        assert corresponding_node is not None

        # Append the corresponding node to the list.
        new_node_list.append(corresponding_node)

    # Re-root the tree object on the original root node.
    t3.set_outgroup(original_outgroup)

    # Return the list of corresponding nodes.
    return new_node_list




def get_nodes_of_interest(tree, type_seq_list):
    """Take an ete3 tree (TreeNode object), a list of type sequence names, and
    return a list of nodes of interest defined by clades (as ete3 TreeNode
    objects).
    """
    # Initiate list of orthogroup nodes of interest.
    orthogroup_nodes = []

    # Make a copy of the input tree.
    t2 = tree.copy()

    # For each "type" sequence, traverse all nodes and find the node with
    # the largest number of child nodes that are leaf (terminal) nodes,
    # containing the "type" sequence of interest, but not containing any of
    # the other "type" sequences.
    ts_num = 0
    first_type_seq_node_name = None
    for ts in type_seq_list:
        ts_num += 1

        if ts_num == 1:
            first_type_seq_node_name = ts
            # Root on another "type" sequence for the first type sequence in
            # the list to get whole clade, then root the tree on the ancestor
            # node of that first clade.

            # Get a node name for a node corresponding to a different
            # "type" sequence.
            other_type_seq_node_name = None
            for i in type_seq_list:
                if i != ts:
                    other_type_seq_node_name = i
                    break
            assert other_type_seq_node_name is not None

            # Get node corresponding to a different "type" sequence.
            other_type_seq_node = None
            for node in t2.traverse():
                if node.name == other_type_seq_node_name:
                    other_type_seq_node = node
                    break
            assert other_type_seq_node is not None, """Could not identify
            a type seq node with name: %s""" % other_type_seq_node_name

            # Root on the other "type" sequence node.
            t2.set_outgroup(other_type_seq_node)

            #print('\n\n\nTree rooted on a type sequence other than the first type sequence.')
            #print(t2)

        elif ts_num == 2:
            # Root on the first "type" sequence node for all subsequent
            # clades.
            first_type_seq_node = None
            for node in t2.traverse():
                leaf_list = []
                for leaf in node.get_leaves():
                    if leaf.name == first_type_seq_node_name:
                        first_type_seq_node = node
                        break
            t2.set_outgroup(first_type_seq_node)
            #print('\n\n\nTree re-rooted on first type sequence:')
            #print(t2)


        # Make a copy of the tree topology to work with for each run
        # through this loop.
        t3 = t2.copy()

        # Make a list of nodes that contain type seq, but not any others. #
        # SEPARATE
        nodes_of_interest = []
        for node in t3.traverse():
            # Search in nodes that contain the type sequence.
            if node.search_nodes(name=ts):
                # Search in nodes that don't contain other type sequences.
                contains_other_type_seqs = False
                for ts2 in type_seq_list:
                    if not ts2 == ts:
                        if node.search_nodes(name=ts2):
                            contains_other_type_seqs = True
                if not contains_other_type_seqs:
                    # Add nodes of interest to list.
                    nodes_of_interest.append(node)

        # Check that a node of interest was found.
        assert len(nodes_of_interest) >= 2

        # find the node with the most child leaf nodes.
        node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
                len(x.get_leaves()), reverse=True)[0]
        #node_w_most_leaves.name = 'X'

        # Add node to list of nodes of interest.
        orthogroup_nodes.append([ts, node_w_most_leaves])

    # Return list of lists containing defining sequence name and node object.
    return orthogroup_nodes


def get_list_of_leaf_names_for_node(node):
    """Given an ete3 TreeNode object, return a list of terminal/leaf node
    names.
    """
    leaf_names = []
    for n in node.traverse():
        if n.is_leaf():
            leaf_names.append(n.name)
    return leaf_names


class TaxonomicInfo():
    """Class for storing taxonomic info.
    """
    def __init__(self, superbranch, supergroup, group, species):
        self.superbranch = superbranch
        self.supergroup = supergroup
        self.group= group
        self.species = species


def get_taxonomic_info(species_name):
    """Given a species name that appears in the genome info csv file in the
    genomes folder, return taxonomic information from that csv file.
    """
    # Define path to info file.
    info_csv = settings.db_info_csv

    # Parse info file with pandas.
    df = pd.read_csv(info_csv)

    # Iterate through rows until one is found with the same species name.
    taxonomic_info_obj = None
    found_row = False
    species = None
    for index, row in df.iterrows():
        species = row['Species (if applicable)']
        for x in range(0, 5):
            parsed_species_name = species_name.strip().split(' ', x)[-1]
            if species.strip().strip('.') == parsed_species_name:
                found_row = True
                superbranch = row['Superbranch']
                supergroup = row['Supergroup']
                group = row['Group']
                taxonomic_info_obj = TaxonomicInfo(superbranch, supergroup, group, species)
            if found_row:
                break
        if found_row:
            break


    assert taxonomic_info_obj is not None, """Could not identify taxonomic
    information for species %s""" % species_name

    return taxonomic_info_obj


def define_nodestyles_dict_for_colourcoding():
    """Define node styles for different nodes of interest that need to be
    colour-coded.
    """
    # Initiate dict.
    nodestyle_dict = {}

    # General default style for tree (apply using the node.set_style(style)
    # function).
    style = NodeStyle()
    style["fgcolor"] = "#0f0f0f"
    style["size"] = 0
    style["vt_line_color"] = "#ff0000"
    style["hz_line_color"] = "#ff0000"
    style["vt_line_width"] = 8
    style["hz_line_width"] = 8
    style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0

    # dashed blue lines
    style1 = NodeStyle()
    style1["fgcolor"] = "#000000"
    style1["shape"] = "circle"
    style1["vt_line_color"] = "#0000aa"
    style1["hz_line_color"] = "#0000aa"
    style1["vt_line_width"] = 2
    style1["hz_line_width"] = 2
    style1["vt_line_type"] = 1 # 0 solid, 1 dashed, 2 dotted
    style1["hz_line_type"] = 1

    # stem_node
    stem_node = NodeStyle()
    stem_node["hz_line_width"] = 5

    # bold lines
    style2 = NodeStyle()
    #style2["vt_line_width"] = 8
    style2["hz_line_width"] = 5

    # bold dashed lines
    style3 = NodeStyle()
    style3["hz_line_width"] = 5
    style3["hz_line_type"] = 1


    # Set dashed blue lines in all leaves
    #style3 = NodeStyle()
    #style3["fgcolor"] = "#000000"
    #style3["shape"] = "circle"
    #style3["vt_line_color"] = "#0000aa"
    #style3["hz_line_color"] = "#0000aa"
    #style3["vt_line_width"] = 2
    #style3["hz_line_width"] = 2
    #style3["vt_line_type"] = 1 # 0 solid, 1 dashed, 2 dotted
    #style3["hz_line_type"] = 1


    # Add styles to dictionary.
    nodestyle_dict['default'] = style1
    nodestyle_dict['stem_node'] = stem_node
    nodestyle_dict['shortest branch'] = style2 # bold.
    nodestyle_dict['longest branch'] = style3 # dashed.
    nodestyle_dict['default'] = style1

    # Return dict.
    return nodestyle_dict


def define_textface_for_labeling_stem(clade_name):
    """Define ete3 textface styles for different nodes of interest that need to be
    labeled.
    """
    # Define TextFace objects for clade stem nodes/branches.
    stem_clade = TextFace(clade_name)

    # Set some attributes
    #stem_clade.margin_top = 10
    #stem_clade.margin_right = 10
    #stem_clade.margin_left = 10
    #stem_clade.margin_bottom = 10
    #stem_clade.opacity = 0.5 # from 0 to 1
    #stem_clade.inner_border.width = 1 # 1 pixel border
    #stem_clade.inner_border.type = 1  # dashed line
    #stem_clade.border.width = 1
    #stem_clade.background.color = "LightGreen"

    # Return dictionary.
    return stem_clade


def get_corresponding_node(node, tree):
    """Return the ete3 TreeNode object in a tree that corresponds to
    a given node.
    """
    # Search for a node with child leaf nodes that have the same set of names.
    corresponding_node = None
    node_leaves = [x.name for x in node.get_leaves()]
    for n in tree.traverse():
        if set([y.name for y in n.get_leaves()]) == set(node_leaves):
            corresponding_node = n
            break
    # If necessary, re-root input tree on midpoint for parsing.
    if corresponding_node is None:
        # Re-root on midpoint.
        tree.set_outgroup(tree.get_midpoint_outgroup())
        # Search for node of interest again.
        for n in tree.traverse():
            if set([y.name for y in n.get_leaves()]) == set(node_leaves):
                corresponding_node = n
                break

    # Try rooting on random leaf nodes until the corresponding node can be
    # found.
    if corresponding_node is None:
        # Loop through random leaf nodes in tree, rooting on them, and then
        # searching for the corresponding node.
        for l in tree.iter_leaves():
            tree.set_outgroup(l)
            # Search for node of interest again.
            found = False
            for n in tree.traverse():
                if set([y.name for y in n.get_leaves()]) == set(node_leaves):
                    corresponding_node = n
                    found = True
                    break
            if found:
                break

    # Check that a corresponding node was identified.
    assert corresponding_node is not None, """Could not find corresponding
    node."""

    # Return corresponding node.
    return corresponding_node



def optimize_sequence_selection(file_with_seqs, model_name, essential_taxa):
    """Add/remove homologous sequences from a dataset to optimize for long stem
    lengths and low average branch lengths for clades of interest.

    The essential_taxa argument is a list of names corresponding to taxonomic
    groups of eukaryotes (such as "Amorphea" or "Diaphoretickes") which appear
    in the genome info csv file in the genomes directory.
    """
    # (make this into a recursive function?)

    # Define a dictionary of ete3 NodeStyle objects for different branch types
    # in the tree to be visualized.
    style_dict = define_nodestyles_dict_for_colourcoding()

    # Get info about model that is relevant.
    model_info = ModelInfoFromCSV(model_name)
    alignment = model_info.alignment_file
    tree = model_info.tree_topology_file
    subs_model = model_info.subs_model
    type_seqs = model_info.type_seqs_file

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(tree, quoted_node_names=False)

    # Print tree.
    print('Input tree:')
    print(t1)

    # Get a list of nodes of interest (ancestral nodes for clades of interest).
    orthogroup_nodes = []

    # Get list of "type" sequences from input.
    type_seq_list = []
    for i in open(type_seqs):
        type_seq_list.append(i.strip().split(',')[0])

    # Copy tree for colour-coding etc.
    t2 = t1.copy()

    # Get list of nodes of interest.
    orthogroup_nodes = get_nodes_of_interest(t1, type_seq_list)

    # Determine characteristics of clades of interest.
    num = 0
    for seq_node_list in orthogroup_nodes:
        num += 1
        type_seq = seq_node_list[0]
        clade_name = get_clade_name_from_model(type_seq, type_seqs)
        node = seq_node_list[1]

        print('\n\nClade defined by sequence ' + type_seq + ':')
        print(node)

        # Get list of leaf names.
        leaf_names = get_list_of_leaf_names_for_node(node)

        # Get stem length for node.
        stem_len = node.dist
        print('\nStem length: ' + str(stem_len))

        # Get branch lengths for all sequences in the clade.
        print('\nBranch lengths:')
        seq_branch_lengths = []
        for n in node.traverse():
            if n.is_leaf():
                length = node.get_distance(n)
                print('\t' + n.name + ': ' + str(length))
                seq_branch_lengths.append([n.name, length])

        # Get longest branch length.
        max_len = max(seq_branch_lengths, key=lambda x: x[1])
        print('\n\tLongest branch: ' + max_len[0] + ' ' + str(max_len[1]))

        # Get shortest branch length.
        min_len = min(seq_branch_lengths, key=lambda x: x[1])
        print('\n\tShortest branch: ' + min_len[0] + ' ' + str(min_len[1]))

        # Get average branch length relative to ancestral/root node in the
        # clade for each sequence/terminal/leaf node.
        mean_len = statistics.mean([x[1] for x in seq_branch_lengths])
        print('\n\tMean length: ' + str(mean_len))

        # Get median branch length relative to ancestral/root node in the
        # clade for each sequence/terminal/leaf node.
        median_len = statistics.median([x[1] for x in seq_branch_lengths])
        print('\n\tMedian length: ' + str(median_len))

        # Calculate ratio of stem length to average branch length.
        stem_branches_ratio = stem_len / mean_len
        print('\n\tRatio of stem len to avg branch len: ' + str(stem_branches_ratio))


        # Format tree for visualization.

        # Root on first clade of interest.
        if num == 1:
            t2.set_outgroup(get_corresponding_node(node, t2))

        # Make stem branch for clade of interest bold.
        stem_node = get_corresponding_node(node, t2)
        stem_node.set_style(style_dict['stem_node'])

        # Add clade name to stem branch.
        #stem_node.name = clade_name
        stem_node.add_face(define_textface_for_labeling_stem(clade_name), column=0, position = "branch-top")

        # Iterate over branches within clade and customize appearance.
        for leaf in node.iter_leaves():
            # Set general features for all leaf nodes.
            style = NodeStyle()
            style["fgcolor"] = "#0f0f0f"
            style["size"] = 0
            style["vt_line_color"] = "#000000"
            style["hz_line_color"] = "#000000"

            #style["vt_line_width"] = 8
            #style["hz_line_width"] = 8
            #style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
            #style["hz_line_type"] = 0

            # Make shortest branch in clade of interest bold.
            if leaf.name == min_len[0]:
                #shortbranchnode = get_corresponding_node(t1.search_nodes(name=min_len[0])[0], t2)
                #shortbranchnode.set_style(style_dict['shortest branch'])
                style["hz_line_width"] = 5

            # Make longest branch in clade of interest bold.
            if leaf.name == max_len[0]:
                #longbranchnode = get_corresponding_node(t1.search_nodes(name=max_len[0])[0], t2)
                #longbranchnode.set_style(style_dict['longest branch'])
                style["hz_line_width"] = 5
                style["hz_line_type"] = 2
            
            # Customize leaf node appearance.
            # Get species name.
            species_name = leaf.name.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
            print(species_name)
            ti = get_taxonomic_info(species_name)
            print('\t' + ' '.join([str(ti.superbranch), str(ti.supergroup),
                str(ti.group), str(ti.species)]))

            # Set style for node.
            get_corresponding_node(leaf, t2).set_style(style)


    # Show tree with colour-coded branches.

    # Remove underscores from leaf names.
    for leaf in t2.iter_leaves():
        leaf.name = leaf.name.replace('__', ' ').replace('_', ' ')

    # Stretch branches.
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale =  120 # 120 pixels per branch length unit

    # Add title.
    tree_title = "[Tree title here]"
    ts.title.add_face(TextFace(tree_title, fsize=20), column=0)

    # Show.
    t2.show(tree_style=ts)

    # Write tree to file.
    t2.render("/Users/Lael/Desktop/mytree.png", w=183, units="mm", tree_style=ts)

    # Return main output path.
    return '[no output path defined]'


def get_ml_tree_branch_lengths(tree, alignment, subs_model, outputdir):
    """Run an ML tree with a constrained topology given a tree (for
    constraining topology), an alignment file, and a substitution model.
    """
    # Convert nexus alignment to afa.
    temp_fa_1 = os.path.join(outputdir,\
            os.path.basename(alignment).rsplit('.', 1)[0] + '_temp1.afa')
    nex_to_afa(alignment, temp_fa_1) 

    # Code names in alignment.
    outalifp = alignment.rsplit('.', 1)[0] + '.C.afa'
    outtablefp = alignment.rsplit('.', 1)[0] + '.C.table'
    code_names_in_ali(outputdir, temp_fa_1, outalifp, outtablefp) 

    # Convert alignment back to nexus.
    outalifpnex = alignment.rsplit('.', 1)[0] + '.C.nex'
    afa_to_nex(outalifp, outalifpnex)

    # Delete temporary files.
    os.remove(temp_fa_1)
    os.remove(outalifp)

    # Define path for intermediate unrooted tree.
    tree2 = outalifpnex.rsplit('.')[0] + '_unrooted_tree.newick'

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(tree, quoted_node_names=False)

    # Unroot tree.
    t1.unroot()

    # Write unrooted tree to a new file.
    t1.write(outfile=tree2, format=9)

    # Define path for constraint tree file.
    constraint_tree_fp_coded = None
    # Only do the following steps once for all the alignments, because all the
    # constrained topology analyses will use the same constraint trees.

    # Put quotation marks around names in constraint tree. # SEPARATE
    #outtreefp = ali.rsplit('.', 1)[0] + '.constraint_tree0.Q.newick'
    constraint_tree_fp = outalifpnex.rsplit('.')[0] + '_constraint_tree.newick'
    with open(tree2) as intreefh, open(constraint_tree_fp, 'w') as o:
        for t in intreefh:
            o.write(quote_tree(t, outtablefp))
    # Make constraint tree with coded names.
    constraint_tree_fp_coded = constraint_tree_fp.rsplit('.', 1)[0] + '.C.tre'
    code_tree(constraint_tree_fp, constraint_tree_fp_coded, outtablefp)

    # Make a copy of the alignment in phylip format for input to IQtree.
    phy_out = outalifpnex.rsplit('.', 1)[0] + '.phy'
    nex_to_phylip(outalifpnex, phy_out)

    # Do phylogenetic analysis.

    # Make subdir for output.
    subdirp = constraint_tree_fp.rsplit('.', 1)[0] + '_IQ-tree_output'
    os.mkdir(subdirp)

    # Use IQtree to do an ML search
    output_file_prefix = os.path.join(subdirp, 'iqtree')
    iqtree_command_list = ['iqtree', '-s', phy_out, '-m', subs_model, '-g',
        constraint_tree_fp_coded, '-pre', output_file_prefix, '-nt', '4']
    stdout_path = output_file_prefix + '.stdout.txt'
    tree_search_start_time = time.time()
    with open(stdout_path, 'w') as o:
        subprocess.call(iqtree_command_list, stdout=o, stderr=subprocess.STDOUT)
    #tree_search_end_time = time.time()
    #with open(logfile, 'a') as logh:
    #    logh.write(' '.join(iqtree_command_list) + '\n')
    #    tree_search_elapsed = tree_search_end_time - tree_search_start_time
    #    logh.write('Run time: ' + str(datetime.timedelta(seconds=tree_search_elapsed)) + '\n')

    # Define path to output file.
    tree_file_path = output_file_prefix + '.treefile'
    # Check that the output file was actually produced.
    assert os.path.isfile(tree_file_path)

    # Parse output tree.
    # Note: parentheses and commas get replaced with underscores.

    # Uncode tree.
    tree_file_path2 = tree_file_path + '2'
    uncode_tree(tree_file_path, tree_file_path2, outtablefp)

    #with open(tree_file_path2) as infh:
    #    for i in infh:
    #        print(i)


    #with open(tree_file_path) as infh, open(tree_file_path2, 'w') as o:
    #    for t in infh:
    #        # Add quotation marks around node names and write to 2nd file.
    #        print('old_tree_text')
    #        print(t)
    #        new_tree_text = quote_tree(t, outtablefp)
    #        print('new_tree_text')
    #        print(new_tree_text)
    #        o.write(new_tree_text)
            
    # Parse new tree file.
    t2 = Tree(tree_file_path2, format=1, quoted_node_names=True)


    # Return the tree.
    return t2 


def get_branch_length_info(t1, type_seqs_dict, outpath=None):
    """Take an ete3 Tree object and a file with type seqs, and return a
    dictionary with keys as type sequence names and values as branch length
    ratios for the respective clades in the input tree.
    """
    # Define a dictionary of ete3 NodeStyle objects for different branch types
    # in the tree to be visualized.
    style_dict = define_nodestyles_dict_for_colourcoding()

    # Get a list of nodes of interest (ancestral nodes for clades of interest).
    orthogroup_nodes = []

    # Get list of "type" sequences from input.
    #type_seq_list = []
    #for i in open(type_seqs):
    #    type_seq_list.append(i.strip().split(',')[0])
    type_seq_list = type_seqs_dict.values()

    # Copy tree for colour-coding etc.
    t2 = t1.copy()
    
    # Root on the midpoint, because ete3 was having trouble identifying clades
    # in some of the unrooted trees.
    t2.set_outgroup(t2.get_midpoint_outgroup())

    # Get list of nodes of interest.
    orthogroup_nodes = get_nodes_of_interest(t1, type_seq_list)


    # Get lengths for internal branches (excluding those within the specific
    # clades of interest) and add those to the dictionary as a list.

    # Traverse all the nodes in the tree, excluding those in the clades of
    # interest.
    internal_branch_lengths = []
    # Copy the tree object.
    t3 = t2.copy()
    # Record initial number of nodes.
    initial_node_num = len(list(t3.traverse()))

    # Record which node is the original outgroup.???
    original_outgroup = list(t3.traverse())[1]

    # Find nodes of interest in new tree object for detaching.
    orthogroup_nodes_t3 = get_nodes_from_actual_tree_obj(t3, [x[1] for x in orthogroup_nodes])

    # Re-root on one of the nodes of interest.
    t3.set_outgroup(orthogroup_nodes_t3[0])

    #print('\nt3 before detaching nodes of interest:')
    #print(t3)

    # Detach all nodes of interest in the copy of the tree.
    #for i in to_remove_list:
    for i in orthogroup_nodes_t3:
        for j in i.iter_descendants():
            j.detach()
    #print('\ntree with nodes removed (is this right???***):')
    #print(t3)

    # Record (reduced) number of nodes.
    reduced_node_num = len(list(t3.traverse()))

    # Check that at least one node was actually removed.
    assert reduced_node_num < initial_node_num
    # Check that the tree has the minimum number of nodes.
    assert reduced_node_num >= len(orthogroup_nodes)

    # Get internal branch lengths from pruned copy of the tree.
    for node in t3.traverse():
        internal_branch_lengths.append(node.dist)

    # Check that a minimum number of internal nodes were identified.
    #assert len(internal_branch_lengths) >= len(orthogroup_nodes), """Too few
    #internal nodes identified."""


    # Initiate output dictionary.
    branch_length_info_dict = {}

    # Determine characteristics of clades of interest.
    num = 0
    for seq_node_list in orthogroup_nodes:
        num += 1
        type_seq = seq_node_list[0]
        clade_name = get_clade_name_from_model2(type_seq, type_seqs_dict)
        node = seq_node_list[1]

        #print('\n\nClade defined by sequence ' + type_seq + ':')
        #print(node)

        # Get list of leaf names.
        leaf_names = get_list_of_leaf_names_for_node(node)

        # Get stem length for node.
        stem_len = node.dist
        #print('\nStem length: ' + str(stem_len))

        # Get branch lengths for all sequences in the clade.
        #print('\nBranch lengths:')
        seq_branch_lengths = []
        for n in node.traverse():
            if n.is_leaf():
                length = node.get_distance(n)
                #print('\t' + n.name + ': ' + str(length))
                seq_branch_lengths.append([n.name, length])

        # Get longest branch length.
        max_len = max(seq_branch_lengths, key=lambda x: x[1])
        #print('\n\tLongest branch: ' + max_len[0] + ' ' + str(max_len[1]))

        # Get shortest branch length.
        min_len = min(seq_branch_lengths, key=lambda x: x[1])
        #print('\n\tShortest branch: ' + min_len[0] + ' ' + str(min_len[1]))

        # Get average branch length relative to ancestral/root node in the
        # clade for each sequence/terminal/leaf node.
        mean_len = statistics.mean([x[1] for x in seq_branch_lengths])
        #print('\n\tMean length: ' + str(mean_len))

        # Get median branch length relative to ancestral/root node in the
        # clade for each sequence/terminal/leaf node.
        median_len = statistics.median([x[1] for x in seq_branch_lengths])
        #print('\n\tMedian length: ' + str(median_len))

        # Calculate ratio of stem length to average branch length.
        stem_branches_ratio = stem_len / mean_len
        #print('\n\tRatio of stem len to avg branch len: ' + str(stem_branches_ratio))


        # Format tree for visualization.

        # Root on first clade of interest.
        if num == 1:
            #print('\n\n\n\n')
            #print('t2:')
            #print(t2)
            #print(t2.write())
            #print('node:')
            #print(node)
            #print('\n\n\n\n')
            t2.set_outgroup(get_corresponding_node(node,t2))

        # Make stem branch for clade of interest bold.
        stem_node = get_corresponding_node(node, t2)
        stem_node.set_style(style_dict['stem_node'])

        # Add clade name to stem branch.
        #stem_node.name = clade_name
        stem_node.add_face(define_textface_for_labeling_stem(clade_name), column=0, position = "branch-top")

        # Iterate over branches within clade and customize appearance.
        for leaf in node.iter_leaves():
            # Set general features for all leaf nodes.
            style = NodeStyle()
            style["fgcolor"] = "#0f0f0f"
            style["size"] = 0
            style["vt_line_color"] = "#000000"
            style["hz_line_color"] = "#000000"

            #style["vt_line_width"] = 8
            #style["hz_line_width"] = 8
            #style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
            #style["hz_line_type"] = 0

            # Make shortest branch in clade of interest bold.
            if leaf.name == min_len[0]:
                #shortbranchnode = get_corresponding_node(t1.search_nodes(name=min_len[0])[0], t2)
                #shortbranchnode.set_style(style_dict['shortest branch'])
                style["hz_line_width"] = 5

            # Make longest branch in clade of interest bold.
            if leaf.name == max_len[0]:
                #longbranchnode = get_corresponding_node(t1.search_nodes(name=max_len[0])[0], t2)
                #longbranchnode.set_style(style_dict['longest branch'])
                style["hz_line_width"] = 5
                style["hz_line_type"] = 2
            
            # Customize leaf node appearance.
            # Get species name.
            species_name = leaf.name.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
            #print(species_name)

            # Check that the name does not end with '.copy'.
            assert not species_name.endswith('.copy'), """A sequence with an
            identical name to a sequence already in the alignment/tree was
            apparently added: %s""" % leaf.name

            ti = get_taxonomic_info(species_name)
            #print('\t' + ' '.join([str(ti.superbranch), str(ti.supergroup),
            #    str(ti.group), str(ti.species)]))

            # Set style for node.
            get_corresponding_node(leaf, t2).set_style(style)

        # Sort seq names by descending branch lengths.
        seq_names_by_br_len = [n.name for n in sorted(node.get_leaves(),\
                key=lambda n: node.get_distance(n), reverse=True)]

        # Identify most "basal" sequence in the clade, and the length of the
        # stem for the remaining seqs.
        # The node that is sister to the most basal sequence would be the one
        # that has the most daughter nodes.
        basal_node = None
        node_sister_to_basal_node = None
        for n in node.iter_descendants():
            if node_sister_to_basal_node is None:
                node_sister_to_basal_node = n
            if len(n.get_leaves()) > len(node_sister_to_basal_node.get_leaves()):
                node_sister_to_basal_node = n
        sister_leaf_names = [x.name for x in node_sister_to_basal_node.get_leaves()]
        for n in node:
            in_sister = True
            for x in n.iter_leaves():
                if not x.name in sister_leaf_names: 
                    in_sister = False
                    break
            if not in_sister:
                basal_node = n
                break
        # Check whether the basal node was identified.
        assert basal_node is not None
        # Define name of basal node.
        basal_node_name = 'clade'
        if basal_node.is_leaf():
            basal_node_name = basal_node.name
        # Define stem length of the clade that is sister to the "basal" clade.
        basal_node_depth = node.get_distance(node_sister_to_basal_node)


        # Add info to output dictionary.
        branch_length_info_dict[clade_name] =\
                {'stem length': stem_len,
                 'average branch length': mean_len,
                 'stem/branch ratio': stem_branches_ratio,
                 'seq names by descending length': seq_names_by_br_len,
                 'basal node name': basal_node_name,
                 'basal node depth': basal_node_depth,
                 'internal branch lengths': internal_branch_lengths
                 }




    # Show tree annotations.

    # Remove underscores from leaf names.
    for leaf in t2.iter_leaves():
        leaf.name = leaf.name.replace('__', ' ').replace('_', ' ')

    # Stretch branches.
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale =  120 # 120 pixels per branch length unit

    # Add title.
    tree_title = "[Tree title here]"
    ts.title.add_face(TextFace(tree_title, fsize=20), column=0)

    # Show.
    #t2.show(tree_style=ts)

    # Write tree to file.
    if outpath is not None:
        t2.render(outpath, w=183, units="mm", tree_style=ts)

    # Return main output path.
    return branch_length_info_dict


def reduce_alignment(alignment_file, output_alignment_file, removal_name_list):
    """Take file path for an existing alignment and an output alignment, and
    a list of sequence names for removal. Write the reduced alignment to the
    output path in nexus format. Note: This function assumes that the names in
    the removal lit and the alignment are all either coded or uncoded.
    """
    # Remove extra mesquite lines if present.
    delete_extra_mesquite_lines(alignment_file)

    # Make a copy of the input alignment file with the listed sequences
    # removed.
    coded_seqs_to_not_remove_from_dataset = []
    with open(alignment_file) as infh, open(output_alignment_file, 'w') as o:
        alignment = AlignIO.read(infh, 'nexus')
        #print('number of seqs in original alignment')
        #print(len(alignment))
        alignment2 = MultipleSeqAlignment([])
        #print('number of seqs in alignment2:')
        #print(len(alignment2))
        #print('Removing sequences from alignment.')
        #ali_len = len(alignment)
        inum = -1 
        for i in alignment:
            inum += 1
            if i.id not in removal_name_list:
                alignment2.append(alignment[inum])

        alignment3 =\
        MultipleSeqAlignment(alignment2[:,:])
        #print('number of seqs in alignment3:')
        #print(len(alignment3))

        # Get list of sequences in new alignment.
        ids_in_alignment3 = []
        for i in alignment3:
            ids_in_alignment3.append(i.id)

        # Write reduced alignment.
        #print('Writing modified alignment file:\n\t%s' % o.name)
        AlignIO.write(alignment3, o, 'nexus')


def recursively_remove_seqs(recursion_num,
                            t1,
                            alignment,
                            subs_model,
                            type_seqs_dict,
                            seqs_attempted_to_remove,
                            essential_taxa,
                            branch_length_info_dict,
                            log_file_path,
                            optimize_internal_branches
                            ):
    """Take a tree, and try to remove sequences to improve the ratio of the
    stem lengths of clades to the branch lengths of the branches in the clades,
    and return the final tree when all the non-taxonomically-essential
    sequences have been removed or no sequence can be removed that will result
    in improving the branch length ratios.
    """
    # Define variable for whether any new tree could be produced.
    new_tree_produced = False

    # Define variable for whether any new tree is better or not.
    new_tree_better = False

    # Define variables to store info for next recursion, if necessary.
    new_tree = None
    new_alignment = None
    new_subs_model = subs_model
    new_type_seqs_dict = None

    # Define variables for writing to log file.
    action = None
    annot_tree_path = 'No output tree'

    # Determine average of stem length to branch length ratio over clades of
    # interest.
    #branch_length_info_dict = get_branch_length_info(t1, type_seqs_dict)
    average_over_clades = statistics.mean([branch_length_info_dict[x]['stem/branch ratio']\
            for x in branch_length_info_dict.keys()])

    # Determine what sequence to remove.

    # Determine the clade with the worst ratio.
    clades_by_ascending_ratio = sorted(branch_length_info_dict.keys(),\
            key=lambda x: branch_length_info_dict[x]['stem/branch ratio'])
    #clade_with_lowest_ratio = clades_by_ascending_ratio[0]

    for clade in clades_by_ascending_ratio:
        # Make a list of taxonomic terms currently represented in the clade.
        taxa_represented = []
        for seqname in\
        branch_length_info_dict[clade]['seq names by descending length']:
            # Get species name from sequence name.
            species_name = seqname.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
            # Get taxonomic information for species name from info csv file.
            ti = get_taxonomic_info(species_name)
            # Add taxon names to taxa represented.
            taxa_represented = taxa_represented\
                    + [ti.superbranch, ti.supergroup, ti.group, ti.species]

        # Determine which seqs in the clade are dispensable, considering their
        # taxonomic placement.
        for seqname in\
        branch_length_info_dict[clade]['seq names by descending length']:
            # Get species name from sequence name.
            species_name = seqname.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
            # Get taxonomic information for species name from info csv file.
            ti = get_taxonomic_info(species_name)
            ti_list = [ti.superbranch, ti.supergroup, ti.group, ti.species]
            # Determine whether any of the taxonomy terms are in the essential
            # list.
            seq_is_essential = False
            essential_terms = []
            for x in ti_list:
                if x in essential_taxa:
                    essential_terms.append(x)
            # Check whether any of the terms appear in the taxa_represented list
            # only once (if so then the sequence is essential).
            for term in essential_terms:
                if taxa_represented.count(term) == 1:
                    seq_is_essential = True
            # Determine whether a sequence can be removed, if not then check the
            # next clade.
            if not seq_is_essential:
                if seqname not in seqs_attempted_to_remove:
                    # Try removing the corresponding sequence from the
                    # tree, and see whether it improves the branch length
                    # ratios.
                    print('\tremoving sequence %s' % seqname)


                    # If the sequence is a type sequence, then you have to
                    # choose another sequence in the clade as a replacement
                    # type sequence.
                    new_type_seqs_dict = type_seqs_dict.copy()
                    is_type_seq = False
                    for cladex in new_type_seqs_dict.keys():
                        if seqname == new_type_seqs_dict[cladex]:
                            is_type_seq = True
                            assert clade == cladex
                    if is_type_seq:
                        # Identify an alternative type sequence for the same
                        # (current) clade.
                        for seqname2 in\
                        reversed(branch_length_info_dict[clade]\
                        ['seq names by descending length'].copy()):
                            if seqname2 != seqname:
                                # Change type seq for clade to a different sequence.
                                new_type_seqs_dict[clade] = seqname2


                    # Define new alignment path.
                    alignment2 = alignment.rsplit('_', 1)[0] + '_'\
                    + str(recursion_num) + '.nex'
                    # Remove seq from alignment.
                    reduce_alignment(alignment, alignment2, [seqname])
                    new_alignment = alignment2

                    print('\tremoved seq from alignment')

                    # Code alignment (get table file).

                    # Remove seq from a copy of the parsed tree object, and
                    # use that as a constraint tree.
                    # Copy tree.
                    t2 = t1.copy()
                    # Get a list of all node objects with the same name as the
                    # sequence of interest.
                    seq_nodes = t2.search_nodes(name=seqname)
                    # Check that only one such node was found.
                    assert len(seq_nodes) == 1
                    # Remove the node from the new tree object.
                    node_to_remove = seq_nodes[0]
                    t2_len_before = len(t2.get_leaves())
                    node_to_remove.delete()

                    # Check that a node was actually removed.
                    t2_len_after = len(t2.get_leaves())
                    assert t2_len_before == t2_len_after + 1

                    # Write modified tree to a file.
                    tfp = os.path.join(alignment2.rsplit('.', 1)[0] + '_newick.tre')
                    t2.write(outfile=tfp, format=9)

                    print('\tremoved seq from tree object')

                    # Run ML tree with IQtree, and get branch length info
                    # dictionary.
                    # Define output directory for files.
                    outputdir = os.path.dirname(alignment)
                    # Call function to run IQtree and return a parsed tree obj.
                    t3 = get_ml_tree_branch_lengths(tfp, alignment2, subs_model, outputdir)
                    print('\tran IQtree')
                    # Set t3 as the new tree.
                    new_tree = t3
                    # Define output path for annotated tree image file.
                    annotated_tree_outpath =\
                    os.path.join(alignment2.rsplit('.', 1)[0] + '_annotated_result.png')
                    annot_tree_path = annotated_tree_outpath
                    # Extract info about branch lengths from tree, and write
                    # annotated tree image file.
                    new_branch_length_info_dict = get_branch_length_info(t3,
                                              new_type_seqs_dict,
                                              annotated_tree_outpath)
                    # Calculate average ratio over all clades from branch
                    # length info. 
                    new_average_over_clades =\
                            statistics.mean([new_branch_length_info_dict[x]['stem/branch ratio']\
                            for x in new_branch_length_info_dict.keys()])

                    print('\taverage ratio over clades: %s' %\
                            str(new_average_over_clades))

                    print('\tparsed IQtree output.')

                    # Write a file with type sequences and clade names for this
                    # iteration.
                    new_type_seqs_file=\
                    os.path.join(alignment2.rsplit('.', 1)[0] + '_type_seqs.csv')
                    with open(new_type_seqs_file, 'w') as o:
                        for key in type_seqs_dict.keys():
                            o.write(type_seqs_dict[key] + ',' + key + '\n')

                    # Compare lengths of internal branches.
                    previous_internal_branch_lengths =\
                    branch_length_info_dict[list(branch_length_info_dict.keys())[0]]['internal branch lengths']
                    new_internal_branch_lengths =\
                    new_branch_length_info_dict[list(branch_length_info_dict.keys())[0]]['internal branch lengths']
                    no_new_internal_branches_are_worse = True
                    at_least_one_internal_branch_longer = False
                    shortest_internal_branch_longer = False
                    shortest_internal_branch_shorter = False
                    for i, j in zip(previous_internal_branch_lengths, new_internal_branch_lengths):
                        if j < i:
                            no_new_internal_branches_are_worse = False
                        if j > i:
                            at_least_one_internal_branch_longer = True
                        if i == min(previous_internal_branch_lengths):
                            if j > i:
                                shortest_internal_branch_longer = True
                            if j < i:
                                shortest_internal_branch_shorter = True
                    # Get sums of previous and new internal branch lengths.
                    previous_sum_of_internal_branch_lengths =\
                    sum(previous_internal_branch_lengths)
                    new_sum_of_internal_branch_lengths =\
                    sum(new_internal_branch_lengths)

                    # Indicate that a new tree was produced.
                    new_tree_produced = True

                    # Compare the branch length ratios of the new tree to
                    # the current tree, and apply criteria to decide whether
                    # the new tree is better or not.
                    if not optimize_internal_branches:
                        if new_average_over_clades > average_over_clades:
                            # Make sure that it does not make the branch length ratio
                            # worse for any clade.
                            #made_a_clade_worse = False
                            #for cladex in branch_length_info_dict.keys():
                            #    old_clade_ratio = branch_length_info_dict\
                            #                      [cladex]['stem/branch ratio']
                            #    new_clade_ratio = new_branch_length_info_dict\
                            #                      [cladex]['stem/branch ratio']
                            #    if new_clade_ratio < old_clade_ratio:
                            #        made_a_clade_worse = True
                            #        break
                            #if not made_a_clade_worse:
                            #    print('\tnew tree is better')
                            #    new_tree_better = True

                            # Make sure that it did not make the ratio worse for
                            # the clade of interest (current worst clade).
                            old_clade_ratio = branch_length_info_dict\
                                              [clade]['stem/branch ratio']
                            new_clade_ratio = new_branch_length_info_dict\
                                              [clade]['stem/branch ratio']
                            if new_clade_ratio >= old_clade_ratio:
                                print('\tnew tree is better')
                                new_tree_better = True
                            
                            else:
                                print('\tnew tree is not better')

                    # Apply different criteria if the
                    # optimize_internal_branches option is selected.
                    elif optimize_internal_branches:
                        #if no_new_internal_branches_are_worse:
                        #    if at_least_one_internal_branch_longer:
                        #        print('\tnew tree is better')
                        #        # Temp.
                        #        assert 2 != 2
                        #        new_tree_better = True
                        if new_sum_of_internal_branch_lengths >\
                        previous_sum_of_internal_branch_lengths:
                            if not shortest_internal_branch_shorter:
                                print('\tnew tree is better')
                                new_tree_better = True

                        #if shortest_internal_branch_longer:
                        #    assert 2 != 2
                        #if shortest_internal_branch_shorter:
                        #    assert 2 != 2



                    # Add sequence to list of sequences which have been considered for
                    # removal.
                    seqs_attempted_to_remove.append(seqname)

                else:
                    # See if the next longest branch can be removed.
                    pass
            
            # Break the loop so that the function can be run again.
            if new_tree_produced:
                break
        # Break the loop so that the function can be run again.
        if new_tree_produced:
            break

    # If no new tree could be produced, then return the final tree.
    if not new_tree_produced:
        with open(log_file_path, 'a') as log:
            log.write(','.join([str(recursion_num), 'Removed a sequence',\
                os.path.basename(annot_tree_path)]) + '\n')

        return [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
                seqs_attempted_to_remove, branch_length_info_dict]

    # Otherwise run the function again.
    else:
        if new_tree_better:
            with open(log_file_path, 'a') as log:
                log.write(','.join([str(recursion_num), 'Removed a sequence',\
                    os.path.basename(annot_tree_path)]) + '\n')
            # Run again on new tree, alignment ,and type sequences.
            return recursively_remove_seqs(recursion_num + 1,
                                           new_tree,
                                           new_alignment,
                                           new_subs_model,
                                           new_type_seqs_dict,
                                           seqs_attempted_to_remove,
                                           essential_taxa,
                                           new_branch_length_info_dict,
                                           log_file_path,
                                           True
                                           )
        else:
            with open(log_file_path, 'a') as log:
                log.write(','.join([str(recursion_num), 'Removed a sequence',\
                    os.path.basename(annot_tree_path)]) + '\n')
            # Run again on same input (except with longer list of seqs that
            # were considered for removal).
            return recursively_remove_seqs(recursion_num + 1,
                                           t1,
                                           alignment,
                                           subs_model,
                                           type_seqs_dict,
                                           seqs_attempted_to_remove,
                                           essential_taxa,
                                           branch_length_info_dict,
                                           log_file_path,
                                           True
                                           )


def recursively_add_seqs(recursion_num,
                         t1,
                         alignment,
                         subs_model,
                         type_seqs_dict,
                         seqs_attempted_to_add,
                         essential_taxa,
                         branch_length_info_dict,
                         all_sequences,
                         log_file_path,
                         optimize_internal_branches
                         ):
    """Take a tree, and try to add sequences to improve the ratio of the
    stem lengths of clades to the branch lengths of the branches in the clades,
    and return the final tree when all the sequences have been added to the
    dataset (at least temporarily).
    """
    # Define variable for whether any new tree could be produced.
    new_tree_produced = False

    # Define variable for whether any new tree is better or not.
    new_tree_better = False

    # Define variables to store info for next recursion, if necessary.
    new_tree = None
    new_alignment = None
    new_subs_model = subs_model

    # Define variables for writing to log file.
    action = None
    annot_tree_path = 'No output tree'

    # Determine average of stem length to branch length ratio over clades of
    # interest.
    #branch_length_info_dict = get_branch_length_info(t1, type_seqs_dict)
    average_over_clades = statistics.mean([branch_length_info_dict[x]['stem/branch ratio']\
            for x in branch_length_info_dict.keys()])


    # Determine the clade with the worst ratio.
    clades_by_ascending_ratio = sorted(branch_length_info_dict.keys(),\
            key=lambda x: branch_length_info_dict[x]['stem/branch ratio'])
    clade_with_lowest_ratio = clades_by_ascending_ratio[0]

    # Get a list of all sequence names for tree.
    seq_names_in_tree = [x.name for x in t1.get_leaves()]

    dblundrscr = re.compile(r'__')
    uniform_spaces = lambda x: dblundrscr.sub('_', x.strip().replace(' ',\
        '_').replace(':', '_'))

    # Randomly pick a sequence that has not already been added, and does not
    # have the same name as any of the sequences already in the tree.
    seq_to_add = None
    shuffle(all_sequences)
    all_potential_seqs_to_add = seqs_attempted_to_add + seq_names_in_tree
    for seq in all_sequences:
        if uniform_spaces(seq.description) not in [uniform_spaces(x) for x in\
                all_potential_seqs_to_add]:
            seq_to_add = seq
            break

    if seq_to_add is not None:
        # Add to list of sequences added to to the tree.
        seqs_attempted_to_add.append(seq_to_add.description)

        # Indicate that a new tree was produced.
        new_tree_produced = True

        # Define logfile path to write to.
        logfile = log_file_path

        #ali = os.path.join(outdirpath1, '1.nex') 
        ali_num = recursion_num
        outdirpath1 = os.path.dirname(alignment)
        ali = alignment.rsplit('_', 1)[0] + '_' + str(recursion_num) + '.nex'
        #add_seq_to_alignment2(seq_to_add, alignment, ali)
        add_seq_to_alignment3(seq_to_add, alignment, ali)
        new_alignment = ali

        # Trim the alignment (any positions where only the new sequence has
        # residues). ????

        # Write input tree to a file.
        tree = ali.rsplit('.', 1)[0] + '_original_newick.tre'
        t1.write(outfile=tree, format=1)

        # Convert nexus alignment to afa.
        temp_fa_1 = ali.rsplit('.', 1)[0] + '_temp1.afa'
        nex_to_afa(ali, temp_fa_1) 

        # Code names in alignment.
        outalifp = ali.rsplit('.', 1)[0] + '.C.afa'
        outtablefp = ali.rsplit('.', 1)[0] + '.C.table'
        code_names_in_ali(outdirpath1, temp_fa_1, outalifp, outtablefp) 

        # Convert alignment back to nexus.
        outalifpnex = ali.rsplit('.', 1)[0] + '.C.nex'
        afa_to_nex(outalifp, outalifpnex)

        # Delete temporary files.
        os.remove(temp_fa_1)
        os.remove(outalifp)

        # Make constraint tree.
        # Put quotation marks around names in constraint tree.
        #outtreefp = ali.rsplit('.', 1)[0] + '.constraint_tree0.Q.newick'
        constraint_tree_fp = os.path.join(os.path.dirname(ali),\
                os.path.basename(tree.rsplit('.')[0] + '_constraint_tree.newick'))
        with open(tree) as intreefh, open(constraint_tree_fp, 'w') as o:
            for t in intreefh:
                o.write(quote_tree(t, outtablefp))
        # Make constraint tree with coded names.
        constraint_tree_fp_coded = constraint_tree_fp.rsplit('.', 1)[0] + '.C.tre'
        code_tree(constraint_tree_fp, constraint_tree_fp_coded, outtablefp)


        # Run an ML search with the alignment (with additional sequence) and
        # constraint tree defined above.

        # Make a copy of the alignment in phylip format for input to IQtree.
        phy_out = outalifpnex.rsplit('.', 1)[0] + '.phy'
        nex_to_phylip(outalifpnex, phy_out)

        # Remove nexus format version of alignment, because it is redundant with
        # the phylip version.
        os.remove(outalifpnex)

        # Do phylogenetic analysis. # SEPARATE
        # Make subdir for output.
        subdirp = constraint_tree_fp.rsplit('.', 1)[0] + '_IQ-tree_output'
        os.mkdir(subdirp)

        # ***Use the get_ml_tree_branch_lengths function here to run IQtree
        # instead???

        # Use IQtree to do an ML search
        output_file_prefix = os.path.join(subdirp, 'iqtree')
        stdout_path = output_file_prefix + '.stdout.txt'
        iqtree_command_list = ['iqtree', '-s', phy_out, '-m', subs_model, '-g',
            constraint_tree_fp_coded, '-pre', output_file_prefix, '-nt', '4']
        with open(stdout_path, 'w') as o:
            subprocess.call(iqtree_command_list, stdout=o, stderr=subprocess.STDOUT)

        # Define path to output file with tree.
        tree_file_path = output_file_prefix + '.treefile'
        # Check that the output file was actually produced.
        assert os.path.isfile(tree_file_path)

        # Parse the resulting tree to determine in which clade of interest, if any,
        # the sequence was placed.

        # Parse tree using ete3.
        # Note: parentheses and commas get replaced with underscores.
        t2 = Tree(tree_file_path, quoted_node_names=True)

        # Convert names in tree back to original names.
        uncode_tree_obj(t2, outtablefp)

        # Define new tree object for further recursions of this function.
        new_tree = t2

        ## Define name for tree without branch lengths.
        #simple_tree = tree_file_path.rsplit('_', 1)[0] + '_TEMP2' 

        ## Write simple tree to a new file 
        #t1.write(format=9, outfile=simple_tree)

        ## Parse simple tree.
        #t2 = Tree(simple_tree)

        # Print simple tree.
        #print('ML tree:')
        #print(t2)

        # Make a copy of the tree object.
        t3 = t2.copy()

        # Get list of "type" sequences from input.
        type_seq_list = list(type_seqs_dict.values())

        # Determine whether the additional sequence was placed in one of the
        # clades of interest, or not.

        # For each "type" sequence, traverse all nodes and find the node with
        # the largest number of child nodes that are leaf (terminal) nodes,
        # containing the "type" sequence of interest, but not containing any of
        # the other "type" sequences.
        ts_num = 0
        first_type_seq_node_name = None
        ts_that_additional_seq_was_placed_in = None
        for ts in type_seq_list:
            ts_num += 1

            if ts_num == 1:
                first_type_seq_node_name = ts
                # Root on another "type" sequence for the first type sequence in
                # the list to get whole clade, then root the tree on the ancestor
                # node of that first clade.

                # Get a node name for a node corresponding to a different
                # "type" sequence.
                other_type_seq_node_name = None
                for i in type_seq_list:
                    if i != ts:
                        other_type_seq_node_name = i
                        break

                # Get node corresponding to a different "type" sequence.
                other_type_seq_node = None
                for node in t3.traverse():
                    if node.name == other_type_seq_node_name:
                        other_type_seq_node = node
                        break

                # Root on the other "type" sequence node.
                t3.set_outgroup(other_type_seq_node)

                #print('\n\n\nTree rooted on a type sequence other than the first type sequence.')
                #print(t3)

            elif ts_num == 2:
                # Root on the first "type" sequence node for all subsequent
                # clades.
                first_type_seq_node = None
                for node in t3.traverse():
                    leaf_list = []
                    for leaf in node.get_leaves():
                        if leaf.name == first_type_seq_node_name:
                            first_type_seq_node = node
                            break
                t3.set_outgroup(first_type_seq_node)
                #print('\n\n\nTree re-rooted on first type sequence:')
                #print(t3)


            # Make a copy of the tree topology to work with for each run
            # through this loop.
            t4 = t3.copy()

            # Make a list of nodes that contain type seq, but not any others. #
            # SEPARATE
            nodes_of_interest = []
            for node in t4.traverse():
                # Search in nodes that contain the type sequence.
                if node.search_nodes(name=ts):
                    # Search in nodes that don't contain other type sequences.
                    contains_other_type_seqs = False
                    for ts2 in type_seq_list:
                        if not ts2 == ts:
                            if node.search_nodes(name=ts2):
                                contains_other_type_seqs = True
                    if not contains_other_type_seqs:
                        # Add nodes of interest to list.
                        nodes_of_interest.append(node)

            # find the node with the most child leaf nodes.
            #node_num = 0
            #for node in nodes_of_interest:
            #    node_num += 1 
            #    print('Node ' + str(node_num) + ' Number of leaves:  ' + str(len(node.get_leaves())))
            #    print(node)
            node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
                    len(x.get_leaves()), reverse=True)[0]
            #node_w_most_leaves.name = 'X'
            #print('\n\nClade of interest defined by sequence ' + ts + ':')
            #print(node_w_most_leaves)

            name_of_additional_seq = seq_to_add.description # should be .id ?
            if name_of_additional_seq.strip() in [x.name for x in node_w_most_leaves.get_leaves()]:
                ts_that_additional_seq_was_placed_in = ts
                #print('\n\n\n%s is in this clade!!!!\n\n\n' % name_of_additional_seq)
                break

        # Check that the clade that the additional sequence was placed in was
        # identified.
        #assert ts_that_additional_seq_was_placed_in is not None, """Sequence was
        #not placed in any of the clades of interest: %s""" % record.description

        # Define output path for annotated tree image file.
        annotated_tree_outpath =\
        os.path.join(ali.rsplit('.', 1)[0] + '_annotated_result.png')
        annot_tree_path = annotated_tree_outpath
        # Extract info about branch lengths from tree, and write
        # annotated tree image file.
        new_branch_length_info_dict = get_branch_length_info(t2,
                                      type_seqs_dict,
                                      annotated_tree_outpath
                                      )

        # Write a file with type sequences and clade names for this
        # iteration.
        new_type_seqs_file=\
        os.path.join(ali.rsplit('.', 1)[0] + '_type_seqs.csv')
        with open(new_type_seqs_file, 'w') as o:
            for key in type_seqs_dict.keys():
                o.write(type_seqs_dict[key] + ',' + key + '\n')


        # If the sequence was placed in one of the clades of interest, then
        # proceed with collecting further information. If not, then assume it
        # is of no use, and move on to the next sequence.
        if ts_that_additional_seq_was_placed_in is not None:
            # Get clade name.
            clade_that_seq_was_placed_in = None
            for c in type_seqs_dict.keys():
                if type_seqs_dict[c] == ts_that_additional_seq_was_placed_in:
                    clade_that_seq_was_placed_in = c
            assert clade_that_seq_was_placed_in is not None, """Could not find
            corresponding clade for sequence %s""" % ts_that_additional_seq_was_placed_in

            print('\tSequence placed in clade: ' +\
                    clade_that_seq_was_placed_in)

            # Get branch length info for new tree.

            # Define output path for annotated tree image file.
            annotated_tree_outpath =\
            os.path.join(ali.rsplit('.', 1)[0] + '_annotated_result.png')
            annot_tree_path = annotated_tree_outpath
            # Extract info about branch lengths from tree, and write
            # annotated tree image file.
            new_branch_length_info_dict = get_branch_length_info(t2,
                                          type_seqs_dict,
                                          annotated_tree_outpath
                                          )
            # Calculate average ratio over all clades from branch
            # length info. 
            new_average_over_clades =\
                    statistics.mean([new_branch_length_info_dict[x]['stem/branch ratio']\
                    for x in new_branch_length_info_dict.keys()])

            print('\tNew average ratio over clades: ' +\
                    str(new_average_over_clades))

            # Compare lengths of internal branches.
            previous_internal_branch_lengths =\
            branch_length_info_dict[list(branch_length_info_dict.keys())[0]]['internal branch lengths']
            new_internal_branch_lengths =\
            new_branch_length_info_dict[list(branch_length_info_dict.keys())[0]]['internal branch lengths']
            no_new_internal_branches_are_worse = True
            at_least_one_internal_branch_longer = False
            shortest_internal_branch_longer = False
            shortest_internal_branch_shorter = False
            for i, j in zip(previous_internal_branch_lengths, new_internal_branch_lengths):
                if j < i:
                    no_new_internal_branches_are_worse = False
                if j > i:
                    at_least_one_internal_branch_longer = True
                if i == min(previous_internal_branch_lengths):
                    if j > i:
                        shortest_internal_branch_longer = True
                    if j < i:
                        shortest_internal_branch_shorter = True
            # Get sums of previous and new internal branch lengths.
            previous_sum_of_internal_branch_lengths =\
            sum(previous_internal_branch_lengths)
            new_sum_of_internal_branch_lengths =\
            sum(new_internal_branch_lengths)


            # Compare branch length ratios between the new tree and the old tree.
            # Including overall average, specifically for the worst clade, and
            # specifically for the clade that the new sequence was placed into.
            if not optimize_internal_branches:
                if new_average_over_clades >= average_over_clades:
                    # Compare for the clade that the sequence was added to.
                    old_clade_ratio = branch_length_info_dict\
                                      [clade_that_seq_was_placed_in]['stem/branch ratio']
                    new_clade_ratio = new_branch_length_info_dict\
                                      [clade_that_seq_was_placed_in]['stem/branch ratio']
                    #print('\tnew ratio for clade %s: ' + str(new_clade_ratio) % clade_that_seq_was_placed_in)
                    if new_clade_ratio >= old_clade_ratio:
                        #print('\tnew tree is better')
                        #new_tree_better = True

                        # Compare for the clade that was the worst clade in the
                        # original tree.
                        old_worst_clade_ratio = branch_length_info_dict\
                                          [clade_with_lowest_ratio]['stem/branch ratio']
                        new_worst_clade_ratio = new_branch_length_info_dict\
                                          [clade_with_lowest_ratio]['stem/branch ratio']
                        #print('\tnew ratio for worst clade (%s): ' +\
                        #        str(new_worst_clade_ratio) %\
                        #        clade_with_lowest_ratio)
                        if new_worst_clade_ratio >= old_worst_clade_ratio:
                            print('\tnew tree is better')
                            new_tree_better = True

            # Apply different criteria if optimizing for internal branches.
            elif optimize_internal_branches:
                if new_sum_of_internal_branch_lengths >\
                previous_sum_of_internal_branch_lengths:
                    if not shortest_internal_branch_shorter:
                        print('\tnew tree is better')
                        new_tree_better = True


            # Make a list of taxonomic terms currently represented in the clade.
            taxa_represented = []
            for seqname in\
            new_branch_length_info_dict[clade_that_seq_was_placed_in]['seq names by descending length']:
                # Get species name from sequence name.
                species_name = seqname.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
                # Get taxonomic information for species name from info csv file.
                ti = get_taxonomic_info(species_name)
                # Add taxon names to taxa represented.
                taxa_represented = taxa_represented\
                        + [ti.superbranch, ti.supergroup, ti.group, ti.species]

            # Determine which seqs in the clade are dispensable, considering their
            # taxonomic placement.
            essential_seq_names = []
            for seqname in\
            new_branch_length_info_dict[clade_that_seq_was_placed_in]['seq names by descending length']:
                # Get species name from sequence name.
                species_name = seqname.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
                # Get taxonomic information for species name from info csv file.
                ti = get_taxonomic_info(species_name)
                ti_list = [ti.superbranch, ti.supergroup, ti.group, ti.species]
                # Determine whether any of the taxonomy terms are in the essential
                # list.
                seq_is_essential = False
                essential_terms = []
                for x in ti_list:
                    if x in essential_taxa:
                        essential_terms.append(x)
                # Check whether any of the terms appear in the taxa_represented list
                # only once (if so then the sequence is essential).
                for term in essential_terms:
                    if taxa_represented.count(term) == 1:
                        seq_is_essential = True
                # Add to list of essential seq names.
                if seq_is_essential:
                    essential_seq_names.append(seqname)

            # If the added sequence is in the list of essential taxa, then
            # assume the tree is better.
            if name_of_additional_seq in essential_seq_names:
                new_tree_better = True


            if not new_tree_better:
                print('\tnew tree is not better')

        else:
            # The additional sequence was not placed in a clade of interest
            print('\tSequence was not placed in a clade of interest.')

    # If no new tree could be produced, then return the final tree.
    if not new_tree_produced:
        with open(log_file_path, 'a') as log:
            log.write(','.join([str(recursion_num), 'No sequences to add',\
                os.path.basename(annot_tree_path)]) + '\n')

        return [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
                seqs_attempted_to_add, branch_length_info_dict]

    # Otherwise run the function again.
    else:
        if new_tree_better:
            with open(log_file_path, 'a') as log:
                log.write(','.join([str(recursion_num), 'Added a sequence',\
                    os.path.basename(annot_tree_path)]) + '\n')
            # Run again on new tree and alignment.
            return recursively_add_seqs(recursion_num + 1,
                                           new_tree,
                                           new_alignment,
                                           new_subs_model,
                                           type_seqs_dict,
                                           seqs_attempted_to_add,
                                           essential_taxa,
                                           new_branch_length_info_dict,
                                           all_sequences,
                                           log_file_path,
                                           optimize_internal_branches
                                           )

        else:
            with open(log_file_path, 'a') as log:
                log.write(','.join([str(recursion_num), 'Added a sequence',\
                    os.path.basename(annot_tree_path)]) + '\n')
            # Run again on same input (except with longer list of seqs that
            # were considered for addition).
            return recursively_add_seqs(recursion_num + 1,
                                           t1,
                                           alignment,
                                           subs_model,
                                           type_seqs_dict,
                                           seqs_attempted_to_add,
                                           essential_taxa,
                                           branch_length_info_dict,
                                           all_sequences,
                                           log_file_path,
                                           optimize_internal_branches
                                           )


def optimize_sequence_selection2(file_with_seqs, model_name, essential_taxa,
        timestamp):
    """Add/remove homologous sequences from a dataset to optimize for long stem
    lengths and low average branch lengths for clades of interest.

    The essential_taxa argument is a list of names corresponding to taxonomic
    groups of eukaryotes (such as "Amorphea" or "Diaphoretickes") which appear
    in the genome info csv file in the genomes directory.
    """
    # Increase recursion limit for python. Potentially dangerous, but necessary for
    # large alignments, like for Giantin and lava lamp proteins!
    sys.setrecursionlimit(10000)

    # Parse data from model.
    # Get info about model that is relevant.
    model_info = ModelInfoFromCSV(model_name)
    alignment = model_info.alignment_file
    tree = model_info.tree_topology_file
    subs_model = model_info.subs_model
    type_seqs = model_info.type_seqs_file

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(tree, quoted_node_names=False)

    # Print tree.
    print('Input tree:')
    print(t1)

    # Define output directory path.
    outputdir = os.path.join(os.path.dirname(file_with_seqs),\
            'select_seqs_output_' + timestamp)
    # Make output directory.
    assert not os.path.isdir(outputdir)
    os.mkdir(outputdir)

    # Construct a dictionary with info from type_seqs_file.
    type_seqs_dict = {}
    with open(type_seqs) as infh:
        for i in infh:
            if not i.startswith('\n'):
                spliti = i.strip().split(',')
                seqname = spliti[0]
                cladename = spliti[1]
                type_seqs_dict[cladename] = seqname

    # Copy alignment file to new path.
    new_ali_path = os.path.join(outputdir, os.path.basename(alignment).rsplit('.', 1)[0] + '_0.nex')
    shutil.copyfile(alignment, new_ali_path)
    assert os.path.isfile(new_ali_path)
    alignment = new_ali_path 

    # Run ML search on model tree to get branch lengths.
    t2 = get_ml_tree_branch_lengths(tree, alignment, subs_model, outputdir)

    # Define output path for annotated tree image file.
    annotated_tree_outpath =\
    os.path.join(alignment.rsplit('.', 1)[0] + '_annotated_result.png')
    # Extract info about branch lengths from tree, and write
    # annotated tree image file.
    branch_length_info_dict = get_branch_length_info(t2,
                              type_seqs_dict,
                              annotated_tree_outpath)

    # Compile list of sequences represented in the input initial tree.
    names_of_seqs_in_tree = [x.name for x in t2.get_leaves()]

    # Parse input sequence file, and construct a list of sequence objects for
    # potential inclusion in the tree.
    seqs = []
    for seq in SeqIO.parse(file_with_seqs, 'fasta'):
        seqs.append(seq)

    # Define path to log file (specifically for select_seqs).
    log_file_path = alignment.rsplit('.', 1)[0] + '_select_seqs_log.txt'
    header_line = "Step number, Action, Output file\n"
    with open(log_file_path, 'w') as o:
        o.write(header_line)
        o.write(','.join(['0', 'Initial version',\
            os.path.basename(annotated_tree_outpath)]) + '\n')


    ########
    ## Load information into a recursive or iterative function/loop to remove
    ## sequences.
    seqs_attempted_to_remove = []
    best_tree_info = recursively_remove_seqs(1,
                                             t2,
                                             alignment,
                                             subs_model,
                                             type_seqs_dict,
                                             seqs_attempted_to_remove,
                                             essential_taxa,
                                             branch_length_info_dict,
                                             log_file_path,
                                             False
                                             )
    ##(returns: [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
    ##        seqs_attempted_to_remove, branch_length_info_dict] )

    ## Try to remove even more sequences.
    #seqs_attempted_to_remove = []
    best_tree_info = recursively_remove_seqs(best_tree_info[0],
                                              best_tree_info[1],
                                              best_tree_info[2],
                                              best_tree_info[3],
                                              best_tree_info[4],
                                              seqs_attempted_to_remove,
                                              essential_taxa,
                                              best_tree_info[6],
                                              log_file_path,
                                              True
                                              )
    #(returns: [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
    #        seqs_attempted_to_remove, branch_length_info_dict] )


    ## Now that all the non-essential sequences have been removed, try adding
    ## sequences.
    seqs_attempted_to_add = []
    best_tree_info1 = recursively_add_seqs(best_tree_info[0],
                                           best_tree_info[1],
                                           best_tree_info[2],
                                           best_tree_info[3],
                                           best_tree_info[4],
                                           seqs_attempted_to_add,
                                           essential_taxa,
                                           best_tree_info[6],
                                           seqs,
                                           log_file_path,
                                           False
                                           )
    seqs_attempted_to_add = []
    best_tree_info1 = recursively_add_seqs(best_tree_info1[0],
                                           best_tree_info1[1],
                                           best_tree_info1[2],
                                           best_tree_info1[3],
                                           best_tree_info1[4],
                                           seqs_attempted_to_add,
                                           essential_taxa,
                                           best_tree_info1[6],
                                           seqs,
                                           log_file_path,
                                           True
                                           )

    #seqs_attempted_to_add = []
    #best_tree_info1 = recursively_add_seqs(1,
    #                                       t2,
    #                                       alignment,
    #                                       subs_model,
    #                                       type_seqs_dict,
    #                                       seqs_attempted_to_add,
    #                                       essential_taxa,
    #                                       branch_length_info_dict,
    #                                       seqs,
    #                                       log_file_path
    #                                       )

    #return [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
    #        seqs_attempted_to_add, branch_length_info_dict]

    seqs_attempted_to_remove = []
    best_tree_info2 = recursively_remove_seqs(best_tree_info1[0],
                                             best_tree_info1[1],
                                             best_tree_info1[2],
                                             best_tree_info1[3],
                                             best_tree_info1[4],
                                             seqs_attempted_to_remove,
                                             essential_taxa,
                                             best_tree_info1[6],
                                             log_file_path,
                                             False
                                             )
    ##(returns: [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
    ##        seqs_attempted_to_remove, branch_length_info_dict] )

    best_tree_info3 = recursively_remove_seqs(best_tree_info2[0],
                                              best_tree_info2[1],
                                              best_tree_info2[2],
                                              best_tree_info2[3],
                                              best_tree_info2[4],
                                              seqs_attempted_to_remove,
                                              essential_taxa,
                                              best_tree_info2[6],
                                              log_file_path,
                                              True
                                              )

    return '[main output path here]'








######

    # Define a dictionary of ete3 NodeStyle objects for different branch types
    ## in the tree to be visualized.
    #style_dict = define_nodestyles_dict_for_colourcoding()

    ## Get info about model that is relevant.
    #model_info = ModelInfoFromCSV(model_name)
    #alignment = model_info.alignment_file
    #tree = model_info.tree_topology_file
    #subs_model = model_info.subs_model
    #type_seqs = model_info.type_seqs_file

    ## Parse tree using ete3.
    ## Note: parentheses and commas get replaced with underscores.
    #t1 = Tree(tree, quoted_node_names=False)

    ## Print tree.
    #print('Input tree:')
    #print(t1)

    ## Get a list of nodes of interest (ancestral nodes for clades of interest).
    #orthogroup_nodes = []

    ## Get list of "type" sequences from input.
    #type_seq_list = []
    #for i in open(type_seqs):
    #    type_seq_list.append(i.strip().split(',')[0])

    ## Copy tree for colour-coding etc.
    #t2 = t1.copy()

    ## Get list of nodes of interest.
    #orthogroup_nodes = get_nodes_of_interest(t1, type_seq_list)

    ## Determine characteristics of clades of interest.
    #num = 0
    #for seq_node_list in orthogroup_nodes:
    #    num += 1
    #    type_seq = seq_node_list[0]
    #    clade_name = get_clade_name_from_model(type_seq, type_seqs)
    #    node = seq_node_list[1]

    #    print('\n\nClade defined by sequence ' + type_seq + ':')
    #    print(node)

    #    # Get list of leaf names.
    #    leaf_names = get_list_of_leaf_names_for_node(node)

    #    # Get stem length for node.
    #    stem_len = node.dist
    #    print('\nStem length: ' + str(stem_len))

    #    # Get branch lengths for all sequences in the clade.
    #    print('\nBranch lengths:')
    #    seq_branch_lengths = []
    #    for n in node.traverse():
    #        if n.is_leaf():
    #            length = node.get_distance(n)
    #            print('\t' + n.name + ': ' + str(length))
    #            seq_branch_lengths.append([n.name, length])

    #    # Get longest branch length.
    #    max_len = max(seq_branch_lengths, key=lambda x: x[1])
    #    print('\n\tLongest branch: ' + max_len[0] + ' ' + str(max_len[1]))

    #    # Get shortest branch length.
    #    min_len = min(seq_branch_lengths, key=lambda x: x[1])
    #    print('\n\tShortest branch: ' + min_len[0] + ' ' + str(min_len[1]))

    #    # Get average branch length relative to ancestral/root node in the
    #    # clade for each sequence/terminal/leaf node.
    #    mean_len = statistics.mean([x[1] for x in seq_branch_lengths])
    #    print('\n\tMean length: ' + str(mean_len))

    #    # Get median branch length relative to ancestral/root node in the
    #    # clade for each sequence/terminal/leaf node.
    #    median_len = statistics.median([x[1] for x in seq_branch_lengths])
    #    print('\n\tMedian length: ' + str(median_len))

    #    # Calculate ratio of stem length to average branch length.
    #    stem_branches_ratio = stem_len / mean_len
    #    print('\n\tRatio of stem len to avg branch len: ' + str(stem_branches_ratio))


    #    # Format tree for visualization.

    #    # Root on first clade of interest.
    #    if num == 1:
    #        t2.set_outgroup(get_corresponding_node(node, t2))

    #    # Make stem branch for clade of interest bold.
    #    stem_node = get_corresponding_node(node, t2)
    #    stem_node.set_style(style_dict['stem_node'])

    #    # Add clade name to stem branch.
    #    #stem_node.name = clade_name
    #    stem_node.add_face(define_textface_for_labeling_stem(clade_name), column=0, position = "branch-top")

    #    # Iterate over branches within clade and customize appearance.
    #    for leaf in node.iter_leaves():
    #        # Set general features for all leaf nodes.
    #        style = NodeStyle()
    #        style["fgcolor"] = "#0f0f0f"
    #        style["size"] = 0
    #        style["vt_line_color"] = "#000000"
    #        style["hz_line_color"] = "#000000"

    #        #style["vt_line_width"] = 8
    #        #style["hz_line_width"] = 8
    #        #style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
    #        #style["hz_line_type"] = 0

    #        # Make shortest branch in clade of interest bold.
    #        if leaf.name == min_len[0]:
    #            #shortbranchnode = get_corresponding_node(t1.search_nodes(name=min_len[0])[0], t2)
    #            #shortbranchnode.set_style(style_dict['shortest branch'])
    #            style["hz_line_width"] = 5

    #        # Make longest branch in clade of interest bold.
    #        if leaf.name == max_len[0]:
    #            #longbranchnode = get_corresponding_node(t1.search_nodes(name=max_len[0])[0], t2)
    #            #longbranchnode.set_style(style_dict['longest branch'])
    #            style["hz_line_width"] = 5
    #            style["hz_line_type"] = 2
    #        
    #        # Customize leaf node appearance.
    #        # Get species name.
    #        species_name = leaf.name.replace('__', ' ').replace('_', ' ').split(' ', 1)[1]
    #        print(species_name)
    #        ti = get_taxonomic_info(species_name)
    #        print('\t' + ' '.join([str(ti.superbranch), str(ti.supergroup),
    #            str(ti.group), str(ti.species)]))

    #        # Set style for node.
    #        get_corresponding_node(leaf, t2).set_style(style)


    ## Show tree with colour-coded branches.

    ## Remove underscores from leaf names.
    #for leaf in t2.iter_leaves():
    #    leaf.name = leaf.name.replace('__', ' ').replace('_', ' ')

    ## Stretch branches.
    #ts = TreeStyle()
    #ts.show_leaf_name = True
    #ts.scale =  120 # 120 pixels per branch length unit

    ## Add title.
    #tree_title = "[Tree title here]"
    #ts.title.add_face(TextFace(tree_title, fsize=20), column=0)

    ## Show.
    #t2.show(tree_style=ts)

    ## Write tree to file.
    #t2.render("/Users/Lael/Desktop/mytree.png", w=183, units="mm", tree_style=ts)

    ## Return main output path.
    #return '[no output path defined]'
