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

"""Contains functions for modifying alignment position/column selection for a
given phylogenetic tree dataset.
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

from trim_nex import trim_nex
from paralogue_counter import get_seq_obj_from_srch_res_csv_info
import column_header_lists
from phylo_clas import ModelInfoFromCSV,\
get_clade_name_from_model, code_names_in_ali, quote_tree, code_tree,\
uncode_tree, uncode_tree_obj
from paralogue_counter import add_seq_to_alignment3

# Import functions for working with sequences, alignments, and trees.
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
from afa_to_nex import delete_extra_mesquite_lines, afa_to_nex, nex_to_afa, nex_to_phylip
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

from select_seqs import get_clade_name_from_model2,\
get_nodes_of_interest, get_list_of_leaf_names_for_node, TaxonomicInfo,\
get_taxonomic_info, define_nodestyles_dict_for_colourcoding,\
define_textface_for_labeling_stem, get_corresponding_node,\
get_ml_tree_branch_lengths, get_branch_length_info, reduce_alignment


def trim_one_column_from_alignment(alignment, alignment2, column_index):
    """Take a nexus alignment file path and remove column at a given index
    (first column has index 0).
    """
    with open(alignment) as infh, open(alignment2, 'w') as o:
        # Parse input alignment file.
        alignment_obj = AlignIO.read(infh, 'nexus')
        # Record input alignment length.
        len1 = len(alignment_obj[0].seq)
        # Remove relevant column from alignment object.
        alignment_obj2 = alignment_obj[:, :column_index] + alignment_obj[:, column_index + 1:]
        # Check that the alignment length is shorter by one column.
        len2 = len(alignment_obj2[0].seq)
        assert len1 - len2 == 1
        # Write modified alignment object to output file path.
        AlignIO.write(alignment_obj2, o, 'nexus')



def check_whether_all_cols_unique_seqs(alignment_path):
    """Take a path to an alignment in nexus format, and return True if all the
    columns in the alignment have unique sequences, and False otherwise.
    """
    # Parse input alignment file and compile a list of column sequences.
    columns = []
    with open(alignment_path) as infh:
        alignment = AlignIO.read(infh, 'nexus')

        # Get length of sequences in alignment.
        seq_len = alignment.get_alignment_length()

        # get a list of columns as strings in the original alignment.
        columns = [alignment[:, col] for col in range(seq_len)] 

    assert columns != []
        
    # Determine whether there are redundant column sequences.
    x = len(columns) - len(list(set(columns)))
    if x == 0:
        return True
    else:
        return False


def recursively_remove_positions(recursion_num,
                                 t1,
                                 alignment,
                                 subs_model,
                                 type_seqs_dict,
                                 cols_attempted_to_remove,
                                 branch_length_info_dict,
                                 log_file_path,
                                 optimize_internal_branches
                                 ):
    """Take a tree, and try to remove positions/columns from the alignment to
    improve the ratio of the stem lengths of clades to the branch lengths of
    the branches in the clades, and return the final tree when as many postions
    as possible have been removed or no sequence can be removed that will
    result in improving the branch length ratios.
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

    # Determine the clade with the worst ratio.
    clades_by_ascending_ratio = sorted(branch_length_info_dict.keys(),\
            key=lambda x: branch_length_info_dict[x]['stem/branch ratio'])
    clade_with_lowest_ratio = clades_by_ascending_ratio[0]

    # Identify the clade for which to check for improved branch length ratios
    # for (an alternative would be to iterate over all of them)***
    clade = clade_with_lowest_ratio


    ####################### Key part that is different from the function for
    # removing sequences instead of positions ########################

    # Parse input alignment, and find a column to try to remove.
    with open(alignment) as infh:
        alignment_obj = AlignIO.read(infh, 'nexus')

        # Get length of sequences in alignment.
        seq_len = alignment_obj.get_alignment_length()

        # get a list of columns as strings in the original alignment.
        columns = [alignment_obj[:, col] for col in range(seq_len)] 

        # Iterate over column sequences.
        column_index = -1
        for column in columns:
            column_index += 1
            if not column in cols_attempted_to_remove:
                # Try removing the column from the tree, and see whether it
                # improves the branch length ratios.

                ## If the sequence is a type sequence, then you have to
                ## choose another sequence in the clade as a replacement
                ## type sequence.
                #new_type_seqs_dict = type_seqs_dict.copy()
                #is_type_seq = False
                #for cladex in new_type_seqs_dict.keys():
                #    if seqname == new_type_seqs_dict[cladex]:
                #        is_type_seq = True
                #        assert clade == cladex
                #if is_type_seq:
                #    # Identify an alternative type sequence for the same
                #    # (current) clade.
                #    for seqname2 in\
                #    reversed(branch_length_info_dict[clade]\
                #    ['seq names by descending length'].copy()):
                #        if seqname2 != seqname:
                #            # Change type seq for clade to a different sequence.
                #            new_type_seqs_dict[clade] = seqname2

                new_type_seqs_dict = type_seqs_dict


                # Define new alignment path.
                alignment2 = alignment.rsplit('_', 1)[0] + '_'\
                + str(recursion_num) + '.nex'

                # Remove column from alignment.
                trim_one_column_from_alignment(alignment, alignment2,
                        column_index)
                new_alignment = alignment2

                print('\tremoved column from alignment')

                # Code alignment (get table file).

                ## Remove seq from a copy of the parsed tree object, and
                ## use that as a constraint tree.
                ## Copy tree.
                t2 = t1.copy()
                ## Get a list of all node objects with the same name as the
                ## sequence of interest.
                #seq_nodes = t2.search_nodes(name=seqname)
                ## Check that only one such node was found.
                #assert len(seq_nodes) == 1
                ## Remove the node from the new tree object.
                #node_to_remove = seq_nodes[0]
                #t2_len_before = len(t2.get_leaves())
                #node_to_remove.delete()

                ## Check that a node was actually removed.
                #t2_len_after = len(t2.get_leaves())
                #assert t2_len_before == t2_len_after + 1

                # Write modified tree to a file.
                tfp = os.path.join(alignment2.rsplit('.', 1)[0] + '_newick.tre')
                t2.write(outfile=tfp, format=9)

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


                # Should make a separate function for comparing current and
                # previous trees.

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



                # Add column sequence to list of column sequences which have
                # been considered for removal.
                cols_attempted_to_remove.append(column)

            else:
                # See if the next longest branch can be removed.
                pass
            
            # Break the loop so that the function can be run again.
            if new_tree_produced:
                break
        ## Break the loop so that the function can be run again.
        #if new_tree_produced:
        #    break

    # If no new tree could be produced, then return the final tree.
    if not new_tree_produced:
        with open(log_file_path, 'a') as log:
            log.write(','.join([str(recursion_num), 'Removed a column',\
                os.path.basename(annot_tree_path)]) + '\n')

        return [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
                cols_attempted_to_remove, branch_length_info_dict]

    # Otherwise run the function again.
    else:
        if new_tree_better:
            with open(log_file_path, 'a') as log:
                log.write(','.join([str(recursion_num), 'Removed a column',\
                    os.path.basename(annot_tree_path)]) + '\n')
            # Run again on new tree, alignment ,and type sequences.
            return recursively_remove_positions(recursion_num + 1,
                                           new_tree,
                                           new_alignment,
                                           new_subs_model,
                                           new_type_seqs_dict,
                                           cols_attempted_to_remove,
                                           new_branch_length_info_dict,
                                           log_file_path,
                                           True
                                           )
        else:
            with open(log_file_path, 'a') as log:
                log.write(','.join([str(recursion_num), 'Removed a column',\
                    os.path.basename(annot_tree_path)]) + '\n')
            # Run again on same input (except with longer list of seqs that
            # were considered for removal).
            return recursively_remove_positions(recursion_num + 1,
                                           t1,
                                           alignment,
                                           subs_model,
                                           type_seqs_dict,
                                           cols_attempted_to_remove,
                                           branch_length_info_dict,
                                           log_file_path,
                                           True
                                           )


def optimize_position_selection(model_name, outdirpath, timestamp):
    """Systematically remove position/columns from an alignment to optimize for
    stem to median branch length ratios for clades of interest as a measure of
    resolution.
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
    outputdir = os.path.join(outdirpath,\
            model_name + '_select_positions_output_' + timestamp)
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

    ## Compile list of sequences represented in the input initial tree.
    #names_of_seqs_in_tree = [x.name for x in t2.get_leaves()]

    ## Parse input sequence file, and construct a list of sequence objects for
    ## potential inclusion in the tree.
    #seqs = []
    #for seq in SeqIO.parse(file_with_seqs, 'fasta'):
    #    seqs.append(seq)

    # Define path to log file (specifically for select_positions).
    log_file_path = alignment.rsplit('.', 1)[0] + '_select_positions_log.txt'
    header_line = "Step number, Action, Output file\n"
    with open(log_file_path, 'w') as o:
        o.write(header_line)
        o.write(','.join(['0', 'Initial version',\
            os.path.basename(annotated_tree_outpath)]) + '\n')

    # Check that all the columns in the input alignment have unique sequences
    # (so that they can be identified according to their sequence rather than
    # their order in the alignment).
    all_columns_have_unique_seqs =\
    check_whether_all_cols_unique_seqs(alignment)
    assert all_columns_have_unique_seqs, """Some columns in the input alignment
    have identical sequences."""

    ########
    ## Load information into a recursive or iterative function/loop to remove
    ## sequences.
    cols_attempted_to_remove = []
    best_tree_info = recursively_remove_positions(1,
                                                  t2,
                                                  alignment,
                                                  subs_model,
                                                  type_seqs_dict,
                                                  cols_attempted_to_remove,
                                                  branch_length_info_dict,
                                                  log_file_path,
                                                  False
                                                  )
    ##(returns: [recursion_num, t1, alignment, subs_model, type_seqs_dict,\
    ##        seqs_attempted_to_remove, branch_length_info_dict] )

    assert 2 != 2, """Remember to make a mask of the input alignment based on
    the columns that were selected for removal."""
    # ...

    return '[main output path here]'



