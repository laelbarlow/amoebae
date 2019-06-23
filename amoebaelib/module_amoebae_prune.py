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
"""Module for amoebae script with functions for parsing trees using the ETE3
library.
"""
# Import built-in modules.
#import argparse
#import sys
import os
import subprocess
import re
import statistics
import numpy as np
#import settings
#import shutil
import glob
#import time
#import pandas as pd

# Import modules from installed libraries/packages.
#from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped

#from get_datatype import get_dbtype
#from module_paralogue_counter import get_seq_obj_from_srch_res_csv_info

from module_afa_to_nex import delete_extra_mesquite_lines
from module_amoebae_name_replace import write_afa_with_code_names,\
    write_newick_tree_with_uncoded_names, get_conversion_dict_from_table,\
    write_newick_tree_with_coded_names

# Import functions for working with sequences, alignments, and trees.
from ete3 import Tree, TreeStyle, Tree, TextFace, add_face_to_node
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment

from visualize_trees import translate_int_node_names_to_support, visualize_tree, find_input_file_in_parent_directory

# Define functions to be used in amoebae.


def reformat_combined_supports2(tree_string, one_support='bootstrap'):
    """For converting IQtree output (using both -bb and -alrt) to a format that
    is readable by ETE3.
    
    Currently have to choose one support value to show either 'boostrap' or
    'alrt'.
    """
    #print('Parsing node support values:')
    #cs = re.compile(r'\)\d+(\.\d)?/\d+:')
    cs = re.compile(r'\)\d+.?\d?/\d+:')
    for instance in cs.findall(tree_string):
        #print(instance)
        # Define a reformatted support string.
        supcomb = instance[1:-1]
        #print('supcomb:')
        #print(supcomb)

        alrt = supcomb.split('/')[0]
        #print('alrt:')
        #print(alrt)
        alrt = str(round(float(alrt)))
        alrt = '0'*(3-len(alrt)) + alrt

        boot = supcomb.split('/')[1]
        #print('boot:')
        #print(boot)
        boot = str(round(float(boot)))
        boot = '0'*(3-len(boot)) + boot


        #supcomb2 = boot + '' + alrt

        assert one_support == 'bootstrap' or one_support == 'alrt', """Support
        type must be 'bootstrap' or 'alrt'."""

        if one_support == 'bootstrap':
            supcomb2 = boot
        elif one_support == 'alrt':
            supcomb2 = alrt

        # Replace the instance with a reformatted support string in the
        # tree_string.
        tree_string = tree_string.replace(instance, instance[0] +\
                supcomb2 + instance[-1])

    # Return modified tree string.
    return tree_string


def make_quoted_table(intablepath, outtablepath):
    """Takes a path to a regular table, and writes a file to a given file path
    with quotation marks around the long names.
    """
    with open(intablepath) as infh, open(outtablepath, 'w') as o:
        for i in infh:
            if not i.startswith('ZZ') and not i.startswith('\n'):
                o.write('\"' + i.strip().replace('(', '_').replace(')',
                    '_').replace(',', '_').replace(':', '_') + '\"\n')
            else:
                o.write(i)


def get_conversion_table_dict(table_with_quotation_marks):
    """Take a conversion table, and return a dictionary with coded names as
    values and original names as keys.
    """
    conversion_table_dict = {}
    with open(table_with_quotation_marks) as infh:
        k = []
        v = []
        for i in infh:
            if not i.startswith('\n'):
                if i.startswith('ZZ'):
                    value = i.strip()
                    v.append(value)
                elif not i.startswith('ZZ'):
                    key = i.strip()
                    k.append(key)

        assert len(k) == len(v), "Error: Could not parse table file."

        for i in range(len(k)):
            conversion_table_dict[k[i].strip('\"')] = v[i]

    return conversion_table_dict


def get_tree_string_with_readable_supports(tree_with_long_names):
    """Change node support values to a format that is readable by ETE3.
    Overwrites input file as output file.
    """
    treestring = None
    with open(tree_with_long_names) as infh:
        num = 0
        for i in infh:
            num += 1
            treestring = i.strip()
            if num == 1:
                break
    if treestring is not None:
        with open(tree_with_long_names, 'w') as o:
            o.write(reformat_combined_supports2(treestring))


def root_tree_manually(t, tx):
    """Re-root a given ete3 Tree manually.
    t is the tree object to modify.
    tx is the tree object to view in terminal.
    """
    # Print simple tree.
    print('Input tree:')
    print(tx)
    # Get node for a key clade. 
    # Or, just require to root it before inputting (figtree can export
    # re-rooted newick)?????
    print("\nPlease define an arbitrary clade of interest for rooting the tree topology\n")
    clade_member_1 = input('Clade member 1: ')
    clade_member_2 = input('Clade member 2: ')
    #ancestor1 = t2.get_common_ancestor(clade_member_1, clade_member_2)
    ancestor1 = t.get_common_ancestor(clade_member_1, clade_member_2)
    #print(ancestor1)

    # Re-root tree. 
    t.set_outgroup(ancestor1)


def get_prune_tree_style():
    # Set up tree style for viewing.
    ts = TreeStyle()
    ts.show_leaf_name = False
    def my_layout(node):
        """Function for showing internal node names (numbers)."""
        #F = TextFace(node.name, tight_text=True)
        F = TextFace(node.name)
        add_face_to_node(F, node, column=0, position="branch-right")
    ts.scale =  100
    ts.branch_vertical_margin = 10
    ts.layout_fn = my_layout
    ts.show_branch_support = True
    # Return the tree style.
    return ts


def get_seqs_to_remove_file_from_user(to_remove_list_file, tree_file_name):
    """Prompt user to list nodes and sequence names to remove in a file.
    """
    # Show tree with node names (render as PDF file).
    #t1.show(tree_style=ts)
    #subprocess.call(['open', tree_file_name]) # uncomment later.

    # Prompt user input.
    print('''\nList node numbers and/or sequence accessions for sequences to
    remove in file:\n%s''' % to_remove_list_file)
    editor = os.environ.get('EDITOR','vim') 
    initial_message = '# List node numbers or sequence accessions for sequences to remove below.'
    with open(to_remove_list_file, 'w') as o:
        o.write(initial_message)
        o.flush()
        subprocess.call([editor, o.name])

    # Check that file was produced.
    assert os.path.isfile(to_remove_list_file), """Did not produce file listing
    sequences to remove from dataset."""


def leaves_all_same_species(leaves):
    """Take a list of TreeNodes (leaves) and return True if they represent
    sequences that all came from the same species.
    """
    if len(set([x.name.split('__')[0] for x in leaves])) == 1:
        return True
    else:
        return False


def get_seqs_to_remove_file(to_remove_list_file, t, remove_redun_seqs, min_support,
        max_bl_iqr_above_third_quartile):
    """Prompt user to list nodes and sequence names to remove in a file.
    """
    ## Show tree with node names (render as PDF file).
    #t1.show(tree_style=ts)
    #subprocess.call(['open', tree_file_name]) # uncomment later.

    ## Prompt user input.
    #print('''\nList node numbers and/or sequence accessions for sequences to
    #remove in file:\n%s''' % to_remove_list_file)
    #editor = os.environ.get('EDITOR','vim') 
    #initial_message = '# List node numbers or sequence accessions for sequences to remove below.'
    #with open(to_remove_list_file, 'w') as o:
    #    o.write(initial_message)
    #    o.flush()
    #    subprocess.call([editor, o.name])



    # Initiate a list of leaf names for removal.
    names_of_seqs_to_remove = []

    # Get furthest leaf node on either side of the root.
    far_leaves = []
    leaves_in_midpoint_outgroup = t.get_midpoint_outgroup().get_leaves()
    names_of_leaves_in_midpoint_outgroup = [x.name for x in leaves_in_midpoint_outgroup]
    leaves_on_other_side = []
    for y in t.get_leaves():
        if y.name not in names_of_leaves_in_midpoint_outgroup:
            leaves_on_other_side.append(y)
    furthest_leaf_in_midpoint_outgroup = sorted(leaves_in_midpoint_outgroup, key=lambda x: t.get_distance(x))[-1]
    far_leaves.append(furthest_leaf_in_midpoint_outgroup)
    furthest_leaf_on_other_side = sorted(leaves_on_other_side, key=lambda x: t.get_distance(x))[-1]
    far_leaves.append(furthest_leaf_on_other_side)
    
    # For each of the leaf nodes, re-root the tree on it, and find redundant
    # and long branches (leaves) to add to list.
    for far_leaf in far_leaves + [t.get_midpoint_outgroup()]:
        # Re-root tree.
        t.set_outgroup(far_leaf)

        #print('\n\nRe-rooted tree:')
        #print(t)

        # Calculate interquartile range of all branch lengths.
        all_leaf_distances = []
        for leaf in t.iter_leaves():
            leaf_distance = t.get_distance(leaf)
            all_leaf_distances.append(leaf_distance)
        quartile3, quartile1 = np.percentile(all_leaf_distances, [75 ,25])
        interquartile_range = quartile3 - quartile1


        # Determine the median branch length.
        distances_for_leaves = []
        for leaf in t.iter_leaves():
            distances_for_leaves.append(t.get_distance(leaf))
        median_branch_length = statistics.median(distances_for_leaves) 

        # Loop over all nodes.
        for node in t.traverse():
            if remove_redun_seqs:
                # Identify all nodes with leaves that are all from the same
                # species, and remove all but the shortest branch in the clade.
                leaves = node.get_leaves()
                if leaves_all_same_species(leaves):
                    # Determine whether support is high enough.
                    # Note: Support values are as modified by the
                    # reformat_combined_supports2 function.
                    if node.support >= float(min_support): 
                        # Identify the longest of the branches.
                        all_but_shortest_leaf = sorted(leaves, key=lambda x: node.get_distance(x))[1:]
                        # List all but the shortest branch for removal.
                        for l in all_but_shortest_leaf:
                            if l.name not in names_of_seqs_to_remove:
                                names_of_seqs_to_remove.append(l.name)

            # Identify leaves with branch lengths (measured from the root node
            # to the leaf) that exceed the maximum number of interquartile
            # ranges above the third quartile of all such lengths in the tree.
            if node.is_leaf():
                leaf_distance = t.get_distance(node)
                #print('\nleaf_distance')
                #print(leaf_distance)
                # Compare to threshold.
                if leaf_distance > float(quartile3 +\
                        (float(max_bl_iqr_above_third_quartile) * interquartile_range)):
                    if node.name not in names_of_seqs_to_remove:
                        names_of_seqs_to_remove.append(node.name)


    # Open output file and write list of sequences to remove.
    with open(to_remove_list_file, 'w') as o:
        for name in names_of_seqs_to_remove:
            o.write(name + '\n')

    # Check that file was produced.
    assert os.path.isfile(to_remove_list_file), """Did not produce file listing
    sequences to remove from dataset."""


def get_key_nodes_list(to_remove_list_file):
    """Get list of names of nodes to remove from input text file.
    """
    #key_nodes_list = key_node_numbers.split(' ')
    key_nodes_list = []
    with open(to_remove_list_file) as infh:
        for i in infh:
            if not i.startswith('#') and not i.startswith('\n'):
                key_nodes_list.append(i.strip())

    # Check that some nodes were found.
    #assert len(key_nodes_list) > 0, """No nodes identified for removal."""

    # Return list of key nodes.
    return key_nodes_list


def get_seqs_to_remove_from_dataset(key_nodes_list, t):
    """Take a list of nodes and an ete3 Tree object, and return a list of
    sequences to remove from dataset.
    """
    seqs_to_remove_from_dataset = []
    for node in t.traverse():
        if node.name in key_nodes_list:
            # Get leaves.
            seqs_to_remove_from_dataset = seqs_to_remove_from_dataset + [x.name for x in node.get_leaves()]
        else:
            for x in key_nodes_list:
                if node.name.split(' ')[0] == x or node.name.split('__')[0] == x:
                    seqs_to_remove_from_dataset = seqs_to_remove_from_dataset + [node.name]
                    break
    # Return list of sequences to remove from dataset.
    return seqs_to_remove_from_dataset


def write_reduced_alignment(alignment_file,
                            coded_seqs_to_remove_from_dataset,
                            conversion_table_dict,
                            output_alignment_file):
    """Take an alignment and a list of sequences to remove, and write a new
    alignment file with those sequences removed.
    """
    # Remove extra mesquite lines if present.
    delete_extra_mesquite_lines(alignment_file)

    coded_seqs_to_not_remove_from_dataset = []
    with open(alignment_file) as infh, open(output_alignment_file, 'w') as o:
        # Parse input nexus alignment.
        alignment = AlignIO.read(infh, 'nexus')

        # Report number of sequences.
        #print('Number of seqs in original alignment:')
        #print(len(alignment))
        #print('\nIDs of seqs in original alignment:')
        #for i in alignment:
        #    #print(i.id + ' (' + conversion_table_dict[i.id] + ')')
        #    print(i.id)

        # Report sequences to remove.
        #print('\nSeqs to remove:')
        #for i in coded_seqs_to_remove_from_dataset:
        #    print(i + ' (' + conversion_table_dict[i] + ')') 

        # Initiate a new alignment object.
        alignment2 = MultipleSeqAlignment([])

        #print('number of seqs in alignment2:')
        #print(len(alignment2))
        #print('Removing sequences from alignment.')
        #ali_len = len(alignment)

        # Loop over sequences in input alignment.
        inum = -1 
        for i in alignment:
            inum += 1
            #if i.id not in coded_seqs_to_remove_from_dataset:
            #print(i.id)
            # Find corresponding key in dict, if present.




            # *****************************************


            # Determine key to use.
            #print(i.id)
            key_to_use = None
            for key in conversion_table_dict.keys():
                #print(key)
                #if key.startswith('ConsensusfromContig11299'):
                #    print('\tkey: ' + key)
                if i.id.strip() == conversion_table_dict[key].strip():
                    key_to_use = key
                    break
                elif ' ' in i.id:
                    if conversion_table_dict[key].startswith(i.id.split(' ')[0]):
                        key_to_use = key
                        break
                    if len(i.id.split(' ')) > 2:
                        if conversion_table_dict[key].startswith(i.id.split(' ')[0] + '_' + i.id.split(' ')[1]):
                            key_to_use = key
                            break
                if key_to_use is None and '_' in key:
                    # ***This is potentially problematic! ***
                    if i.id.startswith(conversion_table_dict[key].split('_')[0]):
                        key_to_use = key


            assert key_to_use is not None, """Could not identify name in
            conversion table that corresponds to name in alignment: %s"""\
            % i.id

            #if conversion_table_dict[key_to_use] not in coded_seqs_to_remove_from_dataset:
            if key_to_use not in coded_seqs_to_remove_from_dataset:
                alignment2.append(alignment[inum])
                coded_seqs_to_not_remove_from_dataset.append(conversion_table_dict[key_to_use])


            ## This code works if the input alignment has coded instead of
            ## original names.

            ## If not in sequences identified for removal, add to new alignment.
            #if i.id not in coded_seqs_to_remove_from_dataset:
            #    alignment2.append(alignment[inum])
            #    # Add to list.
            #    coded_seqs_to_not_remove_from_dataset.append(i.id)
            #else:
            #    print('removed ' + i.id)


            # *****************************************




        # RESOLVED: PROBLEM: the sequences written to the output alignment file are wrong...
        #print('\nNumber of seqs in new alignment:')
        #print(len(alignment2))
        #print('\nIDs of seqs in new alignment:')
        #for i in alignment2:
        #    print(i.id)# + ' (' + conversion_table_dict[i.id] + ')')
        #print('\nNum seqs to remove:')
        #print(len(set(coded_seqs_to_remove_from_dataset)))
        #print('expected number of seqs in new alignment:')
        #print(len(alignment) - len(set(coded_seqs_to_remove_from_dataset)))

        # Check that expectation is met.
        assert len(alignment) - len(set(coded_seqs_to_remove_from_dataset))\
        == len(alignment2), """Number of sequences in new alignment does not
        match expectation."""

        #alignment3 =\
        #MultipleSeqAlignment(alignment2[len(coded_seqs_to_remove_from_dataset) -1:,:])
        alignment3 =\
        MultipleSeqAlignment(alignment2[:,:])
        #print('number of seqs in new alignment:')
        #print(len(alignment3))

        # Get list of sequences in new alignment.
        ids_in_alignment3 = []
        for i in alignment3:
            ids_in_alignment3.append(i.id)

        assert set(ids_in_alignment3) ==\
        set(coded_seqs_to_not_remove_from_dataset), """Sequences in alignment
        do not have the expected set of IDs."""

        ## Decode names in alignment. (actually codes names)
        #print('Decoding names in alignment.') 
        #for i in alignment3:
        #    # Decode name.
        #    for j in conversion_table_dict.keys():
        #        if conversion_table_dict[j] == i.id:
        #            i.id = j
        #            break

        # Write reduced alignment.
        #print('Writing modified alignment file:\n\t%s' % o.name)
        AlignIO.write(alignment3, o, 'nexus')


def manually_select_nodes_and_remove_seqs(input_tree_one, alignment_file,
        output_file=None, name_replace_table=None):
    """Take a tree file (newick) and the corresponding alignment file (nexus)
    and write a new alignment without the selected nodes/sequences.

    ToDo: Automate selection of long branches for removal. Allow manual
    selection of leaf nodes for removal.

    Problem: Interactive tree viewer is not working.
    """
    ## Check that necessary inputs are specified.
    #assert args.input_tree_one != None, "Error: Missing input."
    #assert args.name_replace_table != None, "Error: Missing input."

    # Define output alignment name.
    output_alignment_file = None
    if output_file is None:
        output_alignment_file = alignment_file.rsplit('.', 1)[0] + '_reduced.nex'
    else:
        output_alignment_file = output_file


    # Make a new table file for putting quotation marks around taxon names in
    # newick tree.
    #table_with_quotation_marks = name_replace_table + 'TEMP1'
    #make_quoted_table(name_replace_table, table_with_quotation_marks) # untested 

    # Generate a dictionary from the new table file for conversion of names.
    #conversion_table_dict = get_conversion_table_dict(table_with_quotation_marks)
    #conversion_table_dict = get_conversion_table_dict(name_replace_table)
    conversion_table_dict = get_conversion_dict_from_table(name_replace_table)

    # Run name_replace.pl using new table so that long taxon names can be
    # visualized with ete3.
    tree_with_long_names = input_tree_one + '_TEMP1'
    #subprocess.call(['name_replace.pl', '-b', input_tree_one,
    #    tree_with_long_names, table_with_quotation_marks])
    write_newick_tree_with_uncoded_names(input_tree_one, tree_with_long_names,
            name_replace_table, True)


    # Change node support values to a format that is readable by ETE3.
    treestring = get_tree_string_with_readable_supports(tree_with_long_names)

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(tree_with_long_names, quoted_node_names=True)


    # Define name for tree without branch lengths.
    simple_tree = tree_with_long_names.rsplit('_', 1)[0] + '_TEMP2' 

    # Write simple tree to a new file 
    t1.write(format=9, outfile=simple_tree)
    #t1.write(format=0, outfile=simple_tree)

    # Parse simple tree.
    t2 = Tree(simple_tree)

    # Print simple tree.
    print('Input tree:')
    print(t2)

    # Root tree.
    root_tree_manually(t1, t2)

    # Show new tree rooting.
    print('New tree rooting:')
    t1.write(format=9, outfile=simple_tree)
    t2a = Tree(simple_tree)
    print(t2a)

    # Generate an ete3 tree style.
    ts = get_prune_tree_style()

    # Give all internal nodes unique names, and then visualize with the
    # interactive tree viewer.
    edge = 0
    for node in t1.traverse():
       if not node.is_leaf():
          node.name = str(edge)
          edge += 1

    # Render tree as a pdf.
    tree_file_name = input_tree_one + '_ete3_render.pdf'
    print('Rendering tree.')
    t1.render(tree_file_name, w=183, units="mm", tree_style=ts) #Temp. commented out.

    # Define path of file to list sequences to remove.
    to_remove_list_file = alignment_file + '_to_remove_list.txt'

    # Prompt user to list nodes and sequence names to remove in a file.
    get_seqs_to_remove_file_from_user(to_remove_list_file, tree_file_name)

    # Get list of nodes from file.
    key_nodes_list = get_key_nodes_list(to_remove_list_file)

    # Compile a list of sequences to remove from dataset.
    seqs_to_remove_from_dataset =\
    get_seqs_to_remove_from_dataset(key_nodes_list, t1)

    # Record list of sequences to remove.
    with open(to_remove_list_file, 'w') as o:
        #print('\nSequences to remove:')
        for i in seqs_to_remove_from_dataset:
            print('\t' + i)
            o.write(i + '\n')

    # Convert listed node names back to coded names.
    coded_seqs_to_remove_from_dataset = []
    for name in seqs_to_remove_from_dataset:

        # This code is for if the alignment were coded:
        for key in conversion_table_dict.keys():
            if conversion_table_dict[key] == name:
                coded_name = key
                coded_seqs_to_remove_from_dataset.append(coded_name.strip())

        ## This code is for if the alignment were not coded???:
        #coded_name = conversion_table_dict[name]
        #coded_seqs_to_remove_from_dataset.append(coded_name.strip())

    #print('Coded names:')
    #print(coded_seqs_to_remove_from_dataset)

    # Make a copy of the input alignment file with the listed sequences
    # removed.
    write_reduced_alignment(alignment_file,
                            coded_seqs_to_remove_from_dataset,
                            conversion_table_dict,
                            output_alignment_file)

    # Delete intermediate files.
    os.remove(simple_tree)
    os.remove(tree_with_long_names)
    #os.remove(table_with_quotation_marks)

    # Return main output path.
    return output_alignment_file


def automatically_select_nodes_and_remove_seqs(input_tree_one,
                                               alignment_file,
                                               max_bl_iqr_above_third_quartile,
                                               remove_redun_seqs,
                                               remove_redun_seqs_threshold,
                                               output_file=None,
                                               name_replace_table=None):
    """Take a tree file (newick) and the corresponding alignment file (nexus)
    and write a new alignment without the selected nodes/sequences.

    ToDo: Automate selection of long branches for removal. Allow manual
    selection of leaf nodes for removal.

    Problem: Interactive tree viewer is not working.
    """
    ## Check that necessary inputs are specified.
    #assert args.input_tree_one != None, "Error: Missing input."
    #assert args.name_replace_table != None, "Error: Missing input."

    # Define output alignment name.
    output_alignment_file = None
    if output_file is None:
        output_alignment_file = alignment_file.rsplit('.', 1)[0] + '_reduced.nex'
    else:
        output_alignment_file = output_file


    # Make a new table file for putting quotation marks around taxon names in
    # newick tree.
    #table_with_quotation_marks = name_replace_table + 'TEMP1'
    #make_quoted_table(name_replace_table, table_with_quotation_marks) # untested 

    # Generate a dictionary from the new table file for conversion of names.
    #conversion_table_dict = get_conversion_table_dict(table_with_quotation_marks)
    #conversion_table_dict = get_conversion_table_dict(name_replace_table)
    conversion_table_dict = get_conversion_dict_from_table(name_replace_table)

    # Run name_replace.pl using new table so that long taxon names can be
    # visualized with ete3.
    tree_with_long_names = input_tree_one + '_TEMP1'
    #subprocess.call(['name_replace.pl', '-b', input_tree_one,
    #    tree_with_long_names, table_with_quotation_marks])
    write_newick_tree_with_uncoded_names(input_tree_one, tree_with_long_names,
            name_replace_table, True)

    #*********
    # Change node support values to a format that is readable by ETE3.
    #treestring = get_tree_string_with_readable_supports(tree_with_long_names)

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    #t1 = Tree(tree_with_long_names, quoted_node_names=True)

    # Parse tree using ete3, with branch supports as internal node names.
    t1 = Tree(tree_with_long_names, format=1)

    # Translate node names (with aLRT/bootstrap percent support) to branch
    # support (bootstrap proportion).
    translate_int_node_names_to_support(t1)
    #*********

    # Define name for tree without branch lengths.
    simple_tree = tree_with_long_names.rsplit('_', 1)[0] + '_TEMP2' 

    # Write simple tree to a new file 
    t1.write(format=9, outfile=simple_tree)
    #t1.write(format=0, outfile=simple_tree)

    # Parse simple tree.
    t2 = Tree(simple_tree)

    # Print simple tree.
    #print('Input tree:')
    #print(t2)



    # Begin: Differs from manually_select_nodes_and_remove_seqs**************************************************



    # Root tree on midpoint.
    #root_tree_manually(t1, t2)
    t1.set_outgroup(t1.get_midpoint_outgroup())

    # Show new tree rooting.
    #print('New tree rooting:')
    #t1.write(format=9, outfile=simple_tree)
    #t2a = Tree(simple_tree)
    #print(t2a)

    # Generate an ete3 tree style.
    #ts = get_prune_tree_style()

    # Give all internal nodes unique names, and then visualize with the
    # interactive tree viewer.
    edge = 0
    for node in t1.traverse():
       if not node.is_leaf():
          node.name = str(edge)
          edge += 1

    ## Render tree as a pdf.
    #tree_file_name = input_tree_one + '_ete3_render.pdf'
    #print('Rendering tree.')
    #t1.render(tree_file_name, w=183, units="mm", tree_style=ts) #Temp. commented out.

    # Define path of file to list sequences to remove.
    to_remove_list_file = alignment_file + '_to_remove_list.txt'

    # Identify and list nodes and sequence names to remove in a file.
    #get_seqs_to_remove_file_from_user(to_remove_list_file, tree_file_name)
    get_seqs_to_remove_file(to_remove_list_file, t1,
                            remove_redun_seqs, remove_redun_seqs_threshold,
                            max_bl_iqr_above_third_quartile)



    # End: Differs from manually_select_nodes_and_remove_seqs**************************************************




    # Get list of nodes from file.
    key_nodes_list = get_key_nodes_list(to_remove_list_file)

    # Don't do anything more if no nodes are identified for removal.
    if len(key_nodes_list) == 0:
        print('No nodes identified for removal in alignment %s' % alignment_file)
        print('Not proceeding.')
        return 'None'


    # Compile a list of sequences to remove from dataset.
    seqs_to_remove_from_dataset =\
    get_seqs_to_remove_from_dataset(key_nodes_list, t1)

    # Record list of sequences to remove.
    with open(to_remove_list_file, 'w') as o:
        #print('\nSequences to remove:')
        for i in seqs_to_remove_from_dataset:
            #print('\t' + i)
            o.write(i + '\n')

    # Convert listed node names back to coded names.
    uncoded_leaf_names_to_remove = []
    coded_seqs_to_remove_from_dataset = []
    for name in seqs_to_remove_from_dataset:
        # Also add uncoded name to list.
        uncoded_leaf_names_to_remove.append(name)

        # This code is for if the alignment were coded:
        for key in conversion_table_dict.keys():
            if conversion_table_dict[key] == name:
                coded_name = key
                coded_seqs_to_remove_from_dataset.append(coded_name.strip())

        ## This code is for if the alignment were not coded???:
        #coded_name = conversion_table_dict[name]
        #coded_seqs_to_remove_from_dataset.append(coded_name.strip())

    #print('Coded names:')
    #print(coded_seqs_to_remove_from_dataset)

    # Make a copy of the input alignment file with the listed sequences
    # removed.
    write_reduced_alignment(alignment_file,
                            coded_seqs_to_remove_from_dataset,
                            conversion_table_dict,
                            output_alignment_file)

    # Delete intermediate files.
    os.remove(simple_tree)
    os.remove(tree_with_long_names)
    #os.remove(table_with_quotation_marks)

    # Render a tree as pdf figures.
    method = 'iqtree'
    timestamp = 'N/A'
    taxa_to_root_on = []
    highlight_paralogues = False
    highlight_for_removal = uncoded_leaf_names_to_remove
    file_with_subs_model_name = None
    pdf_files = visualize_tree(method,
                               timestamp,
                               taxa_to_root_on,
                               highlight_paralogues,
                               highlight_for_removal,
                               input_tree_one,
                               name_replace_table,
                               file_with_subs_model_name
                               )

    # Return main output path.
    return output_alignment_file


def automatically_select_nodes_and_remove_seqs_in_dir(indp,
                                                      max_bl_iqr_above_third_quartile,
                                                      remove_redun_seqs,
                                                      remove_redun_seqs_threshold,
                                                      output_file=None):
    """Identify appropriate files and call function to use those to remove
    seqs.
    """
    # Determine tree method ('iqtree', 'mb', etc.).
    # ...

    # Identify input tree file (output from a phylogenetic analysis program,
    # such as IQ-tree).
    input_tree_one = os.path.join(indp, 'output.treefile') # Assumes IQtree. 

    # Check that the tree file was found.
    assert os.path.isfile(input_tree_one), """Could not find tree file at path
    %s.""" % input_tree_one

    # Identify input alignment file (the one with sequences to be removed) with
    # original sequence names.
    disqualifying_strings = ['mask', 'trim'] # Not always necessary? May be for srch_ali_space.
    #disqualifying_strings = []
    alignment_file = find_input_file_in_parent_directory(indp, 'nex',
            disqualifying_strings)
    print('\nUsing this alignment file:\n%s\n' % alignment_file)

    # Identify name replace table to use for decoding names in input tree file
    # so that they match those in the input alignment.
    name_replace_table = find_input_file_in_parent_directory(indp, 'table', [])
    print('\nUsing this file for decoding names:\n%s\n' % name_replace_table)

    # Call function on appropriate files.
    main_out_path = automatically_select_nodes_and_remove_seqs(input_tree_one,
                                               alignment_file,
                                               max_bl_iqr_above_third_quartile,
                                               remove_redun_seqs,
                                               remove_redun_seqs_threshold,
                                               output_file,
                                               name_replace_table)

    # Return path of main output file.
    return main_out_path 



def manually_select_nodes_to_constrain(input_tree_one, name_replace_table,
        output_tree_file=None):
    """Take a tree file (newick) and the corresponding alignment file (nexus)
    and write a new alignment without the selected nodes/sequences.

    ToDo: Automate selection of long branches for removal. Allow manual
    selection of leaf nodes for removal.
    """
    ## Check that necessary inputs are specified.
    #assert args.input_tree_one != None, "Error: Missing input."
    #assert args.name_replace_table != None, "Error: Missing input."

    # Define output alignment name.
    if output_tree_file is None:
        output_tree_file = input_tree_one + '_constraint_tree.txt'

    # Make a new table file for putting quotation marks around taxon names in
    # newick tree.
    table_with_quotation_marks = name_replace_table + 'TEMP1'
    make_quoted_table(name_replace_table, table_with_quotation_marks) # untested 

    # Generate a dictionary from the new table file for conversion of names.
    conversion_table_dict = get_conversion_table_dict(table_with_quotation_marks)

    ## Run name_replace.pl using new table so that long taxon names can be
    ## visualized with ete3.
    #tree_with_long_names = input_tree_one + '_TEMP1'
    #subprocess.call(['name_replace.pl', '-b', input_tree_one,
    #    tree_with_long_names, table_with_quotation_marks])

    # Convert names using new table, so that tree can be parsed using ete3.
    tree_with_long_names = input_tree_one + '_TEMP1'
    write_newick_tree_with_coded_names(input_tree_one, tree_with_long_names,
            table_with_quotation_marks)

    # Change node support values to a format that is readable by ETE3.
    treestring = None
    with open(tree_with_long_names) as infh:
        num = 0
        for i in infh:
            num += 1
            treestring = i.strip()
            if num == 1:
                break
    if treestring is not None:
        with open(tree_with_long_names, 'w') as o:
            o.write(reformat_combined_supports2(treestring))

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(tree_with_long_names, quoted_node_names=True)

    # Define name for tree without branch lengths.
    simple_tree = tree_with_long_names.rsplit('_', 1)[0] + '_TEMP2' 

    # Write simple tree to a new file 
    t1.write(format=9, outfile=simple_tree)
    #t1.write(format=0, outfile=simple_tree)

    # Parse simple tree.
    t2 = Tree(simple_tree)

    # Print simple tree.
    print('Input tree:')
    print(t2)

    # Get node for a key clade. 
    # Or, just require to root it before inputting (figtree can export
    # re-rooted newick)?????
    print("\nPlease define an arbitrary clade of interest for rooting the tree topology\n")
    clade_member_1 = input('Clade member 1: ')
    clade_member_2 = input('Clade member 2: ')
    #ancestor1 = t2.get_common_ancestor(clade_member_1, clade_member_2)
    ancestor1 = t1.get_common_ancestor(clade_member_1, clade_member_2)
    #print(ancestor1)

    # Re-root tree. 
    t1.set_outgroup(ancestor1)

    # Show new tree rooting.
    print('New tree rooting:')
    t1.write(format=9, outfile=simple_tree)
    t2a = Tree(simple_tree)
    print(t2a)

    # Set up tree style for viewing.
    ts = TreeStyle()
    ts.show_leaf_name = False
    def my_layout(node):
        """Function for showing internal node names (numbers)."""
        #F = TextFace(node.name, tight_text=True)
        F = TextFace(node.name)
        add_face_to_node(F, node, column=0, position="branch-right")
    ts.scale =  100
    ts.branch_vertical_margin = 10
    ts.layout_fn = my_layout
    ts.show_branch_support = True

    # Give all internal nodes unique names, and then visualize with the
    # interactive tree viewer.
    edge = 0
    for node in t1.traverse():
       if not node.is_leaf():
          node.name = str(edge)
          edge += 1

    # Show tree with node names (render as PDF file).
    #t1.show(tree_style=ts)
    tree_file_name = input_tree_one + '_ete3_render.pdf'
    t1.render(tree_file_name, w=183, units="mm", tree_style=ts)
    subprocess.call(['open', tree_file_name])

    # Prompt user to list nodes and sequence names to remove in a file.
    to_remove_list_file = input_tree_one + '_seqs_to_constrain_into_clade.txt'
    print('''\nList node numbers and/or sequence accessions for sequences to
    constrain into one clade in output tree file:\n%s''' % to_remove_list_file)
    editor = os.environ.get('EDITOR','vim') 
    initial_message = '# List node numbers or sequence accessions for sequences constrain into a clade below.'
    with open(to_remove_list_file, 'w') as o:
        o.write(initial_message)
        o.flush()
        subprocess.call([editor, o.name])

    # Get list of nodes from file.
    #key_nodes_list = key_node_numbers.split(' ')
    key_nodes_list = []
    with open(to_remove_list_file) as infh:
        for i in infh:
            if not i.startswith('#') and not i.startswith('\n'):
                key_nodes_list.append(i.strip())


    # Compile a list of sequences to remove from dataset, and a list of all
    # other sequences.
    seqs_to_remove_from_dataset = []
    for node in t1.traverse():
        if node.name in key_nodes_list:
            # Get leaves.
            seqs_to_remove_from_dataset = seqs_to_remove_from_dataset + [x.name for x in node.get_leaves()]
        else:
            for x in key_nodes_list:
                if node.name.split(' ')[0] == x:
                    seqs_to_remove_from_dataset = seqs_to_remove_from_dataset + [node.name]
                    break
    other_seqs = []
    for node in t1.get_leaves():
        if node.name not in seqs_to_remove_from_dataset:
            other_seqs = other_seqs + [node.name]

    # Record list of sequences to remove.
    with open(to_remove_list_file, 'w') as o:
        print('\nSequences to constrain into a single clade:')
        for i in seqs_to_remove_from_dataset:
            print('\t' + i)
            o.write(i + '\n')

    # Convert listed node names back to coded names.
    coded_seqs_to_remove_from_dataset = []
    for name in seqs_to_remove_from_dataset:
        coded_name = conversion_table_dict[name]
        coded_seqs_to_remove_from_dataset.append(coded_name.strip())

    # Convert other node names back to coded names.
    other_coded_seqs_to_remove_from_dataset = []
    for name in other_seqs:
        coded_name = conversion_table_dict[name]
        other_coded_seqs_to_remove_from_dataset.append(coded_name.strip())

    # Write newick tree file as constraint.
    with open(output_tree_file, 'w') as o:
        o.write('((')
        o.write(','.join(coded_seqs_to_remove_from_dataset))
        o.write('),')
        o.write(','.join(other_coded_seqs_to_remove_from_dataset))
        o.write(');')

    # Delete intermediate files.
    os.remove(simple_tree)
    os.remove(tree_with_long_names)
    os.remove(table_with_quotation_marks)

    # Return main output path.
    return output_tree_file


def remove_nodes_to_match_ali(alignment, tree_file, output_file):
    """Write tree file with only sequences represented in the alignment.
    """
    # Set output file path.
    if output_file is None:
        output_file = tree_file.rsplit('.', 1)[0] + '_reduced.tre'

    # Remove extra mesquite lines if present.
    delete_extra_mesquite_lines(alignment)

    # Parse alignment file.
    ali = None
    with open(alignment) as infh:
        ali = AlignIO.read(infh, "nexus") 

    # Construct a list of all sequence names.
    ali_seq_names = []
    for i in ali:
        ali_seq_names.append(i.id)

    # Parse tree file.
    t1 = Tree(tree_file)

    # Make a list of terminal/leaf nodes in the tree that are not represented
    # in the alignment, and those that are.
    in_tree_in_ali = []
    in_tree_not_in_ali = []
    for node in t1.traverse():
        if node.is_leaf():
            in_ali = False
            for i in ali_seq_names:
                if node.name.strip('\'') == i or node.name.strip('\'').replace(' ', '_') == i.replace(' ', '_'):
                    in_tree_in_ali.append(node.name)
                    in_ali = True
            if not in_ali:
                in_tree_not_in_ali.append(node.name)


    # Check that all the sequences in the alignment are in the tree.
    assert len(in_tree_in_ali) >= len(ali_seq_names), """Missing seqs in input
    tree."""

    # Make a copy of the parsed tree object for manipulation.
    t2 = t1.copy()

    # Remove all the leaf nodes that do not have corresponding seqs in the
    # alignment.
    for node in t2.traverse():
        if node.is_leaf():
            if node.name in in_tree_not_in_ali:
                #print(node.name)
                node.delete()

    # Write the modified tree object to newick format at the given path.
    t2.write(outfile=output_file, format=9)

    # Return main output path.
    return output_file

