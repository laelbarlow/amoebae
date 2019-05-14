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

"""Contains functions for classifying sequences using existing phylogenetic models.
"""
import settings
import sys
import os
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped
from module_afa_to_nex import delete_extra_mesquite_lines, afa_to_nex, nex_to_afa, nex_to_phylip
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
from module_amoebae_trim_nex import trim_nex
from module_paralogue_counter import get_seq_obj_from_srch_res_csv_info
import pandas as pd
import module_amoebae_column_header_lists
from module_amoebae import mask_nex2
from module_paralogue_counter import add_seq_to_alignment3,\
modify_seq_descr_for_tree
from module_amoebae_name_replace import write_afa_with_code_names






# Need to automatically align a given sequence to an existing alignment (high
# gap extension penalty to prevent wonky alignments?) and trim to get a final
# alignment. (there should be some checking to make sure the right domain was
# aligned, for example in the case of proteins that contain two homologous
# domains such as Qbc SNAREs).


def add_seq_to_alignment(inseqobj, innexpath, outdp):
    """Takes a sequence and adds it to a trimmed alignment for a backbone tree
    (which will be used to classify the sequence).
    """
    # Write input sequence object to a file temporarily for input to MUSCLE.
    fa_temp_1 = os.path.join(outdp, os.path.basename(innexpath) + '_temp1.fa')
    SeqIO.write(inseqobj, fa_temp_1, 'fasta')

    # Temporarily save a copy of the input alignment (nexus format) to fasta
    # file.
    inali = AlignIO.read(innexpath, 'nexus')
    fa_temp_2 = fa_temp_1.replace('_temp1.fa', '_temp2.afa')
    AlignIO.write(inali, fa_temp_2, 'fasta')

    # Check that the input sequence does not have the same name as any of the
    # sequences in the input alignment.
    # Get list of sequence names from alignment.
    existing_seq_names = []
    with open(fa_temp_2) as infh:
        for record in AlignIO.read(infh, 'fasta'):
            extra_text = re.compile(r'\.copy$')
            rec_id_no_extra = extra_text.sub('', record.description) 
            existing_seq_names.append(rec_id_no_extra)
            #existing_seq_names.append(record.id)
    # Get additional sequence name.
    additional_seq_name = inseqobj.description
    # Determine if additional sequence name is in existing names.
    for i in existing_seq_names:
        assert i != additional_seq_name, "Error: Additional sequence with name\
 " + additional_seq_name + " has the same name as a sequence in the\
 input alignment. Make sure that you give the input sequences unique names."

    # Align sequence to alignment (make sure that gap penalties are optimized).
    fa_temp_3 = fa_temp_2.replace('_temp2.afa', '_temp3.afa')
    subprocess.call(["muscle", "-profile", "-in1", fa_temp_2, "-in2",
        fa_temp_1, "-out", fa_temp_3])

    # Convert alignment to nexus format.
    nex_out_1 = fa_temp_3.rsplit('_', 1)[0] + '_plus_' + inseqobj.id.split(' ', 1)[0] + '.nex'
    afa_to_nex(fa_temp_3, nex_out_1) 

    # Mask alignment (mask out any positions in which the last sequence, the
    # additional sequence, has a residue but none of the other sequences have
    # residues at that position).
    nex_out_2_mask = nex_out_1.rsplit('.', 1)[0] + '_mask.nex'
    mask_nex2(nex_out_1, nex_out_2_mask) 

    # Trim alignment.
    nex_out_3_trim = nex_out_2_mask.rsplit('_', 1)[0] + '_trim.nex'
    trim_nex(nex_out_2_mask, nex_out_3_trim)

    # Delete temporary files.
    os.remove(fa_temp_1)
    os.remove(fa_temp_2)
    os.remove(fa_temp_3)
    os.remove(nex_out_1)

    # Return file path of trimmed nexus file for downstream processes.
    return nex_out_3_trim
    

def code_names_in_ali(indp, inalifp, outalifp, outtablefp):
    """Code names using name_replace.pl.
    """
    # Code for doing this with the old name_replace.pl script:
    #subprocess.call(['name_replace.pl', '-f', inalifp, outalifp,\
    #outtablefp], cwd=indp)
    #subprocess.call(['name_replace.pl', '-f', inalifp, outalifp,\
    #outtablefp])

    # Call function from module_amoebae_name_replace.py.
    write_afa_with_code_names(inalifp, outalifp, outtablefp)


def get_conversion_dict(tablefile):
    """Take a Name_Replace-style taxon name conversion table file, and return a
    dictionary for converting between coded names and long names, with long
    names as keys and codes as values.
    """
    conv_dict = {}
    with open(tablefile) as tablefh:
        for i in tablefh:
            if i.startswith('ZZ'):
                cur_code = i.strip()
            else:
                conv_dict[i.strip().replace(':', '_').strip('_')] = cur_code
    return conv_dict


def quote_tree(intree, table):
    """Add quotation marks around long names in tree.
    """
    outtree = intree
    conv_dict = get_conversion_dict(table)

    for key in conv_dict.keys():
        outtree = outtree.replace(key, "\"" + key + "\"")

    #for value in conv_dict.values():
    #    outtree = outtree.replace(value, "\"" + key + "\"")

    return outtree


def code_tree(intreefp, outtreefp, table):
    """Code names in a tree, given a conversion table.
    """
    # Construct a name conversion dict from input table file.
    conv_dict = get_conversion_dict(table)

    # Parse input tree file.
    intree = Tree(intreefp, quoted_node_names=True)

    ## Check that the number of keys in the conversion dict is the same as the
    ## number of leaf nodes in the tree.
    #assert len(conv_dict.keys()) == len(intree.get_leaves()), """Different
    #number of code names than leaf nodes in tree."""

    ## Check that each of the node names is a key in the dictionary.
    #for node in intree.iter_leaves():
    #    assert node.name in conv_dict.keys(), """Node name not in name
    #    conversion dictionary: %s""" % node.name

    # Convert node names in tree.
    for node in intree.iter_leaves():
        node.name = conv_dict[node.name.strip('\"')]

    # Write modified tree to output file.
    intree.write(format=9, outfile=outtreefp)


def uncode_tree(intreefp, outtreefp, table):
    """Code names in a tree file, given a conversion table. Without parsing
    with ete3.
    """
    # Construct a name conversion dict from input table file.
    conv_dict = get_conversion_dict(table)

    with open(intreefp) as infh, open(outtreefp, 'w') as o:
        for t in infh:
            new_text = t
            for longname in conv_dict.keys():
                new_text = new_text.replace(conv_dict[longname], '\"' +\
                        longname + '\"')
            o.write(new_text)


def uncode_tree_obj(intreeobj, table):
    """Uncode names in a tree, given a conversion table.
    """
    # Construct a name conversion dict from input table file.
    conv_dict = get_conversion_dict(table)

    # Convert node names in tree.
    for node in intreeobj.iter_leaves():
        uncoded_name = None
        for key in conv_dict.keys():
            if conv_dict[key] == node.name:
                uncoded_name = key
                break

        # Check that an uncoded name was found.
        assert uncoded_name is not None

        # Change name of node to uncoded name.
        node.name = uncoded_name


def classify_seq_with_constrained_tree(alignment, tree, subs_model, type_seqs, fasta):
    """Classify sequences in a given fasta file using a given backbone tree.
    """

    # Define timestamp.
    cur_time = time.strftime("%Y_%m_%d_%H_%M_%S")
    print('\nTimestamp: ' + cur_time)
    
    # Make output directory.
    outdirpath1 = alignment.rsplit('.', 1)[0] + '_classify_seq_' + cur_time
    os.mkdir(outdirpath1)
    
    ## Take a specific approach, depending on the options specified at runtime.
    #
    ## If bootstrapping option not specified, then use default topology-testing
    ## approach.
    #if not args.boot:
    #
    #    # Call function for doing default analysis.
    
    # Write initial summary csv file.
    summary_csv_fp = os.path.join(outdirpath1, '0_' +\
            os.path.basename(alignment).rsplit('.',\
        1)[0] + '_Summary.csv')
    summary_header_line = ','.join(['Original alignment file', #0
                                    'Original tree file', #1
                                    'Substitution model', #2
                                    'Sequence added to tree', #3
                                    'Hypothesis (tree/clade) number', #4
                                    'Type sequence of clade the additional sequence was constrained into', #5
                                    'Clade name', #6
                                    'Log-likelihood of ML tree',  #7
                                    'bp-RELL: bootstrap proportion using RELL method (Kishino et al. 1990)',  #8
                                    'p-KH: p-value of one sided Kishino-Hasegawa test (1989)',  #8
                                    'p-SH: p-value of Shimodaira-Hasegawa test (2000)',  #8
                                    'c-ELW: Expected Likelihood Weight (Strimmer & Rambaut 2002)',  #8
                                    '"p-AU: p-value of approximately unbiased (AU) test (Shimodaira, 2002)"',  #8
                                    'Accept (+) or reject (-) hypothesis based on consensus of topology tests',  #9
                                    'Masked alignment file', #10
                                    'Sequence name conversion table file', #11
                                    'Trimmed input alignment file', #12
                                    'Constraint tree file', #13
                                    'Starting tree file', #14
                                    'ML tree file', #15
                                    '\n'
                                    ])
    
    with open(summary_csv_fp, 'w') as o:
        o.write(summary_header_line)

    # Align sequences to alignment and trim.
    for record in SeqIO.parse(fasta, 'fasta'):
        add_seq_to_alignment(record, alignment, outdirpath1)
        
    # Take tree topology and alignment (regular names), convert names in alignment
    # generating a conversion table similar to for name_replace.pl, then use that
    # conversion table to code the names in the tree topology as well.
    ali_num = 0
    for ali in glob.glob(os.path.join(outdirpath1, '*trim.nex')):
        ali_num += 1
        
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
    
        # Only do the following steps once for all the alignments, because all the
        # constrained topology analyses will use the same constraint trees.
        if ali_num == 1:
    
            # Put quotation marks around names in constraint tree. # SEPARATE
            outtreefp = ali.rsplit('.', 1)[0] + '.constraint_tree0.Q.newick'
            with open(tree) as intreefh, open(outtreefp, 'w') as o:
                for t in intreefh:
                    o.write(quote_tree(t, outtablefp))
    
    
            # Make topology constraint trees.   
    
            # Parse backbone topology and root on a clade rather than having the root
            # inside a clade.
    
            ######################
            # Parse tree using ete3.
            # Note: parentheses and commas get replaced with underscores.
            t1 = Tree(outtreefp, quoted_node_names=True)
    
            # Define name for tree without branch lengths.
            simple_tree = outtreefp.rsplit('_', 1)[0] + '_TEMP2' 
    
            # Write simple tree to a new file 
            t1.write(format=9, outfile=simple_tree)
    
            # Parse simple tree.
            t2 = Tree(simple_tree)
    
            # Print simple tree.
            print('Input tree:')
            print(t2)
    
            ## Get node for a key clade. 
            ## Or, just require to root it before inputting (figtree can export
            ## re-rooted newick)?????
            #print("\nPlease define an arbitrary clade of interest for rooting the tree topology\n")
            #clade_member_1 = input('Clade member 1: ')
            #clade_member_2 = input('Clade member 2: ')
            #ancestor1 = t2.get_common_ancestor(clade_member_1, clade_member_2)
            ##print(ancestor1)
    
            ## Re-root tree on a clade (as a polytomy) rather than inside a clade, so
            ## that key branches can be more easily identified.
            #t2.set_outgroup(ancestor1)
    
            ######################
    
    
            # Get to topologies sorted out...
    
            # Get constraint topologies for each clade of interest to test
            # alternative topology hypotheses.
    
            # Get list of "type" sequences from input.
            type_seq_list = []
            for i in open(type_seqs):
                type_seq_list.append(i.strip().split(',')[0])
    
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
    
                    # Get node corresponding to a different "type" sequence.
                    other_type_seq_node = None
                    for node in t2.traverse():
                        if node.name == other_type_seq_node_name:
                            other_type_seq_node = node
                            break
    
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
    
                # find the node with the most child leaf nodes.
                #node_num = 0
                #for node in nodes_of_interest:
                #    node_num += 1 
                #    print('Node ' + str(node_num) + ' Number of leaves:  ' + str(len(node.get_leaves())))
                #    print(node)
                node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
                        len(x.get_leaves()), reverse=True)[0]
                node_w_most_leaves.name = 'X'
                #print('\n\nClade defined by sequence ' + ts + ':')
                #print(node_w_most_leaves)
    
    
                # Copy tree and remove node/clades in constraint trees and replace with
                # polytomies composed of the same sequences.
                starting_tree = t3.copy()
    
                node_of_interest = None
                for node in starting_tree.traverse():
                    if node.name == 'X':
                        node_of_interest = node
                        break
    
                # Add node for additional sequence as most basal taxon.
                name_for_additional_node = 'AddSeq'
                node_of_interest.add_child(name=name_for_additional_node)
                node_of_interest.set_outgroup('AddSeq')
                print('node with additional sequence:')
                print(node_of_interest)
    
                # Unroot tree, because ML trees are all unrooted.
                starting_tree.unroot()
    
                # Write modified tree.
                outtreename1 = os.path.join(outdirpath1,
                        'clade_' + str(ts_num) + '_starting_tree_def_seq_' + ts.replace(' ', '_')\
                        + '.tre')
                starting_tree.write(outfile=outtreename1, format=9)
    
    
                # Copy tree and remove node/clades in constraint trees and replace with
                # polytomies composed of the same sequences.
                clade_constr_tree = starting_tree.copy()
    
                # Remove internal child nodes in the clade so that it becomes a
                # polytomy. # SEPARATE
                for node in clade_constr_tree.traverse():
                    if node.name == 'X':
                        nodes_to_delete = []
                        for cn in node.iter_descendants():
                            if not cn.is_leaf():
                                nodes_to_delete.append(cn)
                        for cn in nodes_to_delete:
                            cn.delete()
    
                # Show modified node.
                #print('\nClade as polytomy:')
                for node in clade_constr_tree.traverse():
                    if node.name == 'X':
                        #print(node)
                        pass
    
                # Unroot tree, because ML trees are all unrooted.
                clade_constr_tree.unroot()
                
                # Write modified tree.
                outtreename2 = os.path.join(outdirpath1,
                        'clade_' + str(ts_num) + '_constraint_tree_def_seq_' + ts.replace(' ', '_')\
                        + '.tre')
                #print('outtreename:')
                #print(outtreename)
                #with open(outtreename, 'w') as outfh:
                clade_constr_tree.write(outfile=outtreename2, format=9)
    
            # Remove temporary files.
            os.remove(outtreefp)
            os.remove(simple_tree)
    
    
    # Loop through alignments and make copies of the starting trees and constraint
    # trees with renamed "AddSeq" taxon, as well as copies with coded sequence
    # names for input to IQtree.
    ali_num = 0
    for ali in glob.glob(os.path.join(outdirpath1, '*trim.C.nex')):
        ali_num += 1
    
        # Remove unnecessary uncoded trim file.
        trim_file_to_remove = ali.rsplit('.', 2)[0] + '.nex'
        os.remove(trim_file_to_remove)
    
        # Define corresponding conversion table file path.
        tablefp = ali.rsplit('.', 1)[0] + '.table'
    
        # Get name of additional sequence from table.
        additional_seq_orig_name = None
        with open(tablefp) as tablefh:
            lines = tablefh.readlines()
            additional_seq_orig_name = lines[-1].strip()
    
        # Loop through trees, parse them, change the "AddSeq" node name to the name
        # of the additional sequence for this particular alignment, then write to a
        # file. # SEPARATE
        constraint_tree_list = []
        starting_tree_list = []
        print('\n\n\n')
        for tf in sorted(glob.glob(os.path.join(os.path.dirname(ali), 'clade_*.tre'))):
            print(tf)
            # Parse tree.
            tx = Tree(tf, quoted_node_names=True)
            # Change node name in tree.
            for node in tx.traverse():
                if node.name == 'AddSeq':
                    node.name = additional_seq_orig_name
                    break
            # Write modified tree.
            outtreename = ali.rsplit('.', 2)[0] + '_' + os.path.basename(tf)
            tx.write(outfile=outtreename, format=9)
            # Make copy of tree file with coded names.
            outtreename2 = outtreename.rsplit('.', 1)[0] + '.C.tre'
            code_tree(outtreename, outtreename2, tablefp)
            # Remove copy with uncoded names.
            os.remove(outtreename)
            # Add file paths to lists as appropriate.
            if os.path.basename(tf).split('_')[2] == 'constraint':
                constraint_tree_list.append(outtreename2)
            elif os.path.basename(tf).split('_')[2] == 'starting':
                starting_tree_list.append(outtreename2)
    
    
        # Make a copy of the alignment in phylip format for input to IQtree.
        phy_out = ali.rsplit('.', 1)[0] + '.phy'
        nex_to_phylip(ali, phy_out)
    
        # Remove nexus format version of alignment, because it is redundant with
        # the phylip version.
        os.remove(ali)
    
        # Print to see whether the tree files are in the right order.
        print('\n\n\n')
        for ctf, stf in zip(constraint_tree_list, starting_tree_list):
            print(ctf)
            print(stf)
            print('\n')
        print('\n\n\n')
    
        # Do phylogenetic analysis. # SEPARATE
        outtree_files_list = []
        for ctf, stf in zip(constraint_tree_list, starting_tree_list):
            # Make subdir for output.
            subdirp = ctf.rsplit('.', 1)[0] + '_IQ-tree_output'
            os.mkdir(subdirp)
    
            # Use IQtree to do an ML search
            output_file_prefix = os.path.join(subdirp, 'iqtree')
            subprocess.call(['iqtree', '-s', phy_out, '-m', subs_model, '-g',
                ctf, '-t', stf, '-pre', output_file_prefix, '-nt', '4'])
            # Add output tree file path to list.
            tree_file_path = output_file_prefix + '.treefile'
            outtree_files_list.append(tree_file_path)
            # Check that the output file was actually produced.
            assert os.path.isfile(tree_file_path)
    
        # Do topology tests using yielded trees.
    
        # Make output directory for topology test analysis.
        topo_test_subdirp = phy_out.rsplit('.', 1)[0] + '_topology_tests'
        os.mkdir(topo_test_subdirp)
        concat_tree_file = os.path.join(topo_test_subdirp, 'concat_trees.tre')
        concat_tree_name_file = os.path.join(topo_test_subdirp, 'concat_tree_names.txt')
        f_num = 0
        alt_tree_files = []
        for f in outtree_files_list:
            f_num += 1
            alt_tree_files.append(f)
            with open(f) as infh, open(concat_tree_file, 'a') as o1,\
            open(concat_tree_name_file, 'a') as o2:
                for i in infh:
                    o1.write(i)
                    o2.write(str(f_num) + ': ' + f)
    
        # Use IQ-tree to perform topology tests.
        topo_test_output_prefix = concat_tree_file.rsplit('.', 1)[0] + '_topo_test' 
        topo_test_output_fp = concat_tree_file.rsplit('.', 1)[0] + '_topo_test.iqtree' 
        subprocess.call(['iqtree', '-s', phy_out, '-m', subs_model, '-n', '0',
            '-z', concat_tree_file, '-zb', '10000', '-au', '-pre',
            topo_test_output_prefix])
    
        # Parse IQ-tree output to get results of topology test. # SEPARATE
        data_line = re.compile(r'\d+ +-\d+\.\d+ +\d+\.\d+ +')
        space_char = re.compile(r' +')
        topo_test_res_dict = {}
        with open(topo_test_output_fp) as infh:
            for i in infh:
                if data_line.search(i.strip()):
                    # Parse info from line in IQ-tree output file.
                    parsed_tree_info_list = space_char.split(i.strip())
                    tree_num = int(parsed_tree_info_list[0])
                    original_tree_fp = alt_tree_files[tree_num -1]
                    logL = float(parsed_tree_info_list[1])
                    bp_RELL = float(parsed_tree_info_list[3])
                    p_KH = float(parsed_tree_info_list[5])
                    p_SH = float(parsed_tree_info_list[7])
                    c_ELW = float(parsed_tree_info_list[9])
                    p_AU = float(parsed_tree_info_list[11])
                    #accept_reject = parsed_tree_info_list[12]
                    # Add info to dict.
                    topo_test_res_dict[tree_num] = (original_tree_fp, logL,\
                            bp_RELL, p_KH, p_SH, c_ELW, p_AU)
    
        # Write info to summary spreadsheet file. # SEPARATE
        with open(summary_csv_fp, 'a') as o:
            num_keys = str(len(topo_test_res_dict.keys()))
            key_num = 0
            for key in sorted(topo_test_res_dict.keys(), key=lambda x:\
                    topo_test_res_dict[x][1], reverse=True):
                key_num += 1
                info_tuple = topo_test_res_dict[key]
                # Assemble info to write to spreadsheet.
                relative_tree_path = info_tuple[0].replace(os.path.dirname(os.path.dirname(info_tuple[0])), '')
                orig_align_file = os.path.basename(alignment)
                orig_tree_file = os.path.basename(tree)
                logL = str(info_tuple[1])
                bp_RELL = str(info_tuple[2])
                p_KH = str(info_tuple[3])
                p_SH = str(info_tuple[4])
                c_ELW = str(info_tuple[5])
                p_AU = str(info_tuple[6])
    
                # Apply a consensus approach to determine whether to accept or
                # reject topologies.
                accept_reject = '+'
                consensus_num = 0
                if float(bp_RELL) <= 0.05:
                    consensus_num += 1
                if float(p_KH) <= 0.05:
                    consensus_num += 1
                if float(p_SH) <= 0.05:
                    consensus_num += 1
                if float(c_ELW) <= 0.05:
                    consensus_num += 1
                if float(p_AU) <= 0.05:
                    consensus_num += 1
                if consensus_num >= 4:
                    accept_reject = '-'
    
                hypothesis_num = str(key_num) + ' of ' + num_keys
                relative_table_path = tablefp.replace(os.path.dirname(os.path.dirname(tablefp)), '')
                type_seq_name =\
                os.path.dirname(info_tuple[0]).split('_def_seq_')[1].split('.C_')[0]
    
                masked_ali_file = os.path.join(outdirpath1,'*mask.C.nex').replace(os.path.dirname(outdirpath1), '')
    
                # Write info line to spreadsheet.
                info_line = ','.join([orig_align_file, #0
                                      orig_tree_file, #1
                                      subs_model, #2
                                      additional_seq_orig_name, #3
                                      hypothesis_num, #4
                                      type_seq_name, #5
                                      'Clade ?', #6
                                      logL,  #7
                                      bp_RELL, #8
                                      p_KH, #8A
                                      p_SH, #8B
                                      c_ELW, #8C
                                      p_AU, #8D
                                      accept_reject,
                                      masked_ali_file, #10
                                      relative_table_path, #11
                                      'Trimmed input alignment file', #12
                                      'Constraint tree file', #13
                                      'Starting tree file', #14
                                      relative_tree_path, #15
                                      '\n'
                                      ])
                o.write(info_line)

def get_phylo_class_csv(outdirpath):
    """Given a phylo_class output directory, return the path to the full
    summary file contained in the directory.
    """
    #summary_csv_fp = os.path.join(outdirpath1, '0_' +\
    #        os.path.basename(file_with_seqs).rsplit('.',\
    #    1)[0] + '_Summary.csv')
    summary_csv_fp = os.path.join(outdirpath, '0_full_summary.csv')
    return summary_csv_fp


def classify_seq_with_constrained_tree2(file_with_seqs, model=None,
        subseq=False, abbrev=True, place=False):
    """Classify sequences in a given fasta file using a given backbone tree.
    """
    # Start time.
    start_time = time.time()

    # Define timestamp.
    cur_time = time.strftime("%Y%m%d%H%M%S")
    print('\nTimestamp: ' + cur_time)
    
    # Make output directory.
    #outdirpath1 = alignment.rsplit('.', 1)[0] + '_classify_seq_' + cur_time
    outdirpath1 = file_with_seqs.rsplit('.', 1)[0] + '_classify_seq_' + cur_time
    os.mkdir(outdirpath1)
    
    ## Take a specific approach, depending on the options specified at runtime.
    #
    ## If bootstrapping option not specified, then use default topology-testing
    ## approach.
    #if not args.boot:
    #
    #    # Call function for doing default analysis.
    
    # Write initial summary csv file.
    summary_csv_fp = get_phylo_class_csv(outdirpath1) 

    # Define summary header line string.
    summary_header_line = ''
    if not place:
        summary_header_line = ','.join(['Model/backbone tree name',
                                        'Original alignment file', #0
                                        'Original tree file', #1
                                        'Substitution model', #2
                                        'Sequence added to tree', #3
                                        'Hypothesis (tree/clade) number', #4
                                        'Type sequence of clade the additional sequence was constrained into', #5
                                        'Clade name', #6
                                        'Log-likelihood of ML tree',  #7
                                        'bp-RELL: bootstrap proportion using RELL method (Kishino et al. 1990)',  #8
                                        'p-KH: p-value of one sided Kishino-Hasegawa test (1989)',  #8
                                        'p-SH: p-value of Shimodaira-Hasegawa test (2000)',  #8
                                        'c-ELW: Expected Likelihood Weight (Strimmer & Rambaut 2002)',  #8
                                        '"p-AU: p-value of approximately unbiased (AU) test (Shimodaira, 2002)"',  #8
                                        'Accept (+) or reject (-) hypothesis based on consensus of topology tests',  #9
                                        'Masked alignment file', #10
                                        'Sequence name conversion table file', #11
                                        'Trimmed input alignment file', #12
                                        'Constraint tree file', #13
                                        'Starting tree file', #14
                                        'ML tree file', #15
                                        '\n'
                                        ])
    elif place:
        summary_header_line = ','.join(['Model/backbone tree name',
                                        'Original alignment file', #0
                                        'Original tree file', #1
                                        'Substitution model', #2
                                        'ID of sequence added to tree', #3
                                        'Description of sequence added to tree', #3
                                        'Type sequence of clade the additional sequence was placed in by ML search', #5
                                        'Name of clade that additional sequence was placed in', #6
                                        'Masked alignment file', #10
                                        'Sequence name conversion table file', #11
                                        'Trimmed input alignment file', #12
                                        'Constraint tree file', #13
                                        'ML tree file', #15
                                        '\n'
                                        ])
    
    with open(summary_csv_fp, 'w') as o:
        o.write(summary_header_line)

    # Take tree topology and alignment (regular names), convert names in alignment
    # generating a conversion table similar to for name_replace.pl, then use that
    # conversion table to code the names in the tree topology as well.

    # Determine whether input sequence file is a fasta file or not.
    fasta = True
    if file_with_seqs.endswith('.csv'):
        fasta = False
    
    # Align sequences to alignment and trim.
    if fasta:
        assert model is not None, """A model must be specified if the input
        sequence file is a FASTA file."""
        ## Get info.
        #model_info = ModelInfoFromCSV(model)

        # Classify each sequence.
        ali_num = 0
        for record in SeqIO.parse(file_with_seqs, 'fasta'):
            ali_num += 1
            if not place:
                classify_one_seq_record(ali_num, record, model,
                    outdirpath1, summary_csv_fp, summary_header_line)
            else:
                place_one_seq_record(ali_num, record, model,
                    outdirpath1, summary_csv_fp, summary_header_line)
    else:
        # Parse csv file.
        df = pd.read_csv(file_with_seqs)

        # Get column header for column with Yes/No data.
        yes_no_col = None
        for header in df.columns:
            if header.startswith('Represents an identifiably unique paralogue'):
                yes_no_col = header
                break
        if yes_no_col is None:
            for header in df.columns:
                if header.startswith('Collective interpretation of reverse search results'):
                    yes_no_col = header
                    break
        if yes_no_col is None:
            for header in df.columns:
                if header.startswith('Positive/redundant (+) or negative (-) hit based on E-value criterion'):
                    yes_no_col = header
                    break

        # Iterate over rows in dataframe.
        ali_num = 0
        ali_used_num = 0
        for index, row in df.iterrows():
            ali_num += 1

            # Get a sequence record object and model from info in the row.
            record = None
            model = None
            result = row[yes_no_col]
            if result == 'Yes' or result == '+': 
                ali_used_num += 1
                acc = row['Forward hit accession']
                description = None
                seq = None
                if subseq or row['Forward search method'].startswith('tblastn'):
                    description = row['Forward hit description of subsequence(s) that align(s) to query']
                    seq = row['Forward hit subsequence(s) that align(s) to query'] 
                else:
                    description = row['Forward hit description']
                    seq = row['Forward hit sequence'] 
                db_file = row['Subject database species (if applicable)']

                # Get name of model/tree to use from the Models directory
                # (specified in the settings module).
                model = row['Model/backbone tree name']
                assert model != '-'
                assert model is not None

                # Use the species name if writing the file for input to
                # phylogenetics software.
                if abbrev:
                    description = db_file.rsplit('.', 1)[0]
                
                # Instantiate a sequence object for the sequence of interest to
                # be written to output.
                record = get_seq_obj_from_srch_res_csv_info(acc, description, seq)

                ## Call function to retrieve info about model tree to use.
                #model_info = ModelInfoFromCSV(model)

                # Call function to classify the sequence and add corresponding rows
                # to the (intermediate) summary spreadsheet.
                if not place:
                    classify_one_seq_record(ali_used_num, record, model, 
                        outdirpath1, summary_csv_fp, summary_header_line)
                else:
                    place_one_seq_record(ali_num, record, model, 
                        outdirpath1, summary_csv_fp, summary_header_line)

    # End time.
    end_time = time.time()
    # Record time elapsed.
    elapsed = end_time - start_time
    logfile = get_phylo_class_logfile_path(outdirpath1)
    with open(logfile, 'a') as logh:
        logh.write('Total run time: ' + str(datetime.timedelta(seconds=elapsed)) + '\n')

    # Return main output path.
    return summary_csv_fp


class ModelInfoFromCSV:
    """Class for retrieving information about a named model (backbone tree).
    """
    def __init__(self, model_name):
        self.model_name = model_name

        # Get info from model info csv file using pandas.
        alignmentfp = None
        subs_model = None
        tree_topologyfp = None
        type_seqsfp = None
        df = pd.read_csv(settings.model_info_csv)
        found_model_with_name = False
        for index, row in df.iterrows():
            if row['Model name'] == model_name:
                found_model_with_name = True
                alignmentfp = os.path.join(settings.model_dir_path,
                        row['Alignment file'])
                subs_model = row['Substitution model']
                type_seqsfp = os.path.join(settings.model_dir_path, row['Type sequences and names of the clades that they define'])
                tree_topologyfp = os.path.join(settings.model_dir_path,
                        row['Tree topology file'])
                break
        assert found_model_with_name, """Could not find a model named %s in
        info spreadsheet %s.""" % (model_name, settings.model_info_csv)
        assert alignmentfp is not None and alignmentfp != ''
        assert subs_model is not None and subs_model != ''
        assert tree_topologyfp is not None and tree_topologyfp != ''
        assert type_seqsfp is not None and type_seqsfp != ''

        # Store this information as attributes of the object.
        self.alignment_file = alignmentfp
        self.subs_model = subs_model
        self.tree_topology_file = tree_topologyfp
        self.type_seqs_file = type_seqsfp


def get_clade_name_from_model(type_seq_name, type_seqs):
    """Return clade name for given type sequence name and type seqs info file
    (from models directory).
    """
    # Get clade name.
    clade_name = '?'
    with open(type_seqs) as infh:
        for i in infh:
            spliti = i.strip().split(',')
            n = spliti[0]
            c = spliti[1]
            if n == type_seq_name:
                clade_name = c
                break
            elif n.replace(' ', '_') == type_seq_name:
                clade_name = c
                break
    assert clade_name != '?', """Could not identify clade name for type
    sequence name %s.""" % type_seq_name

    return clade_name
        

def get_phylo_class_logfile_path(phylo_class_outdir):
    """Return a log file path given an output directory.
    """
    logfile = os.path.join(phylo_class_outdir, '0_tree_search_log.txt')
    return logfile


def classify_one_seq_record(ali_num, record, model,
        outdirpath1, summary_csv_fp, summary_header_line):
    """Classify one sequence.
    """
    #assert 1 == 2, """Trying not to use this function.""" # Temp.

    # Check that the input sequence record is not just '-'.
    assert len(record) > 1, """Input sequence is only one residue long."""

    # Get info about model that is relevant.
    model_info = ModelInfoFromCSV(model)
    alignment = model_info.alignment_file
    tree = model_info.tree_topology_file
    subs_model = model_info.subs_model
    type_seqs = model_info.type_seqs_file

    # Define logfile path to write to.
    logfile = get_phylo_class_logfile_path(outdirpath1)

    #ali = os.path.join(outdirpath1, '1.nex') 
    #ali = add_seq_to_alignment3(record, alignment, outdirpath1)
    ali = os.path.join(outdirpath1, str(ali_num) + '.nex') 
    add_seq_to_alignment3(record, alignment, ali)

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


    # Only do the following steps once for all the alignments, because all the
    # constrained topology analyses will use the same constraint trees.
    if ali_num == 1:

        # Put quotation marks around names in constraint tree. # SEPARATE
        outtreefp = ali.rsplit('.', 1)[0] + '.constraint_tree0.Q.newick'
        with open(tree) as intreefh, open(outtreefp, 'w') as o:
            for t in intreefh:
                o.write(quote_tree(t, outtablefp))



        # Make topology constraint trees.   

        # Parse backbone topology and root on a clade rather than having the root
        # inside a clade.

        ######################
        # Parse tree using ete3.
        # Note: parentheses and commas get replaced with underscores.
        t1 = Tree(outtreefp, quoted_node_names=True)

        # Define name for tree without branch lengths.
        simple_tree = outtreefp.rsplit('_', 1)[0] + '_TEMP2' 

        # Write simple tree to a new file 
        t1.write(format=9, outfile=simple_tree)

        # Parse simple tree.
        t2 = Tree(simple_tree)

        # Print simple tree.
        print('Input tree:')
        print(t2)

        ## Get node for a key clade. 
        ## Or, just require to root it before inputting (figtree can export
        ## re-rooted newick)?????
        #print("\nPlease define an arbitrary clade of interest for rooting the tree topology\n")
        #clade_member_1 = input('Clade member 1: ')
        #clade_member_2 = input('Clade member 2: ')
        #ancestor1 = t2.get_common_ancestor(clade_member_1, clade_member_2)
        ##print(ancestor1)

        ## Re-root tree on a clade (as a polytomy) rather than inside a clade, so
        ## that key branches can be more easily identified.
        #t2.set_outgroup(ancestor1)

        ######################


        # Get to topologies sorted out...

        # Get constraint topologies for each clade of interest to test
        # alternative topology hypotheses.

        # Get list of "type" sequences from input.
        type_seq_list = []
        for i in open(type_seqs):
            type_seq_list.append(i.strip().split(',')[0])

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

                # Get node corresponding to a different "type" sequence.
                other_type_seq_node = None
                #print(other_type_seq_node_name)
                for node in t2.traverse():
                    #print(node.name)
                    if node.name == other_type_seq_node_name:
                        other_type_seq_node = node
                        break
                assert other_type_seq_node is not None

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

            # find the node with the most child leaf nodes.
            #node_num = 0
            #for node in nodes_of_interest:
            #    node_num += 1 
            #    print('Node ' + str(node_num) + ' Number of leaves:  ' + str(len(node.get_leaves())))
            #    print(node)
            node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
                    len(x.get_leaves()), reverse=True)[0]
            node_w_most_leaves.name = 'X'
            #print('\n\nClade defined by sequence ' + ts + ':')
            #print(node_w_most_leaves)


            # Copy tree and remove node/clades in constraint trees and replace with
            # polytomies composed of the same sequences.
            starting_tree = t3.copy()

            node_of_interest = None
            for node in starting_tree.traverse():
                if node.name == 'X':
                    node_of_interest = node
                    break

            # Add node for additional sequence as most basal taxon.
            name_for_additional_node = 'AddSeq'
            node_of_interest.add_child(name=name_for_additional_node)
            node_of_interest.set_outgroup('AddSeq')
            print('node with additional sequence:')
            print(node_of_interest)

            # Unroot tree, because ML trees are all unrooted.
            starting_tree.unroot()

            # Write modified tree.
            outtreename1 = os.path.join(outdirpath1,
                    'clade_' + str(ts_num) + '_starting_tree_def_seq_' + ts.replace(' ', '_')\
                    + '.tre')
            starting_tree.write(outfile=outtreename1, format=9)


            # Copy tree and remove node/clades in constraint trees and replace with
            # polytomies composed of the same sequences.
            clade_constr_tree = starting_tree.copy()

            # Remove internal child nodes in the clade so that it becomes a
            # polytomy. # SEPARATE
            for node in clade_constr_tree.traverse():
                if node.name == 'X':
                    nodes_to_delete = []
                    for cn in node.iter_descendants():
                        if not cn.is_leaf():
                            nodes_to_delete.append(cn)
                    for cn in nodes_to_delete:
                        cn.delete()

            # Show modified node.
            #print('\nClade as polytomy:')
            for node in clade_constr_tree.traverse():
                if node.name == 'X':
                    #print(node)
                    pass

            # Unroot tree, because ML trees are all unrooted.
            clade_constr_tree.unroot()
            
            # Write modified tree.
            outtreename2 = os.path.join(outdirpath1,
                    'clade_' + str(ts_num) + '_constraint_tree_def_seq_' + ts.replace(' ', '_')\
                    + '.tre')
            #print('outtreename:')
            #print(outtreename)
            #with open(outtreename, 'w') as outfh:
            clade_constr_tree.write(outfile=outtreename2, format=9)

        # Remove temporary files.
        os.remove(outtreefp)
        os.remove(simple_tree)


    ## Loop through alignments and make copies of the starting trees and constraint
    ## trees with renamed "AddSeq" taxon, as well as copies with coded sequence
    ## names for input to IQtree.
    #ali_num = 0
    #for ali in glob.glob(os.path.join(outdirpath1, '*trim.C.nex')):
    #    ali_num += 1

    # Remove unnecessary uncoded trim file.
    trim_file_to_remove = outalifpnex.rsplit('.', 2)[0] + '.nex'
    os.remove(trim_file_to_remove)

    # Define corresponding conversion table file path.
    tablefp = outalifpnex.rsplit('.', 1)[0] + '.table'

    # Get name of additional sequence from table.
    additional_seq_orig_name = None
    with open(tablefp) as tablefh:
        lines = tablefh.readlines()
        additional_seq_orig_name = lines[-1].strip()

    # Check that the output folder is a directory.
    assert os.path.isdir(os.path.dirname(outalifpnex)), """Directory does not exist: %s""" % os.path.dirname(outalifpnex)

    # Loop through trees, parse them, change the "AddSeq" node name to the name
    # of the additional sequence for this particular alignment, then write to a
    # file. # SEPARATE
    constraint_tree_list = []
    starting_tree_list = []
    print('\n\n\n')
    print(os.path.dirname(outalifpnex))
    trees_list = sorted(glob.glob(os.path.join(os.path.dirname(outalifpnex), 'clade_*.tre')))
    assert len(trees_list) > 0, """No constraint or starting trees identified."""
    for tf in trees_list: 
        print(tf)
        # Parse tree.
        tx = Tree(tf, quoted_node_names=True)
        # Change node name in tree.
        for node in tx.traverse():
            if node.name == 'AddSeq':
                node.name = additional_seq_orig_name
                break
        # Write modified tree.
        outtreename = outalifpnex.rsplit('.', 2)[0] + '_' + os.path.basename(tf)
        tx.write(outfile=outtreename, format=9)
        # Make copy of tree file with coded names.
        outtreename2 = outtreename.rsplit('.', 1)[0] + '.C.tre'
        code_tree(outtreename, outtreename2, tablefp)
        # Remove copy with uncoded names.
        os.remove(outtreename)
        # Add file paths to lists as appropriate.
        if os.path.basename(tf).split('_')[2] == 'constraint':
            constraint_tree_list.append(outtreename2)
        elif os.path.basename(tf).split('_')[2] == 'starting':
            starting_tree_list.append(outtreename2)

    # Check that constraint trees were found.
    assert len(constraint_tree_list) > 0, """No constraint trees in list."""
    assert len(starting_tree_list) > 0, """No starting trees in list."""


    # Make a copy of the alignment in phylip format for input to IQtree.
    phy_out = outalifpnex.rsplit('.', 1)[0] + '.phy'
    nex_to_phylip(outalifpnex, phy_out)

    # Remove nexus format version of alignment, because it is redundant with
    # the phylip version.
    os.remove(outalifpnex)

    # Print to see whether the tree files are in the right order.
    print('\n\n\n')
    for ctf, stf in zip(constraint_tree_list, starting_tree_list):
        print('ctf')
        print(ctf)
        # Check that the constraint tree is a file.
        assert os.path.isfile(ctf), """Could not find constraint tree file: %s""" % ctf
        print('stf')
        print(stf)
        # Check that the starting tree is a file.
        assert os.path.isfile(ctf), """Could not find starting tree file: %s""" % stf
        print('\n')
    print('\n\n\n')


    # Do phylogenetic analysis. # SEPARATE
    outtree_files_list = []
    for ctf, stf in zip(constraint_tree_list, starting_tree_list):
        # Make subdir for output.
        subdirp = ctf.rsplit('.', 1)[0] + '_IQ-tree_output'
        os.mkdir(subdirp)

        # Use IQtree to do an ML search
        output_file_prefix = os.path.join(subdirp, 'iqtree')
        iqtree_command_list = ['iqtree', '-s', phy_out, '-m', subs_model, '-g',
            ctf, '-t', stf, '-pre', output_file_prefix, '-nt', '4']
        tree_search_start_time = time.time()
        subprocess.call(iqtree_command_list)
        tree_search_end_time = time.time()
        with open(logfile, 'a') as logh:
            logh.write(' '.join(iqtree_command_list) + '\n')
            tree_search_elapsed = tree_search_end_time - tree_search_start_time
            logh.write('Run time: ' + str(datetime.timedelta(seconds=tree_search_elapsed)) + '\n')
        # Add output tree file path to list.
        tree_file_path = output_file_prefix + '.treefile'
        outtree_files_list.append(tree_file_path)
        # Check that the output file was actually produced.
        assert os.path.isfile(tree_file_path)

    # Do topology tests using yielded trees.

    # Make output directory for topology test analysis.
    topo_test_subdirp = phy_out.rsplit('.', 1)[0] + '_topology_tests'
    os.mkdir(topo_test_subdirp)
    concat_tree_file = os.path.join(topo_test_subdirp, 'concat_trees.tre')
    concat_tree_name_file = os.path.join(topo_test_subdirp, 'concat_tree_names.txt')
    f_num = 0
    alt_tree_files = []
    for f in outtree_files_list:
        f_num += 1
        alt_tree_files.append(f)
        with open(f) as infh, open(concat_tree_file, 'a') as o1,\
        open(concat_tree_name_file, 'a') as o2:
            for i in infh:
                o1.write(i)
                o2.write(str(f_num) + ': ' + f)

    # Use IQ-tree to perform topology tests.
    topo_test_output_prefix = concat_tree_file.rsplit('.', 1)[0] + '_topo_test' 
    topo_test_output_fp = concat_tree_file.rsplit('.', 1)[0] + '_topo_test.iqtree' 
    subprocess.call(['iqtree', '-s', phy_out, '-m', subs_model, '-n', '0',
        '-z', concat_tree_file, '-zb', '10000', '-au', '-pre',
        topo_test_output_prefix])

    # Parse IQ-tree output to get results of topology test. # SEPARATE
    data_line = re.compile(r'\d+ +-\d+\.\d+ +\d+\.\d+ +')
    space_char = re.compile(r' +')
    topo_test_res_dict = {}
    with open(topo_test_output_fp) as infh:
        for i in infh:
            if data_line.search(i.strip()):
                # Parse info from line in IQ-tree output file.
                parsed_tree_info_list = space_char.split(i.strip())
                tree_num = int(parsed_tree_info_list[0])
                original_tree_fp = alt_tree_files[tree_num -1]
                logL = float(parsed_tree_info_list[1])
                bp_RELL = float(parsed_tree_info_list[3])
                p_KH = float(parsed_tree_info_list[5])
                p_SH = float(parsed_tree_info_list[7])
                c_ELW = float(parsed_tree_info_list[9])
                p_AU = float(parsed_tree_info_list[11])
                #accept_reject = parsed_tree_info_list[12]
                # Add info to dict.
                topo_test_res_dict[tree_num] = (original_tree_fp, logL,\
                        bp_RELL, p_KH, p_SH, c_ELW, p_AU)

    # Write info to summary spreadsheet file. # SEPARATE
    with open(summary_csv_fp, 'a') as o:
        num_keys = str(len(topo_test_res_dict.keys()))
        key_num = 0
        for key in sorted(topo_test_res_dict.keys(), key=lambda x:\
                topo_test_res_dict[x][1], reverse=True):
            key_num += 1
            info_tuple = topo_test_res_dict[key]
            # Assemble info to write to spreadsheet.
            relative_tree_path = info_tuple[0].replace(os.path.dirname(os.path.dirname(info_tuple[0])), '')
            orig_align_file = os.path.basename(alignment)
            orig_tree_file = os.path.basename(tree)
            logL = str(info_tuple[1])
            bp_RELL = str(info_tuple[2])
            p_KH = str(info_tuple[3])
            p_SH = str(info_tuple[4])
            c_ELW = str(info_tuple[5])
            p_AU = str(info_tuple[6])

            # Apply a consensus approach to determine whether to accept or
            # reject topologies.
            accept_reject = '+'
            consensus_num = 0
            if float(bp_RELL) <= 0.05:
                consensus_num += 1
            if float(p_KH) <= 0.05:
                consensus_num += 1
            if float(p_SH) <= 0.05:
                consensus_num += 1
            if float(c_ELW) <= 0.05:
                consensus_num += 1
            if float(p_AU) <= 0.05:
                consensus_num += 1
            if consensus_num >= 4:
                accept_reject = '-'

            hypothesis_num = str(key_num) + ' of ' + num_keys
            relative_table_path = tablefp.replace(os.path.dirname(os.path.dirname(tablefp)), '')
            type_seq_name =\
            os.path.dirname(info_tuple[0]).split('_def_seq_')[1].split('.C_')[0]

            # Get clade name.
            clade_name = get_clade_name_from_model(type_seq_name, type_seqs)

            masked_ali_file = os.path.join(outdirpath1,'*mask.C.nex').replace(os.path.dirname(outdirpath1), '')

            # Write info line to spreadsheet.
            info_line = ','.join([model,
                                  orig_align_file, #0
                                  orig_tree_file, #1
                                  subs_model, #2
                                  additional_seq_orig_name, #3
                                  hypothesis_num, #4
                                  type_seq_name, #5
                                  clade_name, #6
                                  logL,  #7
                                  bp_RELL, #8
                                  p_KH, #8A
                                  p_SH, #8B
                                  c_ELW, #8C
                                  p_AU, #8D
                                  accept_reject,
                                  masked_ali_file, #10
                                  relative_table_path, #11
                                  'Trimmed input alignment file', #12
                                  'Constraint tree file', #13
                                  'Starting tree file', #14
                                  relative_tree_path, #15
                                  '\n'
                                  ])
            o.write(info_line)


def place_one_seq_record(ali_num, record, model,
        outdirpath1, summary_csv_fp, summary_header_line):
    """Classify one sequence.
    """
    # Get info about model that is relevant.
    model_info = ModelInfoFromCSV(model)
    alignment = model_info.alignment_file
    tree = model_info.tree_topology_file
    subs_model = model_info.subs_model
    type_seqs = model_info.type_seqs_file

    # Define logfile path to write to.
    logfile = get_phylo_class_logfile_path(outdirpath1)

    # Modify sequence header for easy parsing.
    additional_seq_description = record.description
    modify_seq_descr_for_tree(record)

    # Add sequence to alignment.
    #ali = os.path.join(outdirpath1, '1.nex') 
    ali = os.path.join(outdirpath1, os.path.basename(alignment).rsplit('.', 1)[0] + '_' + str(record.id) + '.nex')
    masked_ali_file = ali.rsplit('.', 1)[0] + '_mask.nex'
    add_seq_to_alignment3(record, alignment, ali)

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

    # Define path for constraint tree file.
    constraint_tree_fp_coded = None
    # Only do the following steps once for all the alignments, because all the
    # constrained topology analyses will use the same constraint trees.
    if ali_num == 1:

        # Put quotation marks around names in constraint tree. # SEPARATE
        #outtreefp = ali.rsplit('.', 1)[0] + '.constraint_tree0.Q.newick'
        constraint_tree_fp_rooted = os.path.join(os.path.dirname(ali),\
                os.path.basename(tree.rsplit('.')[0] + '_constraint_tree_rooted.newick'))
        constraint_tree_fp = os.path.join(os.path.dirname(ali),\
                os.path.basename(tree.rsplit('.')[0] + '_constraint_tree.newick'))
        with open(tree) as intreefh, open(constraint_tree_fp_rooted, 'w') as o:
            for t in intreefh:
                o.write(quote_tree(t, outtablefp))
        # Unroot constraint tree.
        tobj = Tree(constraint_tree_fp_rooted, quoted_node_names=True)
        tobj.unroot()
        tobj.write(outfile=constraint_tree_fp, format=9)
        assert os.path.isfile(constraint_tree_fp)
        
        # Make constraint tree with coded names.
        constraint_tree_fp_coded = constraint_tree_fp.rsplit('.', 1)[0] + '.C.newick'
        code_tree(constraint_tree_fp, constraint_tree_fp_coded, outtablefp)
        assert os.path.isfile(constraint_tree_fp_coded)
        #print(constraint_tree_fp_coded)
        #assert 2 != 2

    else:
        constraint_tree_fp_coded = glob.glob(os.path.join(os.path.dirname(ali),\
        '*_constraint_tree.C.newick'))[0]

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
    #subdirp = constraint_tree_fp_coded.rsplit('.', 1)[0] + '_IQ-tree_output'
    subdirp = outalifpnex.rsplit('.', 1)[0] + '_IQ-tree_output'
    os.mkdir(subdirp)

    # Use IQtree to do an ML search
    output_file_prefix = os.path.join(subdirp, 'iqtree')
    iqtree_command_list = ['iqtree', '-s', phy_out, '-m', subs_model, '-g',
        constraint_tree_fp_coded, '-pre', output_file_prefix, '-nt', '4']
    tree_search_start_time = time.time()
    subprocess.call(iqtree_command_list)
    tree_search_end_time = time.time()
    with open(logfile, 'a') as logh:
        logh.write(' '.join(iqtree_command_list) + '\n')
        tree_search_elapsed = tree_search_end_time - tree_search_start_time
        logh.write('Run time: ' + str(datetime.timedelta(seconds=tree_search_elapsed)) + '\n')
    # Add output tree file path to list.
    tree_file_path = output_file_prefix + '.treefile'
    # Check that the output file was actually produced.
    assert os.path.isfile(tree_file_path)


    # Parse the resulting tree to determine in which clade of interest, if any,
    # the sequence was placed.

    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(tree_file_path, quoted_node_names=True)

    # Convert names in tree back to original names.
    uncode_tree_obj(t1, outtablefp)

    ## Define name for tree without branch lengths.
    #simple_tree = tree_file_path.rsplit('_', 1)[0] + '_TEMP2' 

    ## Write simple tree to a new file 
    #t1.write(format=9, outfile=simple_tree)

    ## Parse simple tree.
    #t2 = Tree(simple_tree)

    # Print simple tree.
    print('ML tree:')
    print(t1)

    # Make a copy of the tree object.
    t2 = t1.copy()

    # Get list of "type" sequences from input.
    type_seq_list = []
    for i in open(type_seqs):
        type_seq_list.append(i.strip().split(',')[0])

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
            for node in t2.traverse():
                if node.name == other_type_seq_node_name:
                    other_type_seq_node = node
                    break

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

        # find the node with the most child leaf nodes.
        #node_num = 0
        #for node in nodes_of_interest:
        #    node_num += 1 
        #    print('Node ' + str(node_num) + ' Number of leaves:  ' + str(len(node.get_leaves())))
        #    print(node)
        node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
                len(x.get_leaves()), reverse=True)[0]
        node_w_most_leaves.name = 'X'
        print('\n\nClade defined by sequence ' + ts + ':')
        print(node_w_most_leaves)

        #name_of_additional_seq = record.description
        name_of_additional_seq = record.id # + ' ' + record.description

        #print('\n\n\n')
        #print(name_of_additional_seq)
        #print('\n\n\n')
        #for i in [x.name for x in node_w_most_leaves.get_leaves()]:
        #    print(i)
        #print('\n\n\n')

        if name_of_additional_seq.strip() in [x.name for x in node_w_most_leaves.get_leaves()]:
            ts_that_additional_seq_was_placed_in = ts
            print('%s is in this clade.' % name_of_additional_seq)
            break

    # If the sequence could not be placed in any of the clades of interest,
    # then set the name of the type sequence to a string 'None'.
    if ts_that_additional_seq_was_placed_in is None:
        ts_that_additional_seq_was_placed_in = 'None'

    # Check that the clade that the additional sequence was placed in was
    # identified.
    assert ts_that_additional_seq_was_placed_in is not None, """Sequence was
    not placed in any of the clades of interest: %s""" % record.description


    # Write info to summary spreadsheet file. # SEPARATE
    orig_align_file = os.path.basename(alignment) 
    orig_tree_file = os.path.basename(tree) 
    #additional_seq_orig_name = name_of_additional_seq 
    additional_seq_id = name_of_additional_seq.split('__')[0]
    type_seq_name = ts_that_additional_seq_was_placed_in
    clade_name = ts_that_additional_seq_was_placed_in
    if ts_that_additional_seq_was_placed_in is not 'None':
        clade_name = get_clade_name_from_model(type_seq_name, type_seqs)

    masked_ali_file = os.path.basename(masked_ali_file)
    relative_table_path = ''
    relative_tree_path = ''
    with open(summary_csv_fp, 'a') as o:
        # Write info line to spreadsheet.
        info_line = ','.join([model,
                              orig_align_file,
                              orig_tree_file,
                              subs_model,
                              additional_seq_id,
                              additional_seq_description,
                              type_seq_name,
                              clade_name,
                              masked_ali_file,
                              relative_table_path,
                              'Trimmed input alignment file',
                              'Constraint tree file',
                              'Starting tree file',
                              relative_tree_path,
                              '\n'
                              ])
        o.write(info_line)










    ######################### OLD STUFF:

        # Parse backbone topology and root on a clade rather than having the root
        # inside a clade.

    #    ######################
    #    # Parse tree using ete3.
    #    # Note: parentheses and commas get replaced with underscores.
    #    t1 = Tree(outtreefp, quoted_node_names=True)

    #    # Define name for tree without branch lengths.
    #    simple_tree = outtreefp.rsplit('_', 1)[0] + '_TEMP2' 

    #    # Write simple tree to a new file 
    #    t1.write(format=9, outfile=simple_tree)

    #    # Parse simple tree.
    #    t2 = Tree(simple_tree)

    #    # Print simple tree.
    #    print('Input tree:')
    #    print(t2)

    #    ## Get node for a key clade. 
    #    ## Or, just require to root it before inputting (figtree can export
    #    ## re-rooted newick)?????
    #    #print("\nPlease define an arbitrary clade of interest for rooting the tree topology\n")
    #    #clade_member_1 = input('Clade member 1: ')
    #    #clade_member_2 = input('Clade member 2: ')
    #    #ancestor1 = t2.get_common_ancestor(clade_member_1, clade_member_2)
    #    ##print(ancestor1)

    #    ## Re-root tree on a clade (as a polytomy) rather than inside a clade, so
    #    ## that key branches can be more easily identified.
    #    #t2.set_outgroup(ancestor1)

    #    ######################


    #    # Get to topologies sorted out...

    #    # Get constraint topologies for each clade of interest to test
    #    # alternative topology hypotheses.

    #    # Get list of "type" sequences from input.
    #    type_seq_list = []
    #    for i in open(type_seqs):
    #        type_seq_list.append(i.strip().split(',')[0])

    #    # For each "type" sequence, traverse all nodes and find the node with
    #    # the largest number of child nodes that are leaf (terminal) nodes,
    #    # containing the "type" sequence of interest, but not containing any of
    #    # the other "type" sequences.
    #    #ts_num = 0
    #    first_type_seq_node_name = None
    #    #for ts in type_seq_list:
    #    #    ts_num += 1

    #    if ts_num == 1:
    #        first_type_seq_node_name = ts
    #        # Root on another "type" sequence for the first type sequence in
    #        # the list to get whole clade, then root the tree on the ancestor
    #        # node of that first clade.

    #        # Get a node name for a node corresponding to a different
    #        # "type" sequence.
    #        other_type_seq_node_name = None
    #        for i in type_seq_list:
    #            if i != ts:
    #                other_type_seq_node_name = i
    #                break

    #        # Get node corresponding to a different "type" sequence.
    #        other_type_seq_node = None
    #        for node in t2.traverse():
    #            if node.name == other_type_seq_node_name:
    #                other_type_seq_node = node
    #                break

    #        # Root on the other "type" sequence node.
    #        t2.set_outgroup(other_type_seq_node)

    #        #print('\n\n\nTree rooted on a type sequence other than the first type sequence.')
    #        #print(t2)

    #    elif ts_num == 2:
    #        # Root on the first "type" sequence node for all subsequent
    #        # clades.
    #        first_type_seq_node = None
    #        for node in t2.traverse():
    #            leaf_list = []
    #            for leaf in node.get_leaves():
    #                if leaf.name == first_type_seq_node_name:
    #                    first_type_seq_node = node
    #                    break
    #        t2.set_outgroup(first_type_seq_node)
    #        #print('\n\n\nTree re-rooted on first type sequence:')
    #        #print(t2)


    #    # Make a copy of the tree topology to work with for each run
    #    # through this loop.
    #    t3 = t2.copy()

    #    # Make a list of nodes that contain type seq, but not any others. #
    #    # SEPARATE
    #    nodes_of_interest = []
    #    for node in t3.traverse():
    #        # Search in nodes that contain the type sequence.
    #        if node.search_nodes(name=ts):
    #            # Search in nodes that don't contain other type sequences.
    #            contains_other_type_seqs = False
    #            for ts2 in type_seq_list:
    #                if not ts2 == ts:
    #                    if node.search_nodes(name=ts2):
    #                        contains_other_type_seqs = True
    #            if not contains_other_type_seqs:
    #                # Add nodes of interest to list.
    #                nodes_of_interest.append(node)

    #    # find the node with the most child leaf nodes.
    #    #node_num = 0
    #    #for node in nodes_of_interest:
    #    #    node_num += 1 
    #    #    print('Node ' + str(node_num) + ' Number of leaves:  ' + str(len(node.get_leaves())))
    #    #    print(node)
    #    node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
    #            len(x.get_leaves()), reverse=True)[0]
    #    node_w_most_leaves.name = 'X'
    #    #print('\n\nClade defined by sequence ' + ts + ':')
    #    #print(node_w_most_leaves)


    #    # Copy tree and remove node/clades in constraint trees and replace with
    #    # polytomies composed of the same sequences.
    #    starting_tree = t3.copy()

    #    node_of_interest = None
    #    for node in starting_tree.traverse():
    #        if node.name == 'X':
    #            node_of_interest = node
    #            break


    #    ## Add node for additional sequence as most basal taxon.
    #    #name_for_additional_node = 'AddSeq'
    #    #node_of_interest.add_child(name=name_for_additional_node)
    #    #node_of_interest.set_outgroup('AddSeq')
    #    #print('node with additional sequence:')
    #    #print(node_of_interest)

    #    # Unroot tree, because ML trees are all unrooted.
    #    starting_tree.unroot()

    #    # Write modified tree.
    #    outtreename1 = os.path.join(outdirpath1,
    #            'clade_' + str(ts_num) + '_starting_tree_def_seq_' + ts.replace(' ', '_')\
    #            + '.tre')
    #    starting_tree.write(outfile=outtreename1, format=9)


    #    # Copy tree and remove node/clades in constraint trees and replace with
    #    # polytomies composed of the same sequences.
    #    clade_constr_tree = starting_tree.copy()

    #    ## Remove internal child nodes in the clade so that it becomes a
    #    ## polytomy. # SEPARATE
    #    #for node in clade_constr_tree.traverse():
    #    #    if node.name == 'X':
    #    #        nodes_to_delete = []
    #    #        for cn in node.iter_descendants():
    #    #            if not cn.is_leaf():
    #    #                nodes_to_delete.append(cn)
    #    #        for cn in nodes_to_delete:
    #    #            cn.delete()

    #    ## Show modified node.
    #    ##print('\nClade as polytomy:')
    #    #for node in clade_constr_tree.traverse():
    #    #    if node.name == 'X':
    #    #        #print(node)
    #    #        pass

    #    # Unroot tree, because ML trees are all unrooted.
    #    clade_constr_tree.unroot()
    #    
    #    # Write modified tree.
    #    outtreename2 = os.path.join(outdirpath1,
    #            'clade_' + str(ts_num) + '_constraint_tree_def_seq_' + ts.replace(' ', '_')\
    #            + '.tre')
    #    #print('outtreename:')
    #    #print(outtreename)
    #    #with open(outtreename, 'w') as outfh:
    #    clade_constr_tree.write(outfile=outtreename2, format=9)

    #    # Remove temporary files.
    #    os.remove(outtreefp)
    #    os.remove(simple_tree)


    ## Remove unnecessary uncoded trim file.
    #trim_file_to_remove = outalifpnex.rsplit('.', 2)[0] + '.nex'
    #os.remove(trim_file_to_remove)

    ## Define corresponding conversion table file path.
    #tablefp = outalifpnex.rsplit('.', 1)[0] + '.table'

    ## Get name of additional sequence from table.
    #additional_seq_orig_name = None
    #with open(tablefp) as tablefh:
    #    lines = tablefh.readlines()
    #    additional_seq_orig_name = lines[-1].strip()

    ## Loop through trees, parse them, change the "AddSeq" node name to the name
    ## of the additional sequence for this particular alignment, then write to a
    ## file. # SEPARATE
    #constraint_tree_list = []
    #starting_tree_list = []
    #print('\n\n\n')
    ##for tf in sorted(glob.glob(os.path.join(os.path.dirname(outalifpnex), 'clade_*.tre'))):
    ##    print(tf)

    #tf = glob.glob(os.path.join(os.path.dirname(outalifpnex), 'clade_*.tre'))[0]
    #print(tf)

    ## Parse tree.
    #tx = Tree(tf, quoted_node_names=True)
    ## Change node name in tree.
    #for node in tx.traverse():
    #    if node.name == 'AddSeq':
    #        node.name = additional_seq_orig_name
    #        break
    ## Write modified tree.
    #outtreename = outalifpnex.rsplit('.', 2)[0] + '_' + os.path.basename(tf)
    #tx.write(outfile=outtreename, format=9)
    ## Make copy of tree file with coded names.
    #outtreename2 = outtreename.rsplit('.', 1)[0] + '.C.tre'
    #code_tree(outtreename, outtreename2, tablefp)
    ## Remove copy with uncoded names.
    #os.remove(outtreename)
    ## Add file paths to lists as appropriate.
    #if os.path.basename(tf).split('_')[2] == 'constraint':
    #    constraint_tree_list.append(outtreename2)
    #elif os.path.basename(tf).split('_')[2] == 'starting':
    #    starting_tree_list.append(outtreename2)


    ## Make a copy of the alignment in phylip format for input to IQtree.
    #phy_out = outalifpnex.rsplit('.', 1)[0] + '.phy'
    #nex_to_phylip(outalifpnex, phy_out)

    ## Remove nexus format version of alignment, because it is redundant with
    ## the phylip version.
    #os.remove(outalifpnex)

    ## Print to see whether the tree files are in the right order.
    #print('\n\n\n')
    #for ctf, stf in zip(constraint_tree_list, starting_tree_list):
    #    print(ctf)
    #    print(stf)
    #    print('\n')
    #print('\n\n\n')

    ## Do phylogenetic analysis. # SEPARATE
    #outtree_files_list = []
    #for ctf, stf in zip(constraint_tree_list, starting_tree_list):
    #    # Make subdir for output.
    #    subdirp = ctf.rsplit('.', 1)[0] + '_IQ-tree_output'
    #    os.mkdir(subdirp)

    #    # Use IQtree to do an ML search
    #    output_file_prefix = os.path.join(subdirp, 'iqtree')
    #    iqtree_command_list = ['iqtree', '-s', phy_out, '-m', subs_model, '-g',
    #        ctf, '-t', stf, '-pre', output_file_prefix, '-nt', '4']
    #    tree_search_start_time = time.time()
    #    subprocess.call(iqtree_command_list)
    #    tree_search_end_time = time.time()
    #    with open(logfile, 'a') as logh:
    #        logh.write(' '.join(iqtree_command_list) + '\n')
    #        tree_search_elapsed = tree_search_end_time - tree_search_start_time
    #        logh.write('Run time: ' + str(datetime.timedelta(seconds=tree_search_elapsed)) + '\n')
    #    # Add output tree file path to list.
    #    tree_file_path = output_file_prefix + '.treefile'
    #    outtree_files_list.append(tree_file_path)
    #    # Check that the output file was actually produced.
    #    assert os.path.isfile(tree_file_path)

    ## Do topology tests using yielded trees.

    ## Make output directory for topology test analysis.
    #topo_test_subdirp = phy_out.rsplit('.', 1)[0] + '_topology_tests'
    #os.mkdir(topo_test_subdirp)
    #concat_tree_file = os.path.join(topo_test_subdirp, 'concat_trees.tre')
    #concat_tree_name_file = os.path.join(topo_test_subdirp, 'concat_tree_names.txt')
    #f_num = 0
    #alt_tree_files = []
    #for f in outtree_files_list:
    #    f_num += 1
    #    alt_tree_files.append(f)
    #    with open(f) as infh, open(concat_tree_file, 'a') as o1,\
    #    open(concat_tree_name_file, 'a') as o2:
    #        for i in infh:
    #            o1.write(i)
    #            o2.write(str(f_num) + ': ' + f)

    ## Use IQ-tree to perform topology tests.
    #topo_test_output_prefix = concat_tree_file.rsplit('.', 1)[0] + '_topo_test' 
    #topo_test_output_fp = concat_tree_file.rsplit('.', 1)[0] + '_topo_test.iqtree' 
    #subprocess.call(['iqtree', '-s', phy_out, '-m', subs_model, '-n', '0',
    #    '-z', concat_tree_file, '-zb', '10000', '-au', '-pre',
    #    topo_test_output_prefix])

    ## Parse IQ-tree output to get results of topology test. # SEPARATE
    #data_line = re.compile(r'\d+ +-\d+\.\d+ +\d+\.\d+ +')
    #space_char = re.compile(r' +')
    #topo_test_res_dict = {}
    #with open(topo_test_output_fp) as infh:
    #    for i in infh:
    #        if data_line.search(i.strip()):
    #            # Parse info from line in IQ-tree output file.
    #            parsed_tree_info_list = space_char.split(i.strip())
    #            tree_num = int(parsed_tree_info_list[0])
    #            original_tree_fp = alt_tree_files[tree_num -1]
    #            logL = float(parsed_tree_info_list[1])
    #            bp_RELL = float(parsed_tree_info_list[3])
    #            p_KH = float(parsed_tree_info_list[5])
    #            p_SH = float(parsed_tree_info_list[7])
    #            c_ELW = float(parsed_tree_info_list[9])
    #            p_AU = float(parsed_tree_info_list[11])
    #            #accept_reject = parsed_tree_info_list[12]
    #            # Add info to dict.
    #            topo_test_res_dict[tree_num] = (original_tree_fp, logL,\
    #                    bp_RELL, p_KH, p_SH, c_ELW, p_AU)

    ## Write info to summary spreadsheet file. # SEPARATE
    #with open(summary_csv_fp, 'a') as o:
    #    num_keys = str(len(topo_test_res_dict.keys()))
    #    key_num = 0
    #    for key in sorted(topo_test_res_dict.keys(), key=lambda x:\
    #            topo_test_res_dict[x][1], reverse=True):
    #        key_num += 1
    #        info_tuple = topo_test_res_dict[key]
    #        # Assemble info to write to spreadsheet.
    #        relative_tree_path = info_tuple[0].replace(os.path.dirname(os.path.dirname(info_tuple[0])), '')
    #        orig_align_file = os.path.basename(alignment)
    #        orig_tree_file = os.path.basename(tree)
    #        logL = str(info_tuple[1])
    #        bp_RELL = str(info_tuple[2])
    #        p_KH = str(info_tuple[3])
    #        p_SH = str(info_tuple[4])
    #        c_ELW = str(info_tuple[5])
    #        p_AU = str(info_tuple[6])

    #        # Apply a consensus approach to determine whether to accept or
    #        # reject topologies.
    #        accept_reject = '+'
    #        consensus_num = 0
    #        if float(bp_RELL) <= 0.05:
    #            consensus_num += 1
    #        if float(p_KH) <= 0.05:
    #            consensus_num += 1
    #        if float(p_SH) <= 0.05:
    #            consensus_num += 1
    #        if float(c_ELW) <= 0.05:
    #            consensus_num += 1
    #        if float(p_AU) <= 0.05:
    #            consensus_num += 1
    #        if consensus_num >= 4:
    #            accept_reject = '-'

    #        hypothesis_num = str(key_num) + ' of ' + num_keys
    #        relative_table_path = tablefp.replace(os.path.dirname(os.path.dirname(tablefp)), '')
    #        type_seq_name =\
    #        os.path.dirname(info_tuple[0]).split('_def_seq_')[1].split('.C_')[0]

    #        # Get clade name.
    #        clade_name = get_clade_name_from_model(type_seq_name, type_seqs)

    #        masked_ali_file = os.path.join(outdirpath1,'*mask.C.nex').replace(os.path.dirname(outdirpath1), '')

    #        # Write info line to spreadsheet.
    #        info_line = ','.join([model,
    #                              orig_align_file, #0
    #                              orig_tree_file, #1
    #                              subs_model, #2
    #                              additional_seq_orig_name, #3
    #                              hypothesis_num, #4
    #                              type_seq_name, #5
    #                              clade_name, #6
    #                              logL,  #7
    #                              bp_RELL, #8
    #                              p_KH, #8A
    #                              p_SH, #8B
    #                              c_ELW, #8C
    #                              p_AU, #8D
    #                              accept_reject,
    #                              masked_ali_file, #10
    #                              relative_table_path, #11
    #                              'Trimmed input alignment file', #12
    #                              'Constraint tree file', #13
    #                              'Starting tree file', #14
    #                              relative_tree_path, #15
    #                              '\n'
    #                              ])
    #        o.write(info_line)


def write_phylo_class_to_csv(phylo_class_id, outdir,
        csv_file, max_pvalue, timestamp, output_csv_path=None):
    """Parse output of a forward search (from running the fwd_srch command of
    amoebae) and append rows to input csv with information for interpreting
    the forward results. 
    """
    # Get path to full summary csv file.
    full_summary_csv_path = get_phylo_class_csv(outdir)

    # Get path to output spreadsheet.
    if output_csv_path == None:
        output_csv_path = csv_file.rsplit('.', 1)[0] + '_' + phylo_class_id + '.csv'
    
    # Get same list of header titles that the sum_fwd_srch command makes (for
    # parsing).
    fwd_column_label_list = module_amoebae_column_header_lists.fwd_column_label_list 

    # Get data from input file.
    print('\tReading input csv file into a pandas dataframe.')
    df = pd.read_csv(csv_file, index_col=False)

    # Manage parsing and summarization of the intermediate summary file
    # differently depending on whether the --place command was used or not.

    # Determine whether the place command was used or not.
    place = False
    #if 'Accept (+) or reject (-) hypothesis based on consensus of topology tests' not in df.loc[0]:
    #if 'Accept (+) or reject (-) hypothesis based on consensus of topology tests' not in open(full_summary_csv_path).read():
    with open(full_summary_csv_path, 'r') as x:
        if 'Accept (+) or reject (-) hypothesis based on consensus of topology tests' not in x.readline():
            place = True

    if not place:
        # Define a list of new column headers to be appended.
        new_column_label_list = module_amoebae_column_header_lists.phylo_class_column_label_list

        # Initiate new dataframe with columns to be appended/joined to spreadsheet.
        #num_rows = len(list(df.index))
        df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

        # Set default value of all fields to '-'.
        for index, row in df_new_cols.iterrows():
            df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list) 

        # Join constructed dataframe to input dataframe.
        df = df.join(df_new_cols)

        #num_rev_srchs = df['Query title'].count()
        num_potential_seqs = df['Query title'].count()

        # Iterate over rows in input spreadsheet with forward search results, and
        # determine which phylo_class output data needs to be found, and update the
        # row by adding this information.
        row_num = 0
        for index, row in df.iterrows():
            row_num += 1

            # Check that a classification needed to be done for this row.
            if row['Positive/redundant (+) or negative (-) hit based on E-value criterion'] == '-':
                pass

            # Prevent analyzing reverse searches when the forward hit rank is over
            # a certain number.
            elif row['Forward hit rank'] > settings.max_num_reverse_searches_per_database:
                pass

            else:
                # Get info needed to extract relevant info from the full summary
                # csv.
                acc = row['Forward hit accession']
                
                # Parse summary file.
                dff = pd.read_csv(full_summary_csv_path, index_col=False)

                # Get info from corresponding rows in full phylo_class summary csv file.
                info_found = False
                model_name = None
                top_index = None
                au_pvalue = None
                consensus_list = []
                top_clade = None
                second_clade = None
                for indexf, rowf in dff.iterrows():
                    if rowf['Sequence added to tree'].split(' ')[0] == acc:
                        info_found = True
                        consensus_list.append(rowf['Accept (+) or reject (-) hypothesis based on consensus of topology tests'])
                        if rowf['Hypothesis (tree/clade) number'].startswith('1 of '):
                            top_index = indexf
                            model_name = rowf['Model/backbone tree name']
                            top_clade = rowf['Clade name']

                    if top_index is not None and indexf == top_index + 1:
                        # Get info for second best classification clade.
                        au_pvalue = rowf['p-AU: p-value of approximately unbiased (AU) test (Shimodaira, 2002)']
                        second_clade = rowf['Clade name']

                # If a row with the information for the sequence with the currrent
                # Forward hit accession was found in the summary spreadsheet, then
                # write it to the dataframe.
                if info_found:
                    # Check that the info was retrieved.
                    assert model_name is not None
                    assert top_index is not None
                    assert au_pvalue is not None
                    assert top_clade is not None
                    assert second_clade is not None

                    # Determine whether all but the top classification/clade was
                    # rejected by consensus of hypothesis tests.
                    all_but_top_rejected = 'No'
                    if consensus_list.count('+') == 1:
                        all_but_top_rejected = 'Yes'

                    # Add info to output df.
                    #row['Model/backbone tree name'] = model_name
                    row['Classification'] = top_clade
                    row['Second most likely classification'] = second_clade
                    row['AU topology test p-value for comparison with next most likely classification/topology'] = au_pvalue
                    row['All but top classification rejected?'] = all_but_top_rejected
                
                    # Update row in dataframe with new information.
                    df.loc[index] = row

                # Otherwise, write '-' to the dataframe.
                else:
                    pass
    elif place:
        # Define a list of new column headers to be appended.
        new_column_label_list =\
        module_amoebae_column_header_lists.phylo_class_place_column_label_list

        # Initiate new dataframe with columns to be appended/joined to spreadsheet.
        #num_rows = len(list(df.index))
        df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

        # Set default value of all fields to '-'.
        for index, row in df_new_cols.iterrows():
            df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list) 

        # Join constructed dataframe to input dataframe.
        df = df.join(df_new_cols)

        #num_rev_srchs = df['Query title'].count()
        num_potential_seqs = df['Query title'].count()

        # Iterate over rows in input spreadsheet with forward search results, and
        # determine which phylo_class output data needs to be found, and update the
        # row by adding this information.
        row_num = 0
        for index, row in df.iterrows():
            row_num += 1

            # Check that a classification needed to be done for this row.
            if row['Positive/redundant (+) or negative (-) hit based on E-value criterion'] == '-':
                pass

            # Prevent analyzing reverse searches when the forward hit rank is over
            # a certain number.
            elif row['Forward hit rank'] > settings.max_num_reverse_searches_per_database:
                pass

            else:
                # Get info needed to extract relevant info from the full summary
                # csv.
                acc = row['Forward hit accession']
                
                # Parse summary file.
                dff = pd.read_csv(full_summary_csv_path, index_col=False)

                # Get info from corresponding rows in full phylo_class summary csv file.
                info_found = False
                model_name = None
                top_clade = None
                for indexf, rowf in dff.iterrows():
                    if rowf['ID of sequence added to tree'] == acc:
                        info_found = True
                        model_name = rowf['Model/backbone tree name']
                        top_clade = rowf['Name of clade that additional sequence was placed in']

                # If a row with the information for the sequence with the currrent
                # Forward hit accession was found in the summary spreadsheet, then
                # write it to the dataframe.
                if info_found:
                    # Check that the info was retrieved.
                    assert model_name is not None
                    assert top_clade is not None

                    # Add info to output df.
                    row['Model/backbone tree name'] = model_name
                    row['Classification'] = top_clade
                    #row['Second most likely classification'] = second_clade
                    #row['AU topology test p-value for comparison with next most likely classification/topology'] = au_pvalue
                    #row['All but top classification rejected?'] = all_but_top_rejected
                
                    # Update row in dataframe with new information.
                    df.loc[index] = row

                # Otherwise, write '-' to the dataframe.
                else:
                    pass



    # Append phylo_class id to headers (to make unique).
    df.columns = list(df.columns[0:-len(new_column_label_list)]) + [x + ' (' + phylo_class_id + ')' for x in new_column_label_list]

    # Write joined to output path.
    df.to_csv(output_csv_path, index=False)

    # Return main output path.
    return output_csv_path


def polyt_model_backbone(model_name, outfilepath=None):
    """Take a tree and make internal branches/nodes outside specified clades of
    interest into a polytomy.
    """
    # Get info from model.
    model_info = ModelInfoFromCSV(model_name)
    ali = model_info.alignment_file
    tree = model_info.tree_topology_file
    subs_model = model_info.subs_model
    type_seqs = model_info.type_seqs_file

    # Parse topology from tree for model info (need to decode names? similar to
    # previous).

    # identify the clades of interest (using type sequences as before), 
    
    # take those clades/nodes and construct a new topology (with ete3) that is
    # just a polytomy of the clades.

    # Write that new tree to the output file.

    # Get output directory path.
    outdirpath1 = os.path.dirname(ali)

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


    # Get outfilepath.
    if outfilepath is None:
        outfilepath = ali.rsplit('.', 1)[0] + '_polytomy_backbone.newick.tre'

    # Define output file path.
    #if outfilepath is None:
    #    outtreefp = ali.rsplit('.', 1)[0] + '_quoted_names.newick'
    #else:
    #    outtreefp = outfilepath
    outtreefp = ali.rsplit('.', 1)[0] + '_quoted_names.newick'

    # Put quotation marks around names in input tree. # SEPARATE
    with open(tree) as intreefh, open(outtreefp, 'w') as o:
        for t in intreefh:
            o.write(quote_tree(t, outtablefp))


    # Make topology constraint trees.   

    # Parse backbone topology and root on a clade rather than having the root
    # inside a clade.

    ######################
    # Parse tree using ete3.
    # Note: parentheses and commas get replaced with underscores.
    t1 = Tree(outtreefp, quoted_node_names=True)

    # Define name for tree without branch lengths.
    simple_tree = outtreefp.rsplit('_', 1)[0] + '_TEMP2' 

    # Write simple tree to a new file 
    t1.write(format=9, outfile=simple_tree)

    # Parse simple tree.
    t2 = Tree(simple_tree)

    # Print simple tree.
    print('Input tree:')
    print(t2)

    ### Get node for a key clade. 
    ### Or, just require to root it before inputting (figtree can export
    ### re-rooted newick)?????
    ##print("\nPlease define an arbitrary clade of interest for rooting the tree topology\n")
    ##clade_member_1 = input('Clade member 1: ')
    ##clade_member_2 = input('Clade member 2: ')
    ##ancestor1 = t2.get_common_ancestor(clade_member_1, clade_member_2)
    ###print(ancestor1)

    ### Re-root tree on a clade (as a polytomy) rather than inside a clade, so
    ### that key branches can be more easily identified.
    ##t2.set_outgroup(ancestor1)

    #######################


    ## Get to topologies sorted out...

    ## Get constraint topologies for each clade of interest to test
    ## alternative topology hypotheses.

    # Get list of "type" sequences from input.
    type_seq_list = []
    for i in open(type_seqs):
        type_seq_list.append(i.strip().split(',')[0])

    # For each "type" sequence, traverse all nodes and find the node with
    # the largest number of child nodes that are leaf (terminal) nodes,
    # containing the "type" sequence of interest, but not containing any of
    # the other "type" sequences.

    # Make a list of nodes of interest.
    nodes_of_interest_for_polytomy = []
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

            # Get node corresponding to a different "type" sequence.
            other_type_seq_node = None
            for node in t2.traverse():
                if node.name == other_type_seq_node_name:
                    other_type_seq_node = node
                    break

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

        # find the node with the most child leaf nodes.
        #node_num = 0
        #for node in nodes_of_interest:
        #    node_num += 1 
        #    print('Node ' + str(node_num) + ' Number of leaves:  ' + str(len(node.get_leaves())))
        #    print(node)
        node_w_most_leaves = sorted(nodes_of_interest, key=lambda x:\
                len(x.get_leaves()), reverse=True)[0]
        node_w_most_leaves.name = 'X'

        # Get clade name for type sequence.
        clade_name = get_clade_name_from_model(ts, type_seqs)

        # Print clade.
        print('\n\n' + clade_name + ' clade defined by sequence ' + ts + ':')
        print(node_w_most_leaves)

        nodes_of_interest_for_polytomy.append(node_w_most_leaves)


    # Construct a polytomy of the nodes of interest.

    # Count number of nodes.
    numnodes = len(nodes_of_interest_for_polytomy)

    # Make a backbone polytomy. 
    subtrees = []
    for n in nodes_of_interest_for_polytomy:
        # Append modified node to list.
        #subtrees.append(n.write(format=9, quoted_node_names=True).rstrip(';'))
        subtrees.append(n.write(format=9).rstrip(';'))

    newick_backbone = '(' + ','.join(subtrees) + ')'

    print('newick_backbone:')
    print(newick_backbone)
    with open(outfilepath, 'w') as o:
        o.write(newick_backbone)

    # Delete temporary files.
    os.remove(simple_tree)
    os.remove(outtreefp)
    os.remove(outtablefp)
    os.remove(outalifpnex)
