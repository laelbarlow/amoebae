#!/usr/bin/env python3
"""PyTest tests for the select_seqs.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from select_seqs import \
get_clade_name_from_model2, \
get_nodes_from_actual_tree_obj, \
get_nodes_of_interest, \
get_list_of_leaf_names_for_node, \
get_taxonomic_info, \
define_nodestyles_dict_for_colourcoding, \
define_textface_for_labeling_stem, \
get_corresponding_node, \
optimize_sequence_selection, \
get_ml_tree_branch_lengths, \
get_branch_length_info, \
reduce_alignment, \
recursively_remove_seqs, \
recursively_add_seqs, \
optimize_sequence_selection2


def test_get_clade_name_from_model2():  # ***Incomplete test
    """Test the get_clade_name_from_model2 function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    type_seq_name = "type_seq_name"
    type_seqs_dict = "type_seqs_dict"

    ##########################
    # Act.
    #x = get_clade_name_from_model2(type_seq_name,
    #		type_seqs_dict)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_nodes_from_actual_tree_obj():  # ***Incomplete test
    """Test the get_nodes_from_actual_tree_obj function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    t3 = "t3"
    nodes_of_interest = "nodes_of_interest"

    ##########################
    # Act.
    #x = get_nodes_from_actual_tree_obj(t3,
    #		nodes_of_interest)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_nodes_of_interest():  # ***Incomplete test
    """Test the get_nodes_of_interest function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    tree = "tree"
    type_seq_list = "type_seq_list"

    ##########################
    # Act.
    #x = get_nodes_of_interest(tree,
    #		type_seq_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_list_of_leaf_names_for_node():  # ***Incomplete test
    """Test the get_list_of_leaf_names_for_node function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    node = "node"

    ##########################
    # Act.
    #x = get_list_of_leaf_names_for_node(node)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_taxonomic_info():  # ***Incomplete test
    """Test the get_taxonomic_info function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    species_name = "species_name"

    ##########################
    # Act.
    #x = get_taxonomic_info(species_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_define_nodestyles_dict_for_colourcoding():  # ***Incomplete test
    """Test the define_nodestyles_dict_for_colourcoding function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    

    ##########################
    # Act.
    #x = define_nodestyles_dict_for_colourcoding()

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_define_textface_for_labeling_stem():  # ***Incomplete test
    """Test the define_textface_for_labeling_stem function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    clade_name = "clade_name"

    ##########################
    # Act.
    #x = define_textface_for_labeling_stem(clade_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_corresponding_node():  # ***Incomplete test
    """Test the get_corresponding_node function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    node = "node"
    tree = "tree"

    ##########################
    # Act.
    #x = get_corresponding_node(node,
    #		tree)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_optimize_sequence_selection():  # ***Incomplete test
    """Test the optimize_sequence_selection function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    file_with_seqs = "file_with_seqs"
    model_name = "model_name"
    essential_taxa = "essential_taxa"

    ##########################
    # Act.
    #x = optimize_sequence_selection(file_with_seqs,
    #		model_name,
    #		essential_taxa)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_ml_tree_branch_lengths():  # ***Incomplete test
    """Test the get_ml_tree_branch_lengths function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    tree = "tree"
    alignment = "alignment"
    subs_model = "subs_model"
    outputdir = "outputdir"

    ##########################
    # Act.
    #x = get_ml_tree_branch_lengths(tree,
    #		alignment,
    #		subs_model,
    #		outputdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_branch_length_info():  # ***Incomplete test
    """Test the get_branch_length_info function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    t1 = "t1"
    type_seqs_dict = "type_seqs_dict"
    outpath = "outpath=None"

    ##########################
    # Act.
    #x = get_branch_length_info(t1,
    #		type_seqs_dict,
    #		outpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_alignment():  # ***Incomplete test
    """Test the reduce_alignment function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    alignment_file = "alignment_file"
    output_alignment_file = "output_alignment_file"
    removal_name_list = "removal_name_list"

    ##########################
    # Act.
    #x = reduce_alignment(alignment_file,
    #		output_alignment_file,
    #		removal_name_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_recursively_remove_seqs():  # ***Incomplete test
    """Test the recursively_remove_seqs function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    recursion_num = "recursion_num"

    ##########################
    # Act.
    #x = recursively_remove_seqs(recursion_num)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_recursively_add_seqs():  # ***Incomplete test
    """Test the recursively_add_seqs function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    recursion_num = "recursion_num"

    ##########################
    # Act.
    #x = recursively_add_seqs(recursion_num)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_optimize_sequence_selection2():  # ***Incomplete test
    """Test the optimize_sequence_selection2 function in the select_seqs.py file.
    """
    ##########################
    # Arrange.
    file_with_seqs = "file_with_seqs"
    model_name = "model_name"
    essential_taxa = "essential_taxa"

    ##########################
    # Act.
    #x = optimize_sequence_selection2(file_with_seqs,
    #		model_name,
    #		essential_taxa)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


