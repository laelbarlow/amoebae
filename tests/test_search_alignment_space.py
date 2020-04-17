#!/usr/bin/env python3
"""PyTest tests for the search_alignment_space.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from search_alignment_space import \
trim_one_column_from_alignment, \
get_taxa_represented_in_clade, \
get_species_name_from_seq_id, \
get_names_of_dispensable_seqs_in_clade, \
new_seq_in_clade_dispensable, \
get_new_type_seqs_dict, \
reduce_alignment, \
write_constraint_tree_without_extra_parentheses, \
get_seq_placement_in_tree, \
modify_alignment_in_x_way, \
get_taxonomic_info, \
define_textface_for_labeling_stem, \
get_corresponding_node, \
get_list_of_leaf_names_for_node, \
get_clade_name_from_model2, \
get_nodes_from_actual_tree_obj, \
get_nodes_of_interest, \
define_nodestyles_dict_for_colourcoding, \
get_abayes_support_for_node_from_another_tree, \
get_ml_tree_info_dict, \
get_constraint_tree_fp, \
run_tree_for_branch_lengths_and_supports_for_topology, \
get_type_seqs_dict, \
get_essential_taxa_list, \
get_y_measure_of_support, \
get_ali_length, \
new_tree_better, \
search_alignment_space


def test_trim_one_column_from_alignment():  # ***Incomplete test
    """Test the trim_one_column_from_alignment function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"
    alignment2 = "alignment2"
    column_index = "column_index"

    ##########################
    # Act.
    #x = trim_one_column_from_alignment(alignment,
    #		alignment2,
    #		column_index)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_taxa_represented_in_clade():  # ***Incomplete test
    """Test the get_taxa_represented_in_clade function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    clade = "clade"
    ml_tree_info_dict = "ml_tree_info_dict"

    ##########################
    # Act.
    #x = get_taxa_represented_in_clade(clade,
    #		ml_tree_info_dict)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_species_name_from_seq_id():  # ***Incomplete test
    """Test the get_species_name_from_seq_id function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    seq_id = "seq_id"

    ##########################
    # Act.
    #x = get_species_name_from_seq_id(seq_id)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_names_of_dispensable_seqs_in_clade():  # ***Incomplete test
    """Test the get_names_of_dispensable_seqs_in_clade function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    clade = "clade"

    ##########################
    # Act.
    #x = get_names_of_dispensable_seqs_in_clade(clade)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_new_seq_in_clade_dispensable():  # ***Incomplete test
    """Test the new_seq_in_clade_dispensable function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    essential_taxa = "essential_taxa"

    ##########################
    # Act.
    #x = new_seq_in_clade_dispensable(essential_taxa)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_new_type_seqs_dict():  # ***Incomplete test
    """Test the get_new_type_seqs_dict function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    seqname = "seqname"
    ml_tree_info_dict = "ml_tree_info_dict"
    type_seqs_dict = "type_seqs_dict"

    ##########################
    # Act.
    #x = get_new_type_seqs_dict(seqname,
    #		ml_tree_info_dict,
    #		type_seqs_dict)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_alignment():  # ***Incomplete test
    """Test the reduce_alignment function in the search_alignment_space.py file.
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



def test_write_constraint_tree_without_extra_parentheses():  # ***Incomplete test
    """Test the write_constraint_tree_without_extra_parentheses function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    ete3_tree_obj = "ete3_tree_obj"
    tree_outpath = "tree_outpath"

    ##########################
    # Act.
    #x = write_constraint_tree_without_extra_parentheses(ete3_tree_obj,
    #		tree_outpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seq_placement_in_tree():  # ***Incomplete test
    """Test the get_seq_placement_in_tree function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    ml_placement_tree_path = "ml_placement_tree_path"

    ##########################
    # Act.
    #x = get_seq_placement_in_tree(ml_placement_tree_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_modify_alignment_in_x_way():  # ***Incomplete test
    """Test the modify_alignment_in_x_way function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    previous_ali_tree_tuple = "previous_ali_tree_tuple"
    mod_type = "mod_type"

    ##########################
    # Act.
    #x = modify_alignment_in_x_way(previous_ali_tree_tuple,
    #		mod_type)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_taxonomic_info():  # ***Incomplete test
    """Test the get_taxonomic_info function in the search_alignment_space.py file.
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



def test_define_textface_for_labeling_stem():  # ***Incomplete test
    """Test the define_textface_for_labeling_stem function in the search_alignment_space.py file.
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
    """Test the get_corresponding_node function in the search_alignment_space.py file.
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



def test_get_list_of_leaf_names_for_node():  # ***Incomplete test
    """Test the get_list_of_leaf_names_for_node function in the search_alignment_space.py file.
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



def test_get_clade_name_from_model2():  # ***Incomplete test
    """Test the get_clade_name_from_model2 function in the search_alignment_space.py file.
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
    """Test the get_nodes_from_actual_tree_obj function in the search_alignment_space.py file.
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
    """Test the get_nodes_of_interest function in the search_alignment_space.py file.
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



def test_define_nodestyles_dict_for_colourcoding():  # ***Incomplete test
    """Test the define_nodestyles_dict_for_colourcoding function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    

    ##########################
    # Act.
    #x = define_nodestyles_dict_for_colourcoding()

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_abayes_support_for_node_from_another_tree():  # ***Incomplete test
    """Test the get_abayes_support_for_node_from_another_tree function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    node = "node"

    ##########################
    # Act.
    #x = get_abayes_support_for_node_from_another_tree(node)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_ml_tree_info_dict():  # ***Incomplete test
    """Test the get_ml_tree_info_dict function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    ml_tree_path = "ml_tree_path"

    ##########################
    # Act.
    #x = get_ml_tree_info_dict(ml_tree_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_constraint_tree_fp():  # ***Incomplete test
    """Test the get_constraint_tree_fp function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    outalifpnex = "outalifpnex"

    ##########################
    # Act.
    #x = get_constraint_tree_fp(outalifpnex)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_run_tree_for_branch_lengths_and_supports_for_topology():  # ***Incomplete test
    """Test the run_tree_for_branch_lengths_and_supports_for_topology function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    tree = "tree"

    ##########################
    # Act.
    #x = run_tree_for_branch_lengths_and_supports_for_topology(tree)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_type_seqs_dict():  # ***Incomplete test
    """Test the get_type_seqs_dict function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    type_seqs = "type_seqs"

    ##########################
    # Act.
    #x = get_type_seqs_dict(type_seqs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_essential_taxa_list():  # ***Incomplete test
    """Test the get_essential_taxa_list function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    essential_taxa_file = "essential_taxa_file"

    ##########################
    # Act.
    #x = get_essential_taxa_list(essential_taxa_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_y_measure_of_support():  # ***Incomplete test
    """Test the get_y_measure_of_support function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    previous_ali_tree_tuple = "previous_ali_tree_tuple"

    ##########################
    # Act.
    #x = get_y_measure_of_support(previous_ali_tree_tuple)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_ali_length():  # ***Incomplete test
    """Test the get_ali_length function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"

    ##########################
    # Act.
    #x = get_ali_length(alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_new_tree_better():  # ***Incomplete test
    """Test the new_tree_better function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    prev_tree_measure = "prev_tree_measure"

    ##########################
    # Act.
    #x = new_tree_better(prev_tree_measure)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_search_alignment_space():  # ***Incomplete test
    """Test the search_alignment_space function in the search_alignment_space.py file.
    """
    ##########################
    # Arrange.
    model_name = "model_name"

    ##########################
    # Act.
    #x = search_alignment_space(model_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


