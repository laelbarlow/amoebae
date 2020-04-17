#!/usr/bin/env python3
"""PyTest tests for the prune.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from prune import *



def test_reformat_combined_supports2():
    """Test the reformat_combined_supports2 function in the prune.py file.
    """
    ##########################
    # Arrange.
    tree_string = 'tree_string'
    one_support = 'one_support='bootstrap''

    ##########################
    # Act.
    #x = reformat_combined_supports2(tree_string,
    #		one_support)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_make_quoted_table():
    """Test the make_quoted_table function in the prune.py file.
    """
    ##########################
    # Arrange.
    intablepath = 'intablepath'
    outtablepath = 'outtablepath'

    ##########################
    # Act.
    #x = make_quoted_table(intablepath,
    #		outtablepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_conversion_table_dict():
    """Test the get_conversion_table_dict function in the prune.py file.
    """
    ##########################
    # Arrange.
    table_with_quotation_marks = 'table_with_quotation_marks'

    ##########################
    # Act.
    #x = get_conversion_table_dict(table_with_quotation_marks)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_tree_string_with_readable_supports():
    """Test the get_tree_string_with_readable_supports function in the prune.py file.
    """
    ##########################
    # Arrange.
    tree_with_long_names = 'tree_with_long_names'

    ##########################
    # Act.
    #x = get_tree_string_with_readable_supports(tree_with_long_names)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_root_tree_manually():
    """Test the root_tree_manually function in the prune.py file.
    """
    ##########################
    # Arrange.
    t = 't'
    tx = 'tx'

    ##########################
    # Act.
    #x = root_tree_manually(t,
    #		tx)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_prune_tree_style():
    """Test the get_prune_tree_style function in the prune.py file.
    """
    ##########################
    # Arrange.
    

    ##########################
    # Act.
    #x = get_prune_tree_style()

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seqs_to_remove_file_from_user():
    """Test the get_seqs_to_remove_file_from_user function in the prune.py file.
    """
    ##########################
    # Arrange.
    to_remove_list_file = 'to_remove_list_file'
    tree_file_name = 'tree_file_name'

    ##########################
    # Act.
    #x = get_seqs_to_remove_file_from_user(to_remove_list_file,
    #		tree_file_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_leaves_all_same_species():
    """Test the leaves_all_same_species function in the prune.py file.
    """
    ##########################
    # Arrange.
    leaves = 'leaves'

    ##########################
    # Act.
    #x = leaves_all_same_species(leaves)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seqs_to_remove_file():
    """Test the get_seqs_to_remove_file function in the prune.py file.
    """
    ##########################
    # Arrange.
    to_remove_list_file = 'to_remove_list_file'
    t = 't'
    remove_redun_seqs = 'remove_redun_seqs'
    min_support = 'min_support'

    ##########################
    # Act.
    #x = get_seqs_to_remove_file(to_remove_list_file,
    #		t,
    #		remove_redun_seqs,
    #		min_support)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_key_nodes_list():
    """Test the get_key_nodes_list function in the prune.py file.
    """
    ##########################
    # Arrange.
    to_remove_list_file = 'to_remove_list_file'

    ##########################
    # Act.
    #x = get_key_nodes_list(to_remove_list_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seqs_to_remove_from_dataset():
    """Test the get_seqs_to_remove_from_dataset function in the prune.py file.
    """
    ##########################
    # Arrange.
    key_nodes_list = 'key_nodes_list'
    t = 't'

    ##########################
    # Act.
    #x = get_seqs_to_remove_from_dataset(key_nodes_list,
    #		t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_reduced_alignment():
    """Test the write_reduced_alignment function in the prune.py file.
    """
    ##########################
    # Arrange.
    alignment_file = 'alignment_file'

    ##########################
    # Act.
    #x = write_reduced_alignment(alignment_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_manually_select_nodes_and_remove_seqs():
    """Test the manually_select_nodes_and_remove_seqs function in the prune.py file.
    """
    ##########################
    # Arrange.
    input_tree_one = 'input_tree_one'
    alignment_file = 'alignment_file'

    ##########################
    # Act.
    #x = manually_select_nodes_and_remove_seqs(input_tree_one,
    #		alignment_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_automatically_select_nodes_and_remove_seqs():
    """Test the automatically_select_nodes_and_remove_seqs function in the prune.py file.
    """
    ##########################
    # Arrange.
    input_tree_one = 'input_tree_one'

    ##########################
    # Act.
    #x = automatically_select_nodes_and_remove_seqs(input_tree_one)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_automatically_select_nodes_and_remove_seqs_in_dir():
    """Test the automatically_select_nodes_and_remove_seqs_in_dir function in the prune.py file.
    """
    ##########################
    # Arrange.
    indp = 'indp'

    ##########################
    # Act.
    #x = automatically_select_nodes_and_remove_seqs_in_dir(indp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_manually_select_nodes_to_constrain():
    """Test the manually_select_nodes_to_constrain function in the prune.py file.
    """
    ##########################
    # Arrange.
    input_tree_one = 'input_tree_one'
    name_replace_table = 'name_replace_table'

    ##########################
    # Act.
    #x = manually_select_nodes_to_constrain(input_tree_one,
    #		name_replace_table)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_remove_nodes_to_match_ali():
    """Test the remove_nodes_to_match_ali function in the prune.py file.
    """
    ##########################
    # Arrange.
    alignment = 'alignment'
    tree_file = 'tree_file'
    output_file = 'output_file'

    ##########################
    # Act.
    #x = remove_nodes_to_match_ali(alignment,
    #		tree_file,
    #		output_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


