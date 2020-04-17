#!/usr/bin/env python3
"""PyTest tests for the visualize_trees.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from visualize_trees import *



def test_plotImage():
    """Test the plotImage function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    image_file = 'image_file'

    ##########################
    # Act.
    #x = plotImage(image_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_first_leaf_from_species():
    """Test the get_first_leaf_from_species function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    tree_obj = 'tree_obj'
    species_start_string = 'species_start_string'

    ##########################
    # Act.
    #x = get_first_leaf_from_species(tree_obj,
    #		species_start_string)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_largest_subtree():
    """Test the get_largest_subtree function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'

    ##########################
    # Act.
    #x = get_largest_subtree(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_first_subtrees():
    """Test the get_first_subtrees function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'

    ##########################
    # Act.
    #x = get_first_subtrees(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_remove_nodes_that_are_subtrees_of_others():
    """Test the remove_nodes_that_are_subtrees_of_others function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    nodes_with_paralogues = 'nodes_with_paralogues'

    ##########################
    # Act.
    #x = remove_nodes_that_are_subtrees_of_others(nodes_with_paralogues)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_determine_whether_first_subtrees_same_sp_overlap():
    """Test the determine_whether_first_subtrees_same_sp_overlap function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    i = 'i'
    j = 'j'

    ##########################
    # Act.
    #x = determine_whether_first_subtrees_same_sp_overlap(i,
    #		j)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_nodes_with_paralogues():
    """Test the get_nodes_with_paralogues function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'

    ##########################
    # Act.
    #x = get_nodes_with_paralogues(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_label_pdf():
    """Test the label_pdf function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    in_tree_file = 'in_tree_file'
    out_tree_file = 'out_tree_file'
    method = 'method'
    subs_model = 'subs_model'
    timestamp = 'timestamp'

    ##########################
    # Act.
    #x = label_pdf(in_tree_file,
    #		out_tree_file,
    #		method,
    #		subs_model,
    #		timestamp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_add_borders_to_pdf_with_latex():
    """Test the add_borders_to_pdf_with_latex function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    in_tree_file = 'in_tree_file'
    out_tree_file = 'out_tree_file'

    ##########################
    # Act.
    #x = add_borders_to_pdf_with_latex(in_tree_file,
    #		out_tree_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_subs_model_from_iqtree_file():
    """Test the get_subs_model_from_iqtree_file function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    file_with_subs_model_name = 'file_with_subs_model_name'

    ##########################
    # Act.
    #x = get_subs_model_from_iqtree_file(file_with_subs_model_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_subs_model_from_mb_file():
    """Test the get_subs_model_from_mb_file function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    file_with_subs_model_name = 'file_with_subs_model_name'

    ##########################
    # Act.
    #x = get_subs_model_from_mb_file(file_with_subs_model_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_customize_node_styles_for_visualization():
    """Test the customize_node_styles_for_visualization function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'

    ##########################
    # Act.
    #x = customize_node_styles_for_visualization(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_customize_node_styles_for_paralogue_clades():
    """Test the customize_node_styles_for_paralogue_clades function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'
    nodes_with_paralogues = 'nodes_with_paralogues'

    ##########################
    # Act.
    #x = customize_node_styles_for_paralogue_clades(t,
    #		nodes_with_paralogues)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_customize_node_styles_for_clades_to_remove():
    """Test the customize_node_styles_for_clades_to_remove function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'
    highlight_for_removal = 'highlight_for_removal'

    ##########################
    # Act.
    #x = customize_node_styles_for_clades_to_remove(t,
    #		highlight_for_removal)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_lighten_color():
    """Test the lighten_color function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    color = 'color'
    amount = 'amount=0.5'

    ##########################
    # Act.
    #x = lighten_color(color,
    #		amount)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_branch_support_from_node_name():
    """Test the get_branch_support_from_node_name function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    node_name = 'node_name'

    ##########################
    # Act.
    #x = get_branch_support_from_node_name(node_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_translate_int_node_names_to_support():
    """Test the translate_int_node_names_to_support function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'

    ##########################
    # Act.
    #x = translate_int_node_names_to_support(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_translate_int_node_support_to_prob():
    """Test the translate_int_node_support_to_prob function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = 't'

    ##########################
    # Act.
    #x = translate_int_node_support_to_prob(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_visualize_tree():
    """Test the visualize_tree function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    method = 'method'

    ##########################
    # Act.
    #x = visualize_tree(method)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_visualize_tree_in_dir():
    """Test the visualize_tree_in_dir function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    indp = 'indp'
    method = 'method'
    timestamp = 'timestamp'
    taxa_to_root_on = 'taxa_to_root_on'

    ##########################
    # Act.
    #x = visualize_tree_in_dir(indp,
    #		method,
    #		timestamp,
    #		taxa_to_root_on)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


