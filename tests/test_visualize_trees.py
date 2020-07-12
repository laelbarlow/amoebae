#!/usr/bin/env python3
"""PyTest tests for the visualize_trees.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.
from ete3 import Tree

from visualize_trees import \
plotImage, \
get_first_leaf_from_species, \
get_largest_subtree, \
get_first_subtrees, \
remove_nodes_that_are_subtrees_of_others, \
determine_whether_first_subtrees_same_sp_overlap, \
get_nodes_with_paralogues, \
label_pdf, \
add_borders_to_pdf_with_latex, \
get_subs_model_from_iqtree_file, \
get_subs_model_from_mb_file, \
customize_node_styles_for_visualization, \
customize_node_styles_for_paralogue_clades, \
customize_node_styles_for_clades_to_remove, \
lighten_color, \
get_branch_support_from_node_name, \
translate_int_node_names_to_support, \
translate_int_node_support_to_prob, \
visualize_tree, \
visualize_tree_in_dir


####################################################
# Define functions for use within test functions.

def get_output_leaf_name_sets(output_nodes):
    """Take a list of ETE3 TreeNode objects, and return a list of sets of leaf
    (terminal node) names for each TreeNode object.
    """
    output_leaf_name_sets = []
    for i in output_nodes:
        output_leaf_name_sets.append(frozenset([x.name for x in i.get_leaves()]))

    return set(output_leaf_name_sets)


####################################################
# Define test functions.


def test_plotImage():  # ***Incomplete test
    """Test the plotImage function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    image_file = "image_file"

    ##########################
    # Act.
    #x = plotImage(image_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_first_leaf_from_species():  # ***Incomplete test
    """Test the get_first_leaf_from_species function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    tree_obj = "tree_obj"
    species_start_string = "species_start_string"

    ##########################
    # Act.
    #x = get_first_leaf_from_species(tree_obj,
    #		species_start_string)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_largest_subtree():  # ***Incomplete test
    """Test the get_largest_subtree function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"

    ##########################
    # Act.
    #x = get_largest_subtree(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_first_subtrees():  # ***Incomplete test
    """Test the get_first_subtrees function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"

    ##########################
    # Act.
    #x = get_first_subtrees(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_remove_nodes_that_are_subtrees_of_others():  # ***Incomplete test
    """Test the remove_nodes_that_are_subtrees_of_others function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    nodes_with_paralogues = "nodes_with_paralogues"

    ##########################
    # Act.
    #x = remove_nodes_that_are_subtrees_of_others(nodes_with_paralogues)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_determine_whether_first_subtrees_same_sp_overlap():  # ***Incomplete test
    """Test the determine_whether_first_subtrees_same_sp_overlap function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    i = "i"
    j = "j"

    ##########################
    # Act.
    #x = determine_whether_first_subtrees_same_sp_overlap(i,
    #		j)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_nodes_with_paralogues():  # ***Incomplete test
    """Test the get_nodes_with_paralogues function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.

    # Define a simple tree with two clades that should both be identified
    # as containing paralogues.
    input_tree_string_1 = '((A__1, B__1),(A__2, B__2));' 

    # Parse as an ete3 TreeNode object for input to the
    # get_nodes_with_paralogues function.
    input_tree_obj_1 = Tree(input_tree_string_1) 

    # Define a simple tree with two clades that should both be identified
    # as containing paralogues.
    input_tree_string_2 = '((X__1,(A__1, B__1)),(Y__1,(A__2, B__2)));' 

    # Parse as an ete3 TreeNode object for input to the
    # get_nodes_with_paralogues function.
    input_tree_obj_2 = Tree(input_tree_string_2) 

    # Define a simple tree with two clades that should both be identified
    # as containing paralogues.
    input_tree_string_3 = '((A__1,(X__1, B__1)),(B__2,(A__2, Y__1)));' 

    # Parse as an ete3 TreeNode object for input to the
    # get_nodes_with_paralogues function.
    input_tree_obj_3 = Tree(input_tree_string_3) 

    # Define a simple tree with two clades that should both be identified
    # as containing paralogues.
    input_tree_string_4 = '((C__1,(A__1, B__1)),(C__2,(A__2, B__2)));' 

    # Parse as an ete3 TreeNode object for input to the
    # get_nodes_with_paralogues function.
    input_tree_obj_4 = Tree(input_tree_string_4) 

    # Define a simple tree with two clades that should both be identified
    # as containing paralogues.
    input_tree_string_5 = '((A__1, B__1),(A__2, (B__2, ((C__1, D__1),((C__2, D__2),(C__3, D__3))))));' 

    # Parse as an ete3 TreeNode object for input to the
    # get_nodes_with_paralogues function.
    input_tree_obj_5 = Tree(input_tree_string_5) 

    ##########################
    # Act.

    # Get output from the get_nodes_with_paralogues function.
    output_nodes_1 = get_nodes_with_paralogues(input_tree_obj_1)

    # Get the set of leaf names for each output TreeNode object.
    output_leaf_name_sets_1 = get_output_leaf_name_sets(output_nodes_1) 
    # Define the sets of leaf names that should the output TreeNode objects
    # should have.
    correct_sets_1 = set([frozenset(['A__1', 'B__1']), frozenset(['A__2', 'B__2'])])

    # Get output from the get_nodes_with_paralogues function.
    output_nodes_2 = get_nodes_with_paralogues(input_tree_obj_2)

    # Get the set of leaf names for each output TreeNode object.
    output_leaf_name_sets_2 = get_output_leaf_name_sets(output_nodes_2)

    # Define the sets of leaf names that should the output TreeNode objects
    # should have.
    correct_sets_2 = set([frozenset(['A__1', 'B__1']), frozenset(['A__2', 'B__2'])])

    # Get output from the get_nodes_with_paralogues function.
    output_nodes_3 = get_nodes_with_paralogues(input_tree_obj_3)

    # Get the set of leaf names for each output TreeNode object.
    output_leaf_name_sets_3 = get_output_leaf_name_sets(output_nodes_3)

    # Define the sets of leaf names that should the output TreeNode objects
    # should have.
    correct_sets_3 = set([frozenset(['A__1', 'B__1', 'X__1']),
                        frozenset(['A__2', 'B__2', 'Y__1'])])

    # Get output from the get_nodes_with_paralogues function.
    output_nodes_4 = get_nodes_with_paralogues(input_tree_obj_4)

    # Get the set of leaf names for each output TreeNode object.
    output_leaf_name_sets_4 = get_output_leaf_name_sets(output_nodes_4)

    # Define the sets of leaf names that should the output TreeNode objects
    # should have.
    correct_sets_4 = set([frozenset(['C__1', 'A__1', 'B__1']),
                        frozenset(['C__2', 'A__2', 'B__2'])])

    # Get output from the get_nodes_with_paralogues function.
    output_nodes_5 = get_nodes_with_paralogues(input_tree_obj_5)

    # Get the set of leaf names for each output TreeNode object.
    output_leaf_name_sets_5 = get_output_leaf_name_sets(output_nodes_5)

    # Define the sets of leaf names that should the output TreeNode objects
    # should have.
    correct_sets_5 = set([frozenset(['A__1', 'B__1']),
                        frozenset(['A__2', 'B__2', 'C__1', 'D__1', 'C__2', 'D__2', 'C__3', 'D__3']),
                        frozenset(['C__1', 'D__1']),
                        frozenset(['C__2', 'D__2']),
                        frozenset(['C__3', 'D__3'])
                        ])

    ##########################
    # Assert.

    # Check that the output set is the expected set (set of sets).

    # Check that clades with paralogues are identified.
    assert output_leaf_name_sets_1 == correct_sets_1

    # Test that early branching nonparalogous sequences are excluded.
    assert output_leaf_name_sets_2 == correct_sets_2
    
    # Test that nested nonparalogous sequences are included.
    assert output_leaf_name_sets_3 == correct_sets_3

    # Test that redundant subclades are not identified.
    assert output_leaf_name_sets_4 == correct_sets_4

    # Test that relevant subtrees are not ignored.
    assert output_leaf_name_sets_5 == correct_sets_5


def test_label_pdf():  # ***Incomplete test
    """Test the label_pdf function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    in_tree_file = "in_tree_file"
    out_tree_file = "out_tree_file"
    method = "method"
    subs_model = "subs_model"
    timestamp = "timestamp"

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



def test_add_borders_to_pdf_with_latex():  # ***Incomplete test
    """Test the add_borders_to_pdf_with_latex function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    in_tree_file = "in_tree_file"
    out_tree_file = "out_tree_file"

    ##########################
    # Act.
    #x = add_borders_to_pdf_with_latex(in_tree_file,
    #		out_tree_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_subs_model_from_iqtree_file():  # ***Incomplete test
    """Test the get_subs_model_from_iqtree_file function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    file_with_subs_model_name = "file_with_subs_model_name"

    ##########################
    # Act.
    #x = get_subs_model_from_iqtree_file(file_with_subs_model_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_subs_model_from_mb_file():  # ***Incomplete test
    """Test the get_subs_model_from_mb_file function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    file_with_subs_model_name = "file_with_subs_model_name"

    ##########################
    # Act.
    #x = get_subs_model_from_mb_file(file_with_subs_model_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_customize_node_styles_for_visualization():  # ***Incomplete test
    """Test the customize_node_styles_for_visualization function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"

    ##########################
    # Act.
    #x = customize_node_styles_for_visualization(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_customize_node_styles_for_paralogue_clades():  # ***Incomplete test
    """Test the customize_node_styles_for_paralogue_clades function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"
    nodes_with_paralogues = "nodes_with_paralogues"

    ##########################
    # Act.
    #x = customize_node_styles_for_paralogue_clades(t,
    #		nodes_with_paralogues)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_customize_node_styles_for_clades_to_remove():  # ***Incomplete test
    """Test the customize_node_styles_for_clades_to_remove function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"
    highlight_for_removal = "highlight_for_removal"

    ##########################
    # Act.
    #x = customize_node_styles_for_clades_to_remove(t,
    #		highlight_for_removal)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_lighten_color():  # ***Incomplete test
    """Test the lighten_color function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    color = "color"
    amount = "amount=0.5"

    ##########################
    # Act.
    #x = lighten_color(color,
    #		amount)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_branch_support_from_node_name():  # ***Incomplete test
    """Test the get_branch_support_from_node_name function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    node_name = "node_name"

    ##########################
    # Act.
    #x = get_branch_support_from_node_name(node_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_translate_int_node_names_to_support():  # ***Incomplete test
    """Test the translate_int_node_names_to_support function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"

    ##########################
    # Act.
    #x = translate_int_node_names_to_support(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_translate_int_node_support_to_prob():  # ***Incomplete test
    """Test the translate_int_node_support_to_prob function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    t = "t"

    ##########################
    # Act.
    #x = translate_int_node_support_to_prob(t)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_visualize_tree():  # ***Incomplete test
    """Test the visualize_tree function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    method = "method"

    ##########################
    # Act.
    #x = visualize_tree(method)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_visualize_tree_in_dir():  # ***Incomplete test
    """Test the visualize_tree_in_dir function in the visualize_trees.py file.
    """
    ##########################
    # Arrange.
    indp = "indp"
    method = "method"
    timestamp = "timestamp"
    taxa_to_root_on = "taxa_to_root_on"

    ##########################
    # Act.
    #x = visualize_tree_in_dir(indp,
    #		method,
    #		timestamp,
    #		taxa_to_root_on)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


