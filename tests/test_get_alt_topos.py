#!/usr/bin/env python3
"""PyTest tests for the get_alt_topos.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.
from ete3 import Tree

from get_alt_topos import \
get_polytomy_for_treenode, \
determine_if_polytomy, \
get_partitions_of_set, \
get_newick_tree_string_from_set, \
trees, \
get_total_num_topos_for_n_taxa, \
get_leaf_names_to_represent_n_taxa, \
get_all_unique_unrooted_bifurcating_topologies_for_n_taxa, \
get_all_alt_topologies


def test_get_polytomy_for_treenode():
    """Test the get_polytomy_for_treenode function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.

    # Generate a TreeNode object.
    input_treenode_1 = Tree('((A,B),(C,D));')

    ##########################
    # Act.

    # Call function to get polytomy TreeNode object.
    output_treenode_1 = get_polytomy_for_treenode(input_treenode_1)

    ##########################
    # Assert.

    # Check that the output TreeNode is actually a polytomy.
    assert output_treenode_1.write(format=9) == '(A,B,C,D);'


def test_determine_if_polytomy():  # ***Incomplete test
    """Test the determine_if_polytomy function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    subtree = "subtree"

    ##########################
    # Act.
    #x = determine_if_polytomy(subtree)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_partitions_of_set():  # ***Incomplete test
    """Test the get_partitions_of_set function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    s = "s"

    ##########################
    # Act.
    #x = get_partitions_of_set(s)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_newick_tree_string_from_set():  # ***Incomplete test
    """Test the get_newick_tree_string_from_set function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    s = "s"

    ##########################
    # Act.
    #x = get_newick_tree_string_from_set(s)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_trees():
    """Test the trees function in the get_alt_topos.py file.  This function
    recursively generates all alternative trees of sets of elements of the input
    list of strings.
    """
    ##########################
    # Arrange.
    leaves_1 = ['A', 'B', 'C']
    leaves_2 = ['A', 'B', 'C', 'D']
    leaves_3 = ['A', 'B', 'C', 'D', 'E']
    leaves_4 = ['A', 'B', 'C', 'D', 'E', 'F']
    leaves_5 = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

    ##########################
    # Act.
    all_trees_1 = frozenset(
            {frozenset({tree, leaves_1[0]}) for tree in trees(set(leaves_1[1:]))})
    all_trees_2 = frozenset(
            {frozenset({tree, leaves_2[0]}) for tree in trees(set(leaves_2[1:]))})
    all_trees_3 = frozenset(
            {frozenset({tree, leaves_3[0]}) for tree in trees(set(leaves_3[1:]))})
    all_trees_4 = frozenset(
            {frozenset({tree, leaves_4[0]}) for tree in trees(set(leaves_4[1:]))})
    all_trees_5 = frozenset(
            {frozenset({tree, leaves_5[0]}) for tree in trees(set(leaves_5[1:]))})

    ##########################
    # Assert.
    assert len(list(all_trees_1)) == get_total_num_topos_for_n_taxa(len(leaves_1))
    assert len(list(all_trees_2)) == get_total_num_topos_for_n_taxa(len(leaves_2))
    assert len(list(all_trees_3)) == get_total_num_topos_for_n_taxa(len(leaves_3))
    assert len(list(all_trees_4)) == get_total_num_topos_for_n_taxa(len(leaves_4))
    assert len(list(all_trees_5)) == get_total_num_topos_for_n_taxa(len(leaves_5))


def test_get_total_num_topos_for_n_taxa():
    """Test the get_total_num_topos_for_n_taxa function in the get_alt_topos.py
    file.  This function is for determining how many unique unrooted topologies
    should be produced in total.
    """
    ##########################
    # Arrange.
    n3 = 3
    n4 = 4
    n5 = 5
    n6 = 6

    ##########################
    # Act.
    x3 = get_total_num_topos_for_n_taxa(n3)
    x4 = get_total_num_topos_for_n_taxa(n4)
    x5 = get_total_num_topos_for_n_taxa(n5)
    x6 = get_total_num_topos_for_n_taxa(n6)

    ##########################
    # Assert.
    assert x3 == 1
    assert x4 == 3
    assert x5 == 15
    assert x6 == 105


def test_get_leaf_names_to_represent_n_taxa():  # ***Incomplete test
    """Test the get_leaf_names_to_represent_n_taxa function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    n = "n"

    ##########################
    # Act.
    #x = get_leaf_names_to_represent_n_taxa(n)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_all_unique_unrooted_bifurcating_topologies_for_n_taxa():  # ***Incomplete test
    """Test the get_all_unique_unrooted_bifurcating_topologies_for_n_taxa function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    n = "n"

    ##########################
    # Act.
    #x = get_all_unique_unrooted_bifurcating_topologies_for_n_taxa(n)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_all_alt_topologies():
    """Test the get_all_alt_topologies function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.

    # Generate lists of subtrees (ete3 TreeNode objects) as input to
    # function for finding alternative topologies.

    subtrees_1 = [] 
    for i in range(1, 4):
        string = '(C__' + str(i) + ',(A__' + str(i) + ', B__' + str(i) + '));'
        subtrees_1.append(Tree(string))

    subtrees_2 = [] 
    for i in range(1, 5):
        string = '(C__' + str(i) + ',(A__' + str(i) + ', B__' + str(i) + '));'
        subtrees_2.append(Tree(string))

    subtrees_3 = [] 
    for i in range(1, 6):
        string = '(C__' + str(i) + ',(A__' + str(i) + ', B__' + str(i) + '));'
        subtrees_3.append(Tree(string))

    ##########################
    # Act.

    # Call the function to get all alternative tree topologies.
    result_1 = get_all_alt_topologies(subtrees_1)
    result_2 = get_all_alt_topologies(subtrees_2)
    result_3 = get_all_alt_topologies(subtrees_3)

    ##########################
    # Assert.

    # Loop over the lists of alternative topologies (TreeNode objects).
    for alt_topos in [result_1, result_2, result_3]:
        leaf_name_lists = []
        ascii_topos = []
        for alt_topo in alt_topos:
            # Add leaf names to list.
            #print(alt_topo[0])
            leaf_names = str([leaf.name for leaf in Tree(alt_topo[0]).get_leaves()].sort())
            leaf_name_lists.append(leaf_names)

            # Add ASCII string of topology to list.
            t = Tree(alt_topo[0])
            ascii_str = t.write()
            ascii_topos.append(ascii_str)

        # Check that the lists of leaf names are all the same (i.e., that each
        # alternative topology has the same leaves).
        assert len(list(set(leaf_name_lists))) == 1

        # Check that the ETE3 ASCII tree strings are different.
        assert len(ascii_topos) == len(list(set(ascii_topos)))



