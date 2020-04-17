#!/usr/bin/env python3
"""PyTest tests for the get_alt_topos.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

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


def test_get_polytomy_for_treenode():  # ***Incomplete test
    """Test the get_polytomy_for_treenode function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    treenode = "treenode"

    ##########################
    # Act.
    #x = get_polytomy_for_treenode(treenode)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



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



def test_trees():  # ***Incomplete test
    """Test the trees function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    leaves = "leaves"
    polytomies = "polytomies=False"

    ##########################
    # Act.
    #x = trees(leaves,
    #		polytomies)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_total_num_topos_for_n_taxa():  # ***Incomplete test
    """Test the get_total_num_topos_for_n_taxa function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    n = "n"

    ##########################
    # Act.
    #x = get_total_num_topos_for_n_taxa(n)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



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



def test_get_all_alt_topologies():  # ***Incomplete test
    """Test the get_all_alt_topologies function in the get_alt_topos.py file.
    """
    ##########################
    # Arrange.
    subtrees = "subtrees"

    ##########################
    # Act.
    #x = get_all_alt_topologies(subtrees)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


