#!/usr/bin/env python3
"""PyTest tests for the replace_nodes_in_tree.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from replace_nodes_in_tree import \
replace_nodes_of_interest_in_tree


def test_replace_nodes_of_interest_in_tree():  # ***Incomplete test
    """Test the replace_nodes_of_interest_in_tree function in the replace_nodes_in_tree.py file.
    """
    ##########################
    # Arrange.
    t = "t"
    nodes = "nodes"
    replacement_nodes = "replacement_nodes"

    ##########################
    # Act.
    #x = replace_nodes_of_interest_in_tree(t,
    #		nodes,
    #		replacement_nodes)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


