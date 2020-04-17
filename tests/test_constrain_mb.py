#!/usr/bin/env python3
"""PyTest tests for the constrain_mb.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from constrain_mb import \
get_taxon_number_dict, \
constrain_mb_with_tree


def test_get_taxon_number_dict():  # ***Incomplete test
    """Test the get_taxon_number_dict function in the constrain_mb.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"

    ##########################
    # Act.
    #x = get_taxon_number_dict(alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_constrain_mb_with_tree():  # ***Incomplete test
    """Test the constrain_mb_with_tree function in the constrain_mb.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"
    tree = "tree"
    out_alignment = "out_alignment=None"

    ##########################
    # Act.
    #x = constrain_mb_with_tree(alignment,
    #		tree,
    #		out_alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


