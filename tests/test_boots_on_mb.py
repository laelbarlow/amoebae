#!/usr/bin/env python3
"""PyTest tests for the boots_on_mb.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from boots_on_mb import *



def test_reformat_combined_supports():
    """Test the reformat_combined_supports function in the boots_on_mb.py file.
    """
    ##########################
    # Arrange.
    tree_string = 'tree_string'

    ##########################
    # Act.
    #x = reformat_combined_supports(tree_string)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_combine_supports():
    """Test the combine_supports function in the boots_on_mb.py file.
    """
    ##########################
    # Arrange.
    boot_newick = 'boot_newick'
    prob_newick = 'prob_newick'
    combined_figtree_newick = 'combined_figtree_newick'

    ##########################
    # Act.
    #x = combine_supports(boot_newick,
    #		prob_newick,
    #		combined_figtree_newick)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_mbcontre_to_newick_w_probs():
    """Test the mbcontre_to_newick_w_probs function in the boots_on_mb.py file.
    """
    ##########################
    # Arrange.
    intreepath = 'intreepath'
    outtreepath = 'outtreepath'

    ##########################
    # Act.
    #x = mbcontre_to_newick_w_probs(intreepath,
    #		outtreepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_contre_to_newick():
    """Test the contre_to_newick function in the boots_on_mb.py file.
    """
    ##########################
    # Arrange.
    intreepath = 'intreepath'
    outtreepath = 'outtreepath'

    ##########################
    # Act.
    #x = contre_to_newick(intreepath,
    #		outtreepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


