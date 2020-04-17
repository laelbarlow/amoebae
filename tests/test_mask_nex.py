#!/usr/bin/env python3
"""PyTest tests for the mask_nex.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from mask_nex import \
apply_mask_criteria, \
mask_alignment, \
mask_nex


def test_apply_mask_criteria():  # ***Incomplete test
    """Test the apply_mask_criteria function in the mask_nex.py file.
    """
    ##########################
    # Arrange.
    column = "column"

    ##########################
    # Act.
    #x = apply_mask_criteria(column)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_mask_alignment():  # ***Incomplete test
    """Test the mask_alignment function in the mask_nex.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"

    ##########################
    # Act.
    #x = mask_alignment(alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_mask_nex():  # ***Incomplete test
    """Test the mask_nex function in the mask_nex.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"
    outfilepath = "outfilepath=None"

    ##########################
    # Act.
    #x = mask_nex(infilepath,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


