#!/usr/bin/env python3
"""PyTest tests for the trim_nex.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from trim_nex import \
trim_alignment, \
trim_nex


def test_trim_alignment():  # ***Incomplete test
    """Test the trim_alignment function in the trim_nex.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"

    ##########################
    # Act.
    #x = trim_alignment(alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_trim_nex():  # ***Incomplete test
    """Test the trim_nex function in the trim_nex.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"
    outfilename = "outfilename"

    ##########################
    # Act.
    #x = trim_nex(infilepath,
    #		outfilename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


