#!/usr/bin/env python3
"""PyTest tests for the nex_to_hmm.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from nex_to_hmm import *



def test_nex_to_hmm():
    """Test the nex_to_hmm function in the nex_to_hmm.py file.
    """
    ##########################
    # Arrange.
    infilepath = 'infilepath'
    outfilepath = 'outfilepath=None'

    ##########################
    # Act.
    #x = nex_to_hmm(infilepath,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_afa_to_hmm():
    """Test the afa_to_hmm function in the nex_to_hmm.py file.
    """
    ##########################
    # Arrange.
    infilepath = 'infilepath'
    outfilepath = 'outfilepath=None'

    ##########################
    # Act.
    #x = afa_to_hmm(infilepath,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


