#!/usr/bin/env python3
"""PyTest tests for the align_to_profile_iter.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from align_to_profile_iter import *



def test_split_fasta():
    """Test the split_fasta function in the align_to_profile_iter.py file.
    """
    ##########################
    # Arrange.
    infilepath = 'infilepath'
    temp_subdirname = 'temp_subdirname'

    ##########################
    # Act.
    #x = split_fasta(infilepath,
    #		temp_subdirname)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_align_iteratively():
    """Test the align_iteratively function in the align_to_profile_iter.py file.
    """
    ##########################
    # Arrange.
    inalignfilepath = 'inalignfilepath'
    temp_subdirname = 'temp_subdirname'
    outfilepath = 'outfilepath'

    ##########################
    # Act.
    #x = align_iteratively(inalignfilepath,
    #		temp_subdirname,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_do_align_iteratively():
    """Test the do_align_iteratively function in the align_to_profile_iter.py file.
    """
    ##########################
    # Arrange.
    innexpath = 'innexpath'
    infapath = 'infapath'
    outnexpath = 'outnexpath=None'

    ##########################
    # Act.
    #x = do_align_iteratively(innexpath,
    #		infapath,
    #		outnexpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


