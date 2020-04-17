#!/usr/bin/env python3
"""PyTest tests for the replace_seqs.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from replace_seqs import *



def test_afa_to_fa():
    """Test the afa_to_fa function in the replace_seqs.py file.
    """
    ##########################
    # Arrange.
    infilepath = 'infilepath'
    outfilepath = 'outfilepath'

    ##########################
    # Act.
    #x = afa_to_fa(infilepath,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_replace_seqs_in_alignment_with_seqs_from_fasta():
    """Test the replace_seqs_in_alignment_with_seqs_from_fasta function in the replace_seqs.py file.
    """
    ##########################
    # Arrange.
    alignment = 'alignment'
    fasta = 'fasta=None'

    ##########################
    # Act.
    #x = replace_seqs_in_alignment_with_seqs_from_fasta(alignment,
    #		fasta)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


