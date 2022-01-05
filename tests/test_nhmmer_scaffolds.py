#!/usr/bin/env python3
"""PyTest tests for the nhmmer_scaffolds.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))

from nhmmer_scaffolds import \
frameify_seq


def test_frameify_seq():
    """Test the frameify_seq function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    inseq1 = 'ATGATGATG'
    hmmfrom1 = '1'
    expect1 = 'ATGATGATG'

    inseq2 = 'TGATGATG'
    hmmfrom2 = '2'
    expect2 = 'ATGATG'

    inseq3 = 'GATGATG'
    hmmfrom3 = '3'
    expect3 = 'ATGATG'

    ##########################
    # Act.
    result1 = frameify_seq(inseq1, hmmfrom1)
    result2 = frameify_seq(inseq2, hmmfrom2)
    result3 = frameify_seq(inseq3, hmmfrom3)

    ##########################
    # Assert.
    assert result1 == expect1 
    assert result2 == expect2
    assert result3 == expect3


