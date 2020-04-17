#!/usr/bin/env python3
"""PyTest tests for the run_exonerate.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from run_exonerate import *



def test_get_subseq_from_nucl():
    """Test the get_subseq_from_nucl function in the run_exonerate.py file.
    """
    ##########################
    # Arrange.
    subj_seq_id = 'subj_seq_id'

    ##########################
    # Act.
    #x = get_subseq_from_nucl(subj_seq_id)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_run_exonerate_as_subprocess():
    """Test the run_exonerate_as_subprocess function in the run_exonerate.py file.
    """
    ##########################
    # Arrange.
    query_prot_faa = 'query_prot_faa'

    ##########################
    # Act.
    #x = run_exonerate_as_subprocess(query_prot_faa)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


