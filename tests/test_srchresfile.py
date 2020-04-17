#!/usr/bin/env python3
"""PyTest tests for the srchresfile.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from srchresfile import *



def test_get_srch_file_info():
    """Test the get_srch_file_info function in the srchresfile.py file.
    """
    ##########################
    # Arrange.
    search_result_path = 'search_result_path'

    ##########################
    # Act.
    #x = get_srch_file_info(search_result_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hmmer_hit_seq_coord():
    """Test the get_hmmer_hit_seq_coord function in the srchresfile.py file.
    """
    ##########################
    # Arrange.
    searchio_hit_obj = 'searchio_hit_obj'
    db_file = 'db_file'
    extra = 'extra=0'

    ##########################
    # Act.
    #x = get_hmmer_hit_seq_coord(searchio_hit_obj,
    #		db_file,
    #		extra)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


