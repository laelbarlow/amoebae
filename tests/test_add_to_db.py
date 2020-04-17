#!/usr/bin/env python3
"""PyTest tests for the add_to_db.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from add_to_db import *



def test_make_easel_index():
    """Test the make_easel_index function in the add_to_db.py file.
    """
    ##########################
    # Arrange.
    outfp = 'outfp'

    ##########################
    # Act.
    #x = make_easel_index(outfp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_convert_headers():
    """Test the convert_headers function in the add_to_db.py file.
    """
    ##########################
    # Arrange.
    infp = 'infp'
    outfp = 'outfp'
    split_char = 'split_char=' ''
    split_pos = 'split_pos='0''

    ##########################
    # Act.
    #x = convert_headers(infp,
    #		outfp,
    #		split_char,
    #		split_pos)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_make_blast_db():
    """Test the make_blast_db function in the add_to_db.py file.
    """
    ##########################
    # Arrange.
    infp = 'infp'

    ##########################
    # Act.
    #x = make_blast_db(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_corr_fasta_exten():
    """Test the get_corr_fasta_exten function in the add_to_db.py file.
    """
    ##########################
    # Arrange.
    infp = 'infp'

    ##########################
    # Act.
    #x = get_corr_fasta_exten(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_dbtype_from_file_exten():
    """Test the get_dbtype_from_file_exten function in the add_to_db.py file.
    """
    ##########################
    # Arrange.
    infp = 'infp'

    ##########################
    # Act.
    #x = get_dbtype_from_file_exten(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_update_csv():
    """Test the update_csv function in the add_to_db.py file.
    """
    ##########################
    # Arrange.
    outfp = 'outfp'
    csv_file = 'csv_file'

    ##########################
    # Act.
    #x = update_csv(outfp,
    #		csv_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


