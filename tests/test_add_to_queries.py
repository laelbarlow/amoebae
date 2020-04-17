#!/usr/bin/env python3
"""PyTest tests for the add_to_queries.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from add_to_queries import \
get_query_file_type, \
get_hmm_datatype, \
get_mod_query_path, \
is_single_fasta, \
fasta_seqs_all_same_len, \
update_query_csv


def test_get_query_file_type():  # ***Incomplete test
    """Test the get_query_file_type function in the add_to_queries.py file.
    """
    ##########################
    # Arrange.
    query_file = "query_file"

    ##########################
    # Act.
    #x = get_query_file_type(query_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hmm_datatype():  # ***Incomplete test
    """Test the get_hmm_datatype function in the add_to_queries.py file.
    """
    ##########################
    # Arrange.
    query_file = "query_file"

    ##########################
    # Act.
    #x = get_hmm_datatype(query_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_mod_query_path():  # ***Incomplete test
    """Test the get_mod_query_path function in the add_to_queries.py file.
    """
    ##########################
    # Arrange.
    query_file = "query_file"
    filetype = "filetype"
    datatype = "datatype"
    query_dir = "query_dir"

    ##########################
    # Act.
    #x = get_mod_query_path(query_file,
    #		filetype,
    #		datatype,
    #		query_dir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_is_single_fasta():  # ***Incomplete test
    """Test the is_single_fasta function in the add_to_queries.py file.
    """
    ##########################
    # Arrange.
    query_file = "query_file"

    ##########################
    # Act.
    #x = is_single_fasta(query_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_fasta_seqs_all_same_len():  # ***Incomplete test
    """Test the fasta_seqs_all_same_len function in the add_to_queries.py file.
    """
    ##########################
    # Arrange.
    query_file = "query_file"

    ##########################
    # Act.
    #x = fasta_seqs_all_same_len(query_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_update_query_csv():  # ***Incomplete test
    """Test the update_query_csv function in the add_to_queries.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    mod_query_path = "mod_query_path"
    datatype = "datatype"

    ##########################
    # Act.
    #x = update_query_csv(csv_file,
    #		mod_query_path,
    #		datatype)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


