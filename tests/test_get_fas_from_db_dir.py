#!/usr/bin/env python3
"""PyTest tests for the get_fas_from_db_dir.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from get_fas_from_db_dir import \
get_seq_obj_from_db_fasta, \
get_fas_from_db_dir


def test_get_seq_obj_from_db_fasta():  # ***Incomplete test
    """Test the get_seq_obj_from_db_fasta function in the get_fas_from_db_dir.py file.
    """
    ##########################
    # Arrange.
    acc_list = "acc_list"
    fa_path = "fa_path"

    ##########################
    # Act.
    #x = get_seq_obj_from_db_fasta(acc_list,
    #		fa_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_fas_from_db_dir():  # ***Incomplete test
    """Test the get_fas_from_db_dir function in the get_fas_from_db_dir.py file.
    """
    ##########################
    # Arrange.
    db_name = "db_name"
    acc_list = "acc_list"
    dbdirpath = "dbdirpath"
    prot_name = "prot_name=None"

    ##########################
    # Act.
    #x = get_fas_from_db_dir(db_name,
    #		acc_list,
    #		dbdirpath,
    #		prot_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


