#!/usr/bin/env python3
"""PyTest tests for the find_pos_hmmer_hits.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from find_pos_hmmer_hits import \
get_red_acc_list, \
write_pos_seqs, \
get_csv_name, \
find_pos_hmmer_hitsx


def test_get_red_acc_list():  # ***Incomplete test
    """Test the get_red_acc_list function in the find_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"

    ##########################
    # Act.
    #x = get_red_acc_list(infilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_pos_seqs():  # ***Incomplete test
    """Test the write_pos_seqs function in the find_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"
    dbdirpath = "dbdirpath"
    outfilepath = "outfilepath"
    prot_name = "prot_name=None"

    ##########################
    # Act.
    #x = write_pos_seqs(infilepath,
    #		dbdirpath,
    #		outfilepath,
    #		prot_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_csv_name():  # ***Incomplete test
    """Test the get_csv_name function in the find_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"
    cur_time = "cur_time"

    ##########################
    # Act.
    #x = get_csv_name(infilepath,
    #		cur_time)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_find_pos_hmmer_hitsx():  # ***Incomplete test
    """Test the find_pos_hmmer_hitsx function in the find_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    infilepath1 = "infilepath1"
    infilepath2 = "infilepath2"
    fwdeval = "fwdeval"
    reveval = "reveval"

    ##########################
    # Act.
    #x = find_pos_hmmer_hitsx(infilepath1,
    #		infilepath2,
    #		fwdeval,
    #		reveval)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


