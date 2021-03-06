#!/usr/bin/env python3
"""PyTest tests for the hmmscan.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from hmmscan import \
all_nex_to_hmm, \
get_hmm_paths, \
cat_hmm, \
call_hmmpress, \
make_hmm_db, \
determine_if_is_hmmdb


def test_all_nex_to_hmm():  # ***Incomplete test
    """Test the all_nex_to_hmm function in the hmmscan.py file.
    """
    ##########################
    # Arrange.
    indirpath = "indirpath"

    ##########################
    # Act.
    #x = all_nex_to_hmm(indirpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hmm_paths():  # ***Incomplete test
    """Test the get_hmm_paths function in the hmmscan.py file.
    """
    ##########################
    # Arrange.
    indirpath = "indirpath"

    ##########################
    # Act.
    #x = get_hmm_paths(indirpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_cat_hmm():  # ***Incomplete test
    """Test the cat_hmm function in the hmmscan.py file.
    """
    ##########################
    # Arrange.
    hmm_paths = "hmm_paths"
    outfilepath = "outfilepath"

    ##########################
    # Act.
    #x = cat_hmm(hmm_paths,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_call_hmmpress():  # ***Incomplete test
    """Test the call_hmmpress function in the hmmscan.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"

    ##########################
    # Act.
    #x = call_hmmpress(infilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_make_hmm_db():  # ***Incomplete test
    """Test the make_hmm_db function in the hmmscan.py file.
    """
    ##########################
    # Arrange.
    indirpath = "indirpath"
    hmmdbname = "hmmdbname=None"
    outdirpath = "outdirpath=None"

    ##########################
    # Act.
    #x = make_hmm_db(indirpath,
    #		hmmdbname,
    #		outdirpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_determine_if_is_hmmdb():  # ***Incomplete test
    """Test the determine_if_is_hmmdb function in the hmmscan.py file.
    """
    ##########################
    # Arrange.
    infp = "infp"

    ##########################
    # Act.
    #x = determine_if_is_hmmdb(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


