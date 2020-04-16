#!/usr/bin/env python3
"""PyTest tests for the test_dacksify_pos_hmmer_hits.py module.
"""
from test_dacksify_pos_hmmer_hits.py import *



def test_get_rblast_string():
    """Test the get_rblast_string function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    relev_blast_res = 'relev_blast_res'
    red_acc_file = 'red_acc_file'

    ##########################
    # Act.
    #x = get_rblast_string(relev_blast_res,
    #		red_acc_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_find_relev_blast():
    """Test the find_relev_blast function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    blast_res = 'blast_res'
    red_acc_file = 'red_acc_file'

    ##########################
    # Act.
    #x = find_relev_blast(blast_res,
    #		red_acc_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_acc_for_prot_name():
    """Test the get_acc_for_prot_name function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    prot_name = 'prot_name'
    key_file = 'key_file'

    ##########################
    # Act.
    #x = get_acc_for_prot_name(prot_name,
    #		key_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_rblast_list():
    """Test the get_rblast_list function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    linestring = 'linestring'

    ##########################
    # Act.
    #x = get_rblast_list(linestring)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_is_redun_acc():
    """Test the is_redun_acc function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    redun_acc_file = 'redun_acc_file'
    acc = 'acc'

    ##########################
    # Act.
    #x = is_redun_acc(redun_acc_file,
    #		acc)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_find_unabbrev_name():
    """Test the find_unabbrev_name function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    abbrev = 'abbrev'
    tablefilepath = 'tablefilepath'

    ##########################
    # Act.
    #x = find_unabbrev_name(abbrev,
    #		tablefilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_red_acc_listx():
    """Test the get_red_acc_listx function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    infilepath = 'infilepath'

    ##########################
    # Act.
    #x = get_red_acc_listx(infilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_query_retrieved():
    """Test the query_retrieved function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    hmmer_string = 'hmmer_string'
    red_acc_file = 'red_acc_file'

    ##########################
    # Act.
    #x = query_retrieved(hmmer_string,
    #		red_acc_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_dacksify():
    """Test the dacksify function in the dacksify_pos_hmmer_hits.py file.
    """
    ##########################
    # Arrange.
    infp = 'infp'
    outfp = 'outfp'
    tablefilepath = 'tablefilepath'
    red_acc_file = 'red_acc_file'
    evaldiffset = 'evaldiffset'

    ##########################
    # Act.
    #x = dacksify(infp,
    #		outfp,
    #		tablefilepath,
    #		red_acc_file,
    #		evaldiffset)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


