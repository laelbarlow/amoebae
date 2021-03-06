#!/usr/bin/env python3
"""PyTest tests for the run_searches.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from run_searches import \
get_query_subdir, \
get_query_list_from_file, \
get_db_list_from_file, \
get_out_bash_path, \
get_out_query_list_path, \
get_out_db_list_path, \
determine_search_method, \
search_result_filepath, \
get_out_hmm_path, \
run_any_search, \
run_all_searches


def test_get_query_subdir():  # ***Incomplete test
    """Test the get_query_subdir function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = get_query_subdir(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_query_list_from_file():  # ***Incomplete test
    """Test the get_query_list_from_file function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    infp = "infp"

    ##########################
    # Act.
    #x = get_query_list_from_file(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_db_list_from_file():  # ***Incomplete test
    """Test the get_db_list_from_file function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    infp = "infp"

    ##########################
    # Act.
    #x = get_db_list_from_file(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_out_bash_path():  # ***Incomplete test
    """Test the get_out_bash_path function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = get_out_bash_path(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_out_query_list_path():  # ***Incomplete test
    """Test the get_out_query_list_path function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = get_out_query_list_path(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_out_db_list_path():  # ***Incomplete test
    """Test the get_out_db_list_path function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = get_out_db_list_path(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_determine_search_method():  # ***Incomplete test
    """Test the determine_search_method function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    query_exten = "query_exten"
    db_exten = "db_exten"

    ##########################
    # Act.
    #x = determine_search_method(query_exten,
    #		db_exten)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_search_result_filepath():  # ***Incomplete test
    """Test the search_result_filepath function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    query_filename = "query_filename"
    db_filename = "db_filename"
    dirpath = "dirpath"

    ##########################
    # Act.
    #x = search_result_filepath(query_filename,
    #		db_filename,
    #		dirpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_out_hmm_path():  # ***Incomplete test
    """Test the get_out_hmm_path function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    new_afa_path = "new_afa_path"

    ##########################
    # Act.
    #x = get_out_hmm_path(new_afa_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_run_any_search():  # ***Incomplete test
    """Test the run_any_search function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    queryfile = "queryfile"

    ##########################
    # Act.
    #x = run_any_search(queryfile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_run_all_searches():  # ***Incomplete test
    """Test the run_all_searches function in the run_searches.py file.
    """
    ##########################
    # Arrange.
    query_file_list = "query_file_list"

    ##########################
    # Act.
    #x = run_all_searches(query_file_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


