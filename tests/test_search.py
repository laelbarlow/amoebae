#!/usr/bin/env python3
"""PyTest tests for the search.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from search import \
intersect, \
get_search_software_version, \
get_search_output_file_type, \
fwd_search_results_to_csv, \
get_corr_evalue, \
get_redun_hits_in_dbs, \
get_query_len, \
get_subsequences_matching_query, \
get_hit_obj_for_hsp_cluster, \
get_rows_for_fwd_srch_df, \
write_fwd_srch_res_to_csv, \
get_rev_query_path_from_fwd_srch_row, \
get_rev_queries, \
get_csv_with_rev_path, \
get_redun_hit_dict, \
get_evaldiff, \
write_rev_srch_res_to_csv, \
write_interp_csv, \
write_fwd_srch_interp_csv, \
write_redun_hit_interp_csv


def test_intersect():  # ***Incomplete test
    """Test the intersect function in the search.py file.
    """
    ##########################
    # Arrange.
    a = "a"
    b = "b"

    ##########################
    # Act.
    #x = intersect(a,
    #		b)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_search_software_version():  # ***Incomplete test
    """Test the get_search_software_version function in the search.py file.
    """
    ##########################
    # Arrange.
    method = "method"

    ##########################
    # Act.
    #x = get_search_software_version(method)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_search_output_file_type():  # ***Incomplete test
    """Test the get_search_output_file_type function in the search.py file.
    """
    ##########################
    # Arrange.
    filepath = "filepath"

    ##########################
    # Act.
    #x = get_search_output_file_type(filepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_fwd_search_results_to_csv():  # ***Incomplete test
    """Test the fwd_search_results_to_csv function in the search.py file.
    """
    ##########################
    # Arrange.
    search_outdir = "search_outdir"
    csv_file = "csv_file"

    ##########################
    # Act.
    #x = fwd_search_results_to_csv(search_outdir,
    #		csv_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_corr_evalue():  # ***Incomplete test
    """Test the get_corr_evalue function in the search.py file.
    """
    ##########################
    # Arrange.
    evalue = "evalue"

    ##########################
    # Act.
    #x = get_corr_evalue(evalue)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_redun_hits_in_dbs():  # ***Incomplete test
    """Test the get_redun_hits_in_dbs function in the search.py file.
    """
    ##########################
    # Arrange.
    query_title = "query_title"

    ##########################
    # Act.
    #x = get_redun_hits_in_dbs(query_title)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_query_len():  # ***Incomplete test
    """Test the get_query_len function in the search.py file.
    """
    ##########################
    # Arrange.
    query_filename = "query_filename"

    ##########################
    # Act.
    #x = get_query_len(query_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_subsequences_matching_query():  # ***Incomplete test
    """Test the get_subsequences_matching_query function in the search.py file.
    """
    ##########################
    # Arrange.
    hit_obj = "hit_obj"

    ##########################
    # Act.
    #x = get_subsequences_matching_query(hit_obj)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_obj_for_hsp_cluster():  # ***Incomplete test
    """Test the get_hit_obj_for_hsp_cluster function in the search.py file.
    """
    ##########################
    # Arrange.
    query_res_obj = "query_res_obj"
    cluster = "cluster"

    ##########################
    # Act.
    #x = get_hit_obj_for_hsp_cluster(query_res_obj,
    #		cluster)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_rows_for_fwd_srch_df():  # ***Incomplete test
    """Test the get_rows_for_fwd_srch_df function in the search.py file.
    """
    ##########################
    # Arrange.
    df = "df"

    ##########################
    # Act.
    #x = get_rows_for_fwd_srch_df(df)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_fwd_srch_res_to_csv():  # ***Incomplete test
    """Test the write_fwd_srch_res_to_csv function in the search.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = write_fwd_srch_res_to_csv(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_rev_query_path_from_fwd_srch_row():  # ***Incomplete test
    """Test the get_rev_query_path_from_fwd_srch_row function in the search.py file.
    """
    ##########################
    # Arrange.
    query_subdir = "query_subdir"
    df_row = "df_row"
    aasubseq = "aasubseq=False"

    ##########################
    # Act.
    #x = get_rev_query_path_from_fwd_srch_row(query_subdir,
    #		df_row,
    #		aasubseq)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_rev_queries():  # ***Incomplete test
    """Test the get_rev_queries function in the search.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    query_subdir = "query_subdir"
    aasubseq = "aasubseq"
    nafullseq = "nafullseq"

    ##########################
    # Act.
    #x = get_rev_queries(csv_file,
    #		query_subdir,
    #		aasubseq,
    #		nafullseq)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_csv_with_rev_path():  # ***Incomplete test
    """Test the get_csv_with_rev_path function in the search.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"

    ##########################
    # Act.
    #x = get_csv_with_rev_path(csv_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_redun_hit_dict():  # ***Incomplete test
    """Test the get_redun_hit_dict function in the search.py file.
    """
    ##########################
    # Arrange.
    redun_hit_csv = "redun_hit_csv"

    ##########################
    # Act.
    #x = get_redun_hit_dict(redun_hit_csv)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_evaldiff():  # ***Incomplete test
    """Test the get_evaldiff function in the search.py file.
    """
    ##########################
    # Arrange.
    evalue1 = "evalue1"
    evalue2 = "evalue2"

    ##########################
    # Act.
    #x = get_evaldiff(evalue1,
    #		evalue2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_rev_srch_res_to_csv():  # ***Incomplete test
    """Test the write_rev_srch_res_to_csv function in the search.py file.
    """
    ##########################
    # Arrange.
    rev_srch_id = "rev_srch_id"

    ##########################
    # Act.
    #x = write_rev_srch_res_to_csv(rev_srch_id)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_interp_csv():  # ***Incomplete test
    """Test the write_interp_csv function in the search.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    outfp = "outfp"
    fwd_evalue_cutoff = "fwd_evalue_cutoff"
    rev_evalue_cutoff = "rev_evalue_cutoff"

    ##########################
    # Act.
    #x = write_interp_csv(csv_file,
    #		outfp,
    #		fwd_evalue_cutoff,
    #		rev_evalue_cutoff)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_fwd_srch_interp_csv():  # ***Incomplete test
    """Test the write_fwd_srch_interp_csv function in the search.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    outfp = "outfp"
    score_cutoff = "score_cutoff"

    ##########################
    # Act.
    #x = write_fwd_srch_interp_csv(csv_file,
    #		outfp,
    #		score_cutoff)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_redun_hit_interp_csv():  # ***Incomplete test
    """Test the write_redun_hit_interp_csv function in the search.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    outfp = "outfp"

    ##########################
    # Act.
    #x = write_redun_hit_interp_csv(csv_file,
    #		outfp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


