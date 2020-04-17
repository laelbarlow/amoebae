#!/usr/bin/env python3
"""PyTest tests for the search_scaffolds.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from search_scaffolds import *



def test_check_if_two_hsp_ranges_overlap():
    """Test the check_if_two_hsp_ranges_overlap function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    lists = 'lists'

    ##########################
    # Act.
    #x = check_if_two_hsp_ranges_overlap(lists)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_if_two_hsps_overlap():
    """Test the check_if_two_hsps_overlap function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_obj1 = 'hsp_obj1'
    hsp_obj2 = 'hsp_obj2'

    ##########################
    # Act.
    #x = check_if_two_hsps_overlap(hsp_obj1,
    #		hsp_obj2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_if_hsps_overlapping():
    """Test the check_if_hsps_overlapping function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = 'hsp_objs'

    ##########################
    # Act.
    #x = check_if_hsps_overlapping(hsp_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_strands():
    """Test the check_strands function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = 'hsp_objs'
    fwd_strand = 'fwd_strand'

    ##########################
    # Act.
    #x = check_strands(hsp_objs,
    #		fwd_strand)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_whether_fwd_strand():
    """Test the whether_fwd_strand function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp = 'hsp'

    ##########################
    # Act.
    #x = whether_fwd_strand(hsp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_for_duplicate_ranges():
    """Test the check_for_duplicate_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list = 'range_list'

    ##########################
    # Act.
    #x = check_for_duplicate_ranges(range_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_sort_left_ranges():
    """Test the sort_left_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    left_ranges = 'left_ranges'

    ##########################
    # Act.
    #x = sort_left_ranges(left_ranges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_left_ranges():
    """Test the reduce_left_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sorted_left_ranges = 'sorted_left_ranges'
    top_hsp_range = 'top_hsp_range'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = reduce_left_ranges(sorted_left_ranges,
    #		top_hsp_range,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_sort_right_ranges():
    """Test the sort_right_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    right_ranges = 'right_ranges'

    ##########################
    # Act.
    #x = sort_right_ranges(right_ranges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_right_ranges():
    """Test the reduce_right_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sorted_right_ranges = 'sorted_right_ranges'
    top_hsp_range = 'top_hsp_range'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = reduce_right_ranges(sorted_right_ranges,
    #		top_hsp_range,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_proximate_hsp_ranges():
    """Test the get_proximate_hsp_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list = 'range_list'
    top_hsp_range = 'top_hsp_range'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_proximate_hsp_ranges(range_list,
    #		top_hsp_range,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_remove_hsps_with_redun_ranges():
    """Test the remove_hsps_with_redun_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list = 'range_list'
    hsp_objs = 'hsp_objs'

    ##########################
    # Act.
    #x = remove_hsps_with_redun_ranges(range_list,
    #		hsp_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_to_best_hsps():
    """Test the reduce_to_best_hsps function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sorted_hsps = 'sorted_hsps'
    max_gap = 'max_gap'
    recursion_num = 'recursion_num'
    name = 'name'
    leftright = 'leftright'

    ##########################
    # Act.
    #x = reduce_to_best_hsps(sorted_hsps,
    #		max_gap,
    #		recursion_num,
    #		name,
    #		leftright)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_find_missed_hsp_ranges():
    """Test the find_missed_hsp_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    complete_ranges = 'complete_ranges'

    ##########################
    # Act.
    #x = find_missed_hsp_ranges(complete_ranges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_best_proximate_hsps():
    """Test the get_best_proximate_hsps function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = 'hsp_objs'
    top_hsp = 'top_hsp'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_best_proximate_hsps(hsp_objs,
    #		top_hsp,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_proximate_hsp_objs():
    """Test the get_proximate_hsp_objs function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = 'hsp_objs'
    top_hsp = 'top_hsp'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_proximate_hsp_objs(hsp_objs,
    #		top_hsp,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_proximate_hsp_objs2():
    """Test the get_proximate_hsp_objs2 function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = 'hsp_objs'
    top_hsp = 'top_hsp'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_proximate_hsp_objs2(hsp_objs,
    #		top_hsp,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_seq_record_and_coord():
    """Test the get_hit_seq_record_and_coord function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = 'hit'
    proximate_hsp_objs = 'proximate_hsp_objs'

    ##########################
    # Act.
    #x = get_hit_seq_record_and_coord(hit,
    #		proximate_hsp_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_seq_record_and_coord2():
    """Test the get_hit_seq_record_and_coord2 function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    search_result_path = 'search_result_path'

    ##########################
    # Act.
    #x = get_hit_seq_record_and_coord2(search_result_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_seq_obj():
    """Test the get_hit_seq_obj function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = 'hit'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_hit_seq_obj(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_tblastn_hit_seq_obj_and_coord():
    """Test the get_tblastn_hit_seq_obj_and_coord function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = 'hit'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_tblastn_hit_seq_obj_and_coord(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_recursively_add_clusters():
    """Test the recursively_add_clusters function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_list = 'hsp_list'
    hsp_clusters = 'hsp_clusters'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = recursively_add_clusters(hsp_list,
    #		hsp_clusters,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hsp_clusters():
    """Test the get_hsp_clusters function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = 'hit'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_hsp_clusters(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_strand_word():
    """Test the get_strand_word function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_hit_frame = 'hsp_hit_frame'

    ##########################
    # Act.
    #x = get_strand_word(hsp_hit_frame)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_cluster_range():
    """Test the get_cluster_range function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_list = 'hsp_list'

    ##########################
    # Act.
    #x = get_cluster_range(hsp_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_clusters_overlap():
    """Test the clusters_overlap function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    cluster1 = 'cluster1'
    cluster2 = 'cluster2'

    ##########################
    # Act.
    #x = clusters_overlap(cluster1,
    #		cluster2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_print_cluster():
    """Test the print_cluster function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    clusterplus = 'clusterplus'
    hit_num = 'hit_num'
    cluster_num = 'cluster_num'
    num_dots = 'num_dots'
    startend = 'startend=None'

    ##########################
    # Act.
    #x = print_cluster(clusterplus,
    #		hit_num,
    #		cluster_num,
    #		num_dots,
    #		startend)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_split_tblastn_hits_into_separate_genes():
    """Test the split_tblastn_hits_into_separate_genes function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    query_res_obj = 'query_res_obj'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = split_tblastn_hits_into_separate_genes(query_res_obj,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_blastp_hit_seq_obj_and_coord():
    """Test the get_blastp_hit_seq_obj_and_coord function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = 'hit'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = get_blastp_hit_seq_obj_and_coord(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_parse_tblastn():
    """Test the parse_tblastn function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    tblastn_out_path = 'tblastn_out_path'
    seq_out_path = 'seq_out_path'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = parse_tblastn(tblastn_out_path,
    #		seq_out_path,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_search_scaffolds():
    """Test the search_scaffolds function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    query_file = 'query_file'
    subject_db = 'subject_db'
    max_gap = 'max_gap'

    ##########################
    # Act.
    #x = search_scaffolds(query_file,
    #		subject_db,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


