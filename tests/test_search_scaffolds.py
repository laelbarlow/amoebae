#!/usr/bin/env python3
"""PyTest tests for the search_scaffolds.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.
from Bio import SearchIO

from search_scaffolds import \
check_if_two_hsp_ranges_overlap, \
check_if_two_hsps_overlap, \
check_if_hsps_overlapping, \
check_strands, \
whether_fwd_strand, \
check_for_duplicate_ranges, \
sort_left_ranges, \
reduce_left_ranges, \
sort_right_ranges, \
reduce_right_ranges, \
get_proximate_hsp_ranges, \
remove_hsps_with_redun_ranges, \
reduce_to_best_hsps, \
find_missed_hsp_ranges, \
get_best_proximate_hsps, \
get_proximate_hsp_objs, \
get_proximate_hsp_objs2, \
get_hit_seq_record_and_coord, \
get_hit_seq_record_and_coord2, \
get_hit_seq_obj, \
get_tblastn_hit_seq_obj_and_coord, \
recursively_add_clusters, \
get_hsp_clusters, \
get_strand_word, \
get_cluster_range, \
clusters_overlap, \
print_cluster, \
split_tblastn_hits_into_separate_genes, \
get_blastp_hit_seq_obj_and_coord, \
parse_tblastn, \
search_scaffolds


def test_bio_searchio_coordinates():
    """Test whether the Biopython Bio.SearchIO module returns HSP range
    coordinates as python string slices or as the actual positions of the first
    and last residues in the HSP (of the subject sequence). This is important
    for the search_scaffolds module, specifically the
    check_if_hsps_overlapping function (which may be obsolete). 
    """
    ##########################
    # Arrange.

    # Example BLASTP input file is at this path:
    #   data/blastp_query_1_database_1_output.xml
    #
    # The text format output for the same search is here:
    #   data/blastp_query_1_database_1_output.txt

    ##########################
    # Act.

    # Initiate variable to store a SearchIO HSP object.
    hsp_1 = None

    # Open a BLASTP output file (XML format).
    with open('tests/data/blastp_query_1_database_1_output.xml') as infh:

        # Parse the BLASTP output file using SearchIO (from Biopython).
        search_results = SearchIO.parse(infh, 'blast-xml')

        # Get the first HSP from the first hit in the first record.
        for record in search_results:
            for hit in record:
                for hsp in hit:
                    hsp_1 = hsp
                    break
                break
            break

    ##########################
    # Assert.

    # The top hit is the original query, which is 288 residues long. So, the
    # first (only) HSP for the top/first hit should be the length of the
    # original sequence from residue 1 to residue 288. Using the Python string
    # slicing conventions (counting from zero, and not including the last
    # number in the slice), the HSP coordinates would be (0, 288) which
    # corresponds to residues 0 to 287 (still 288 residues). 
    assert hsp_1.hit_range == (0, 288)


def test_check_if_two_hsp_ranges_overlap():  
    """Test the check_if_two_hsp_ranges_overlap function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    input_value1 = [[1, 5], [3, 7]]
    input_value2 = [[1, 5], [7, 9]]
    input_value3 = [[1, 5], [1, 5]]
    input_value4 = [[1, 5], [5, 7]]

    ##########################
    # Act.
    result_value1 = check_if_two_hsp_ranges_overlap(input_value1)
    result_value2 = check_if_two_hsp_ranges_overlap(input_value2)
    result_value3 = check_if_two_hsp_ranges_overlap(input_value3)
    result_value4 = check_if_two_hsp_ranges_overlap(input_value4)

    ##########################
    # Assert.
    assert result_value1 == True
    assert result_value2 == False
    assert result_value3 == True
    assert result_value4 == False





def test_check_for_duplicate_ranges():
    """Test the check_for_duplicate_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    tuple_list1 = [(10,11),(10,11)]
    tuple_list2 = [(0,11),(10,11)]

    ##########################
    # Act.
    result_1 = check_for_duplicate_ranges(tuple_list1)
    result_2 = check_for_duplicate_ranges(tuple_list2)

    ##########################
    # Assert.
    assert result_1 is True
    assert result_2 is False



def test_sort_left_ranges():
    """Test the sort_left_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list_1 = [(1,10),(11,20)]
    expect_1 = [(11,20),(1,10)]

    ##########################
    # Act.
    result_1 = sort_left_ranges(range_list_1)

    ##########################
    # Assert.
    assert result_1 == expect_1



def test_reduce_left_ranges():
    """Test the reduce_left_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.

    # Test for inclusion of ranges within maximum gap.
    sorted_left_ranges_1 = [(51,70),(31,40),(11,20)]
    top_hsp_range_1 = (81,100)
    max_gap_1 = 15
    expect_1 = [(51,70),(31,40),(11,20)]        

    # Test for exclusion of ranges outside maximum gap.
    sorted_left_ranges_2 = [(51,70),(31,40),(11,20)]
    top_hsp_range_2 = (81,100)
    max_gap_2 = 5
    expect_2 = []

    # Test for exclusion of overlapping ranges.
    sorted_left_ranges_3 = [(51,82),(31,52),(11,32)]
    top_hsp_range_3 = (81,100)
    max_gap_3 = 15
    expect_3 = []

    # Test for correct measurement of maximum gap length.
    sorted_left_ranges_4 = [(1,2)]
    top_hsp_range_4 = (3,4)
    max_gap_4 = 0
    expect_4 = [(1,2)]

    ##########################
    # Act.
    result_1 = reduce_left_ranges(sorted_left_ranges_1,
                                  top_hsp_range_1,
                                  max_gap_1)
    result_2 = reduce_left_ranges(sorted_left_ranges_2,
                                  top_hsp_range_2,
                                  max_gap_2)
    result_3 = reduce_left_ranges(sorted_left_ranges_3,
                                  top_hsp_range_3,
                                  max_gap_3)
    result_4 = reduce_left_ranges(sorted_left_ranges_4,
                                  top_hsp_range_4,
                                  max_gap_4)

    ##########################
    # Assert.
    assert result_1 == expect_1
    assert result_2 == expect_2
    assert result_3 == expect_3
    assert result_4 == expect_4



def test_sort_right_ranges():
    """Test the sort_right_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list_1 = [(11,20),(1,10)]
    expect_1 = [(1,10),(11,20)]

    ##########################
    # Act.
    result_1 = sort_right_ranges(range_list_1)

    ##########################
    # Assert.
    assert result_1 == expect_1



def test_reduce_right_ranges():
    """Test the reduce_right_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.

    # Test for inclusion of ranges within maximum gap.
    sorted_right_ranges_1 = [(21,30),(41,60),(71,90)]        
    top_hsp_range_1 = (1,10)
    max_gap_1 = 15
    expect_1 = [(21,30),(41,60),(71,90)]        

    # Test for exclusion of ranges outside maximum gap.
    sorted_right_ranges_2 = [(21,30),(41,60),(71,90)]        
    top_hsp_range_2 = (1,10)
    max_gap_2 = 5
    expect_2 = []

    # Test for exclusion of overlapping ranges.
    sorted_right_ranges_3 = [(21,42),(41,72),(71,90)]        
    top_hsp_range_3 = (1,22)
    max_gap_3 = 15
    expect_3 = []

    # Test for correct measurement of maximum gap length. 
    sorted_right_ranges_4 = [(3,4)]        
    top_hsp_range_4 = (1,2)
    max_gap_4 = 0
    expect_4 = [(3,4)]

    ##########################
    # Act.
    result_1 = reduce_right_ranges(sorted_right_ranges_1,
                                   top_hsp_range_1,
                                   max_gap_1)
    result_2 = reduce_right_ranges(sorted_right_ranges_2,
                                   top_hsp_range_2,
                                   max_gap_2)
    result_3 = reduce_right_ranges(sorted_right_ranges_3,
                                   top_hsp_range_3,
                                   max_gap_3)
    result_4 = reduce_right_ranges(sorted_right_ranges_4,
                                   top_hsp_range_4,
                                   max_gap_4)

    ##########################
    # Assert.
    assert result_1 == expect_1
    assert result_2 == expect_2
    assert result_3 == expect_3
    assert result_4 == expect_4



def test_get_proximate_hsp_ranges():  # ***Incomplete test
    """Test the get_proximate_hsp_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list = "range_list"
    top_hsp_range = "top_hsp_range"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_proximate_hsp_ranges(range_list,
    #		top_hsp_range,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_remove_hsps_with_redun_ranges():  # ***Incomplete test
    """Test the remove_hsps_with_redun_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list = "range_list"
    hsp_objs = "hsp_objs"

    ##########################
    # Act.
    #x = remove_hsps_with_redun_ranges(range_list,
    #		hsp_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_to_best_hsps():  # ***Incomplete test
    """Test the reduce_to_best_hsps function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sorted_hsps = "sorted_hsps"
    max_gap = "max_gap"
    recursion_num = "recursion_num"
    name = "name"
    leftright = "leftright"

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
    a1 =  [(1,10), (20,30), (40,50), (60,70), (80,90), (19,21), (29,41), (49,61), (69,71)]
    b1 =  [        (20,30), (40,50), (60,70)         ]
    c1 =  [   (19,21), (29,41), (49,61), (69,71)     ]

    a2 =  [(1,10), (20,30), (40,50), (60,70), (80,90), (19,21), (29,41), (49,61), (69,71)]
    b2 =  [        (20,30), (40,50), (60,70)         ]
    c2 =  [   (19,21),          (49,61), (69,71)     ]

    a3 =  [(1,10), (20,30), (40,50), (60,70), (80,90), (19,21), (29,41), (49,61), (69,71), (32,35)]
    b3 =  [        (20,30), (40,50), (60,70)         ]
    c3 =  [   (19,21), (29,41), (49,61), (69,71)     ]

    a4 =  [(8466018, 8466147), (8466353, 8466488), (8466585, 8466801),
            (8468034, 8468208), (8468456, 8468645), (8468717, 8468861),
            (8469013, 8469343)]
    b4 =  [(8466018, 8466147), (8466353, 8466488), (8466585, 8466801),
                                (8468456, 8468645), (8468717, 8468861),
            (8469013, 8469343)]
    c4 =  [ ]

    ##########################
    # Act.
    result_value1 = find_missed_hsp_ranges(a1,b1,c1)
    result_value2 = find_missed_hsp_ranges(a2,b2,c2)
    result_value3 = find_missed_hsp_ranges(a3,b3,c3)
    result_value4 = find_missed_hsp_ranges(a4,b4,c4)

    ##########################
    # Assert.
    assert result_value1 == False
    assert result_value2 == True
    assert result_value3 == True
    assert result_value4 == True




