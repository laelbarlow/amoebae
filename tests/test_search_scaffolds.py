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


def test_check_if_two_hsps_overlap(): # ***Incomplete test
    """Test the check_if_two_hsps_overlap function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    lists = "lists"

    ##########################
    # Act.
    #x = check_if_two_hsps_overlap(lists)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


def test_check_if_hsps_overlapping():  # ***Incomplete test
    """Test the check_if_hsps_overlapping function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = "hsp_objs"

    ##########################
    # Act.
    #x = check_if_hsps_overlapping(hsp_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_strands():  # ***Incomplete test
    """Test the check_strands function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = "hsp_objs"
    fwd_strand = "fwd_strand"

    ##########################
    # Act.
    #x = check_strands(hsp_objs,
    #		fwd_strand)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_whether_fwd_strand():  # ***Incomplete test
    """Test the whether_fwd_strand function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp = "hsp"

    ##########################
    # Act.
    #x = whether_fwd_strand(hsp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_for_duplicate_ranges():  # ***Incomplete test
    """Test the check_for_duplicate_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    range_list = "range_list"

    ##########################
    # Act.
    #x = check_for_duplicate_ranges(range_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_sort_left_ranges():  # ***Incomplete test
    """Test the sort_left_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    left_ranges = "left_ranges"

    ##########################
    # Act.
    #x = sort_left_ranges(left_ranges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_left_ranges():  # ***Incomplete test
    """Test the reduce_left_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sorted_left_ranges = "sorted_left_ranges"
    top_hsp_range = "top_hsp_range"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = reduce_left_ranges(sorted_left_ranges,
    #		top_hsp_range,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_sort_right_ranges():  # ***Incomplete test
    """Test the sort_right_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    right_ranges = "right_ranges"

    ##########################
    # Act.
    #x = sort_right_ranges(right_ranges)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_reduce_right_ranges():  # ***Incomplete test
    """Test the reduce_right_ranges function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sorted_right_ranges = "sorted_right_ranges"
    top_hsp_range = "top_hsp_range"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = reduce_right_ranges(sorted_right_ranges,
    #		top_hsp_range,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



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


def test_get_best_proximate_hsps():  # ***Incomplete test
    """Test the get_best_proximate_hsps function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = "hsp_objs"
    top_hsp = "top_hsp"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_best_proximate_hsps(hsp_objs,
    #		top_hsp,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_proximate_hsp_objs():  # ***Incomplete test
    """Test the get_proximate_hsp_objs function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = "hsp_objs"
    top_hsp = "top_hsp"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_proximate_hsp_objs(hsp_objs,
    #		top_hsp,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_proximate_hsp_objs2():  # ***Incomplete test
    """Test the get_proximate_hsp_objs2 function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_objs = "hsp_objs"
    top_hsp = "top_hsp"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_proximate_hsp_objs2(hsp_objs,
    #		top_hsp,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_seq_record_and_coord():  # ***Incomplete test
    """Test the get_hit_seq_record_and_coord function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = "hit"
    proximate_hsp_objs = "proximate_hsp_objs"

    ##########################
    # Act.
    #x = get_hit_seq_record_and_coord(hit,
    #		proximate_hsp_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_seq_record_and_coord2():  # ***Incomplete test
    """Test the get_hit_seq_record_and_coord2 function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    search_result_path = "search_result_path"

    ##########################
    # Act.
    #x = get_hit_seq_record_and_coord2(search_result_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_seq_obj():  # ***Incomplete test
    """Test the get_hit_seq_obj function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = "hit"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_hit_seq_obj(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_tblastn_hit_seq_obj_and_coord():  # ***Incomplete test
    """Test the get_tblastn_hit_seq_obj_and_coord function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = "hit"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_tblastn_hit_seq_obj_and_coord(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_recursively_add_clusters():  # ***Incomplete test
    """Test the recursively_add_clusters function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_list = "hsp_list"
    hsp_clusters = "hsp_clusters"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = recursively_add_clusters(hsp_list,
    #		hsp_clusters,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hsp_clusters():  # ***Incomplete test
    """Test the get_hsp_clusters function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = "hit"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_hsp_clusters(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_strand_word():  # ***Incomplete test
    """Test the get_strand_word function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_hit_frame = "hsp_hit_frame"

    ##########################
    # Act.
    #x = get_strand_word(hsp_hit_frame)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_cluster_range():  # ***Incomplete test
    """Test the get_cluster_range function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hsp_list = "hsp_list"

    ##########################
    # Act.
    #x = get_cluster_range(hsp_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_clusters_overlap():  # ***Incomplete test
    """Test the clusters_overlap function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    cluster1 = "cluster1"
    cluster2 = "cluster2"

    ##########################
    # Act.
    #x = clusters_overlap(cluster1,
    #		cluster2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_print_cluster():  # ***Incomplete test
    """Test the print_cluster function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    clusterplus = "clusterplus"
    hit_num = "hit_num"
    cluster_num = "cluster_num"
    num_dots = "num_dots"
    startend = "startend=None"

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



def test_split_tblastn_hits_into_separate_genes():  # ***Incomplete test
    """Test the split_tblastn_hits_into_separate_genes function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    query_res_obj = "query_res_obj"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = split_tblastn_hits_into_separate_genes(query_res_obj,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_blastp_hit_seq_obj_and_coord():  # ***Incomplete test
    """Test the get_blastp_hit_seq_obj_and_coord function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit = "hit"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = get_blastp_hit_seq_obj_and_coord(hit,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_parse_tblastn():  # ***Incomplete test
    """Test the parse_tblastn function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    tblastn_out_path = "tblastn_out_path"
    seq_out_path = "seq_out_path"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = parse_tblastn(tblastn_out_path,
    #		seq_out_path,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_search_scaffolds():  # ***Incomplete test
    """Test the search_scaffolds function in the search_scaffolds.py file.
    """
    ##########################
    # Arrange.
    query_file = "query_file"
    subject_db = "subject_db"
    max_gap = "max_gap"

    ##########################
    # Act.
    #x = search_scaffolds(query_file,
    #		subject_db,
    #		max_gap)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


