#!/usr/bin/env python3
"""PyTest tests for the srchresfile.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from srchresfile import \
get_srch_file_info, \
get_hmmer_hit_seq_coord


def test_srchresfile_class_methods(): # ***Incomplete test
    """Test the methods in the SrchResFile class that are used for parsing
    search result files.
    """
    pass


def test_get_srch_file_info():  # ***Incomplete test
    """Test the get_srch_file_info function in the srchresfile.py file.
    """
    ##########################
    # Arrange.
    search_result_path = "search_result_path"

    ##########################
    # Act.
    #x = get_srch_file_info(search_result_path)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


#command_line_list = sys.argv
##blastpfp = str(command_line_list[1])
##hmmsearchfp = str(command_line_list[1])
#
#num = 0
#for i in command_line_list[1:]:
#    num += 1
#    hit_rank = 0
#    print('\nSearch file ' + str(num))
#    parsed_file_obj = SrchResFile(i)
#    print(os.path.basename(parsed_file_obj.filepath))
#    print('program')
#    print(parsed_file_obj.program)
#    print('version')
#    print(parsed_file_obj.version)
#    print('hit rank')
#    print(hit_rank)
#    print('hit id')
#    print(parsed_file_obj.hit_id(hit_rank))
#    print('hit descr')
#    print(parsed_file_obj.hit_description(hit_rank))
#    print('hit evalue')
#    print(str(parsed_file_obj.hit_evalue(hit_rank)))
#    print('hit score')
#    print(str(parsed_file_obj.hit_score(hit_rank)))
#    #print('hit seq')
#    #print(str(parsed_file_obj.hit_sequence(0).seq))
#    print('hit subseq description')
#    print(str(parsed_file_obj.hit_subsequence_and_coord(0)[0].description))
#    print('hit subseq')
#    print(str(parsed_file_obj.hit_subsequence_and_coord(0)[0].seq))
#    print('hit subseq coord')
#    print(str(parsed_file_obj.hit_subsequence_and_coord(0)[1]))

# query length?
# fwd_hit_score?
# query length?



def test_get_hmmer_hit_seq_coord():  # ***Incomplete test
    """Test the get_hmmer_hit_seq_coord function in the srchresfile.py file.
    """
    ##########################
    # Arrange.
    searchio_hit_obj = "searchio_hit_obj"
    db_file = "db_file"
    extra = "extra=0"

    ##########################
    # Act.
    #x = get_hmmer_hit_seq_coord(searchio_hit_obj,
    #		db_file,
    #		extra)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


