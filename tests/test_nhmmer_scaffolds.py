#!/usr/bin/env python3
"""PyTest tests for the test_nhmmer_scaffolds.py module.
"""
from test_nhmmer_scaffolds.py import *



def test_frameify_seq():
    """Test the frameify_seq function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    inseq = 'inseq'
    hmmfrom = 'hmmfrom'

    ##########################
    # Act.
    #x = frameify_seq(inseq,
    #		hmmfrom)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_if_hits_overlapping():
    """Test the check_if_hits_overlapping function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    hit_objs = 'hit_objs'

    ##########################
    # Act.
    #x = check_if_hits_overlapping(hit_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_db_fasta_filepath():
    """Test the get_db_fasta_filepath function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    dir_string = 'dir_string'
    dbname = 'dbname'

    ##########################
    # Act.
    #x = get_db_fasta_filepath(dir_string,
    #		dbname)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_run_nhmmer_search():
    """Test the run_nhmmer_search function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    query_path = 'query_path'
    db_path = 'db_path'
    nhmmer_outfilepath = 'nhmmer_outfilepath'

    ##########################
    # Act.
    #x = run_nhmmer_search(query_path,
    #		db_path,
    #		nhmmer_outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_run_hmmsearch():
    """Test the run_hmmsearch function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    query_path = 'query_path'
    db_path = 'db_path'
    hmmsearch_outfilepath = 'hmmsearch_outfilepath'

    ##########################
    # Act.
    #x = run_hmmsearch(query_path,
    #		db_path,
    #		hmmsearch_outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_parse_nhmmer_out():
    """Test the parse_nhmmer_out function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    nhmmer_outfilepath = 'nhmmer_outfilepath'

    ##########################
    # Act.
    #x = parse_nhmmer_out(nhmmer_outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_nhmmer_hit_fasta():
    """Test the write_nhmmer_hit_fasta function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    nhmmer_hit_objs = 'nhmmer_hit_objs'
    fasta_filepath = 'fasta_filepath'

    ##########################
    # Act.
    #x = write_nhmmer_hit_fasta(nhmmer_hit_objs,
    #		fasta_filepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_cluster_ali_recur():
    """Test the cluster_ali_recur function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    sublists = 'sublists'
    origlist = 'origlist'

    ##########################
    # Act.
    #x = cluster_ali_recur(sublists,
    #		origlist)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_concat_seq():
    """Test the get_concat_seq function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    nhmmer_hit_objs = 'nhmmer_hit_objs'
    dtype = 'dtype='nucl''

    ##########################
    # Act.
    #x = get_concat_seq(nhmmer_hit_objs,
    #		dtype)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_gene_region():
    """Test the get_gene_region function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    nhmmer_hit_objs = 'nhmmer_hit_objs'

    ##########################
    # Act.
    #x = get_gene_region(nhmmer_hit_objs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_display_dict():
    """Test the display_dict function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    nhmmer_hit_objs_dict2 = 'nhmmer_hit_objs_dict2'

    ##########################
    # Act.
    #x = display_dict(nhmmer_hit_objs_dict2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_three_frame_transl_string():
    """Test the get_three_frame_transl_string function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    inputs = 'inputs'

    ##########################
    # Act.
    #x = get_three_frame_transl_string(inputs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_three_frame_transl():
    """Test the get_three_frame_transl function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    seq_filepath = 'seq_filepath'
    three_frame_translation_dirpath = 'three_frame_translation_dirpath'

    ##########################
    # Act.
    #x = get_three_frame_transl(seq_filepath,
    #		three_frame_translation_dirpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_nhmmer_hit_gene_region_seqs():
    """Test the get_nhmmer_hit_gene_region_seqs function in the nhmmer_scaffolds.py file.
    """
    ##########################
    # Arrange.
    nhmmer_hit_objs = 'nhmmer_hit_objs'
    seq_filepath = 'seq_filepath'

    ##########################
    # Act.
    #x = get_nhmmer_hit_gene_region_seqs(nhmmer_hit_objs,
    #		seq_filepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


