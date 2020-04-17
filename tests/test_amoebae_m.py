#!/usr/bin/env python3
"""PyTest tests for the amoebae_m.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from amoebae_m import \
get_corr_fasta_exten, \
get_dbtype_from_file_exten, \
get_seqs_from_fasta_db, \
get_subseq_from_fasta_db, \
get_query_taxon_from_filename, \
get_query_title_from_filename, \
get_query_title_from_csv, \
get_query_taxon_from_csv, \
get_species_from_db_csv, \
get_db_filename_for_query_from_db_csv, \
get_species_for_db_filename, \
write_seqs_to_fasta, \
record_amoebae_info_in_log_file, \
mask_alignment2, \
mask_nex2, \
get_seq_obj_from_srch_res_csv_info, \
get_hit_range_from_hsp_ranges, \
apply_mask_criteria2, \
disqualifying_string_in_file_basename, \
find_input_file_in_parent_directory


def test_get_corr_fasta_exten():  # ***Incomplete test
    """Test the get_corr_fasta_exten function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    infp = "infp"

    ##########################
    # Act.
    #x = get_corr_fasta_exten(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_dbtype_from_file_exten():  # ***Incomplete test
    """Test the get_dbtype_from_file_exten function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    infp = "infp"

    ##########################
    # Act.
    #x = get_dbtype_from_file_exten(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seqs_from_fasta_db():  # ***Incomplete test
    """Test the get_seqs_from_fasta_db function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    db_name = "db_name"
    accs = "accs"
    slow = "slow=False"

    ##########################
    # Act.
    #x = get_seqs_from_fasta_db(db_name,
    #		accs,
    #		slow)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_subseq_from_fasta_db():  # ***Incomplete test
    """Test the get_subseq_from_fasta_db function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    db_name = "db_name"
    acc = "acc"
    subseq_coord = "subseq_coord"

    ##########################
    # Act.
    #x = get_subseq_from_fasta_db(db_name,
    #		acc,
    #		subseq_coord)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_query_taxon_from_filename():  # ***Incomplete test
    """Test the get_query_taxon_from_filename function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    query_filename = "query_filename"

    ##########################
    # Act.
    #x = get_query_taxon_from_filename(query_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_query_title_from_filename():  # ***Incomplete test
    """Test the get_query_title_from_filename function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    query_filename = "query_filename"

    ##########################
    # Act.
    #x = get_query_title_from_filename(query_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_query_title_from_csv():  # ***Incomplete test
    """Test the get_query_title_from_csv function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    query_filename = "query_filename"

    ##########################
    # Act.
    #x = get_query_title_from_csv(query_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_query_taxon_from_csv():  # ***Incomplete test
    """Test the get_query_taxon_from_csv function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    query_filename = "query_filename"

    ##########################
    # Act.
    #x = get_query_taxon_from_csv(query_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_species_from_db_csv():  # ***Incomplete test
    """Test the get_species_from_db_csv function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    taxon = "taxon"

    ##########################
    # Act.
    #x = get_species_from_db_csv(taxon)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_db_filename_for_query_from_db_csv():  # ***Incomplete test
    """Test the get_db_filename_for_query_from_db_csv function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    taxon = "taxon"

    ##########################
    # Act.
    #x = get_db_filename_for_query_from_db_csv(taxon)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_species_for_db_filename():  # ***Incomplete test
    """Test the get_species_for_db_filename function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    db_filename = "db_filename"

    ##########################
    # Act.
    #x = get_species_for_db_filename(db_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_seqs_to_fasta():  # ***Incomplete test
    """Test the write_seqs_to_fasta function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    output_dir = "output_dir"
    abbrev = "abbrev=False"

    ##########################
    # Act.
    #x = write_seqs_to_fasta(csv_file,
    #		output_dir,
    #		abbrev)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_record_amoebae_info_in_log_file():  # ***Incomplete test
    """Test the record_amoebae_info_in_log_file function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    commandline = "commandline"
    outdir = "outdir"
    start_time = "start_time"
    end_time = "end_time"

    ##########################
    # Act.
    #x = record_amoebae_info_in_log_file(commandline,
    #		outdir,
    #		start_time,
    #		end_time)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_mask_alignment2():  # ***Incomplete test
    """Test the mask_alignment2 function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"

    ##########################
    # Act.
    #x = mask_alignment2(alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_mask_nex2():  # ***Incomplete test
    """Test the mask_nex2 function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    infilepath = "infilepath"
    outfilepath = "outfilepath"

    ##########################
    # Act.
    #x = mask_nex2(infilepath,
    #		outfilepath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seq_obj_from_srch_res_csv_info():  # ***Incomplete test
    """Test the get_seq_obj_from_srch_res_csv_info function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    acc = "acc"
    description = "description"
    seq = "seq"
    abbrev = "abbrev=False"

    ##########################
    # Act.
    #x = get_seq_obj_from_srch_res_csv_info(acc,
    #		description,
    #		seq,
    #		abbrev)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_hit_range_from_hsp_ranges():  # ***Incomplete test
    """Test the get_hit_range_from_hsp_ranges function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    subseq_coord = "subseq_coord"

    ##########################
    # Act.
    #x = get_hit_range_from_hsp_ranges(subseq_coord)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_apply_mask_criteria2():  # ***Incomplete test
    """Test the apply_mask_criteria2 function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    column = "column"

    ##########################
    # Act.
    #x = apply_mask_criteria2(column)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_disqualifying_string_in_file_basename():  # ***Incomplete test
    """Test the disqualifying_string_in_file_basename function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    f = "f"
    disqualifying_strings = "disqualifying_strings"

    ##########################
    # Act.
    #x = disqualifying_string_in_file_basename(f,
    #		disqualifying_strings)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_find_input_file_in_parent_directory():  # ***Incomplete test
    """Test the find_input_file_in_parent_directory function in the amoebae_m.py file.
    """
    ##########################
    # Arrange.
    indp = "indp"
    extension = "extension"
    disqualifying_strings = "disqualifying_strings"

    ##########################
    # Act.
    #x = find_input_file_in_parent_directory(indp,
    #		extension,
    #		disqualifying_strings)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


