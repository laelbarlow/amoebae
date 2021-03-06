#!/usr/bin/env python3
"""PyTest tests for the paralogue_counter.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from paralogue_counter import \
remove_irrelevant_data_from_alignment, \
check_for_redundant_seq_names, \
rank_seqs_by_hmm, \
get_sig_overlap, \
get_percent_identity, \
add_seq_to_alignment2, \
modify_seq_descr_for_tree, \
add_seq_to_alignment3, \
get_sequence_comparison_dict, \
get_all_comparison_output_filepath, \
find_redun_model_recursively, \
count_paralogues, \
get_seq_obj_x, \
get_seq_obj_x2, \
get_cur_ali_num, \
count_paralogues2, \
get_overlapping_genes_from_gff3, \
two_proteins_same_gene_in_gff3, \
cds_record_matches_prot_id, \
two_proteins_same_gene_in_gff3_2, \
two_proteins_same_gene_in_gff3_3, \
annot_id_matches_prot_id, \
get_relevant_annot_sql_file, \
count_paralogues3, \
add_alignment_column, \
add_model_column


def test_remove_irrelevant_data_from_alignment():  # ***Incomplete test
    """Test the remove_irrelevant_data_from_alignment function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    ali_plus_top_seq_plus_1 = "ali_plus_top_seq_plus_1"

    ##########################
    # Act.
    #x = remove_irrelevant_data_from_alignment(ali_plus_top_seq_plus_1)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_check_for_redundant_seq_names():  # ***Incomplete test
    """Test the check_for_redundant_seq_names function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"

    ##########################
    # Act.
    #x = check_for_redundant_seq_names(alignment)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_rank_seqs_by_hmm():  # ***Incomplete test
    """Test the rank_seqs_by_hmm function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inali = "inali"
    inseqs = "inseqs"

    ##########################
    # Act.
    #x = rank_seqs_by_hmm(inali,
    #		inseqs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_sig_overlap():  # ***Incomplete test
    """Test the get_sig_overlap function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inali = "inali"

    ##########################
    # Act.
    #x = get_sig_overlap(inali)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_percent_identity():  # ***Incomplete test
    """Test the get_percent_identity function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inali = "inali"

    ##########################
    # Act.
    #x = get_percent_identity(inali)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_add_seq_to_alignment2():  # ***Incomplete test
    """Test the add_seq_to_alignment2 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inseqobj = "inseqobj"
    innexpath = "innexpath"
    outnexpath = "outnexpath"

    ##########################
    # Act.
    #x = add_seq_to_alignment2(inseqobj,
    #		innexpath,
    #		outnexpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_modify_seq_descr_for_tree():  # ***Incomplete test
    """Test the modify_seq_descr_for_tree function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inseqobj = "inseqobj"

    ##########################
    # Act.
    #x = modify_seq_descr_for_tree(inseqobj)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_add_seq_to_alignment3():  # ***Incomplete test
    """Test the add_seq_to_alignment3 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inseqobj = "inseqobj"
    innexpath = "innexpath"
    outnexpath = "outnexpath"

    ##########################
    # Act.
    #x = add_seq_to_alignment3(inseqobj,
    #		innexpath,
    #		outnexpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_sequence_comparison_dict():  # ***Incomplete test
    """Test the get_sequence_comparison_dict function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    inali = "inali"
    score_df = "score_df"

    ##########################
    # Act.
    #x = get_sequence_comparison_dict(inali,
    #		score_df)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_all_comparison_output_filepath():  # ***Incomplete test
    """Test the get_all_comparison_output_filepath function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = get_all_comparison_output_filepath(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_find_redun_model_recursively():  # ***Incomplete test
    """Test the find_redun_model_recursively function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    outdir = "outdir"

    ##########################
    # Act.
    #x = find_redun_model_recursively(outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_count_paralogues():  # ***Incomplete test
    """Test the count_paralogues function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"
    fasta = "fasta"
    fwdeval = "fwdeval"
    identity = "identity"

    ##########################
    # Act.
    #x = count_paralogues(alignment,
    #		fasta,
    #		fwdeval,
    #		identity)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seq_obj_x():  # ***Incomplete test
    """Test the get_seq_obj_x function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    acc_list = "acc_list"
    fastadir = "fastadir"

    ##########################
    # Act.
    #x = get_seq_obj_x(acc_list,
    #		fastadir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_seq_obj_x2():  # ***Incomplete test
    """Test the get_seq_obj_x2 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    acc_list = "acc_list"

    ##########################
    # Act.
    #x = get_seq_obj_x2(acc_list)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_cur_ali_num():  # ***Incomplete test
    """Test the get_cur_ali_num function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    alignmentdir = "alignmentdir"

    ##########################
    # Act.
    #x = get_cur_ali_num(alignmentdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_count_paralogues2():  # ***Incomplete test
    """Test the count_paralogues2 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"
    alignmentdir = "alignmentdir"
    fastadir = "fastadir"
    fwdeval = "fwdeval"
    metric_name = "metric_name"

    ##########################
    # Act.
    #x = count_paralogues2(csv_file,
    #		alignmentdir,
    #		fastadir,
    #		fwdeval,
    #		metric_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_overlapping_genes_from_gff3():  # ***Incomplete test
    """Test the get_overlapping_genes_from_gff3 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    sql_database = "sql_database"

    ##########################
    # Act.
    #x = get_overlapping_genes_from_gff3(sql_database)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_two_proteins_same_gene_in_gff3():  # ***Incomplete test
    """Test the two_proteins_same_gene_in_gff3 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    sql_database = "sql_database"
    prot_id_1 = "prot_id_1"
    prot_id_2 = "prot_id_2"

    ##########################
    # Act.
    #x = two_proteins_same_gene_in_gff3(sql_database,
    #		prot_id_1,
    #		prot_id_2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_cds_record_matches_prot_id():  # ***Incomplete test
    """Test the cds_record_matches_prot_id function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    cds_obj = "cds_obj"
    prot_id = "prot_id"

    ##########################
    # Act.
    #x = cds_record_matches_prot_id(cds_obj,
    #		prot_id)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_two_proteins_same_gene_in_gff3_2():  # ***Incomplete test
    """Test the two_proteins_same_gene_in_gff3_2 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    sql_database = "sql_database"
    prot_id_1 = "prot_id_1"
    prot_id_2 = "prot_id_2"

    ##########################
    # Act.
    #x = two_proteins_same_gene_in_gff3_2(sql_database,
    #		prot_id_1,
    #		prot_id_2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_two_proteins_same_gene_in_gff3_3():  # ***Incomplete test
    """Test the two_proteins_same_gene_in_gff3_3 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    sql_database = "sql_database"
    prot_id_1 = "prot_id_1"
    prot_id_2 = "prot_id_2"

    ##########################
    # Act.
    #x = two_proteins_same_gene_in_gff3_3(sql_database,
    #		prot_id_1,
    #		prot_id_2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_annot_id_matches_prot_id():  # ***Incomplete test
    """Test the annot_id_matches_prot_id function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    annot_id = "annot_id"

    ##########################
    # Act.
    #x = annot_id_matches_prot_id(annot_id)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_relevant_annot_sql_file():  # ***Incomplete test
    """Test the get_relevant_annot_sql_file function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    sequence_filename = "sequence_filename"

    ##########################
    # Act.
    #x = get_relevant_annot_sql_file(sequence_filename)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_count_paralogues3():  # ***Incomplete test
    """Test the count_paralogues3 function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    csv_file = "csv_file"

    ##########################
    # Act.
    #x = count_paralogues3(csv_file)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_add_alignment_column():  # ***Incomplete test
    """Test the add_alignment_column function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    incsv = "incsv"
    outcsv = "outcsv"

    ##########################
    # Act.
    #x = add_alignment_column(incsv,
    #		outcsv)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_add_model_column():  # ***Incomplete test
    """Test the add_model_column function in the paralogue_counter.py file.
    """
    ##########################
    # Arrange.
    incsv = "incsv"
    outcsv = "outcsv"

    ##########################
    # Act.
    #x = add_model_column(incsv,
    #		outcsv)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


