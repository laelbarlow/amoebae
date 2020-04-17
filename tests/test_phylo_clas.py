#!/usr/bin/env python3
"""PyTest tests for the phylo_clas.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from phylo_clas import \
add_seq_to_alignment, \
code_names_in_ali, \
get_conversion_dict, \
quote_tree, \
code_tree, \
uncode_tree, \
uncode_tree_obj, \
classify_seq_with_constrained_tree, \
get_phylo_class_csv, \
classify_seq_with_constrained_tree2, \
get_clade_name_from_model, \
get_phylo_class_logfile_path, \
classify_one_seq_record, \
place_one_seq_record, \
write_phylo_class_to_csv, \
get_all_alt_model_backbones


def test_add_seq_to_alignment():  # ***Incomplete test
    """Test the add_seq_to_alignment function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    inseqobj = "inseqobj"
    innexpath = "innexpath"
    outdp = "outdp"

    ##########################
    # Act.
    #x = add_seq_to_alignment(inseqobj,
    #		innexpath,
    #		outdp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_code_names_in_ali():  # ***Incomplete test
    """Test the code_names_in_ali function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    indp = "indp"
    inalifp = "inalifp"
    outalifp = "outalifp"
    outtablefp = "outtablefp"

    ##########################
    # Act.
    #x = code_names_in_ali(indp,
    #		inalifp,
    #		outalifp,
    #		outtablefp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_conversion_dict():  # ***Incomplete test
    """Test the get_conversion_dict function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    tablefile = "tablefile"

    ##########################
    # Act.
    #x = get_conversion_dict(tablefile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_quote_tree():  # ***Incomplete test
    """Test the quote_tree function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    intree = "intree"
    table = "table"

    ##########################
    # Act.
    #x = quote_tree(intree,
    #		table)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_code_tree():  # ***Incomplete test
    """Test the code_tree function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    intreefp = "intreefp"
    outtreefp = "outtreefp"
    table = "table"

    ##########################
    # Act.
    #x = code_tree(intreefp,
    #		outtreefp,
    #		table)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_uncode_tree():  # ***Incomplete test
    """Test the uncode_tree function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    intreefp = "intreefp"
    outtreefp = "outtreefp"
    table = "table"

    ##########################
    # Act.
    #x = uncode_tree(intreefp,
    #		outtreefp,
    #		table)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_uncode_tree_obj():  # ***Incomplete test
    """Test the uncode_tree_obj function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    intreeobj = "intreeobj"
    table = "table"

    ##########################
    # Act.
    #x = uncode_tree_obj(intreeobj,
    #		table)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_classify_seq_with_constrained_tree():  # ***Incomplete test
    """Test the classify_seq_with_constrained_tree function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    alignment = "alignment"
    tree = "tree"
    subs_model = "subs_model"
    type_seqs = "type_seqs"
    fasta = "fasta"

    ##########################
    # Act.
    #x = classify_seq_with_constrained_tree(alignment,
    #		tree,
    #		subs_model,
    #		type_seqs,
    #		fasta)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_phylo_class_csv():  # ***Incomplete test
    """Test the get_phylo_class_csv function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    outdirpath = "outdirpath"

    ##########################
    # Act.
    #x = get_phylo_class_csv(outdirpath)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_classify_seq_with_constrained_tree2():  # ***Incomplete test
    """Test the classify_seq_with_constrained_tree2 function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    file_with_seqs = "file_with_seqs"
    model = "model=None"

    ##########################
    # Act.
    #x = classify_seq_with_constrained_tree2(file_with_seqs,
    #		model)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_clade_name_from_model():  # ***Incomplete test
    """Test the get_clade_name_from_model function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    type_seq_name = "type_seq_name"
    type_seqs = "type_seqs"

    ##########################
    # Act.
    #x = get_clade_name_from_model(type_seq_name,
    #		type_seqs)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_phylo_class_logfile_path():  # ***Incomplete test
    """Test the get_phylo_class_logfile_path function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    phylo_class_outdir = "phylo_class_outdir"

    ##########################
    # Act.
    #x = get_phylo_class_logfile_path(phylo_class_outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_classify_one_seq_record():  # ***Incomplete test
    """Test the classify_one_seq_record function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    ali_num = "ali_num"
    record = "record"
    model = "model"

    ##########################
    # Act.
    #x = classify_one_seq_record(ali_num,
    #		record,
    #		model)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_place_one_seq_record():  # ***Incomplete test
    """Test the place_one_seq_record function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    ali_num = "ali_num"
    record = "record"
    model = "model"

    ##########################
    # Act.
    #x = place_one_seq_record(ali_num,
    #		record,
    #		model)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_phylo_class_to_csv():  # ***Incomplete test
    """Test the write_phylo_class_to_csv function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    phylo_class_id = "phylo_class_id"
    outdir = "outdir"

    ##########################
    # Act.
    #x = write_phylo_class_to_csv(phylo_class_id,
    #		outdir)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_all_alt_model_backbones():  # ***Incomplete test
    """Test the get_all_alt_model_backbones function in the phylo_clas.py file.
    """
    ##########################
    # Arrange.
    model_name = "model_name"

    ##########################
    # Act.
    #x = get_all_alt_model_backbones(model_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


