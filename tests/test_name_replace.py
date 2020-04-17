#!/usr/bin/env python3
"""PyTest tests for the name_replace.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from name_replace import *



def test_generate_code_name():
    """Test the generate_code_name function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    number = 'number'

    ##########################
    # Act.
    #x = generate_code_name(number)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_afa_with_code_names():
    """Test the write_afa_with_code_names function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    infile = 'infile'
    outfile = 'outfile'
    tablefile = 'tablefile'

    ##########################
    # Act.
    #x = write_afa_with_code_names(infile,
    #		outfile,
    #		tablefile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_nex_with_code_names():
    """Test the write_nex_with_code_names function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    infile = 'infile'
    outfile = 'outfile'
    tablefile = 'tablefile'

    ##########################
    # Act.
    #x = write_nex_with_code_names(infile,
    #		outfile,
    #		tablefile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_conversion_dict_from_table():
    """Test the get_conversion_dict_from_table function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    tablefile = 'tablefile'

    ##########################
    # Act.
    #x = get_conversion_dict_from_table(tablefile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_newick_tree_with_uncoded_names():
    """Test the write_newick_tree_with_uncoded_names function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    infile = 'infile'
    outfile = 'outfile'
    tablefile = 'tablefile'

    ##########################
    # Act.
    #x = write_newick_tree_with_uncoded_names(infile,
    #		outfile,
    #		tablefile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_write_newick_tree_with_coded_names():
    """Test the write_newick_tree_with_coded_names function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    infile = 'infile'
    outfile = 'outfile'
    tablefile = 'tablefile'

    ##########################
    # Act.
    #x = write_newick_tree_with_coded_names(infile,
    #		outfile,
    #		tablefile)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_codenames_nex():
    """Test the codenames_nex function in the name_replace.py file.
    """
    ##########################
    # Arrange.
    infp = 'infp'
    aamodelpr = 'aamodelpr=None'

    ##########################
    # Act.
    #x = codenames_nex(infp,
    #		aamodelpr)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


