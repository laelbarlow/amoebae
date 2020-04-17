#!/usr/bin/env python3
"""PyTest tests for the get_datatype.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from get_datatype import \
get_datatype_for_sequence_string, \
get_dbtype


def test_get_datatype_for_sequence_string():  # ***Incomplete test
    """Test the get_datatype_for_sequence_string function in the get_datatype.py file.
    """
    ##########################
    # Arrange.
    concat_seq = "concat_seq"

    ##########################
    # Act.
    #x = get_datatype_for_sequence_string(concat_seq)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_dbtype():  # ***Incomplete test
    """Test the get_dbtype function in the get_datatype.py file.
    """
    ##########################
    # Arrange.
    f = "f"

    ##########################
    # Act.
    #x = get_dbtype(f)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


