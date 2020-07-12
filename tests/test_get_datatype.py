#!/usr/bin/env python3
"""PyTest tests for the get_datatype.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from data.datatype_test_seqs import testseq1

from get_datatype import \
get_datatype_for_sequence_string, \
get_dbtype


def test_get_datatype_for_sequence_string():
    """Test the get_datatype_for_sequence_string function in the get_datatype.py file.
    """
    ##########################
    # Arrange.
    x1 = 'ATGC'*3 + 'Y'
    x2 = 'atgc'*3 + 'y'
    x3 = 'ATGCV'
    x4 = 'atgcv'
    x5 = 'tccaaaaaatcgaaTTATYttattccccaccttcttttctcattttttga'
    x6 = testseq1

    ##########################
    # Act.
    dbtype_1 = get_datatype_for_sequence_string(x1)
    dbtype_2 = get_datatype_for_sequence_string(x2)
    dbtype_3 = get_datatype_for_sequence_string(x3)
    dbtype_4 = get_datatype_for_sequence_string(x4)
    dbtype_5 = get_datatype_for_sequence_string(x5)
    dbtype_6 = get_datatype_for_sequence_string(x6)

    ##########################
    # Assert.
    assert dbtype_1 == 'nucl'
    assert dbtype_2 == 'nucl'
    assert dbtype_3 == 'prot'
    assert dbtype_4 == 'prot'
    assert dbtype_5 == 'nucl'
    assert dbtype_6 == 'nucl'


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


