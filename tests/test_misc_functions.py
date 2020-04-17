#!/usr/bin/env python3
"""PyTest tests for the misc_functions.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from misc_functions import \
launch, \
query_yes_no, \
get_fa_record_text_from_obj, \
get_abbrev_fa_record_text_from_obj


def test_launch():  # ***Incomplete test
    """Test the launch function in the misc_functions.py file.
    """
    ##########################
    # Arrange.
    app = "app"
    thing = "thing"

    ##########################
    # Act.
    #x = launch(app,
    #		thing)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_query_yes_no():  # ***Incomplete test
    """Test the query_yes_no function in the misc_functions.py file.
    """
    ##########################
    # Arrange.
    question = "question"
    default = "default='yes'"

    ##########################
    # Act.
    #x = query_yes_no(question,
    #		default)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_fa_record_text_from_obj():  # ***Incomplete test
    """Test the get_fa_record_text_from_obj function in the misc_functions.py file.
    """
    ##########################
    # Arrange.
    obj = "obj"

    ##########################
    # Act.
    #x = get_fa_record_text_from_obj(obj)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_abbrev_fa_record_text_from_obj():  # ***Incomplete test
    """Test the get_abbrev_fa_record_text_from_obj function in the misc_functions.py file.
    """
    ##########################
    # Arrange.
    obj = "obj"
    db_name = "db_name"
    prot_name = "prot_name=None"

    ##########################
    # Act.
    #x = get_abbrev_fa_record_text_from_obj(obj,
    #		db_name,
    #		prot_name)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


