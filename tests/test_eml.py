#!/usr/bin/env python3
"""PyTest tests for the eml.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from eml import *



def test_send_eml():
    """Test the send_eml function in the eml.py file.
    """
    ##########################
    # Arrange.
    subject = 'subject'
    body = 'body'

    ##########################
    # Act.
    #x = send_eml(subject,
    #		body)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


