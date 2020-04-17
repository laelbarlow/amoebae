#!/usr/bin/env python3
"""PyTest tests for the generate_sankey_diagram.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from generate_sankey_diagram import \
generate_sankey


def test_generate_sankey():  # ***Incomplete test
    """Test the generate_sankey function in the generate_sankey_diagram.py file.
    """
    ##########################
    # Arrange.
    title = "title"

    ##########################
    # Act.
    #x = generate_sankey(title)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


