#!/usr/bin/env python3
"""PyTest tests for proper accessibility of dependencies installed in the
environment (singularity container).
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

import subprocess

#def test_blast_plus_install():  # ***Incomplete test
#    """Test the make_easel_index function in the add_to_db.py file.
#    """
#    ##########################
#    # Arrange.
#    outfp = "outfp"
#
#    ##########################
#    # Act.
#    blastp_path = subprocess.check_output(['command', '-v', 'blastp'])
#    expected_blastp_path = "/opt/ncbi-blast-2.9.0+/bin/blastp"
#
#    ##########################
#    # Assert.
#    assert blastp_path == expected_blastp_path



def test_stuff():  # ***Incomplete test
    """Test...
    """
    ##########################
    # Arrange.
    outfp = "outfp"

    ##########################
    # Act.
    #...

    ##########################
    # Assert.
    assert True == True # ***Temporary.
