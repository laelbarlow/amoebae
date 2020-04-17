#!/usr/bin/env python3
"""PyTest tests for the add_to_models.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from add_to_models import \
update_models_csv


def test_update_models_csv():  # ***Incomplete test
    """Test the update_models_csv function in the add_to_models.py file.
    """
    ##########################
    # Arrange.
    model_name = "model_name"
    csv_file = "csv_file"
    alignmentfp = "alignmentfp"
    topologyfp = "topologyfp"
    subs_model = "subs_model"
    type_seqsfp = "type_seqsfp"

    ##########################
    # Act.
    #x = update_models_csv(model_name,
    #		csv_file,
    #		alignmentfp,
    #		topologyfp,
    #		subs_model,
    #		type_seqsfp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


