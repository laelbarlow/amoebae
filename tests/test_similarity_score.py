#!/usr/bin/env python3
"""PyTest tests for the similarity_score.py module.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib')) # Customize.

from similarity_score import \
get_similarity_score, \
get_score_dataframe_from_file, \
get_similarity_scores


def test_get_similarity_score():  # ***Incomplete test
    """Test the get_similarity_score function in the similarity_score.py file.
    """
    ##########################
    # Arrange.
    aa1 = "aa1"
    aa2 = "aa2"
    score_df = "score_df"

    ##########################
    # Act.
    #x = get_similarity_score(aa1,
    #		aa2,
    #		score_df)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_score_dataframe_from_file():  # ***Incomplete test
    """Test the get_score_dataframe_from_file function in the similarity_score.py file.
    """
    ##########################
    # Arrange.
    infp = "infp"

    ##########################
    # Act.
    #x = get_score_dataframe_from_file(infp)

    ##########################
    # Assert.
    assert True == True # ***Temporary.



def test_get_similarity_scores():  # ***Incomplete test
    """Test the get_similarity_scores function in the similarity_score.py file.
    """
    ##########################
    # Arrange.
    aas1 = "aas1"
    aas2 = "aas2"

    ##########################
    # Act.
    #x = get_similarity_scores(aas1,
    #		aas2)

    ##########################
    # Assert.
    assert True == True # ***Temporary.


