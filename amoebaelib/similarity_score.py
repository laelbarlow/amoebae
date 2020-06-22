#!/usr/bin/env python3
# Copyright 2018 Lael D. Barlow
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
"""Module for getting scores of amino acid substitutions from a scoring matrix
(e.g., BLOSUM62 from NCBI).

The main functions to use are get_score_dataframe_from_file and then
get_similarity_score for each pair of amino acids.
"""
import pandas as pd
from random import shuffle

def get_similarity_score(aa1, aa2, score_df):
    """Take two amino acid (residues) and return their similarity score
    according to the blosum62 matrix from NCBI.
    """
    aa1 = aa1.upper()
    aa2 = aa2.upper()
    valid_aa_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*']

    assert aa1 in valid_aa_list and aa2 in valid_aa_list, """Invalid input amino acid
    codes."""

    return int(score_df.loc[aa1][aa2])


def get_score_dataframe_from_file(infp):
    """Take a filepath for a csv file with a scoring matrix, and return a
    pandas dataframe with the data it contains.
    """
    return pd.read_csv(infp, index_col=0, encoding='utf-8')


def get_similarity_scores(aas1, aas2):
    """bla.
    """
    # Get score matrix as a pandas dataframe.
    score_df = get_score_dataframe_from_file('blosum62.csv')

    # Get scores for each pair of residues.
    for aa1, aa2 in zip(aas1, aas2):
        score = get_similarity_score(aa1, aa2, score_df)
        print(aa1 + ' vs. ' + aa2 + ': ' + str(score))
        if score > 0:
            print('positive substitution')
        else:
            print('negative substitution')



if __name__ == '__main__':

    valid_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*']

    aas1 = ['A', 'R', 'N', 'D', 'C']
    shuffle(valid_list)
    aas1 = valid_list.copy()

    aas2 = ['A', 'R', 'N', 'D', 'C']
    shuffle(valid_list)
    aas2 = valid_list.copy()

    get_similarity_scores(aas1, aas2)
