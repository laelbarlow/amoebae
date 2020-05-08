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
"""This test module contains functions for testing the
predict_redun_hits_selections.py module. 
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))

from predict_redun_hit_selection import process_query_title_hit_id_list_dict, compare_two_hit_id_lists


def test_compare_two_hit_id_lists():  
    """Test the compare_two_hit_id_lists function.
    """
    ##########################
    # Arrange.
    ids1 = ['A', 'B', 'C', 'D', 'E', 'F']
    ids2 = ['D', 'E', 'F', 'A', 'B', 'C']

    ids3 = ['A', 'D', 'B', 'C', 'E', 'F']
    ids4 = ['D', 'E', 'F', 'A', 'B', 'C']

    AP3B = \
"""NP_003655.3
XP_005248675.1
NP_001258698.1
XP_016865490.1
XP_005248676.1
NP_004635.2
NP_001265441.1
NP_001265440.1
XP_016879776.1
XP_005257998.1
XP_011522753.1
XP_011522752.1
XP_011522751.1
XP_011522750.1
XP_016879773.1
NP_001273.1
XP_005257995.1
XP_005257994.1
NP_001025177.1
NP_001118.3
NP_663782.2
NP_001365491.1
NP_001159491.1
NP_001365492.1
NP_001365495.1
NP_001365493.1
NP_001365494.1
XP_016879775.1
XP_011522757.1
XP_011522756.1
XP_011522755.1
XP_011522754.1
XP_016879774.1
NP_001335369.1
NP_001335370.1
NP_006585.2
NP_001240781.1
XP_016855582.1
XP_024308191.1
XP_016855578.1
NP_001240782.1
XP_011538825.1
XP_016855579.1
NP_001295241.1
XP_024308190.1
XP_011538827.1
XP_024308203.1
XP_016855580.1
XP_024308209.1
NP_057535.1
NP_001137534.1
NP_001137533.1
NP_001119.3
NP_001025178.1
XP_005268230.1
XP_005268229.1
NP_003908.1
XP_016855581.1""".split('\n')

    AP3Bmod = \
"""NP_003655.3
XP_005248675.1
NP_001258698.1
XP_016865490.1
XP_005248676.1
NP_004635.2
NP_001265441.1
NP_001265440.1""".split('\n')
    
    AP4B = \
"""NP_006585.2
NP_001240781.1
NP_001240782.1
XP_011538825.1
XP_024308190.1
XP_011538827.1
XP_016855579.1
NP_001295241.1
XP_024308191.1
XP_016855582.1
XP_016855578.1
XP_016855581.1
XP_024308203.1
XP_024308209.1
XP_016855580.1
XP_016879776.1
XP_005257998.1
XP_016879773.1
NP_001273.1
XP_005257995.1
XP_005257994.1
NP_001025177.1
XP_011522753.1
XP_011522752.1
XP_011522751.1
XP_011522750.1
NP_001365495.1
NP_001159491.1
NP_001118.3
NP_001365491.1
NP_001365492.1
NP_663782.2
NP_001365494.1
XP_016879775.1
XP_011522757.1
XP_016879774.1
NP_001365493.1
XP_011522756.1
XP_011522755.1
XP_011522754.1
NP_004635.2
NP_001265441.1
XP_005248676.1
NP_001258698.1
XP_016865490.1
XP_005248675.1
NP_003655.3
NP_001265440.1
NP_036265.3
NP_001276962.1
NP_001335369.1
NP_057535.1
NP_001137534.1
NP_001137533.1
XP_005268230.1
XP_005268229.1
NP_003908.1
XP_016877234.1""".split('\n')

    AP4Bmod = \
"""NP_006585.2
NP_001240781.1
NP_001240782.1
XP_011538825.1
XP_024308190.1
XP_011538827.1
XP_016855579.1
NP_001295241.1
XP_024308191.1
XP_016855582.1
XP_016855578.1
XP_016855581.1
XP_024308203.1
XP_024308209.1
XP_016855580.1""".split('\n')



    ##########################
    # Act.
    ids1a, ids2a = compare_two_hit_id_lists(ids1, ids2)

    ids3a, ids4a = compare_two_hit_id_lists(ids3, ids4)

    AP3Ba, AP4Ba = compare_two_hit_id_lists(AP3B, AP4B)


    ##########################
    # Assert.

    assert ids1a == ['A', 'B', 'C']
    assert ids2a == ['D', 'E', 'F']

    ## Should only remove elements if it actually reduces the number of
    ## overlapping elements.
    #assert ids3a == ['A', 'D', 'B', 'C']
    #assert ids4a == ['D', 'E', 'F']

    #assert AP3Ba == AP3Bmod
    assert AP4Ba == AP4Bmod

    # More generally, the lists should not contain overlapping sets of ids.
    #assert len(set(ids1a).intersection(set(ids2a))) == 0


def test_process_query_title_hit_id_list_dict():
    """Test the process_query_title_hit_id_list_dict function.
    """
    ##########################
    # Arrange.
    query_title_id_list_dict = {'qt1': ['A', 'B', 'C', 'D', 'E', 'F'],
                                'qt2': ['D', 'E', 'F', 'A', 'B', 'C'],
                                'qt3': ['G', 'H', 'I', 'J', 'K', 'L'],
                                'qt4': ['J', 'K', 'L', 'G', 'H', 'I']
                                }

    ##########################
    # Act.
    #ids1a, ids2a = compare_two_hit_id_lists(ids1, ids2)
    query_title_id_list_dict2 =\
    process_query_title_hit_id_list_dict(query_title_id_list_dict)

    ##########################
    # Assert.
    #assert ids1a == ['A', 'B', 'C']
    #assert ids2a == ['D', 'E', 'F']
    # More generally, the lists should not contain overlapping sets of ids.
    #assert len(set(ids1a).intersection(set(ids2a))) == 0

    assert query_title_id_list_dict2 ==\
                               {'qt1': ['A', 'B', 'C'],
                                'qt2': ['D', 'E', 'F'],
                                'qt3': ['G', 'H', 'I'],
                                'qt4': ['J', 'K', 'L']
                                }








