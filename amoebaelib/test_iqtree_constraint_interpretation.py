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
"""Tests functions that run simple constrained IQtree tree searches to test how
IQtree interprets constraints.
"""

# Import necessary modules.
import unittest
import os
import time
import shutil
import subprocess
from ete3 import Tree
from iqtree_constraint_interpretation import run_simple_constrained_iqtree_analysis

# Define a class with functions for performing tests.
class TestConstInterp(unittest.TestCase):

    def test_that_constraint_works(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 4 10',
                                      'A  VLVLVLLLLL',
                                      'B  VLVLFLLLLL',
                                      'C  FLFLVLLLLL',
                                      'D  FLFLFLLLLL'
                                      ])

        constraint_tree_string = "(A, C, (B, D));"
        constraint_tree_obj = Tree(constraint_tree_string)

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_obj)

        # Check that output topology is the same as the constraint topology.
        self.assertTrue('(B,D)' in output_topo.replace(' ', '')\
                        or '(A,C)' in output_topo.replace(' ', ''))


    def test_that_polytomy_resolved(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 4 10',
                                      'A  VLVLVLLLLL',
                                      'B  VLVLFLLLLL',
                                      'C  FLFLVLLLLL',
                                      'D  FLFLFLLLLL'
                                      ])

        constraint_tree_string = "(A, B, C, D);"
        constraint_tree_obj = Tree(constraint_tree_string)

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_obj)

        # Check that output topology is the same as the constraint topology.
        self.assertNotEqual(output_topo.replace(' ', ''), constraint_tree_string.replace(' ', ''))

    def test_that_an_outgroup_can_be_placed_on_any_branch(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 6 10',
                                      'A  LLLLLLLLLL',
                                      'B  LLLLLLLLLL',
                                      'C  FFFFFVLLLL',
                                      'D  FFFFFLLLLL',
                                      'X  FFFFFVLLLL',
                                      'Y  FFFFFVLLLL'
                                      ])

        #constraint_tree_string = "(X, (A, B, (C, D)));"
        constraint_tree_string = "(X, Y, (A, B, (C, D)));"
        constraint_tree_obj = Tree(constraint_tree_string)

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_obj)
        print('\n')
        print(output_topo)
        print('\n')

        # Check that output topology is the same as the constraint topology.
        self.assertTrue('(C,D,X)' in output_topo.replace(' ', ''))




if __name__ == '__main__':
    unittest.main()
