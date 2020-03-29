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
"""Tests functions in visualize_trees module.
"""

# Import necessary modules.
import unittest
from ete3 import Tree
from visualize_trees import get_nodes_with_paralogues

# Define a class with functions for performing tests.
class TestGetNodesWithParalogues(unittest.TestCase):

    def get_output_leaf_name_sets(self, output_nodes):

        output_leaf_name_sets = []
        for i in output_nodes:
            output_leaf_name_sets.append(frozenset([x.name for x in i.get_leaves()]))

        return set(output_leaf_name_sets)


    def test_that_clades_with_paralogues_identified(self):
        
        # Define a simple tree with two clades that should both be identified
        # as containing paralogues.
        input_tree_string = '((A__1, B__1),(A__2, B__2));' 

        # Parse as an ete3 TreeNode object for input to the
        # get_nodes_with_paralogues function.
        input_tree_obj = Tree(input_tree_string) 

        # Get output from the get_nodes_with_paralogues function.
        output_nodes = get_nodes_with_paralogues(input_tree_obj)

        # Get the set of leaf names for each output TreeNode object.
        output_leaf_name_sets = self.get_output_leaf_name_sets(output_nodes) 
        # Define the sets of leaf names that should the output TreeNode objects
        # should have.
        correct_sets = set([frozenset(['A__1', 'B__1']), frozenset(['A__2', 'B__2'])])

        # Check that the output set is the expected set (set of sets).
        self.assertEqual(output_leaf_name_sets, correct_sets)


    def test_that_early_branching_nonparalogous_seqs_excluded(self):
        
        # Define a simple tree with two clades that should both be identified
        # as containing paralogues.
        input_tree_string = '((X__1,(A__1, B__1)),(Y__1,(A__2, B__2)));' 

        # Parse as an ete3 TreeNode object for input to the
        # get_nodes_with_paralogues function.
        input_tree_obj = Tree(input_tree_string) 

        # Get output from the get_nodes_with_paralogues function.
        output_nodes = get_nodes_with_paralogues(input_tree_obj)

        # Get the set of leaf names for each output TreeNode object.
        output_leaf_name_sets = self.get_output_leaf_name_sets(output_nodes)

        # Define the sets of leaf names that should the output TreeNode objects
        # should have.
        correct_sets = set([frozenset(['A__1', 'B__1']), frozenset(['A__2', 'B__2'])])

        # Check that the output set is the expected set (set of sets).
        self.assertEqual(output_leaf_name_sets, correct_sets)


    def test_that_nested_nonparalogous_seqs_included(self):
        
        # Define a simple tree with two clades that should both be identified
        # as containing paralogues.
        input_tree_string = '((A__1,(X__1, B__1)),(B__2,(A__2, Y__1)));' 

        # Parse as an ete3 TreeNode object for input to the
        # get_nodes_with_paralogues function.
        input_tree_obj = Tree(input_tree_string) 

        # Get output from the get_nodes_with_paralogues function.
        output_nodes = get_nodes_with_paralogues(input_tree_obj)

        # Get the set of leaf names for each output TreeNode object.
        output_leaf_name_sets = self.get_output_leaf_name_sets(output_nodes)

        # Define the sets of leaf names that should the output TreeNode objects
        # should have.
        correct_sets = set([frozenset(['A__1', 'B__1', 'X__1']),
                            frozenset(['A__2', 'B__2', 'Y__1'])])

        # Check that the output set is the expected set (set of sets).
        self.assertEqual(output_leaf_name_sets, correct_sets)


    def test_that_redundant_subclades_not_identified(self):
        
        # Define a simple tree with two clades that should both be identified
        # as containing paralogues.
        input_tree_string = '((C__1,(A__1, B__1)),(C__2,(A__2, B__2)));' 

        # Parse as an ete3 TreeNode object for input to the
        # get_nodes_with_paralogues function.
        input_tree_obj = Tree(input_tree_string) 

        # Get output from the get_nodes_with_paralogues function.
        output_nodes = get_nodes_with_paralogues(input_tree_obj)

        # Get the set of leaf names for each output TreeNode object.
        output_leaf_name_sets = self.get_output_leaf_name_sets(output_nodes)

        # Define the sets of leaf names that should the output TreeNode objects
        # should have.
        correct_sets = set([frozenset(['C__1', 'A__1', 'B__1']),
                            frozenset(['C__2', 'A__2', 'B__2'])])

        # Check that the output set is the expected set (set of sets).
        self.assertEqual(output_leaf_name_sets, correct_sets)


    # This feature is less useful, so not testing for it.
    #def test_basal_seqs_of_same_sp_set_not_excluded(self):
    #    
    #    # Define a simple tree with two clades that should both be identified
    #    # as containing paralogues.
    #    input_tree_string = '((A__3,(A__1, B__1)),(B__3,(A__2, B__2)));' 

    #    # Parse as an ete3 TreeNode object for input to the
    #    # get_nodes_with_paralogues function.
    #    input_tree_obj = Tree(input_tree_string) 

    #    # Get output from the get_nodes_with_paralogues function.
    #    output_nodes = get_nodes_with_paralogues(input_tree_obj)

    #    # Get the set of leaf names for each output TreeNode object.
    #    output_leaf_name_sets = self.get_output_leaf_name_sets(output_nodes)

    #    # Define the sets of leaf names that should the output TreeNode objects
    #    # should have.
    #    correct_sets = set([frozenset(['A__3', 'A__1', 'B__1']),
    #                        frozenset(['B__3', 'A__2', 'B__2'])])
    #    print('\ncorrect_sets')
    #    print(correct_sets)

    #    print('\noutput_leaf_name_sets')
    #    print(output_leaf_name_sets

    #    # Check that the output set is the expected set (set of sets).
    #    self.assertEqual(output_leaf_name_sets, correct_sets)


    def test_that_relevant_subtrees_not_ignored(self):
        
        # Define a simple tree with two clades that should both be identified
        # as containing paralogues.
        input_tree_string = '((A__1, B__1),(A__2, (B__2, ((C__1, D__1),((C__2, D__2),(C__3, D__3))))));' 

        # Parse as an ete3 TreeNode object for input to the
        # get_nodes_with_paralogues function.
        input_tree_obj = Tree(input_tree_string) 

        # Get output from the get_nodes_with_paralogues function.
        output_nodes = get_nodes_with_paralogues(input_tree_obj)

        # Get the set of leaf names for each output TreeNode object.
        output_leaf_name_sets = self.get_output_leaf_name_sets(output_nodes)

        # Define the sets of leaf names that should the output TreeNode objects
        # should have.
        correct_sets = set([frozenset(['A__1', 'B__1']),
                            frozenset(['A__2', 'B__2', 'C__1', 'D__1', 'C__2', 'D__2', 'C__3', 'D__3']),
                            frozenset(['C__1', 'D__1']),
                            frozenset(['C__2', 'D__2']),
                            frozenset(['C__3', 'D__3'])
                            ])

        # Check that the output set is the expected set (set of sets).
        self.assertEqual(output_leaf_name_sets, correct_sets)





if __name__ == '__main__':
    unittest.main()
