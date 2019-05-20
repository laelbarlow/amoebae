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
"""Tests functions in the get_alt_topos.py module.
"""

# Import necessary modules.
import unittest
from ete3 import Tree
from get_alt_topos import get_total_num_topos_for_n_taxa,\
    get_all_alt_topologies,\
    get_all_unique_unrooted_bifurcating_topologies_for_n_taxa,\
    trees

# Define a class with functions for performing tests.
class TestGetAltTopos(unittest.TestCase):

    # Test the function for determining how many unique unrooted topologies
    # should be produced in total.
    def test_get_total_num_topos_for_n_taxa_1(self):
        n = 3
        x = get_total_num_topos_for_n_taxa(n)
        correct_answer = 1
        self.assertEqual(x, correct_answer)

    def test_get_total_num_topos_for_n_taxa_2(self):
        n = 4
        x = get_total_num_topos_for_n_taxa(n)
        correct_answer = 3
        self.assertEqual(x, correct_answer)

    def test_get_total_num_topos_for_n_taxa_3(self):
        n = 5
        x = get_total_num_topos_for_n_taxa(n)
        correct_answer = 15
        self.assertEqual(x, correct_answer)

    def test_get_total_num_topos_for_n_taxa_4(self):
        n = 6
        x = get_total_num_topos_for_n_taxa(n)
        correct_answer = 105
        self.assertEqual(x, correct_answer)

    # Test function for recursively generating all alternative sets of sets of
    # elements of the input list of strings.
    def test_trees_1(self):
        leaves = ['A', 'B', 'C']
        all_trees = frozenset(
                {frozenset({tree, leaves[0]}) for tree in trees(set(leaves[1:]))})
        self.assertEqual(len(list(all_trees)),
                get_total_num_topos_for_n_taxa(len(leaves)))

    def test_trees_2(self):
        leaves = ['A', 'B', 'C', 'D']
        all_trees = frozenset(
                {frozenset({tree, leaves[0]}) for tree in trees(set(leaves[1:]))})
        self.assertEqual(len(list(all_trees)),
                get_total_num_topos_for_n_taxa(len(leaves)))

    def test_trees_3(self):
        leaves = ['A', 'B', 'C', 'D', 'E']
        all_trees = frozenset(
                {frozenset({tree, leaves[0]}) for tree in trees(set(leaves[1:]))})
        self.assertEqual(len(list(all_trees)),
                get_total_num_topos_for_n_taxa(len(leaves)))

    def test_trees_4(self):
        leaves = ['A', 'B', 'C', 'D', 'E', 'F']
        all_trees = frozenset(
                {frozenset({tree, leaves[0]}) for tree in trees(set(leaves[1:]))})
        self.assertEqual(len(list(all_trees)),
                get_total_num_topos_for_n_taxa(len(leaves)))

    def test_trees_5(self):
        leaves = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        all_trees = frozenset(
                {frozenset({tree, leaves[0]}) for tree in trees(set(leaves[1:]))})
        self.assertEqual(len(list(all_trees)),
                get_total_num_topos_for_n_taxa(len(leaves)))


    # Test function for generating all unique unrooted topologies as newick
    # tree strings (with arbitrary taxon names).
    def test_get_all_unique_unrooted_bifurcating_topologies_for_n_taxa_1(self):
        get_all_unique_unrooted_bifurcating_topologies_for_n_taxa(3)
    def test_get_all_unique_unrooted_bifurcating_topologies_for_n_taxa_2(self):
        get_all_unique_unrooted_bifurcating_topologies_for_n_taxa(4)
    def test_get_all_unique_unrooted_bifurcating_topologies_for_n_taxa_3(self):
        get_all_unique_unrooted_bifurcating_topologies_for_n_taxa(5)


    # Test function for getting all alternative unrooted topologies. 
    def test_get_all_alt_topologies_1(self):
        # Generate a list of subtrees (ete3 TreeNode objects) as input to
        # function for finding alternative topologies.
        subtrees = [] 
        for i in range(1, 4):
            string = '(C__' + str(i) + ',(A__' + str(i) + ', B__' + str(i) + '));'
            subtrees.append(Tree(string))

        x = get_all_alt_topologies(subtrees)

    def test_get_all_alt_topologies_2(self):
        # Generate a list of subtrees (ete3 TreeNode objects) as input to
        # function for finding alternative topologies.
        subtrees = [] 
        for i in range(1, 5):
            string = '(C__' + str(i) + ',(A__' + str(i) + ', B__' + str(i) + '));'
            subtrees.append(Tree(string))

        x = get_all_alt_topologies(subtrees)

    def test_get_all_alt_topologies_3(self):
        # Generate a list of subtrees (ete3 TreeNode objects) as input to
        # function for finding alternative topologies.
        subtrees = [] 
        for i in range(1, 6):
            string = '(C__' + str(i) + ',(A__' + str(i) + ', B__' + str(i) + '));'
            subtrees.append(Tree(string))

        x = get_all_alt_topologies(subtrees)




if __name__ == '__main__':
    unittest.main()
