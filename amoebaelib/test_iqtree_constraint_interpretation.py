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
import os
import time
import shutil
from random import shuffle
import subprocess
import unittest
from ete3 import Tree


# Define functions for running simple constrained IQtree analyses.
def run_simple_constrained_iqtree_analysis(seqs,
                                           constraint_tree_obj):
    """Takes a list of seq objects, aligns them, and runs an IQtree analysis
    that is constained by the given tree. Returns the topologies of the output
    trees.
    """
    # Make a temporary directory to store output files.
    timestamp = time.strftime("%Y%m%d%H%M%S")
    chars = ['A','B','C','D','E','F','G','H','I','G','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    shuffle(chars)
    random_char = chars[0]
    outdirpath = 'TEMP_CONSTRAINT_TREE_TEST_OUTPUT_' + timestamp + random_char
    os.mkdir(outdirpath)

    # Write input aligned sequences to a phylip alignment file.
    alignment = os.path.join(outdirpath, 'alignment.phy')
    with open(alignment, 'w') as o:
        o.write(seqs)

    # Write constraint tree to a file.
    constraint_tree_file = os.path.join(outdirpath, 'constraint_tree.tre')
    constraint_tree_obj.write(format=9, outfile=constraint_tree_file)

    # Run a tree search with IQtree.
    output_prefix = os.path.join(outdirpath, 'output')
    iqtree_command_list = ['iqtree',
                           '-s', alignment,
                           '-m', 'LG',
                           '-g', constraint_tree_file,
                           '-pre', output_prefix
                           ]
    #subprocess.call(iqtree_command_list)
    with open(output_prefix + '.stdout.txt', 'w') as o:
        subprocess.call(iqtree_command_list, stdout=o, stderr=subprocess.STDOUT)

    # Parse the relevant output tree file.
    output_tree = output_prefix + '.treefile'
    assert os.path.isfile(output_tree), """IQtree did not produce a treefile
    with the expected path."""
    t1 = Tree(output_tree)
    topology_as_string = t1.write(format=9)

    # Delete temporary directory.
    shutil.rmtree(outdirpath)

    # Return the output topology.
    return topology_as_string



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

    def test_that_an_outgroup_not_in_constraint_can_go_anywhere(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 5 10',
                                      'A  LLLLLLLLLL',
                                      'B  LLLLLLLLLL',
                                      'C  FFFFFVLLLL',
                                      'D  FFFFFLLLLL',
                                      'X  FFFFFVLLLL'
                                      ])

        constraint_tree_string = "(A, B, (C, D));"
        constraint_tree_obj = Tree(constraint_tree_string)

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_obj)

        # Check that output topology is the same as the constraint topology.
        self.assertTrue('(C,X)' in output_topo.replace(' ', ''))


    def test_constrained_subtrees_are_constrained(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 8 10',
                                      'A  LLLLLLLLLL',
                                      'B  LLLLLLLLAA',
                                      'C  LLLLLLLAAA',
                                      'D  LLLLLLLAAA',
                                      'E  VVVVVVVVVV',
                                      'F  FFVVVVVVVV',
                                      'G  FFFVVVVVVV',
                                      'H  FFFVVVVVVV',
                                      ])

        # Without constraint produces: (A,(B,(C,D)),(E,(F,(G,H))));
        #constraint_tree_string = "((A, B), C, D, E, F, G, H);"
        # produces: (A,B,((C,D),(E,(F,(G,H)))));
        #constraint_tree_string = "(A, B, C, D, (E, F), G, H);"
        # produces: (A,(B,(C,D)),((E,F),(G,H)));
        #constraint_tree_string = "((A, (B, (E, F))), C, D, G, H);"
        # produces: (A,(B,(E,F)),((C,D),(G,H)));
        #constraint_tree_string = "((A, (B, (E, F))), (C, G), (D, H));"
        # produces: (A,(B,(E,F)),((C,G),(D,H)));
        constraint_tree_string = "((A, (B, (E, F))), (C, G), (D, H));"
        constraint_tree_obj = Tree(constraint_tree_string)

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_obj)
        #print('\n')
        #print(output_topo)
        #print('\n')

        # Check that output topology is the same as the constraint topology.
        self.assertEqual('(A,(B,(E,F)),((C,G),(D,H)));', output_topo.replace(' ', ''))


    def test_unrooted_subtrees_are_rooted(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 8 10',
                                      'A  LLLLLLLLLL',
                                      'B  LLLLLLLLAA',
                                      'C  LLLLLLLAAA',
                                      'D  LLLLLLLAAA',
                                      'E  VVVVVVVVVV',
                                      'F  FFVVVVVVVV',
                                      'G  FFFVVVVVVV',
                                      'H  FFFVVVVVVV',
                                      ])

        # Without constraint produces: (A,(B,(C,D)),(E,(F,(G,H))));
        constraint_tree_string = "((A, B, (C, D)), E, F, (G, H));"
        constraint_tree_obj = Tree(constraint_tree_string)

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_obj)
        #print('\n')
        #print(output_topo)
        #print('\n')

        # Check that output topology is the same as the constraint topology.
        self.assertEqual('(A,(B,(C,D)),(E,(F,(G,H))));', output_topo.replace(' ', ''))

        # So, all you need to do to construct a constraint tree for pairwise
        # rooting of family trees is make them in the form:
        # "(unrootedfamilytree1, (unrootedfamilytree2));"
        # Where "unrootedfamilytree1" and "unrootedfamilytree2" do not have the
        # outer set of parentheses that would make this constraint tree rooted.



if __name__ == '__main__':
    unittest.main()
