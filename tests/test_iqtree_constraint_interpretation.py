#!/usr/bin/env python3 -W ignore::DeprecationWarning
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

The only unusual package required is ete3. This and it's dependencies can be
installed using appropriate package managers. See the ETE3 website:
    http://etetoolkit.org/new_download/
"""
import os
import time
import shutil
from random import shuffle
import subprocess
import unittest
from ete3 import Tree


def run_simple_constrained_iqtree_analysis(alignment_str,
                                           constraint_tree_str):
    """Takes a string composing the contents of a multiple sequence alignment
    file in PHYLIP format, writes it to a file, and runs an IQtree analysis
    that is constained by a given tree (given as an ete3 Tree object). Returns
    the topology of the output tree as a string in newick format.
    """
    # Make a temporary directory to store output files.
    timestamp = time.strftime("%Y%m%d%H%M%S")
    chars = ['A','B','C','D','E','F','G','H','I','G','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    shuffle(chars)
    random_char = chars[0]
    outdirpath = 'TEMP_CONSTRAINT_TREE_TEST_OUTPUT_' + timestamp + random_char
    os.mkdir(outdirpath)

    # Write input multiple sequence alignment string to a phylip alignment
    # file.
    alignment = os.path.join(outdirpath, 'alignment.phy')
    with open(alignment, 'w') as o:
        o.write(alignment_str)

    # Write constraint tree to a file.
    constraint_tree_obj = Tree(constraint_tree_str)
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

    # Re-root the tree on a particular terminal node so that the output is
    # consistent between runs of IQ-TREE.
    t1.set_outgroup("A")

    # Define the tree topology as a string in newick format.
    topology_as_string = t1.write(format=9)

    # Delete temporary directory.
    shutil.rmtree(outdirpath)

    # Return the output topology.
    return topology_as_string


def run_special_constrained_iqtree_analysis(seqs,
                                           constraint_tree_obj,
                                           starting_tree_obj):
    """Takes a string composing the contents of a multiple sequence alignment
    file in PHYLIP format, writes it to a file, and runs an IQtree analysis
    that is constained by a given tree (given as an ete3 Tree object), starting
    with a given starting tree topology. Returns the topology of the output
    tree as a string in newick format.  """
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

    # Write starting tree to a file.
    starting_tree_file = os.path.join(outdirpath, 'starting_tree.tre')
    starting_tree_obj.write(format=9, outfile=starting_tree_file)

    # Run a tree search with IQtree.
    output_prefix = os.path.join(outdirpath, 'output')
    iqtree_command_list = ['iqtree',
                           '-s', alignment,
                           '-m', 'LG',
                           '-g', constraint_tree_file,
                           '-t', starting_tree_file,
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



class TestConstInterp(unittest.TestCase):
    """Contains methods for performing unit tests using the unittest module.
    """

    def test_that_constraint_works(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 4 10',
                                      'A  VLVLVLLLLL',
                                      'B  VLVLFLLLLL',
                                      'C  FLFLVLLLLL',
                                      'D  FLFLFLLLLL'
                                      ])

        # Define a constraint tree topology in newick format, and parse it
        # using ete3.
        constraint_tree_string = "(A, C, (B, D));"

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)

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

        # Define a constraint tree topology in newick format, and parse it
        # using ete3.
        constraint_tree_string = "(A, B, C, D);"

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)

        # Check that output topology is the same as the constraint topology.
        self.assertNotEqual(output_topo.replace(' ', ''), constraint_tree_string.replace(' ', ''))

    def test_that_an_outgroup_not_in_constraint_can_go_anywhere(self):
        """This method tests whether a taxon not specified in the constraint
        tree topology can be placed as sister to taxa nested within the
        constraint topology
        """
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 5 10',
                                      'A  LLLLLLLLLL',
                                      'B  LLLLLLLLLL',
                                      'C  FFFFFVLLLL',
                                      'D  FFFFFLLLLL',
                                      'X  FFFFFVLLLL'
                                      ])

        # Given this alignment, the actual topology of the tree should be:
        # "(A,B,((C,X),D));"

        # Define a constraint tree topology in newick format, and parse it
        # using ete3.
        constraint_tree_string = "(A, B, (C, D));"

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)

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

        # Define a constraint tree topology in newick format, and parse it
        # using ete3.

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

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)

        # Check that output topology is the same as the constraint topology.
        self.assertEqual('(A,((B,(E,F)),((C,G),(D,H))));', output_topo.replace(' ', ''))


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

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)
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


    #def test_unrooted_subtrees_can_be_rooted_on_any_branch(self):
    #    # Define a string constituting an alignment in phylip format.
    #    alignment_string = '\n'.join([' 8 10',
    #                                  'A  LLLLLLLLLL',
    #                                  'B  LLLLLLLLAA',
    #                                  'C  LLLLLLAAAA',
    #                                  'D  LLLLLLAAAA',
    #                                  'E  VVVVVVVVVV',
    #                                  'F  FFVVVVVVVV',
    #                                  'G  FFFVVVVVVV',
    #                                  'H  FFFVVVVVVV',
    #                                  ])

    #    # Without constraint produces:  (A,(B,(C,D)),(E,(F,(G,H))));
    #    #constraint_tree_string = "(A, B, C, D, E, F, G, H);"
    #    constraint_tree_string = "((C, D, (A, B)), E, F, (G, H));"
    #    constraint_tree_obj = Tree(constraint_tree_string)

    #    # Call function.
    #    output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
    #                                                         constraint_tree_obj)
    #    print('\n\n')
    #    print(output_topo)
    #    print('\n\n')

    #    # Check that output topology is the same as the constraint topology.
    #    self.assertEqual('(A,(B,(C,D)),(E,(F,(G,H))));', output_topo.replace(' ', ''))

    #    # The failure of this test indicates that unrooted subtrees cannot be
    #    # rooted on any clade of interest.

    def test_constrained_analysis_with_starting_tree(self):
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

        # Define starting tree.

        #starting_tree_string = "(A, B, C, D, E, F, G, H);"
        # produces: (A,B,C,D,E,F,G,H); (Error)

        #starting_tree_string = "((A, B, (C, D)), E, F, (G, H));"
        # produces: (A,B,(C,D),(E,F,(G,H))); (Error)

        # *** Why does the output tree always match the starting tree exactly,
        # even when it is multifurcating??? There's actually an error. It
        # doesn't go beyond the starting tree.

        #starting_tree_string = "((A, (B, (C, D))), E, (F, (G, H)));"
        # produces: (A,(B,(C,D)),(E,(F,(G,H)))); (No error here!)

        #starting_tree_string = "((A, (B, (C, D))), (E, (F, (G, H))));"
        # produces: (A,(B,(C,D)),(E,(F,(G,H)))); (No error)

        #starting_tree_string = "((C, (D, (A, B))), (G, (H, (E, F))));"
        # produces: "ERROR: Initial tree is not compatible with constraint
        # tree"

        # Starting tree that is not identical to the constraint tree, but is
        # compatible?
        starting_tree_string = "((F, (E, (G, H))), (B, (A, (C, D))));"
        # produces: (A,(B,(C,D)),(E,(F,(G,H)))); (Works)

        # Parse starting tree string.
        starting_tree_obj = Tree(starting_tree_string)

        # Call function.
        output_topo = run_special_constrained_iqtree_analysis(alignment_string,
                                                              constraint_tree_obj,
                                                              starting_tree_obj)

        # Check that output topology is the same as the constraint topology.
        self.assertEqual('(A,(B,(C,D)),(E,(F,(G,H))));', output_topo.replace(' ', ''))

        # So, to make a usable starting tree, just need to root each of the two
        # subtrees (so that the starting topology is a fully resolved,
        # bifurcating tree) in a way that is compatible with the constraint
        # tree. Maybe this can be done with a script.


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

        # Produces: (A,(B,(C,D)),(E,(F,(G,H)))); with the following constraint:
        constraint_tree_string = "((A, B, (C, D)), E, F, (G, H));"

        #constraint_tree_string = "((E, F, (C, D)), A, B, (G, H));"

        #constraint_tree_string_1 = "(A, B, (C, D));"
        #constraint_tree_string_2 = "(E, F, (G, H));"

        # Call function.
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)
        #print('\n')
        #print(output_topo)
        #print('\n')

        # Check that output topology is the same as the constraint topology.
        #self.assertEqual('(A,(B,(C,D)),(E,(F,(G,H))));', output_topo.replace(' ', ''))
        self.assertEqual('(A,((B,(C,D)),(E,(F,(G,H)))));', output_topo.replace(' ', ''))


    #def test_unrooted_subtrees_are_rootedXXX(self):
    #    # Define a string constituting an alignment in phylip format.
    #    alignment_string = '\n'.join([' 12 10',
    #                                  'A  LLLLLLLLLL',
    #                                  'B  LLLLLLLLAA',
    #                                  'C  LLLLLLLAAA',
    #                                  'D  LLLLLLLAAA',
    #                                  'E  VVVVVVVVVV',
    #                                  'F  FFVVVVVVVV',
    #                                  'G  FFFVVVVVVV',
    #                                  'H  FFFVVVVVVV',
    #                                  'I  KKKKKKKKKK',
    #                                  'J  DDKKKKKKKK',
    #                                  'K  DDDKKKKKKK',
    #                                  'L  DDDKKKKKKK',
    #                                  ])

    #    # Without constraint produces: (A,(B,(C,D)),(E,(F,(G,H))));

    #    # Produces: (A,(B,(C,D)),(E,(F,(G,H)))); with the following constraint:
    #    constraint_tree_string = "((A, B, (C, D)), E, F, (G, H));"

    #    #constraint_tree_string = "((E, F, (C, D)), A, B, (G, H));"

    #    #constraint_tree_string_1 = "(A, B, (C, D));"
    #    #constraint_tree_string_2 = "(E, F, (G, H));"

    #    # Call function.
    #    output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
    #                                                         constraint_tree_string)
    #    #print('\n')
    #    #print(output_topo)
    #    #print('\n')

    #    # Check that output topology is the same as the constraint topology.
    #    #self.assertEqual('(A,(B,(C,D)),(E,(F,(G,H))));', output_topo.replace(' ', ''))
    #    self.assertEqual('(A,((B,((C,D),(I,(J,(K,L))))),(E,(F,(G,H)))));', output_topo.replace(' ', ''))


    def test_two_constraint_trees(self):
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

        ## No constraints.
        #constraint_tree_string = "(A, B, C, D, E, F, G, H);"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string)
        #print('\n\n#######################\nNo constraints\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)

        ## First subtree constrained.
        constraint_tree_string_1 = "(A,C,(B,D));"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string_1)
        #print('\n\n#######################\nFirst subtree constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string_1)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)

        ## Second subtree constrained.
        constraint_tree_string_2 = "(E,G,(F,H));"
        #output_topo_2 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string_2)
        #print('\n\n#######################\nSecond subtree constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string_2)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_2)
        #print(t1)

        # Both constraints.
        constraint_tree_string = constraint_tree_string_1 + '\n' + constraint_tree_string_2
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)
        #print('\n\n#######################\nBoth subtrees constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo)
        #print(t1)

        ## Print tree to illustrate expected results with both constraints.
        #print('\nExpected topology when both subtrees are constrained:')
        #t2 = Tree("(A,((C,(B,D)),(E,(G,(F,H)))));")
        #print(t2)


        # Check that output topology for clade 1 is the same as the constraint topology.
        self.assertTrue(\
                '(B,D)' in output_topo.replace(' ', '') \
                or '(A,C)' in output_topo.replace(' ', '') #\
                )

        # Check that output topology for clade 2 is the same as the constraint topology.
        self.assertTrue(\
                '(F,H)' in output_topo.replace(' ', '') \
                or '(E,G)' in output_topo.replace(' ', '')
                )


    def test_multiple_constraint_trees(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 12 15',
                                      'A  AAAAAAAAAAAAAAA',
                                      'B  AAAAAAAAAAWWWWW',
                                      'C  AAAAAAAWWWWWWWW',
                                      'D  AAAAAAAWWWWWWWW',
                                      'E  FDDDDDDDDDDDDDD',
                                      'F  FFDDDDDDDDDDDDD',
                                      'G  FFFDDDDDDDDDDDD',
                                      'H  FFFDDDDDDDDDDDD',
                                      'I  KEEEEEEEEEEEEEE',
                                      'J  KKEEEEEEEEEEEEE',
                                      'K  KKKEEEEEEEEEEEE',
                                      'L  KKKEEEEEEEEEEEE',
                                      ])

        ## No constraints.
        #constraint_tree_string = "(A, B, C, D, E, F, G, H, I, J, K, L);"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string)
        #print('\n\n#######################\nNo constraints\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)

        ## First subtree constrained.
        constraint_tree_string_1 = "(A,C,(B,D));"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string_1)
        #print('\n\n#######################\nFirst subtree constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string_1)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)

        ## Second subtree constrained.
        constraint_tree_string_2 = "(E,G,(F,H));"
        #output_topo_2 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string_2)
        #print('\n\n#######################\nSecond subtree constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string_2)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_2)
        #print(t1)

        ## Third subtree constrained.
        constraint_tree_string_3 = "(I,K,(J,L));"
        #output_topo_3 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string_3)
        #print('\n\n#######################\nThird subtree constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string_3)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_3)
        #print(t1)

        # All three constraints.
        constraint_tree_string = constraint_tree_string_1 + '\n' \
                + constraint_tree_string_2 + '\n' \
                + constraint_tree_string_3
        output_topo = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)
        #print('\n\n#######################\nBoth subtrees constrained\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo)
        #print(t1)


        # Check that output topology for clade 1 is the same as the constraint topology.
        self.assertTrue(\
                '(B,D)' in output_topo.replace(' ', '') \
                or '(A,C)' in output_topo.replace(' ', '') #\
                )

        # Check that output topology for clade 2 is the same as the constraint topology.
        self.assertTrue(\
                '(F,H)' in output_topo.replace(' ', '') \
                or '(E,G)' in output_topo.replace(' ', '')
                )

        # Check that output topology for clade 3 is the same as the constraint topology.
        self.assertTrue(\
                '(J,L)' in output_topo.replace(' ', '') \
                or '(I,K)' in output_topo.replace(' ', '')
                )


    def test_possibility_of_different_rooting_of_constrained_subtrees(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 12 15',
                                      'A  AAAAAAAAAAAAAAA',
                                      'B  WAAAAAAAAAAAAAA',
                                      'C  WWAAAAAAAAAAAAA',
                                      'D  WWWAAAAAAAAAAAA',
                                      'E  WWWWAAAAAAAAAAA',
                                      'F  WWWWWAAAAAAAAAA',
                                      'G  WWWWWWAAAAAAAAA',
                                      'H  WWWWWWWAAAAAAAA',
                                      'I  WWWWWWWWAAAAAAA',
                                      'J  WWWWWWWWWAAAAAA',
                                      'K  WWWWWWWWWWAAAAA',
                                      'L  WWWWWWWWWWWAAAA',
                                      ])

        ## No constraints.
        #constraint_tree_string = "(A, B, C, D, E, F, G, H, I, J, K, L);"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string)
        #print('\n\n#######################\nNo constraints\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)

        ## Multiple small constraint trees.
        ##constraint_tree_string = ';\n'.join(["(
        #constraint_tree_string = "((L,J),H,(I,K));"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                       constraint_tree_string)
        #print('\n\n#######################\nMultiple small constraint trees\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)


        ## Define a string constituting an alignment in phylip format.
        #alignment_string = '\n'.join([' 12 15',
        #                              'A  AAAAAAAAAAAAAAA',
        #                              'B  WAAAAAAAAAAAAAA',
        #                              'C  WWAAAAAAAAAAAAA',
        #                              'D  WWWAAAAAAAAAAAA',
        #                              'E  WWWWAAAAAAAAAAA',
        #                              'F  WWWWWAAAAAAAAAA',
        #                              'G  WWWWWWAAAAAAAAA',
        #                              'L  WWWWWWWAAAAAAAA',
        #                              'K  WWWWWWWWAAAAAAA',
        #                              'J  WWWWWWWWWAAAAAA',
        #                              'I  WWWWWWWWWWAAAAA',
        #                              'H  WWWWWWWWWWWAAAA',
        #                              ])

        ## Multiple small constraint trees.
        ##constraint_tree_string = ';\n'.join(["(
        #constraint_tree_string = "((L,J),H,(I,K));"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string)
        #print('\n\n#######################\nMultiple small constraint trees\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)


        self.assertTrue(True)


    def test_alternate_rooting_of_constrained_subtrees(self):
        # Define a string constituting an alignment in phylip format.
        alignment_string = '\n'.join([' 10 15',
                                      'A         AAAAAAAAAAAAAAA',
                                      'GolgiQb   WAAAAAAAAAAAAAA',
                                      'Endo1Qb   WWAAAAAAAAAAAAA',
                                      'Endo2Qb   WWWAAAAAAAAAAAA',
                                      'PMQb      WWWWAAAAAAAAAAA',
                                      'PMQc      WWWWWAAAAAAAAAA',
                                      'Endo2Qc   WWWWWWAAAAAAAAA',
                                      'Endo1Qc   WWWWWWWAAAAAAAA',
                                      'GolgiQc   WWWWWWWWAAAAAAA',
                                      'ERQc      WWWWWWWWWAAAAAA',
                                      ])



        # No constraints.
        constraint_tree_string = "(A,GolgiQb,Endo1Qb,Endo2Qb,PMQb,ERQc,GolgiQc,Endo1Qc,Endo2Qc,PMQc);"
        output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
                                                             constraint_tree_string)
        print('\n\n#######################\nNo constraints\n')
        print('Contents of input constraint tree file:')
        print(constraint_tree_string)
        print('Resulting topology:')
        t1 = Tree(output_topo_1)
        print(t1)

        ## Multiple small constraint trees.
        ##constraint_tree_string = ';\n'.join(["(
        #constraint_tree_string = "((L,J),H,(I,K));"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                       constraint_tree_string)
        #print('\n\n#######################\nMultiple small constraint trees\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)


        ## Define a string constituting an alignment in phylip format.
        #alignment_string = '\n'.join([' 12 15',
        #                              'A  AAAAAAAAAAAAAAA',
        #                              'B  WAAAAAAAAAAAAAA',
        #                              'C  WWAAAAAAAAAAAAA',
        #                              'D  WWWAAAAAAAAAAAA',
        #                              'E  WWWWAAAAAAAAAAA',
        #                              'F  WWWWWAAAAAAAAAA',
        #                              'G  WWWWWWAAAAAAAAA',
        #                              'L  WWWWWWWAAAAAAAA',
        #                              'K  WWWWWWWWAAAAAAA',
        #                              'J  WWWWWWWWWAAAAAA',
        #                              'I  WWWWWWWWWWAAAAA',
        #                              'H  WWWWWWWWWWWAAAA',
        #                              ])

        ## Multiple small constraint trees.
        ##constraint_tree_string = ';\n'.join(["(
        #constraint_tree_string = "((L,J),H,(I,K));"
        #output_topo_1 = run_simple_constrained_iqtree_analysis(alignment_string,
        #                                                     constraint_tree_string)
        #print('\n\n#######################\nMultiple small constraint trees\n')
        #print('Contents of input constraint tree file:')
        #print(constraint_tree_string)
        #print('Resulting topology:')
        #t1 = Tree(output_topo_1)
        #print(t1)


        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
