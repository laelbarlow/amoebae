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

"""Contains functions for replacing nodes of interest in a given ete3 TreeNode
object with alternative nodes.

***This module is not in use yet.
"""
import sys
import os
#from math import factorial
from ete3 import Tree
from search_alignment_space import get_corresponding_node
#from itertools import product

def replace_nodes_of_interest_in_tree(t, nodes, replacement_nodes):
    """Takes a TreeNode object, a list of nodes (TreeNode objects) of interest,
    and a list of replacement nodes, and makes the appropriate node
    replacements.
    """
    t2 = t.copy()
    print(t2)
    for n, replacement_n in zip(nodes, replacement_nodes):
        t2n = get_corresponding_node(n, t2)
        # Replace t2n with replacement_n
        ...XXX...


    return t2


#def get_polytomy_for_treenode(treenode):
#    """Take an ETE3 TreeNode object and return a TreeNode object with all the
#    same leaf nodes, but just as a polytomy.
#    """
#    # Initiate a new TreeNode object.
#    new_treenode = Tree()
#
#    # Add each terminal (leaf) node from the original TreeNode object as a
#    # child to the new TreeNode object. 
#    for leaf in treenode.get_leaves():
#        new_treenode.add_child(leaf)
#
#    # Return polytomy.
#    return new_treenode
#
#
#def determine_if_polytomy(subtree):
#    """Determine if the given subtree (represented as a set, which may contain
#    sets) is a polytomy or not.
#    """
#    is_polytomy = False
#    # Basically, input sets cannot contain more than two elements. Otherwise,
#    # they represent a polytomy.
#    if len(subtree) > 2:
#        is_polytomy = True
#
#    # Return True if it's a polytomy.
#    return is_polytomy
#
#
#def get_partitions_of_set(s):
#    """Recursively partition elements in an input set.
#    """
#    # Add a new singleton subset.
#    if len(s) == 1:
#        yield frozenset(s)
#        return
#    # Extract one element from the set
#    elem, *_ = s
#    rest = frozenset(s - {elem})
#    for partition in get_partitions_of_set(rest):
#        for subset in partition:
#            # Insert the element in the subset
#            try:
#                augmented_subset = frozenset(subset | frozenset({elem}))
#            except TypeError:
#                # subset is actually an atomic element
#                augmented_subset = frozenset({subset} | frozenset({elem}))
#            yield frozenset({augmented_subset}) | (partition - {subset})
#        # Case with the element in its own extra subset
#        yield frozenset({elem}) | partition
#
#
#def get_newick_tree_string_from_set(s):
#    if type(s) not in (set, frozenset):
#        return str(s)
#    return "(" + ",".join(sorted(map(get_newick_tree_string_from_set, s))) + ")"
#
#
#def trees(leaves, polytomies=False):
#    """Generate all alternative trees (in newick format) for a given list of
#    leaf names, by recursively partitioning the input set of leaf names in
#    every possible way. By default polytomies are excluded, but can be included
#    by setting the polytomies option to True.
#    """
#    if type(leaves) not in (set, frozenset):
#        # It actually is a single leaf
#        yield leaves
#        # Don't try to yield any more trees
#        return
#    # Otherwise, we will have to consider all the possible
#    # partitions of the set of leaves, and for each partition,
#    # construct the possible trees for each part
#    for partition in get_partitions_of_set(leaves):
#        # We need to skip the case where the partition
#        # has only one subset (the initial set itself),
#        # otherwise we will try to build an infinite
#        # succession of nodes with just one subtree
#        if len(partition) == 1:
#            part, *_ = partition
#            # Just to be sure the assumption is correct
#            assert part == leaves
#            continue
#        # We recursively apply the trees function to each part and obtain the
#        # possible trees by making the product of the sets of possible
#        # subtrees.
#        for subtree in product(*map(trees, partition)):
#            # Only return trees that are polytomies if the option is used.
#            yield_subtree = False
#            if polytomies:
#                yield_subtree = True
#            else:
#                # Determine whether the subtree is a polytomy.
#                is_polytomy = determine_if_polytomy(subtree)
#                if not is_polytomy:
#                    yield_subtree = True
#            # Now yield the subtree if appropriate.
#            if yield_subtree:
#                yield frozenset(subtree)
#            else:
#                # Ignore topologies that are not interesting.
#                continue
#
#
#def get_total_num_topos_for_n_taxa(n):
#    """Given a number of taxa (integer greater than 2), return an integer
#    representing the theoretical total number of possible topologies (without
#    polytomies) for relationships between those taxa.
#    """
#    assert isinstance(n, int), """Input must be an int."""
#    assert n > 2, """Input number must be greater than 2."""
#    return factorial(2*n - 5) / (factorial(n - 3) * 2**(n - 3))
#
#
#def get_leaf_names_to_represent_n_taxa(n):
#    """Take n integer and return a list of the same number of unique strings
#    that can be easily replaced with other taxon names or subtrees.
#    """
#    longest_leaf_list = [x*10 for x in ["A", "B", "C", "D", "E", "F", "G"]]
#    leaves = longest_leaf_list[:n]
#    return leaves
#
#
#def get_all_unique_unrooted_bifurcating_topologies_for_n_taxa(n):
#    """Given an integer from 3 to 7, generate all unique unrooted bifurcating
#    tree topologies for that given number of taxa.
#    """
#    # Check that input is an integer between 3 and 7.
#    assert isinstance(n, int), """Input must be an int."""
#    assert 2 < n < 8, """Input number must be greater than 2 and less than 8."""
#
#    # Generate a set of leaf names to represent n taxa.
#    leaves = get_leaf_names_to_represent_n_taxa(n)
#
#    # Call recursive function (trees) to generate all unique bifurcating rooted
#    # trees for n taxa (rooted on a particular leaf).
#    all_trees = frozenset(
#            {frozenset({tree, leaves[0]}) for tree in trees(set(leaves[1:]))})
#
#    ## Print resulting set of trees.
#    #print('\nTree topologies for ' + str(n) + ' taxa:')
#    #for tree in all_trees:
#    #    print(get_newick_tree_string_from_set(tree) + ";")
#
#    # Check that the length of the resulting set of topologies is the same as
#    # predicted by the formula.
#    assert len(all_trees) == get_total_num_topos_for_n_taxa(n), """Different
#    number of tree topologies produced than expected."""
#
#    # Return the set of all alternative tree topologies as a list.
#    return list(all_trees)
#
#
#def get_all_alt_topologies(subtrees):
#    """Take a list of ete3 TreeNode objects and return a list of TreeNode
#    objects representing all possible backbone topologies for relationships
#    between the input clades.
#    """
#    # Find all alternative topologies for tree with the same number of taxa as
#    # the number of input subtrees (with taxa represented by arbitrary
#    # strings).
#    all_tree_sets =\
#        get_all_unique_unrooted_bifurcating_topologies_for_n_taxa(len(subtrees))
#
#    # Get all sets as newick tree strings.
#    all_rooted_newick_trees = [get_newick_tree_string_from_set(x) + ';' for x in all_tree_sets]
#
#    ## Print rooted trees.
#    #print('First three rooted trees:')
#    #for i in all_rooted_newick_trees[:3]:
#    #    print(i)
#    #print('')
#
#    # (Re)generate the list of leaf names in the backbone trees. 
#    leaf_names_in_all_trees = get_leaf_names_to_represent_n_taxa(len(subtrees))
#
#    # Unroot all the trees.
#    all_tree_obj = [Tree(x) for x in all_rooted_newick_trees]
#    for t in all_tree_obj:
#        t.unroot()
#    all_unrooted_newick_trees = [t.write(format=9) for t in all_tree_obj]
#
#    ## Print unrooted trees.
#    #print('First three unrooted trees:')
#    #for i in all_unrooted_newick_trees[:3]:
#    #    print(i)
#    #print('')
#
#    # Get list of input subtrees as newick trees (without trailing semicolons).
#    subtree_newick_strings = [subtree.write(format=9) for subtree in subtrees]
#
#    ## Print subtrees.
#    #print('First three subtrees:')
#    #for i in subtree_newick_strings[:3]:
#    #    print(i[:-1])
#    #print('')
#
#    # Get list of alternative topologies by replacing strings in unrooted trees
#    # with strings representing the input subtrees.
#    all_unrooted_trees_with_subtrees = all_unrooted_newick_trees
#    for place_holder_name, subtree_string in\
#            zip(leaf_names_in_all_trees, subtree_newick_strings):
#        #print('Replacing %s with %s in all alternative topologies.' %\
#        #        (place_holder_name, subtree_string[:-1]))
#        all_unrooted_trees_with_subtrees =\
#                [t.replace(place_holder_name, subtree_string[:-1]) for t in\
#                        all_unrooted_trees_with_subtrees] 
#
#    ## Print trees with subtrees.
#    #print('\nFirst three unrooted trees with subtrees:')
#    #for i in all_unrooted_trees_with_subtrees[:3]:
#    #    print(i)
#    #print('')
#
#    # Return list of alternative topologies as newick strings.
#    return all_unrooted_trees_with_subtrees


