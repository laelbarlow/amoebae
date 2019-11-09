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
"""Takes an newick tree file and a list of sequence/taxon names, and removes
all the leaves with the names in the list by using the ETE3 prune function to
keep all the rest.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from ete3 import Tree

command_line_list = sys.argv
infile1 = str(command_line_list[1]) # Newick tree file.
infile2 = str(command_line_list[2]) # Text file listing sequence names on separate lines.
outfile1 = str(command_line_list[3]) # Output tree file.

# Report inputs.
print('\nInput tree file: ' + infile1)
print('\nInput list file: ' + infile2)
print('\nOutput tree file: ' + outfile1)

# Parse the input tree using ETE3.
t1 = Tree(infile1)
print('\nParsed input tree:')
print(t1)

# Count number of leaves in input tree:
num_leaves_1 = len(t1.get_leaves())

# Parse the sequence name list file.
leaf_names_to_remove = []
with open(infile2, 'r') as infh1:
    for i in infh1:
        if not i.startswith('\n'):
            leaf_names_to_remove.append(i.strip())

# Record number of leaf names in input list.
num_leaves_to_remove = len(leaf_names_to_remove)

# Check that at least one leaf name was found.
assert num_leaves_to_remove >= 1, """Could not identify any sequence names
in input list file."""

# Get a list of all the leaf names.
all_leaf_names_in_original_tree = [x.name for x in t1.get_leaves()]

# Check that all the leaf names in the input list are in the input tree.
for leaf in leaf_names_to_remove:
    assert leaf in all_leaf_names_in_original_tree, """Could not identify leaf
    with name %s in input newick tree.""" % leaf

# Get a list of all the leaf names excluding those in the input list. 
reduced_leaf_names = list(set(all_leaf_names_in_original_tree) - set(leaf_names_to_remove))

# Check that it is the expected number of leaf names to include.
assert len(reduced_leaf_names) == len(all_leaf_names_in_original_tree) - len(leaf_names_to_remove)

# Use ETE3 to prune the input tree, only including the specified leaf names.
t2 = t1.copy()
t2.prune(reduced_leaf_names)
print('\nPruned tree:')
print(t2)

# Check that the number of leaf names is reduced to the expected number.
assert len(t2.get_leaves()) == len(reduced_leaf_names)

# Write pruned tree to output path.
t2.write(format=9, outfile=outfile1)





































