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
"""Module with functions to re-root newick trees, which can be tested to see
whether support values in node labels are properly interpreted as branch
attributes during the re-rooting process. This may not be the case, as
discussed by Czech et al. 2017.
"""
import sys
from ete3 import Tree, TreeStyle, TextFace

cmd = sys.argv
treefp = cmd[1]

t1 = Tree(treefp)#, quoted_node_names=True, format=1)
print(t1)
print(t1.write())

ts = TreeStyle()
ts.show_branch_support = True
t1.show(tree_style=ts)

print('rerooted on leaf node X:')

t2 = t1.copy()
t2.set_outgroup('X')
print(t1)
print(t1.write())

ts = TreeStyle()
ts.show_branch_support = True
t2.show(tree_style=ts)

# Appears to be correct to this point.
# Maybe the problem arises because I use node names instead of branch support?
# What exactly did I do in amoebae visualize_tree? I used node names as branch
# support instead of branch support, because branch supports can only be a
# single float in ete3.

print('Showing node names on branches as "support":')

t3 = t2.copy()
for node in t3.traverse():
    if not node.is_leaf():
        print(node.name)
        node.add_face(TextFace(node.name), column=0, position='branch-top')
# Show tree.
ts = TreeStyle()
ts.show_branch_support = True
t3.show(tree_style=ts)

# The problem is that I had to parse the tree with "format=1", so it takes the
# node support values as node names. Perhaps the node names are not translated
# as attributes of branches when rerooting?

print('Parsing tree with format=1, so branch supports should be interpreted as node names')
t4 = Tree(treefp, format=1)
for node in t4.traverse():
    if not node.is_leaf():
        node.add_face(TextFace(node.name), column=0, position='branch-top')

ts = TreeStyle()
ts.show_branch_support = True
t4.show(tree_style=ts)

print('rerooting tree')
t4.set_outgroup('X')
t4.show(tree_style=ts)

# This reproduces the error that I observed before when running amoebae
# visualize_tree command. The branch supports are not translated correctly when
# rerooting, when they are parsed as node names with format=1.

# Possible solutions:
#     1. Modify support values to be a single float value (only one source of
#     support) and parse with Tree(xxx, format=0), then add additional support
#     values back in later, if necessary, and only set up TextFaces at the end
#     after all rerooting has been done.

#     *2. DONE - Parse with format=1,
#         then manually set node support values by parsing out one of the numbers.
#         Only after doing that, then re-root.
#         Then, add branch supports as TextFaces to nodes.
#         Then render.

#         Also, write functions to combine multiple supports onto a single topology as
#         TextFaces before rendering.



