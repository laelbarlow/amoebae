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
"""Defines a function that runs simple constrained IQtree tree searches to test
how IQtree interprets constraints.
"""

# Import necessary modules.
import os
import time
import shutil
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
    outdirpath = 'TEMP_CONSTRAINT_TREE_TEST_OUTPUT_' + timestamp
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
    t1 = Tree(output_tree)
    topology_as_string = t1.write(format=9)

    # Delete temporary directory.
    shutil.rmtree(outdirpath)

    # Return the output topology.
    return topology_as_string


