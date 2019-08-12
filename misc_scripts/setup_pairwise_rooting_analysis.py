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
"""Sets up input files for tree searches (IQtree and MrBayes) for finding the
root of homologous protein family trees.
"""

# Import built-in modules.
import sys
import os
import argparse
import time
import shutil
import subprocess

# Define paths to directories containing AMOEBAE modules.
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))

# Import functions and classes from AMOEBAE modules.
import settings
from module_amoebae_phylo_clas import ModelInfoFromCSV
from search_alignment_space import get_type_seqs_dict, get_nodes_of_interest,\
     get_corresponding_node
from module_afa_to_nex import nex_to_afa, afa_to_nex
from module_amoebae import mask_nex2
from module_amoebae_trim_nex import trim_nex
from module_amoebae_name_replace import codenames_nex

# Import modules from python3 packages.
from ete3 import Tree

from module_amoebae_name_replace import write_newick_tree_with_coded_names

# Handle command line arguments.
parser = argparse.ArgumentParser(
    description='''Sets up pairwise protein family tree rooting input files.''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('model_name1', help='''Name of model/backbone tree
        to use (other info provided in the model info csv file).''')
parser.add_argument('model_name2', help='''Name of another model/backbone tree
        to use (other info provided in the model info csv file).''')
parser.add_argument('out_dir_path', help='''Path to directory in which
        output directory will be written.''')
parser.add_argument('--polytomy_clades', help='''Make subtrees/clades of
        interest polytomies in constraint topology.''',
        action='store_true')
args = parser.parse_args(sys.argv[1:])

# Get current time for timestamp.
timestamp = time.strftime("%Y%m%d%H%M%S")

# Define output directory path.
main_out_path = os.path.join(args.out_dir_path, 'pairwise_rooting_' + timestamp)

# Make output directory.
os.mkdir(main_out_path)

# Unpack paths from first model.
model_info1 = ModelInfoFromCSV(args.model_name1)
alignment1 = model_info1.alignment_file
tree1 = model_info1.tree_topology_file
subs_model1 = model_info1.subs_model
type_seqs_dict1 = get_type_seqs_dict(model_info1.type_seqs_file)

# Unpack paths from first model.
model_info2 = ModelInfoFromCSV(args.model_name2)
alignment2 = model_info2.alignment_file
tree2 = model_info2.tree_topology_file
subs_model2 = model_info2.subs_model
type_seqs_dict2 = get_type_seqs_dict(model_info2.type_seqs_file)

# Write copies of the alignment files for both family trees.
alignment1_copy = os.path.join(main_out_path, os.path.basename(alignment1))
shutil.copyfile(alignment1, alignment1_copy)
alignment2_copy = os.path.join(main_out_path, os.path.basename(alignment2))
shutil.copyfile(alignment2, alignment2_copy)

# Convert alignments from nexus to aligned fasta format.
alignment1_fasta = alignment1_copy.rsplit('.', 1)[0] + '.afaa'
nex_to_afa(alignment1_copy, alignment1_fasta)
alignment2_fasta = alignment2_copy.rsplit('.', 1)[0] + '.afaa'
nex_to_afa(alignment2_copy, alignment2_fasta)

# Align the two files to each other.
#alignment_combined_fasta = alignment1_copy.rsplit('.', 1)[0] + '__'\
#                           + os.path.basename(alignment2_copy).rsplit('.', 1)[0] + '.afaa'
alignment_combined_fasta = os.path.join(main_out_path, 'combined') + '.afaa'
subprocess.call(['muscle',
                 '-profile',
                 '-in1', alignment1_fasta,
                 '-in2', alignment2_fasta,
                 '-out', alignment_combined_fasta
                 ])

# Convert afaa file to nexus format.
alignment_combined_nex = alignment_combined_fasta.rsplit('.', 1)[0] + '.nex'
afa_to_nex(alignment_combined_fasta, alignment_combined_nex)

# May not be necessary:
## Mask (include all positions) and Trim the combined alignment file.
#alignment_combined_mask = alignment_combined_nex.rsplit('.', 1)[0] + '.mask.nex'
#mask_nex2(alignment_combined_nex, alignment_combined_mask)
#alignment_combined_trim = alignment_combined_mask.rsplit('.', 1)[0] + '.trim.nex'
#trim_nex(alignment_combined_mask, alignment_combined_trim)

# Write coded alignment file.
codenames_nex(alignment_combined_nex)

# Define path to table file.
tablefile = alignment_combined_nex.rsplit('.', 1)[0] + '_C.table'

# Copy tree topology files to output directory.
unrooted_topology1 = os.path.join(main_out_path, os.path.basename(tree1))
shutil.copyfile(tree1, unrooted_topology1)
unrooted_topology2 = os.path.join(main_out_path, os.path.basename(tree2))
shutil.copyfile(tree2, unrooted_topology2)

# Generate constraint tree of the form:
# "(familytree1, (familytree2));". Where familytree1 and familytree2 are
# unrooted and do not include the outer set of parentheses that would make this
# constraint tree a rooted tree.
# e.g., "(A, B, (C, D), (E, F, (G, H)));"

# Parse unrooted topology files.
t1 = Tree(unrooted_topology1)
t2 = Tree(unrooted_topology2)

# Initiate variables for storing unrooted topology strings.
unrooted_topology1_string = None
unrooted_topology2_string = None

# Make clades of interest polytomies if specified by the command line option.
if args.polytomy_clades:
    # Make copies of the unrooted tree topologies.
    t1_copy = t1.copy()
    t2_copy = t2.copy()

    # Loop over both trees.
    for t, tsd in zip([t1_copy, t2_copy], [type_seqs_dict1, type_seqs_dict2]):
        print(t)
        # Get a list of nodes of interest.
        orthogroup_nodes = get_nodes_of_interest(t, tsd.values())

        # Iterate through corresponding nodes in the copies of the tree objects.
        for orthogroup_node in [x[1] for x in orthogroup_nodes]:
            print('orthogroup_node:')
            print(orthogroup_node)
            # Make a list of leaf node objects in clade of interest.
            leaf_nodes = orthogroup_node.get_leaves()

            # Remove all child nodes of clade of interest.
            child_nodes_to_remove = []
            child_num = 0
            for child_node in orthogroup_node.traverse():
                child_num += 1
                if child_num > 1:
                    child_nodes_to_remove.append(child_node)
            for child_node in child_nodes_to_remove:
                orthogroup_node.remove_child(child_node)
            print('orthogroup_node without children:')
            print(orthogroup_node)
            ...XXX...
                    


            # Add all original leaf node objects as child nodes to the clade of
            # interest (to make a polytomy).

    # Get topology strings without outer set of parentheses from modified
    # topologies.
    unrooted_topology1_string = t1_copy.write(format=9)[1:-2]
    unrooted_topology2_string = t2_copy.write(format=9)[1:-2]

else:
    # Just get topology strings without outer set of parentheses.
    #unrooted_topology1_string = open(unrooted_topology1, 'r').read().strip()[1:-2]
    #unrooted_topology2_string = open(unrooted_topology2, 'r').read().strip()[1:-2]
    unrooted_topology1_string = t1.write(format=9)[1:-2]
    unrooted_topology2_string = t2.write(format=9)[1:-2]

# Construct constraint topology string.
constraint_tree_string = '(' + unrooted_topology1_string + ',(' +\
                         unrooted_topology2_string + '));'
constraint_tree_file = os.path.join(main_out_path, 'constraint.tre')
with open(constraint_tree_file, 'w') as o:
    o.write(constraint_tree_string)

# Generate starting tree of the form:
# "((familytree1), (familytree2));". Where familytree1 and familytree2 are
# rooted on one of the clades of interest in each (compatible with the
# constraint tree).
# e.g., "((B, (A, (C, D))), (F, (E, (G, H))));"

# Identify one of the clades of interest in the trees.
orthogroup_nodes1 = get_nodes_of_interest(t1, type_seqs_dict1.values())
orthogroup_nodes2 = get_nodes_of_interest(t2, type_seqs_dict2.values())

# Root each tree object on a clade of interest (to get a rooted tree topology
# that is consistent with the constraint tree).
t1.set_outgroup(get_corresponding_node(orthogroup_nodes1[0][1], t1))
t2.set_outgroup(get_corresponding_node(orthogroup_nodes2[0][1], t2))

# Get newick strings for rooted trees.
rooted_topology1_string = t1.write(format=9)[1:-2]
rooted_topology2_string = t2.write(format=9)[1:-2]

# Construct combined starting tree string.
starting_tree_string = '((' + rooted_topology1_string + '),(' +\
                       rooted_topology2_string + '));'

# Write starting tree string to a file.
starting_tree_file = os.path.join(main_out_path, 'starting.tre')
with open(starting_tree_file, 'w') as o:
    o.write(starting_tree_string)

# Code names in constraint tree.
constraint_tree_file_coded = constraint_tree_file.rsplit('.', 1)[0] + '_C.tre'
write_newick_tree_with_coded_names(constraint_tree_file,
        constraint_tree_file_coded, tablefile)

# Code names in starting tree.
starting_tree_file_coded = starting_tree_file.rsplit('.', 1)[0] + '_C.tre'
write_newick_tree_with_coded_names(starting_tree_file,
        starting_tree_file_coded, tablefile)
 

# Write script to perform an IQtree bootstrap analysis using the constraint and
# starting trees.

# Apply constraints and starting tree to MrBayes code block.





