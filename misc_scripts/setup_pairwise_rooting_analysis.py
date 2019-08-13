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
from string import Template

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
from module_amoebae_constrain_mb import constrain_mb_with_tree
from module_amoebae_name_replace import write_newick_tree_with_coded_names

# Import modules from python3 packages.
from ete3 import Tree


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

# Define path to coded phylip alignment file.
phylip_coded = alignment_combined_nex.rsplit('.', 1)[0] + '_C.phy'

# Define path to coded nexus for MrBayes alignment file.
mb_coded = alignment_combined_nex.rsplit('.', 1)[0] + '_C.mb.nex'

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

# Identify one of the clades of interest in the trees.
orthogroup_nodes1 = get_nodes_of_interest(t1, type_seqs_dict1.values())
orthogroup_nodes2 = get_nodes_of_interest(t2, type_seqs_dict2.values())

# Root each tree object on a clade of interest.
t1.set_outgroup(get_corresponding_node(orthogroup_nodes1[0][1], t1))
t2.set_outgroup(get_corresponding_node(orthogroup_nodes2[0][1], t2))

# Unroot each tree object in such a way that it will be compatible with the
# starting tree.
for t in [t1, t2]:
    child_num = 0
    for child_node in t.traverse():
        child_num += 1
        if child_num == 3:
            child_node.delete()

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
        #print(t)
        # Get a list of nodes of interest.
        orthogroup_nodes = get_nodes_of_interest(t, tsd.values())

        # Iterate through corresponding nodes in the copies of the tree objects.
        for orthogroup_node in [x[1] for x in orthogroup_nodes]:
            #print('orthogroup_node:')
            #print(orthogroup_node)

            # Get corresponding node in actual tree object.
            orthogroup_node = get_corresponding_node(orthogroup_node, t)

            # Remove all child nodes that are not leaves from clade of interest
            # (making the clade a polytomy).
            child_nodes_to_remove = []
            child_num = 0
            for child_node in orthogroup_node.traverse():
                child_num += 1
                if child_num > 1:
                    if not child_node.is_leaf():
                        child_nodes_to_remove.append(child_node)
            for child_node in child_nodes_to_remove:
                child_node.delete()

            #print('orthogroup_node as polytomy:')
            #print(orthogroup_node)

        #print(t)
                    
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

# Root each tree object on a clade of interest.
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
 
# Perform model selection using IQtree ModelFinder.
subs_model = None
model_finder_dir = os.path.join(main_out_path, 'ModelFinder_output')
os.mkdir(model_finder_dir)
prefix = os.path.join(model_finder_dir, 'output')
iqtree_command_list = ['iqtree',
                       '-s', phylip_coded,
                       '-m', 'MFP',
                       '-mset', 'LG,WAG,JTT,Dayhoff,VT',
                       '-pre', prefix 
                       ]
#subprocess.call(iqtree_command_list)
#file_with_model = prefix + '.iqtree'
#with open(file_with_model) as fh:
#    for i in fh:
#        if i.startswith('Best-fit model'):
#            subs_model = i.split(' ')[2]
#            break

# Check that a model was found.
subs_model = 'LG'
assert subs_model is not None

# Change substitution model in MrBayes file. 
mbmodelregex = re.compile(r'aamodelpr=fixed(.+)')
temp_mb_file = mb_coded + '_TEMP'
with open(mb_coded,'r') as infh, open(temp_mb_file, 'w') as o:
    contents = infh.read()
assert os.path.isfile(temp_mb_file)
os.remove(mb_coded)
shutil.copyfile(temp_mb_file, mb_coded)
os.remove(temp_mb_file)


# Write script to perform an IQtree bootstrap analysis using the constraint and
# starting trees.
iqtree_script = os.path.join(main_out_path, '0_run_constrained_iqtree.sh')
with open(iqtree_script, 'w') as o:
    script_template = Template("""\
#!/usr/bin/env bash
mkdir IQtree_output
iqtree -s $alignment_var1 -m $subs_model_var -g $constraint_var -t $starting_var -pre $alignment_var2 -nt AUTO -bb 1000 -alrt 1000
""")
    script_text = script_template.substitute(alignment_var1 = os.path.basename(phylip_coded),
                                             alignment_var2 =\
                                             os.path.basename(phylip_coded)\
                                             + '_iqtree_output/output',
                                             subs_model_var = subs_model,
                                             constraint_var =\
                                             os.path.basename(constraint_tree_file_coded),
                                             starting_var =\
                                             os.path.basename(starting_tree_file_coded)
                                             )
    o.write(script_text)

# Make a subdirectory to contain MrBayes output files.
mb_subdir = os.path.join(main_out_path, 'MrBayes_output')
os.mkdir(mb_subdir)

## Copy constraint tree into MrBayes subdir.
#shutil.copyfile(constraint_tree_file_coded, os.path.join(mb_subdir,
#    os.path.basename(constraint_tree_file_coded)))

# Apply constraints to MrBayes code block.
# Note: MrBayes appears not to have a way of defining a starting tree topology.
mb_coded_constrained = os.path.join(mb_subdir,
        os.path.basename(mb_coded).rsplit('.', 1)[0] + '.constrained.nex')
constrain_mb_with_tree(mb_coded,
                       constraint_tree_file_coded,
                       mb_coded_constrained
                       )


# Write script for performing constrained tree search with MrBayes.
mrbayes_script = os.path.join(main_out_path, '0_run_constrained_mrbayes.sh')
with open(mrbayes_script, 'w') as o:
    script_template = Template("""\
#!/usr/bin/env bash
cd $mb_subdir_var
mb $mb_nex_var
"""
    )
    script_text = script_template.substitute(mb_nex_var =\
                     os.path.basename(mb_coded_constrained),
                     mb_subdir_var = os.path.basename(mb_subdir))
    o.write(script_text)



