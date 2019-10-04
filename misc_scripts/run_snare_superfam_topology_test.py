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
"""Runs constrained ML searches and topology tests on the resulting trees using
IQ-tree.
"""
import os
import sys
import glob
import subprocess
from datetime import datetime

# Take coded phylip alignment file as input.
cmdln = sys.argv
phylip = cmdln[1]
model = cmdln[2]
constraint_indir = cmdln[3]

assert os.path.isfile(phylip)
assert os.path.isdir(constraint_indir)

# Make an output directory with a unique name.
now = datetime.now()
timestamp = datetime.timestamp(now)
outputdir = os.path.join(constraint_indir, 'topo_test_output_' +\
        str(timestamp))
os.mkdir(outputdir)

# Get best model with ModelFinder.
# ...
#model = 'LG+F+G4'

# Run constrained tree searches
constraint_tree_files = glob.glob(os.path.join(constraint_indir,
    '*_clean.C.tre'))
for i in constraint_tree_files:
    print('Running ML search with constraint tree ' + i)
    subprocess.call(['iqtree',
                    '-s', phylip,
                    '-m', model,
                    '-g', i,
                    '-t', i.rsplit('.', 1)[0] + '_starting.tre',
                    '-pre', os.path.join(outputdir, os.path.basename(i)),
                    '-allnni',
                    '--runs', '5'
                    ])

# Concatenate treefiles.
output_treefiles = glob.glob(os.path.join(outputdir, '*.treefile'))
concat_path = os.path.join(outputdir, 'concatenated.treefile')
tree_order_file = concat_path + '_tree_order.txt'
with open(concat_path, 'w') as o1, open(tree_order_file, 'w') as o2:
    treenum = 0
    for i in output_treefiles:
        treenum += 1
        # Append tree to concatenated tree file.
        o1.write(open(i, 'r').read().strip() + '\n')
        # Record in tree order file.
        o2.write('Tree ' + str(treenum) + ': ' + os.path.basename(i) + '\n')

# Run topology test.
topo_test_out = concat_path + '_topo_test_output'
subprocess.call(['iqtree',
                 '-s', phylip,
                 '-m', model,
                 '-z', concat_path,
                 '-n', '0',
                 '-zb', '10000',
                 '-au', 
                 '-pre', topo_test_out
                 ])
                 
print('Topology test results written to ' + topo_test_out + '.iqtree')
                 

