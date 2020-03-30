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
"""Takes a directory and runs RAxML map bootstrap values from one file to a
newick tree in another file.

Assumes specific names for input files.

You have to make sure that the best tree and the bootstrap files are both in
the input directory.

Usage:
    - Move conversion table file into input directory.
    - Move file with bootstrap pseudoreplicate trees into the input directory.

    boots_on_best_ml.py <input directory path>

"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import subprocess
import glob
import settings
# from afa_to_nex import nex_to_afa, afa_to_nex, delete_extra_mesquite_lines
from name_replace import write_newick_tree_with_uncoded_names

command_line_list = sys.argv
indp = str(command_line_list[1])
raxmlname = settings.raxmlname

# Names for input files.
inboots = 'RAxML_bootstrap.result'
inbesttree = 'RAxML_bestTree.result'

# Check that the files are present in the input directory.
assert os.path.isfile(os.path.join(indp, inboots)), 'Error: Bootstrap file with filename\
RAxML_bootstrap.result missing from input directory.'
assert os.path.isfile(os.path.join(indp, inbesttree)), 'Error: Newick tree file with filename\
RAxML_bestTree.result missing from input directory.'

# Get tablefilename list.
tablefilenames = glob.glob(indp + '/*.table')

# Make sure one and only one table file path was found.
assert len(tablefilenames) > 0, 'Error: No conversion tables identified in\
 directory ' + indp

assert len(tablefilenames) < 2, 'Error: More than one potential conversion\
 table files identified in directory '+ indp + ' \nFiles: ' +\
str(tablefilenames)

# Get tablefilename
tablefilename = tablefilenames[0]

# Addition to output filename.
outname = 'mapped'


# Call RAxML.
subprocess.call([raxmlname, '-f', 'b', '-t', inbesttree, '-z', inboots,\
'-n', outname, '-m', 'PROTGAMMALG4X'], cwd=indp) 

# Output filename.
outfilename = 'RAxML_bipartitions.' + outname

# name_replace.pl output.
outfilename2 = outfilename + '_DC.tre'

# Call name_replace.pl (have to put the conversion table in the input
# directory).
#subprocess.call(['name_replace.pl', '-b', outfilename, outfilename2,\
#tablefilename], cwd=indp)

# Use module to uncode names.
write_newick_tree_with_uncoded_names(outfilename, outfilename2, tablefilename)



#if __name__ == '__main__':
