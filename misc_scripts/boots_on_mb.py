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

You have to make sure that the MrBayes tree and the bootstrap files are both in
the input directory.

Make sure that you have defined the name of RAxML (e.g., raxmlHPC-AVX2) in your
settings.py file (and added it to your .bash_profile file).

Usage:
    - Move conversion table file into input directory.
    - Move file with bootstrap pseudoreplicate trees into the input directory.

    boots_on_mb.py <input directory path>

"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
import subprocess
import glob
import re
from ete3 import Tree
from settings import raxmlname
from module_amoebae_name_replace import write_newick_tree_with_uncoded_names

from module_boots_on_mb import reformat_combined_supports, combine_supports,\
mbcontre_to_newick_w_probs, contre_to_newick 


command_line_list = sys.argv
indp = str(command_line_list[1])



## Convert MrBayes output .con.tre file to necessary newick files.

# Convert .tre file to newick format.
assert len(glob.glob(os.path.join(indp, '*.con.tre'))) == 1, """More or less
than one .con.tre file identified in the input directory."""
inbesttreecontre = glob.glob(os.path.join(indp, '*.con.tre'))[0] # Assumes only one .tre file.
outbesttreenewick = inbesttreecontre + '.newick'
contre_to_newick(inbesttreecontre, outbesttreenewick)

# Convert .tre file to newick format with simple percent probability supports
# for internal nodes.
outbesttreenewick2 = inbesttreecontre + '_probs.newick'
mbcontre_to_newick_w_probs(inbesttreecontre, outbesttreenewick2)




# Names for input files.
inboots = os.path.join(indp, 'RAxML_bootstrap.result')
if not os.path.isfile(inboots):
    inboots = os.path.join(indp, 'output.ufboot')
if not os.path.isfile(inboots):
    inboots = os.path.join(indp, 'output.boottrees')
#inbesttree = glob.glob(os.path.join(indp, '*.newick'))[0] # Assumes only one .newick file.
inbesttree = outbesttreenewick 

# Check that the files are present in the input directory.
assert os.path.isfile(os.path.join(indp, inboots)), 'Error: Bootstrap file with filename \
RAxML_bootstrap.result missing from input directory.'
assert inbesttree != ['/'], 'Error'
assert os.path.isfile(os.path.join(indp, inbesttree)), 'Error: Newick tree missing from input directory.'

# Get tablefilename list.
tablefilenames = glob.glob(indp + '/*.table')

# Make sure one and only one table file path was found.
assert len(tablefilenames) > 0, 'Error: No conversion tables identified in\
 directory ' + indp

assert len(tablefilenames) < 2, 'Error: More than one potential conversion\
_table files identified in directory (remove all but one)' + indp + ' \nFiles: ' +\
str(tablefilenames)

# Get tablefilename
tablefilename = tablefilenames[0]

# Addition to output filename.
outname = 'mapped'

# Call RAxML (needs to be in $PATH) to map bootstraps onto newick of mb tree.
subprocess.call([raxmlname, '-f', 'b', '-t', inbesttree, '-z', inboots,\
'-n', outname, '-m', 'PROTGAMMALG4X'], cwd=indp) 

# Output filename.
outfilename = os.path.join(indp, 'RAxML_bipartitions.' + outname)



## Decode RAxML tree.

# name_replace.pl output.
outfilename2 = outfilename + '_DC.tre'

# Call name_replace.pl (have to put the conversion table in the input
# directory).
#subprocess.call(['name_replace.pl', '-b', outfilename, outfilename2,\
#tablefilename], cwd=indp)

# Use module to uncode names.
write_newick_tree_with_uncoded_names(outfilename, outfilename2, tablefilename)
assert os.path.isfile(outfilename2)



## Decode taxon names in .con.tre file.

# Make a new conversion table file with no space characters.
tablefilename2 = tablefilename + '.forMB'
with open(tablefilename) as infh, open(tablefilename2, 'w') as o:
    for i in infh:
        o.write(i.replace(' ', '_'))

# Define name for uncoded file.
outcontrefile = inbesttreecontre + '.DC.tre'

# Call name_replace.pl on .con.tre file.
#subprocess.call(['name_replace.pl', '-b', inbesttreecontre, outcontrefile,\
#tablefilename2], cwd=indp)

# Use module to uncode names.
#write_newick_tree_with_uncoded_names(inbesttreecontre, outcontrefile,
#        tablefilename2)
write_newick_tree_with_uncoded_names(outbesttreenewick, outcontrefile,
        tablefilename2)
assert os.path.isfile(outcontrefile)



## Combine RAxML and MrBayes support values onto a single tree.

# Parse the newick trees with bootstraps and percent probabilities to combine
# node labels.
# newick with bootstrap values: outfilename
# newick with posterior probabilities (as a percentage): outbesttreenewick2
combined_figtree_newick = inbesttreecontre + '_mb+raxml.tre'
combine_supports(outfilename, outbesttreenewick2, combined_figtree_newick)

# Call name_replace.pl (have to put the conversion table in the input
# directory).
combined_figtree_newick_dc = combined_figtree_newick.rsplit('.', 1)[0] + '_DC.tre'
#subprocess.call(['name_replace.pl', '-b', combined_figtree_newick,\
#    combined_figtree_newick_dc, tablefilename], cwd=indp)

# Use module to uncode names.
write_newick_tree_with_uncoded_names(combined_figtree_newick,
        combined_figtree_newick_dc,
        tablefilename)




# Remove unnecessary files.
os.remove(combined_figtree_newick)
os.remove(outbesttreenewick2)

