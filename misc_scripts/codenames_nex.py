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
"""Takes an nex alignment file, and generates files for input to phylogenetic
analysis software, including RAxML and MrBayes.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import subprocess
from module_afa_to_nex import nex_to_afa, afa_to_nex, delete_extra_mesquite_lines, nex_to_phylip, nex_to_mbnex
import re
from module_amoebae_name_replace import write_afa_with_code_names

command_line_list = sys.argv
infp = str(command_line_list[1])



# Find the number of taxa in the original file.
ntaxre = re.compile(r"(NTAX|ntax)=\d+")
ntaxa = None
with open(infp) as f:
    for l in f:
        if ntaxre.search(l):
            s = ntaxre.search(l)
            ntaxa = int(s.group().split('=')[1])
            break
assert ntaxa is not None, 'Error: Could not identify number of taxa in input alignment.'

# Make a name for afa file.
tempafa1 = infp.replace('.nex', '.afa_temp')

# Write an afa file from the nex file.
nex_to_afa(infp, tempafa1)

# Make name for temporary afa file with new taxon names.
tempafa2 = tempafa1.replace('.afa_temp', '_C.afa_temp')

# Make a table output file name.
tablefilename = tempafa2.replace('.afa_temp', '.table')
 
# Call name_replace.pl (needs to be in $PATH).
#subprocess.call(['name_replace.pl', '-f', tempafa1, tempafa2, tablefilename]) 

# Code names with module.
write_afa_with_code_names(tempafa1, tempafa2, tablefilename)

# Check that the number of taxa in the input alignment is the same as the
# number listed in the output conversion table file.
num_lines = sum(1 for line in open(tablefilename))
ntaxa_out = int(num_lines / 2)
assert ntaxa_out == ntaxa, 'Error: Different number of taxa in input alignment\
 compared to output conversion table.'

# Check for redundant sequence names in output table file.
lines = [line for line in open(tablefilename)]
old_names = [line.strip() for line in lines[1::2]]
old_name_nonredun = list(set(old_names))
assert len(old_names) <= len(old_name_nonredun), 'Error: Redundant sequences\
 identified in input alignment.'

for name in old_names:
    assert not name.endswith('.copy'), 'Error: duplicate sequences detected.'

# Check that number of sequences is the same in input and in table.
assert len(old_names) == ntaxa, 'Error: Number of taxa in table is different\
 from number of taxa in input alignment.'

# Make name for new nex file with new taxon names.
newnex = tempafa2.replace('.afa_temp', '.nex')
 
# Convert afa file to nex.
afa_to_nex(tempafa2, newnex)

# Convert nex output file to phylip.
newphylip = newnex.rsplit('.', 1)[0] + '.phy'
nex_to_phylip(newnex, newphylip) 


# Convert nex output file to nexus for MrBayes:

# Define name for output file.
mbnex = newnex.rsplit('.', 1)[0] + '.mb.nex'

# Define code blocks for mrbayes input nexus alignment files, using output name.
# Note: There is probably a better way to handle the variables for MrBayes, but
# the code blocks can be manually edited after this file is produced.
mrbayescodeblocks =\
"""
BEGIN mrbayes;

   set autoclose=yes;
   log start filename = mbLogFile.log.txt;
   prset aamodelpr=fixed(LG4X);   lset rates=gamma Ngammacat=4;
   mcmc ngen=10000000 stoprule=Yes stopval=0.01 printfreq=10000 samplefreq=1000 nchains=4 savebrlens=yes filename=%s;
   log stop;
    quit;

END;



[
Begin mrbayes;
log start filename=%s.con.log replace;
sumt filename=%s relburnin=yes burninfrac=0.25 contype=allcompat;
end;
]
""" % (os.path.basename(mbnex), os.path.basename(mbnex.rsplit('.', 1)[0]), os.path.basename(mbnex))

# Use code blocks and output name to write mrbayes nexus file.
nex_to_mbnex(newnex, mbnex, mrbayescodeblocks)

 
# # Delete unnecessary files.
os.remove(tempafa1)
os.remove(tempafa2)
