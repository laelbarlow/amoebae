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
"""This module is for converting aligned fasta files to fasta files.

Bug: This script realigns the input file as well as the output file.

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import shutil
import re
from afa_to_nex import afa_to_nex, nex_to_afa, delete_extra_mesquite_lines
from afa_to_fa import afa_to_fa
from afa_to_nex import align_one_fa
from Bio import SeqIO

def realign_nex(infilepath, outdirpath=None, text_append='_realigned'):
    """Takes a nexus file that has been modified in some way since a previous
    version and re-aligns it with muscle, while performing the necessary file
    format conversions and deleting intermediate afa files.
    """
    # Set the output dir to the input dir by default.
    if outdirpath is None:
        outdirpath = os.path.dirname(infilepath)

    # Remove extra blocks added by Mesquite, if present.  This prevents
    # Biopython from raising an error about an unmatched 'end' in one of the
    # nexus blocks.  The extra blocks are not necessary anyway
    delete_extra_mesquite_lines(infilepath)

    # Convert nex to afa.
    outafapath = os.path.join(outdirpath,
            os.path.basename(infilepath).replace('.nex', '_temp.afa'))
    nex_to_afa(infilepath, outafapath)

    # Convert afa to fa, and remove afa.
    outfapath = outafapath.replace('.afa', '.fa')
    afa_to_fa(outafapath, outfapath)
    os.remove(outafapath)
    
    # Align fa (output afa).
    align_one_fa(outfapath)
    assert os.path.isfile(outfapath), 'Error: No alignment produced.'
    os.remove(outfapath)

    # Convert afa to nex.
    outnexpath = outafapath.replace('_temp.afa', text_append + '.nex') 
    afa_to_nex(outafapath, outnexpath) 
    os.remove(outafapath)

    ## This replaces the original file, so commented out.
    ## Copy nex for mod.
    #outnexpath2 = outnexpath.replace(text_append + '.nex', '.nex')
    #if outnexpath2 != outnexpath:
    #    shutil.copyfile(outnexpath, outnexpath2)
    
if __name__ == '__main__':
    command_line_list = sys.argv
    infp = str(command_line_list[1])

    realign_nex(infp)    
