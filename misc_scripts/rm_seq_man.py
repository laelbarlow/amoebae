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
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import shutil
from parse_mod_num import update_mod_num_alphabetic
from misc_functions import launch, query_yes_no
from realign_nex import realign_nex

command_line_list = sys.argv
infp = str(command_line_list[1])

def rm_seq_man(infilepath):
    """Prompts manual removal of sequences from an appropriately named nexus
    file, and re-aligns it if necessary.
    """
    # Save nex with appropriate mod to modification number
    outfilepath = update_mod_num_alphabetic(infilepath)
    shutil.copyfile(infilepath, outfilepath)
    print('\nSaved as ' + outfilepath)

    # Open the output nexus file in Mesquite.
    launch('Mesquite', outfilepath)

    # Prompt manual removal from new nex file. 
    x = input('\nDone removing sequences?')

    # Determine whether or not to re-align the sequences.
    if query_yes_no('\nRe-align the nexus file?', 'yes'):
        # Re-align nex with call to realign_nex.py (adding '_realigned' should
        # not matter.
        realign_nex(outfilepath, None, '')
        print('\nRe-aligned ' + outfilepath)
       
if __name__ == '__main__':
    rm_seq_man(infp)    
