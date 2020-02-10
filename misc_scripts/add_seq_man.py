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
import time
from module_afa_to_nex import afa_to_nex, nex_to_afa
from afa_to_fa import afa_to_fa
from module_afa_to_nex import align_one_fa
from subprocess import call
from parse_mod_num import update_mod_num_numeric

command_line_list = sys.argv
infp = str(command_line_list[1])

def add_seq_man(infilepath):
    """Prompts manual addition of sequences from an appropriately named nexus
    file, and re-aligns it.
    """
    # Save nex with appropriate mod to modification number
    outfilepath = update_mod_num_numeric(infilepath)
    shutil.copyfile(infilepath, outfilepath)

    # Convert new nex to afa.
    afafilepath = outfilepath.replace('.nex', '.afa')
    nex_to_afa(outfilepath, afafilepath)
    os.remove(outfilepath)

    # Convert afa to fa.
    fafilepath = afafilepath.replace('.afa', '.fa')
    afa_to_fa(afafilepath, fafilepath)
    #os.remove(afafilepath)

    # Take fasta sequences from vim input (may be better if from a file)
    EDITOR = os.environ.get('EDITOR','nano')
    tempfilepath = outfilepath.replace('.nex', '_temporary_file.txt')
    tempfilehandle = open(tempfilepath, 'w')
    initial_message = '' # optional
    tempfilehandle.write(initial_message)
    tempfilehandle.flush() # This is so initial_message is written.
    call([EDITOR, tempfilepath])
    tempfilehandle.close()

    # Append the input sequences to the .fa file
    fafilehandle = open(fafilepath, 'a')
    with open(fafilepath, 'a') as fa, open(tempfilepath) as tf:
        for line in tf:
            fa.write(line)
    os.remove(tempfilepath)

    # Check whether done adding fasta sequences.
    input('\nDone adding sequences to the .fa file (ready to align it)?')

    # Align the fa file.
    afafilepath2 = fafilepath.replace('.fa', '.afa')
    afadirpath = os.path.dirname(afafilepath2)
    align_one_fa(fafilepath, afadirpath)
    os.remove(fafilepath)

    # Convert afa to nex
    afa_to_nex(afafilepath2, outfilepath)
    os.remove(afafilepath2)

    print('\nSaved as ' + outfilepath)
  
if __name__ == '__main__':
    add_seq_man(infp)    
