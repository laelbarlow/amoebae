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
"""This module is for removing illegal characters from fasta files.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import shutil
import time
from afa_to_nex import afa_to_nex, nex_to_afa
from afa_to_fa import afa_to_fa
from afa_to_nex import align_one_fa
from subprocess import call
from parse_mod_num import update_mod_num_numeric

command_line_list = sys.argv
infp = str(command_line_list[1])

def rp_ill_char(infilepath):
    """Writes a new file with illegal characters replaced by X.
    """
    outfilepath = infilepath.replace('.fa', '_cleaned.fa')
    with open(infilepath) as infh, open(outfilepath, 'w') as outfh:
        for i in infh:
            if '>' not in i:
                new = i.replace('+', 'X').replace('#', 'X').replace('*', 'X')
                outfh.write(new)
            else:
                outfh.write(i)

  
if __name__ == '__main__':
    rp_ill_char(infp)    
