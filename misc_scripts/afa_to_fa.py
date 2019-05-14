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
from Bio import SeqIO
from misc_functions import get_fa_record_text_from_obj


def afa_to_fa(infilepath, outfilepath):
    """Takes an aligned fasta file (protein) and writes a .fa file.
    """
    inhandle = open(infilepath)
    outhandle = open(outfilepath, 'w')
    afa_file = SeqIO.parse(inhandle, "fasta")
    for record in afa_file:
        x = str(record.seq)
        record.seq = x.replace('-', '')
        outhandle.write(get_fa_record_text_from_obj(record))
    inhandle.close()
    outhandle.close()

if __name__ == '__main__':
    command_line_list = sys.argv
    infp = str(command_line_list[1])
    afa_to_fa(infp, infp.replace('.afa', '.fa'))    
