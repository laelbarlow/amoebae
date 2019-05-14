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
"""Trim an nex alignment using the last taxon "MASK" sequence as a guide.

Usage:
    trim_nex.py <path to nex alignment with MASK>
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#import os
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
#from module_afa_to_nex import delete_extra_mesquite_lines
from module_amoebae_trim_nex import trim_nex

# Increase recursion limit for python. Potentially dangerous, but necessary for
# large alignments, like for Giantin and lava lamp proteins!
sys.setrecursionlimit(10000)

command_line_list = sys.argv
infp = str(command_line_list[1])


if __name__ == '__main__':
    outfilename = infp.rsplit('.', 1)[0] + '.trim.nex'
    trim_nex(infp, outfilename)
