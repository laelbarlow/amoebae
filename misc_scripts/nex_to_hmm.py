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
"""This module is for making HMMs from alignment files in nexus format by first
converting them to aligned fasta files and then running hmmbuild.
"""


import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#import subprocess
#from afa_to_nex import nex_to_afa, delete_extra_mesquite_lines
from nex_to_hmm import nex_to_hmm 

command_line_list = sys.argv
infp1 = str(command_line_list[1])
#infp2 = str(command_line_list[1])

if __name__ == '__main__':
    nex_to_hmm(infp1)    
