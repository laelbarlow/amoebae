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
"""For converting nex in a directory to fa and afa format.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
from module_afa_to_nex import nex_to_afa, delete_extra_mesquite_lines
from afa_to_fa import afa_to_fa
from nex_to_fa import nex_to_fa

command_line_list = sys.argv
indp = str(command_line_list[1])

#print(glob.glob(os.path.join(indp, '*.nex')))

def all_nex_to_fa_and_afa(indirpath):
    """Converts all nex files in a directory to fa and afa using functions from other files.
    """
    # Loop over files in directory
    for nex in glob.glob(os.path.join(indirpath, '*.nex')):
        # Convert nex to afa
        nex_to_afa(nex, nex.replace('.nex', '.afa'))

        # Convert nex to fa
        #nex_to_fa(nex, nex.replace('.nex', '.fa'))

if __name__ == '__main__':
    all_nex_to_fa_and_afa(indp)

