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
"""For converting nex files to fa.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
from afa_to_nex import nex_to_afa, delete_extra_mesquite_lines
from afa_to_fa import afa_to_fa

def nex_to_fa(infilepath, outfilepath):
    """Converts nex to fa using functions from other files.
    """
    # Name intermediate file.
    intfilepath = infilepath.replace('.nex', '.afa')

    # Convert nex input file to afa intermediate.
    delete_extra_mesquite_lines(infilepath)
    nex_to_afa(infilepath, intfilepath)

    # Convert afa to fa outfile.
    afa_to_fa(intfilepath, outfilepath)

    # remove intermediate file.
    os.remove(intfilepath)

if __name__ == '__main__':
    command_line_list = sys.argv
    infp = str(command_line_list[1])

    nex_to_fa(infp, infp.replace('.nex', '.fa'))
