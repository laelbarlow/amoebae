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

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob

command_line_list = sys.argv
indp = str(command_line_list[1])


def delete_extra_mesquite_lines(infilepath):
    """Remove extra blocks added by Mesquite, if present.  This prevents
    Biopython from raising an error about an unmatched 'end' in one of the
    nexus blocks.  The extra blocks are not necessary anyway.
    """
    inhandle = open(infilepath, 'r')
    lines = inhandle.readlines()
    inhandle.close()
    inhandle = open(infilepath, 'w')
    x = False
    for line in lines:
        if line == 'BEGIN ASSUMPTIONS;\n':
            x = True
        if not x:
            inhandle.write(line)
    inhandle.close()

if __name__ == '__main__':
    for f in glob.glob(os.path.join(indp, '*.nex')):
        delete_extra_mesquite_lines(f)
