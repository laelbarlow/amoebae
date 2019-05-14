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
"""For concatenating csv files in the current working directory.
Does not write lines if they start with 'Protein name'.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob


with open('0_concat_csvs_out.csv', 'w') as o:
    fnum = 0
    for f in glob.glob('*.csv'):
        fnum += 1
        if os.stat(f).st_size != 0:
            with open(f) as infile:
                inum = 0
                for i in infile:
                    inum += 1
                    if fnum == 1 and inum == 1:
                        o.write(i)
                    if not i.startswith('Protein name') and not\
                        i.startswith('\n') and not i.startswith('Forward'):
                        o.write(i) # + '\n')
    
