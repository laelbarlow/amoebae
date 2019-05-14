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
"""Removes lines that start with commas from all csv files in a given
directory.  These commas are added to blank lines by microsoft excel when
editing csv files, but mess with my scripts that attempt to parse lines in
these files.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import re

command_line_list = sys.argv
indp1 = str(command_line_list[1])
outdp1 = indp1.rstrip('/') + '_noextracommas'
os.mkdir(outdp1)

infiles = glob.glob(os.path.join(indp1.rstrip('/'), '*.csv'))

ec = re.compile(r',+\n')

for f in infiles:
    print(f)
    with open(f) as infh, open(os.path.join(outdp1, os.path.basename(f)), 'w') as outfh:
        for i in infh:
            if not i.startswith(',,'):
                outfh.write(ec.sub('\n', i))
