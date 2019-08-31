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
"""Removes unnecessary quotation marks and branch lengths from newick trees.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
from ete3 import Tree

command_line_list = sys.argv
indp = str(command_line_list[1])


relevant_newick_file_paths = glob.glob(os.path.join(indp, '*_cat.tre'))
for newick in relevant_newick_file_paths:
    t = Tree(newick)
    for i in t.get_leaves():
        i.name = i.name.strip('\'')
    outfp = newick.rsplit('.', 1)[0] + '_clean.tre'
    t.write(format=9, outfile=outfp)
