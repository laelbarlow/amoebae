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
"""Takes an newick tree file, and renames the taxa with codes from a given
table.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from name_replace import write_newick_tree_with_coded_names

command_line_list = sys.argv
infile = str(command_line_list[1])
outfile = str(command_line_list[2])
tablefile = str(command_line_list[3])

write_newick_tree_with_coded_names(infile, outfile, tablefile)



