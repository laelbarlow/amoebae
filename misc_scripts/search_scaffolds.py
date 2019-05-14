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
"""Script to search in nuclotide sequences in fasta format using protein
queries in another fasta file.

Usage:
    search_scaffolds.py <fasta with prot queries> <fasta with nucleotide seqs>

"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import settings
import module_search_scaffolds


command_line_list = sys.argv
query_file = str(command_line_list[1])
subject_db = str(command_line_list[2])

max_gap = settings.max_gap


if __name__ == '__main__':
    module_search_scaffolds.search_scaffolds(query_file, subject_db, max_gap)
    #module_search_scaffolds.parse_tblastn(query_file, subject_db, max_gap) # but pass tblastn output and sequence file
