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
"""Takes fasta files in one directory and copies them to subdirectories in
another directory, and aligns them with muscle, and converts the alignments to
nexus format.  This is to speed up the process of using positive hit sequences
from BLAST to run HMMer.  The alignments must be viewed and potentially
modified to remove false-positives from previous steps (in nexus format with
Mesquite).
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import shutil
from afa_to_nex import align_one_fa
from afa_to_nex import afa_to_nex


command_line_list = sys.argv
indp = str(command_line_list[1])
outdp = str(command_line_list[2])

def make_dir_aln(indirpath, outdirpath, name_end=''):
    fa_files = glob.glob(indirpath + '/*' + name_end + '.fa')
    if len(fa_files) == 0:
        print('No input files identified.')
    for f in fa_files:
        basename = os.path.basename(f)
        dir_name = basename.replace('.fa', '')
        dir_path = os.path.join(outdirpath, dir_name)
        os.mkdir(dir_path)
        new_fa_path = os.path.join(dir_path, basename)
        shutil.copyfile(f, new_fa_path)
        align_one_fa(new_fa_path, dir_path)
        new_afa_path = new_fa_path.replace('.fa', '.afa')
        # The extra '.' makes it so that Finder sorts the files in a more
        # logical way.
        new_nex_path = new_afa_path.replace('.afa', '_mod0..nex')
        afa_to_nex(new_afa_path, new_nex_path)
        os.remove(new_fa_path)
        os.remove(new_afa_path)

if __name__ == '__main__':
    if indp == None or outdp == None:
        print('Please specify input and output files in the commandline.') 
    else:
        x = indp
        y = outdp
        make_dir_aln(x, y, '')
 
