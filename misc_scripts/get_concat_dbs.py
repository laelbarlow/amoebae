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
"""For concatenating protein databases in databases directory.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import settings
import time

# command_line_list = sys.argv
# indp = str(command_line_list[1])
# outdb = str(command_line_list[2])

dbdirpath = settings.dbdirpath
cur_time = time.strftime("%Y_%m_%d_%H_%M_%S")

def get_dbname_from_path(path):
    return path.split('/')[-3]

# Get a list of paths for files to concatenate. 
fa_paths = glob.glob(os.path.join(dbdirpath, '*/db_Prot/*.fa'))

# Define a list of db names to remove.
to_remove = ['Bsaltans',
             'Tbrucei',
             'Cintestinalis',
             'Cparadoxa',
             'Ctobin',
             'Ehuxleyi',
             'Gtheta',
             'Monocerc',
             'Ngruberi',
             'Otauri',
             'Psojae',
             'Spombe',
             'Tpseudonana',
             'Ngruberi',
             'Nvectensis',
             'Bnatans']

# Remove certain dbs from the list.
filtered_paths = [x for x in fa_paths if not get_dbname_from_path(x) in to_remove]

# Write a concatenated file.
outfile = os.path.join(dbdirpath, 'Concat_db_' + cur_time + '.fa')
with open(outfile, 'w') as o:
    for f in filtered_paths:
        with open(f) as infile:
            for i in infile:
                o.write(i)


