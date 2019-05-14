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
"""This script is for reading two spreadsheets and writing a new spreadsheet
with rows that appear in the second input spreadsheet but not in the first
input spreadsheet.

Usage:
    find_unique_rows.py <input sheet 1> <input sheet 2>
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))

command_line_list = sys.argv
infilepath1 = str(command_line_list[1])
infilepath2 = command_line_list[2]

def get_csv_name(infp2):
    return infp2.replace('.csv', '_unique.csv')

def find_unique_rows(infp1, infp2):
    # Open first sheet and make a list of lists of db and acc fields for each
    # line.
    i1_list = []
    with open(infilepath1) as i1:
        for i in i1:
            i1_list.append(i)

    # Loop over lines in second sheet and only write the ones that don't
    # correspond to hits that already appear in the first sheet.
    with open(infp2) as i2, open(get_csv_name(infp2), 'w') as o:
        for i in i2:
            if not i in i1_list:
                o.write(i)

if __name__ == '__main__':
    find_unique_rows(infilepath1, infilepath2)

