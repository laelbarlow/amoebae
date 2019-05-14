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
"""Gets taxon numbers from a name conversion table, which should correspond to
their numbers in the corresponding nexus file 

Option 2 is designed to generate text for constraining MrBayes and RAxML
analyses to reduce the search space. Currently it allows subfamily clades to be
easily constrained as monophyletic, leaving the algorithm to search for the
best backbone topology and the topology within the subfamilies.

Usage:
    get_taxon_numbers_from_table.py <conversion table file path> <option number (1/2)>

    *** Modify the output from option 1 (add clade name values in first
    column), and use as input for option 2. If there are sequences that are not
    to be included in a constrained clade, then assign them to the 'Not' clade.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))

command_line_list = sys.argv

infp1 = str(command_line_list[1]) # conversion table as output by name_replace.pl. 
in2 = str(command_line_list[2]) # 

# Option 1: Get a table with useful numbers.
if in2 == '1':
    with open(infp1) as table, open(infp1 + '_taxon_numbers.csv', 'w') as o:
        table_list = []
        for i in table:
            table_list.append(i.strip())

        old_names = table_list[1::2]
        new_names = table_list[::2]
        
        # Write a new file with numbers.
        num = 0
        for i in old_names:
            num += 1
            o.write(',' + str(num) + ',' + i + ',' + new_names[num - 1] + '\n')

# Option 2: Format lists for input to MrBayes topology constraint command.
# Input a modified output file from option 1.
elif in2 == '2':
    with open(infp1) as table,\
         open(infp1.replace('.csv','_MrBayes_constr.txt'), 'w') as o,\
         open(infp1.replace('.csv','_RAxML_constr_tree.txt'), 'w') as o2:
        table_dict = {}
        table_dict2 = {}
        for i in table:
            if not i.startswith(','):
                spliti = i.split(',')
                if spliti[0] not in table_dict.keys():
                    table_dict[spliti[0]] = [spliti[1]]
                    table_dict2[spliti[0]] = [spliti[3].strip()]
                else:
                    table_dict[spliti[0]] = table_dict[spliti[0]] + [spliti[1]]
                    table_dict2[spliti[0]] = table_dict2[spliti[0]] + [spliti[3].strip()]

        # Make MrBayes constraint commands for all clades.
        constraint_strings = []
        for key in table_dict.keys():
            # Do not constrain sequences in the 'Not' clade.
            if key.strip() != 'Not':
                constraint_strings.append('constraint' + ' ' + key +\
                ' -1' + ' = ' + ' '.join(table_dict[key]))

        o.write('MrBayes constraint commands:\n\n')

        for constraint in constraint_strings:
            o.write(constraint + ';\n')

        constraint_command = 'prset topologypr=constraints(' +\
        ', '.join([x for x in table_dict.keys() if x != 'Not']) + ');'

        o.write(constraint_command)

        # Make RAxML constraint tree
        constraint_strings2 = []
        for key in table_dict2.keys():
            # Do not constrain sequences in the 'Not' clade.
            if key.strip() != 'Not':
                newconstraint = '(' + ','.join(table_dict2[key]) + ')'
                constraint_strings2.append(newconstraint)
            else:
                newconstraint = ','.join(table_dict2[key])
                constraint_strings2.append(newconstraint)

        o2.write('(' + ','.join(constraint_strings2) + ');')

        # Make a garli command to constrain all the clades.
        # o.write('\n\n\nGARLI constraint command:\n\n')

        # garli_clades = []
        # for key in table_dict.keys():
        #     garli_clades.append('(' + ','.join(table_dict[key]) + ')')

        # o.write('+(' + ','.join(garli_clades) + ');')

# # Option 3: Make a constraint tree for RAxML.
# elif in2 == '3':
#     with open(infp1) as table, open(infp1.replace('.csv', '_RAxMLconstr.txt'), 'w') as o:
#         table_dict = {}
#         for i in table:
#             spliti = i.split(',')
#             if spliti[0] not in table_dict.keys():
#                 table_dict[spliti[0]] = [spliti[1]]
#             else:
#                 table_dict[spliti[0]] = table_dict[spliti[0]] + [spliti[1]]
#         constraint_strings = []
#         for key in table_dict.keys():
#             constraint_strings.append('constraint' + ' ' + key +\
#             ' -1' + ' = ' + ' '.join(table_dict[key]))
# 
#         #o.write('\n\n\nRAxML constraint tree:\n\n')
        



            

