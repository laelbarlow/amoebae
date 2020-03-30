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
sys.path.append(os.path.dirname(sys.path[0]))
from amoebae_m import get_seqs_from_fasta_db
from datetime import datetime


command_line_list = sys.argv
db_name = str(command_line_list[1])
acc_list = command_line_list[2:]

if __name__ == '__main__':
    startTime = datetime.now()
    #seq_objs = get_seqs_from_fasta_db(db_name, acc_list)
    seq_objs = get_seqs_from_fasta_db(db_name, acc_list, True)
    for s in seq_objs:
        #print(s.id)
        #print('')
        print('>' + s.description)
        #print('>' + s.id)
        print(s.seq)

    #print('Time to execute: ')
    #print(datetime.now() - startTime)

