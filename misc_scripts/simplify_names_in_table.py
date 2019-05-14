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

command_line_list = sys.argv
infp = str(command_line_list[1])

outfp = infp + '_simple.table'
with open(infp) as infh, open(outfp, 'w') as o:
    for i in infh:
        spliti = i.split(' ')
        if len(spliti) == 1:
            o.write(i)
        elif len(spliti) >= 1:
            o.write(spliti[0] + ' ' + spliti[1] + '\n')
        elif i.startswith('\n'):
            o.write(i)



