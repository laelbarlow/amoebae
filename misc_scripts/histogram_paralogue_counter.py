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
"""Make a histogram.

This will need to be updated to work with formatting of current output files.

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
import re
import pylab
import subprocess


command_line_list = sys.argv
infilep = str(command_line_list[1])
#infasta = str(command_line_list[2])

# Get values from column as list.
column_index = 9
values = []
with open(infilep) as infh:
    for i in infh:
        if not i.startswith('ID') and not i.startswith('\n'):
            spliti = i.split(',')
            value = spliti[column_index]
            values.append(value)

# Make histogram of scores.
pylab.hist(values, bins=30) #specify the number of bins for the histogram
pylab.title("Histogram of values")
pylab.xlabel("Value")
pylab.ylabel("Number of values in bin")
#pylab.show() #can do this instead of the savefig method if just want to view
pylab.savefig(infilep.rstrip('.csv') + '_histo.output.pdf') #works for pdf or png
            



