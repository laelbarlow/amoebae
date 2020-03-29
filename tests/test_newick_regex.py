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
"""Test regular expressions on newick trees to see if they can identify
elements.

Probably obsolete.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
#from module_afa_to_nex import delete_extra_mesquite_lines
#import numpy as np
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import collections
import re
from ete3 import Tree


command_line_list = sys.argv
infp = str(command_line_list[1])

re1 = re.compile(r'[,\(][\w-]+')
re2 = re.compile(r'\((\w| )+,')
re3 = re.compile(r'\(\w+-?\w+')
re4 = re.compile(r',\w+-?\w+|\(\w+-?\w+')
re5 = re.compile(r'((,)|(\())\w+-?\w+')


def check_regex(infp):
    """Open a newick tree file and test regular expressions on the newick tree
    string.
    """
    newick_string = open(infp, 'r').read()

    x = re1.findall(newick_string)
    num = 0
    for i in x:
        num += 1
        print(num)
        print(i)

def count_leaves_in_tree(infp):
    """Take a newick tree, and print the number of leaf/terminal nodes that are
    present in it.
    """
    t1 = Tree(infp, quoted_node_names=False)

    node_count = 0
    for node in t1.traverse():
        if node.is_leaf():
            node_count += 1
            print(node.name.strip('\''))

    print(node_count)


if __name__ == '__main__':
    #check_regex(infp)
    count_leaves_in_tree(infp)
