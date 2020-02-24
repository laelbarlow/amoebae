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
"""Mask a nex alignment by adding a new taxon, "MASK" sequence.

This is only a rough mask that is applied based on very simple criteria.

It would be good to eventually add more sophistocated criteria, especially
criteria based on similarity (according to a scoring matrix) rather than just
identity, or number of sequences without gaps.

Usage:
    mask_nex.py <path to nex alignment to MASK>
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
from module_mask_nex import mask_nex


if __name__ == '__main__':
    command_line_list = sys.argv
    infp = str(command_line_list[1])

    mask_nex(infp)
