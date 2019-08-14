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

"""Contains functions for running qt.
"""
import os
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

def run_qt():
    """Runs python code that is dependent on qt.
    """
    # Define temporary output file path.
    temp_image_path = 'temporary_amoebae_test_image.pdf'

    # Define node styles for tree.
    print('\nX1')
    style = NodeStyle()
    style["fgcolor"] = "#0f0f0f"
    style["size"] = 0
    style["vt_line_color"] = "#ff0000"
    style["hz_line_color"] = "#ff0000"
    style["vt_line_width"] = 8
    style["hz_line_width"] = 8
    style["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
    style["hz_line_type"] = 0

    print('\nX2')
    t = Tree()
    t.populate(10, random_branches=True)
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.scale =  120 # 120 pixels per branch length unit

    # Apply node style.
    print('\nX3')
    for node in t.traverse():
        node.set_style(style)

    # Add title to tree.
    print('\nX4')
    tree_title = "[Tree title here]"
    ts.title.add_face(TextFace(tree_title, fsize=20), column=0)

    # Write tree to file.
    print('\nX5')
    t.render(temp_image_path, w=183, units="mm", tree_style=ts)

    # Remove temporary file.
    os.remove(temp_image_path)

