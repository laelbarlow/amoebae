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
"""Module for visualizing phylogenetic trees (parsing newick trees and writing
as formatted pdf files with images).
"""

import sys
import os
#sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#sys.path.append(os.path.dirname(sys.path[0]))
import subprocess
import glob
import re
# from module_afa_to_nex import nex_to_afa, afa_to_nex, delete_extra_mesquite_lines
from module_amoebae_name_replace import write_newick_tree_with_uncoded_names
from ete3 import Tree, NodeStyle, TreeStyle, TextFace
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy.misc import imread
import numpy as np
from string import Template
from PyPDF2 import PdfFileWriter, PdfFileReader
import io
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
import time
import matplotlib.colors as mc
import colorsys
from module_boots_on_mb import mbcontre_to_newick_w_probs
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from search_alignment_space import get_type_seqs_dict, get_nodes_of_interest, get_clade_name_from_model2
from module_amoebae import find_input_file_in_parent_directory

## Get directory path from input.
#command_line_list = sys.argv
#indp = str(command_line_list[1])
#
#
## Set time for timestamp.
#timestamp = time.strftime("%Y%m%d%H%M")



def plotImage(image_file):
    im = imread(image_file).astype(np.float32) / 255
    plt.imshow(im)
    a = plt.gca()
    a.get_xaxis().set_visible(False) # We don't need axis ticks
    a.get_yaxis().set_visible(False)

def get_first_leaf_from_species(tree_obj, species_start_string):
    """Take an ete3 TreeNode object and return a node object with a name that
    starts with a given string.
    """
    first_leaf = None
    for leaf in tree_obj.iter_leaves():
        if leaf.name.startswith(species_start_string):
            first_leaf = leaf
            break
    #print(species_start_string)
    #assert first_leaf is not None
    return first_leaf

def get_largest_subtree(t):
    """Take an ete3 TreeNode object and return the largest subtree.
    """
    nodes = []
    for node in t.traverse():
        nodes.append(node)
    node_with_most_leaves = sorted(nodes, key=lambda x: len(x.get_leaves()), reverse=True)[1]
    return node_with_most_leaves


def get_first_subtrees(t):
    """Take an ete3 TreeNode object and return the clades on either side of the
    root node.
    """
    # Initiate a list of the two clades on either side of the root.
    first_subtrees = []

    # Add largest subtree to list.
    #largest_subtree = get_largest_subtree(t)
    #largest_subtree_leaf_names_set = set([x.name for x in largest_subtree])
    #first_subtrees.append(largest_subtree)

    # Add the other clade to the list.
    #other_nodes = []
    #for node in t.traverse():
    #    node_leaf_names_set = set([x.name for x in node.get_leaves()])
    #    # *****This line seems to be the problem... ********
    #    if node_leaf_names_set.intersection(largest_subtree_leaf_names_set) == 0:
    #        other_nodes.append(node)
    #other_node_with_most_leaves = sorted(other_nodes, key=lambda x: len(x.get_leaves()), reverse=True)[0]
    #first_subtrees.append(other_node_with_most_leaves)

    #print('\n\nt:')
    #print(t)
    # Loop over all the nodes in the TreeNode object.
    num = 0
    for subtree in t.traverse():
        num += 1
        # Look for just the first two subtrees (either side of the root).
        if 1 < num < 4:
            #print(subtree)
            # Add the clade to the list.
            first_subtrees.append(subtree)

    # Return the list of two nodes.
    return first_subtrees


def remove_nodes_that_are_subtrees_of_others(nodes_with_paralogues):
    """Take a list of ete3 TreeNode objects and return only those that are not
    subtrees of others in the list.
    """
    nodes_with_paralogues2 = []
    for i in nodes_with_paralogues:
        nested = False
        for j in nodes_with_paralogues:
            if not i == j:
                i_leaf_names = set([x.name for x in i.get_leaves()])
                j_leaf_names = set([x.name for x in j.get_leaves()])

                # See if they are nested.
                if not len(i_leaf_names.intersection(j_leaf_names)) == 0: 
                    # Include if j is a subtree of i, but not vice versa.
                    if len(j_leaf_names) > len(i_leaf_names):
                        nested = True
        if not nested:
            nodes_with_paralogues2.append(i)

    # Make list of nodes with paralogues non-redundant.
    nodes_with_paralogues2 = list(set(nodes_with_paralogues2))

    # Return reduced list.
    return nodes_with_paralogues2


def determine_whether_first_subtrees_same_sp_overlap(i, j,
        overlapping_species):
    """Take two TreeNode objects and determine whether the first object (i) has
    a top-level subtree that has the same overlap as the two nodes do.
    """
    # Get set of species represented in the other clade.
    j_leaves = j.get_leaves()
    j_leaf_names = set([x.name for x in j_leaves])
    j_sp_names = set([x.split('__')[0] for x in j_leaf_names])

    # Call function to get top-level subtrees of the first node.
    i_first_subtrees = get_first_subtrees(i)
    
    # Get species name set for first subtree.
    i_first_subtree_sp_set1 =\
    set([x.name.split('__')[0] for x in\
        i_first_subtrees[0].get_leaves()])

    # Get species name set for second subtree.
    i_first_subtree_sp_set2 =\
    set([x.name.split('__')[0] for x in\
        i_first_subtrees[1].get_leaves()])

    # Determine whether either set overlaps just as much with the second node
    # as the entire first node does.
    first_subtrees_same_sp_overlap = False
    if i_first_subtree_sp_set1.intersection(j_sp_names) == overlapping_species:
        first_subtrees_same_sp_overlap = True
    elif i_first_subtree_sp_set2.intersection(j_sp_names) == overlapping_species:
        first_subtrees_same_sp_overlap = True

    # Return the result.
    return first_subtrees_same_sp_overlap


def get_nodes_with_paralogues(t):
    """Find nodes that contain paralogues, and return as a list.
    """
    # Initiate a final list of nodes with paralogues (can be nested).  
    established_nodes_with_paralogues = [t] # remember to slice this first element out later.
    #print('established_nodes_with_paralogues')
    #print(len(established_nodes_with_paralogues))

    # Initiate variable to store True until no new nodes with paralogues are
    # identified.
    finding_new_nodes_with_paralogues = True

    new_nodes_with_paralogues = established_nodes_with_paralogues
    num = 0
    while finding_new_nodes_with_paralogues:
        num += 1

        found_new_nodes_this_round = False

        # Initiate list of new nodes with paralogues found this round (of the
        # while loop).
        new_nodes_found_this_round = []

        # Initiate list of nodes with paralogues (to be highlighted).
        nodes_with_paralogues = []

        # Loop over new nodes.
        #print('looping over ' + str(len(new_nodes_with_paralogues)) + ' new nodes with paralogues iteration # ' + str(num))
        for t1 in new_nodes_with_paralogues:

            # Loop over all nodes to identify those with paralogues.
            for i in t1.traverse():

                # Determine whether the clade may contain paralogues of those in
                # another clade in the tree.
                # Find the set of species names.
                i_leaves = i.get_leaves()
                i_leaf_names = set([x.name for x in i_leaves])
                i_sp_names = set([x.split('__')[0] for x in i_leaf_names])
                # Only consider clades with two or more species represented.
                if len(i_sp_names) >= 2:

                    # Identify other clades that may contain paralogues of the
                    # sequences in the clade.
                    for j in t1.traverse():

                        # Get set of species represented in the other clade.
                        j_leaves = j.get_leaves()
                        j_leaf_names = set([x.name for x in j_leaves])
                        j_sp_names = set([x.split('__')[0] for x in j_leaf_names])

                        # Only consider the second clade if it is not nested with the
                        # first.
                        if len(i_leaf_names.intersection(j_leaf_names)) == 0:

                            # Determine whether the second clade has an overlapping set
                            # of species represented compared to the first clade (at
                            # least two of the same species).
                            overlapping_species = i_sp_names.intersection(j_sp_names)
                            if len(overlapping_species) >= 2:

                                # Determine whether the two subtrees on either side of
                                # the root of the first clade also have the same number
                                # of overlapping species names with the second clade.
                                first_subtrees_same_sp_overlap =\
                                determine_whether_first_subtrees_same_sp_overlap(i, j, overlapping_species)
                                if not first_subtrees_same_sp_overlap:

                                    # Add node to the list of nodes with paralogues.
                                    nodes_with_paralogues.append(i)

            # Remove nodes from the list that are subtrees of other nodes in the list.
            nodes_with_paralogues2 =\
            remove_nodes_that_are_subtrees_of_others(nodes_with_paralogues)

            # If any new nodes identified...
            if len(nodes_with_paralogues2) > 0:

                #print('Identified nodes:')
                #for i in nodes_with_paralogues2:
                #    print(i)

                # Add new nodes to list of established nodes with paralogues.
                established_nodes_with_paralogues = established_nodes_with_paralogues + nodes_with_paralogues2

                # Add nodes to list of nodes found this round.
                new_nodes_found_this_round = new_nodes_found_this_round + nodes_with_paralogues2

                # Change variable to True so that the while loop continues.
                found_new_nodes_this_round = True

        # The list of new nodes becomes the list of nodes compiled
        # here.
        new_nodes_with_paralogues = new_nodes_found_this_round
        #print(str(len(new_nodes_found_this_round)) + ' new nodes found this round')

        # Stop looking if no new nodes found.
        if not found_new_nodes_this_round:
            finding_new_nodes_with_paralogues = False

    # Return final non-redundant list of established nodes with paralogues.
    return list(set(established_nodes_with_paralogues[1:]))


def label_pdf(in_tree_file, out_tree_file, method, subs_model, timestamp):
    """Taka pdf file of a phylogenetic tree and label with a given substitution
    model and timestamp.
    """
    # Determine what to report the support as.
    if method == 'iqtree':
        node_support_type = 'Bootstrap proportion'
    elif method == 'mb':
        node_support_type = 'Prior probability'
    packet = io.BytesIO()
    # create a new PDF with Reportlab
    can = canvas.Canvas(packet, pagesize=letter)
    ## Define string and position to add.
    #can.setFont("Helvetica", 40) # Choose your font type and size (Arial is tricky).
    #can.drawString(450, 30, 'Model: ' + subs_model)
    t = can.beginText()
    #t.setCharSpace(3)

    #t.setFont('Helvetica', 12)
    #t.setTextOrigin(250, 45)

    t.setFont('Helvetica', 30)
    t.setTextOrigin(350, 145)

    t.textLines('Model: ' + subs_model + '\nTopology: IQ-tree' + '\nNode support: ' + node_support_type + '\nTimestamp: ' + timestamp)
    can.drawText(t)
    can.save()
    
    # Move to the beginning of the StringIO buffer
    packet.seek(0)
    new_pdf = PdfFileReader(packet)
    # Read your existing PDF
    existing_pdf = PdfFileReader(open(in_tree_file, "rb"))
    output = PdfFileWriter()
    # add the "watermark" (which is the new pdf) on the existing page
    page = existing_pdf.getPage(0)
    page.mergePage(new_pdf.getPage(0))
    output.addPage(page)
    # finally, write "output" to a real file
    outputStream = open(out_tree_file, "wb")
    output.write(outputStream)
    outputStream.close()


def add_borders_to_pdf_with_latex(in_tree_file, out_tree_file):
    """Make a latex file, add the existing pdf as an image, and write to pdf,
    so that there are borders.
    """

    # Define template for latex file text.
    latex_template_string = Template(\
    r"""
    \documentclass[12pt]{article}
    \usepackage{graphicx}
    \graphicspath{ {./} }
    \usepackage[margin=0.5in]{geometry}
    \begin{document}
    \pagestyle{empty}
    %\includegraphics[width=\textwidth]{$treefile}
    %\includeimage[width=\textwidth,height=\textheight,keepaspectration=true]{myimage}
    \includegraphics[width=\textwidth,height=\textheight,keepaspectratio=true]{$treefile}
    \end{document}
    """)
    
    # Define text for latex file.
    latex_file_contents =\
    latex_template_string.substitute(treefile=os.path.basename(in_tree_file))
    
    # Write latex to text file.
    latex_file_path = out_tree_file.rsplit('.', 1)[0] + '.tex'
    with open(latex_file_path, 'w') as o:
        o.write(latex_file_contents)
    assert os.path.isfile(latex_file_path)
    
    # Make pdf file from latex file.
    #os.chdir(indp)
    #os.system("pdflatex " + os.path.basename(latex_file_path))
    tempstdout = os.path.join(os.path.dirname(in_tree_file), 'tempfile.txt')
    with open(tempstdout, 'w') as o:
        # If you want to show latex output for debugging:
        #subprocess.call(["pdflatex", os.path.basename(latex_file_path)], cwd=indp)
    
        subprocess.call(["pdflatex", os.path.basename(latex_file_path)],
                cwd=os.path.dirname(latex_file_path),
                stdout=o)
    
    # Remove temporary files.
    os.remove(latex_file_path)
    os.remove(latex_file_path.rsplit('.', 1)[0] + '.log')
    os.remove(latex_file_path.rsplit('.', 1)[0] + '.aux')
    #os.remove(tempstdout)

    # Check that correct output file path exists.
    assert os.path.isfile(out_tree_file)


def get_subs_model_from_iqtree_file(file_with_subs_model_name):
    """Take the input directory path, and return the substitution model used
    (listed in the output.iqtree file).
    """
    subs_model = None
    # Open the appropriate file.
    with open(file_with_subs_model_name) as infh:
        for i in infh:
            # Find the substitution mode.
            if i.startswith('Model of substitution:'):
                subs_model = i.strip().split(': ')[1]
                break
    # Check that a substitution model was found.
    assert subs_model is not None, """Could not find substitution model."""
    # Return the substitution model.
    return subs_model


def get_subs_model_from_mb_file(file_with_subs_model_name):
    """Take the input directory path, and return the substitution model used
    (listed in the .nex file).
    """
    subs_model = None
    # Open the appropriate file.
    with open(file_with_subs_model_name) as infh:
        for i in infh:
            # Find the substitution mode.
            if 'aamodelpr' in i:
                subs_model = i.split('aamodelpr=fixed(')[1].split(')')[0]
                break
    # Check that a substitution model was found.
    assert subs_model is not None, """Could not find substitution model."""
    # Return the substitution model.
    return subs_model


def customize_node_styles_for_visualization(t):
    """Take an ete3 tree object , and modify the node styles for better
    visualization. 
    """
    # Add a face to internal nodes with branch support string.
    for node in t.traverse():
        if not node.is_leaf():
            #node.add_face(TextFace(node.name, fsize=5), column=0, position='branch-top')
            node.add_face(TextFace(node.support, fsize=5), column=0, position='branch-top')
    
    ## Change node names to branch support values for internal nodes.
    #for node in t.traverse():
    #    if not node.is_leaf():
    #        #node.support = float(node.name.split("/")[1])
    #        if node.name != '':
    #            node.support = float(node.name.split("/")[1])
    #        else:
    #            node.support = float(0.0)

    # Remove blue dots before leaf names.
    # Draws nodes as small red spheres of diameter equal to 10 pixels
    nstyle = NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 0
    # Gray dashed branch lines
    #nstyle["hz_line_type"] = 1
    #nstyle["hz_line_color"] = "#cccccc"
    #nstyle["fgcolor"] = "#0f0f0f"
    nstyle["vt_line_color"] = "#000000"
    nstyle["hz_line_color"] = "#000000"
    nstyle["vt_line_width"] = 1.5
    nstyle["hz_line_width"] = 1.5
    nstyle["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
    nstyle["hz_line_type"] = 0
    # Applies the same static style to all nodes in the tree. Note that,
    # if "nstyle" is modified, changes will affect to all nodes
    for n in t.traverse():
       n.set_style(nstyle)
       if n.is_leaf():
           # Set font size for leaf/taxon labels.
           name_face = TextFace(n.name, fgcolor="black", fsize=5)
           n.add_face(name_face, column=0, position='branch-right')
    

def customize_node_styles_for_paralogue_clades(t, nodes_with_paralogues):
    """Take an ete3 tree object and a list of nodes with paralogues, and modify
    the node styles to emphasize important features. 
    """
    # Check that there are no redundant nodes listed in the input list of
    # TreeNode objects.
    assert len(set(nodes_with_paralogues)) == len(nodes_with_paralogues), """Redundant objects in input list."""

    # Make a list of colours as an iterable.
    # Colour scheme source: https://www.nature.com/articles/nmeth.1618?WT.ec_id=NMETH-201106
    colour_iterable = iter(['#%02x%02x%02x' % (230, 159, 0),#Orange
                            '#%02x%02x%02x' % (86, 180, 233),#Sky blue
                            '#%02x%02x%02x' % (0, 158, 115),#Bluish green
                            '#%02x%02x%02x' % (240, 228, 66),#Yellow
                            '#%02x%02x%02x' % (0, 114, 178),#Blue
                            '#%02x%02x%02x' % (213, 94, 0),#Vermillion
                            '#%02x%02x%02x' % (204, 121, 167)#Reddish purple
                            ]*100)

    # Sort list of nodes with paralogues by descending leaf node count.
    nodes_with_paralogues.sort(key=lambda x: len(x.get_leaves()), reverse=True)
    
    # Highlight nodes with identified paralogues.
    num_nodes = len(nodes_with_paralogues)
    index = -1 
    for i in range(0, num_nodes):
        index += 1
        node = nodes_with_paralogues[index]

        # Get current highlight colour, if any.
        current_bgcolor = node.img_style['bgcolor']

        # Set up a node style for clades with paralogues.
        nstyle2 = NodeStyle()
        nstyle2["shape"] = "sphere"
        nstyle2["size"] = 0
        nstyle2["vt_line_color"] = "#000000"
        nstyle2["hz_line_color"] = "#000000"
        nstyle2["vt_line_width"] = 2
        nstyle2["hz_line_width"] = 2
        nstyle2["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
        nstyle2["hz_line_type"] = 0

        # Make sure that the new colour is different than the old color (so that
        # all highlighted clades are actually visible.
        same_bgcolor = True
        new_bgcolor = None
        while same_bgcolor:
            new_bgcolor = lighten_color(next(colour_iterable), 0.5)
            if new_bgcolor != current_bgcolor: 
                same_bgcolor = False
        assert new_bgcolor != current_bgcolor

        # Set new background color.
        nstyle2["bgcolor"] = new_bgcolor

        ## Temp.
        #print('\n\nNode:')
        #print(node)
        #print('current_bgcolor')
        #print(current_bgcolor)
        #print('new_bgcolor')
        #print(new_bgcolor)

        # Set the node style for node and all subnodes.
        node.set_style(nstyle2)
        for sn in node.traverse():
            sn.set_style(nstyle2)


def customize_node_styles_for_clades_to_remove(t, highlight_for_removal):
    """Take an ete3 tree object and a list of nodes with paralogues, and modify
    the node styles to emphasize important features. 
    """
    # Define highlight colour.
    # Colour scheme source: https://www.nature.com/articles/nmeth.1618?WT.ec_id=NMETH-201106
    highlight_hex_colour = '#%02x%02x%02x' % (240, 228, 66) #Yellow

    # Highlight leaf nodes identified for removal. 
    for node in t.iter_leaves():
        if node.name in highlight_for_removal: 
            # Set up a node style for leaf. Mostly similar to that defiined in
            # customize_node_styles_for_paralogue_clades and
            # customize_node_styles_for_visualization functions.
            nstyle2 = NodeStyle()
            nstyle2["shape"] = "sphere"
            nstyle2["size"] = 0
            nstyle2["vt_line_color"] = "#000000"
            nstyle2["hz_line_color"] = "#000000"
            nstyle2["vt_line_width"] = 2
            nstyle2["hz_line_width"] = 2
            nstyle2["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
            nstyle2["hz_line_type"] = 0
            #nstyle2["bgcolor"] = lighten_color(next(colour_iterable), 0.5)
            nstyle2["bgcolor"] = highlight_hex_colour
            node.set_style(nstyle2)


def lighten_color(color, amount=0.5):
    """Takes a hex colour code and returns a hex code for the same colour but
    lighter by multiplying (1-luminosity) by the given amount (fraction).
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    ct = colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
    return '#%02x%02x%02x' % (int(ct[0]*255), int(ct[1]*255), int(ct[2]*255))


def get_branch_support_from_node_name(node_name):
    """Parse node name of internal nodes resulting from parsing a newick tree
    with ete3 with the format=1 option (branch supports become node names for
    internal nodes). Return a single float corresponding to the bootstrap
    proportion.
    """
    # Get bootstrap proportion from tree generated by IQ-tree with bootstrap
    # and aLRT options specified (first value is aLRT).
    bootstrap_proportion = float(node_name.split('/')[1]) / 100
    # Return bootstrap proportion as a float.
    return float(bootstrap_proportion)


def translate_int_node_names_to_support(t):
    """Parse branch support values out of internal node names, add one of them
    to each node as branch support (a single float), and delete the internal
    node names.  """
    for node in t.traverse():
        if not node.is_leaf():
            if node.name != '':
                # Assign support attribute value to node.
                node.support = get_branch_support_from_node_name(node.name)
                # Delete existing node name attribute for node.
                node.name = ''


def translate_int_node_support_to_prob(t):
    """Parse branch support values and change them to probabilities instead of
    perentages.
    """
    for node in t.traverse():
        if not node.is_leaf():
            #if node.name != '':
            # Change node support value.
            node.support = float(format(node.support * 0.01, '3.2f'))


#def add_clade_names_to_conversion_table(tablefilename, typeseqsfilename):
#    """Take paths to a taxon name conversion file and a file defining clade
#    names corresponding to type seqs (sequences defining clades of interest).
#    And, append corresponding clade names to each taxon name in the conversion
#    table file.
#    """
#    # Parse tree file.
#    t1 = Tree(...XXX...)
#
#    # Get list of nodes of interest.
#    orthogroup_nodes = get_nodes_of_interest(t1, type_seq_list)
#
#    # Construct a dictionary storing the clade name corresponding to each
#    # sequence name in the input tree.
#    seq_clade_name_dict = {}
#    for nl in orthogroup_nodes:
#        clade_name = get_clade_name_from_model2(nl[0], type_seqs_dict)
#        for ln in [x.name for x in nl[1].get_leaves()]:
#            seq_clade_name_dict[ln] = clade_name
#
#    # Define path to original table file.
#    original_table_file =\
#    find_input_file_in_parent_directory(os.path.dirname(ml_tree_path), 'table', [])
#
#    # Make a copy of the original table file.
#    original_table_file_copy =\
#    original_table_file + '_original'
#    shutil.copyfile(original_table_file, original_table_file_copy)
#
#    # Define path to temporary modified file.
#    original_table_file_temp =\
#    original_table_file + '_TEMP'
#
#    # Loop over lines in original table.
#    with open(original_table_file) as tablefh,\
#        open(original_table_file_temp,'w') as o:
#        for i in tablefh:
#            if i.startswith('ZZ') and '_' not in i:
#                o.write(i)
#            elif i.strip() == '':
#                o.write(i)
#            else:
#                o.write(i.strip() + '__' + seq_clade_name_dict[i.strip()].strip() + '\n')
#
#    # Copy over the original table file path with the temporary file contents.
#    os.remove(original_table_file)
#    os.rename(original_table_file_temp, original_table_file)


def visualize_tree(method,
                   timestamp,
                   taxa_to_root_on,
                   highlight_paralogues,
                   highlight_for_removal,
                   tree_file,
                   tablefilename,
                   typeseqsfilename,
                   add_clade_names,
                   file_with_subs_model_name
                   ):
    """Take relevant files and write PDF files with tree figures, and return
    output file paths.
    """
    ## Check that necessary files are input if not just highlighting for
    ## removal.
    #if len(highlight_for_removal) == 0:
    #    assert file_with_subs_model_name is not None, """No file with
    #    substitution model name was input."""

    # If method is mb, then generate a .tre file.
    if method == 'mb':
        # Modify nexus file to use second code block. 
        temp_nexus = tree_file + '_temp'
        with open(tree_file, 'r') as infh, open(temp_nexus, 'w') as o:
            # Get text from file.
            text = infh.read()
            # Parse text.
            delimiters = re.compile(r'(BEGIN \w+;\n|END;\n)', re.IGNORECASE)
            split_on_delimiters = delimiters.split(text)
            inum = -1
            #for i in split_on_delimiters:
            #    inum += 1
            #    print('\n\n\n\n\n')
            #    print('<' + str(inum) + '>')
            #    print(i)
            assert len(split_on_delimiters) == 13
            line1 = split_on_delimiters[0]
            ali_block = split_on_delimiters[1] + split_on_delimiters[2] + split_on_delimiters[3]
            code_block1 = split_on_delimiters[5] + split_on_delimiters[6] + split_on_delimiters[7]
            code_block2 = split_on_delimiters[9] + split_on_delimiters[10] + split_on_delimiters[11]
            # Make new text with brackets around first block instead of second.
            new_text = line1 + ali_block + '\n[\n' + code_block1 + ']\n\n' + code_block2
            #print(new_text)
            # Write new text to temp file.
            o.write(new_text)
        # Replace old file with new file.
        os.remove(tree_file)
        os.rename(temp_nexus, tree_file)

        # Generate line graphs from .p files.
        for f in glob.glob(os.path.join(os.path.dirname(tree_file), '*.p')):
            # Parse the data.
            generations = []
            log_likelihoods = []
            with open(f) as infh:
                for i in infh:
                    if not i.startswith('\n') and not i.startswith('[ID') and not i.startswith('Gen'):
                        # Parse relevant lines add relevant values to data
                        # lists.
                        spliti = i.split('\t')
                        generation = spliti[0]
                        generations.append(int(generation))
                        log_likelihood = spliti[1]
                        log_likelihoods.append(float(log_likelihood))

            # Make data into an array.
            x = np.array(generations)
            y = np.array(log_likelihoods)

            # Use matplotlib to plot the data array.
            fig, ax = plt.subplots()
            ax.plot(x, y)
            ax.set(xlabel='Generation', ylabel='Log likelihood',
                   title='Plot of log likelihood change over ' + str(generations[-1]) + ' MCMC generations')
            ax.grid()

            # Write the plot to a pdf file.
            fig.savefig(f + '_line_graph.pdf')

        # Add end; to .t files if necessary.
        for f in glob.glob(os.path.join(os.path.dirname(tree_file), '*.t')):
            # Check for end;.
            has_end = False
            with open(f, 'r') as infh:
                text = infh.read()
                if text.strip().endswith('end;'):
                    has_end = True
            # If no end; then append an end;.
            if not has_end:
                with open(f, 'a') as infh:
                    infh.write('\nend;')

        # Run mb on the modified nexus alignment file.
        subprocess.call(['mb', os.path.basename(tree_file)],
                cwd=os.path.dirname(tree_file)) 

        # Change "tree_file" to actual tree file (.con.tre).
        tree_file =\
        glob.glob(os.path.join(os.path.dirname(tree_file),'*.con.tre'))[0]

        # Convert nexus .tre file to newick format.
        newick_file_path = tree_file.rsplit('.', 1)[0] + '.newick' 
        #contre_to_newick(tree_file, newick_file_path) # This erases supports!
        mbcontre_to_newick_w_probs(tree_file, newick_file_path)

        # Set tree file path to the newick file instead of the nexus file.
        tree_file = newick_file_path

    # Define name for decoded tree file.
    tree_file2 = os.path.join(os.path.dirname(tree_file), os.path.basename(tree_file).replace('.', '_') + '_DC.tre')
    if len(highlight_for_removal) > 0:
        tree_file2 = os.path.join(os.path.dirname(tree_file),\
                os.path.basename(tree_file).replace('.', '_') +\
                '_DC_for_removal.tre')

    # Use module to uncode names.
    write_newick_tree_with_uncoded_names(tree_file, tree_file2, tablefilename)
    
    # Parse tree from file with ete3.
    t1 = None
    if method == 'iqtree':
        t1 = Tree(tree_file2, format=1)
        # Parse branch support values out of internal node names, add one of
        # them to each node as branch support (a single float), and delete the
        # internal node names.
        translate_int_node_names_to_support(t1)

    elif method == 'mb':
        t1 = Tree(tree_file2, format=0)
        # Translate branch support values by multiplying by 0.01 (to get prior
        # probabilities instead of percentages.
        translate_int_node_support_to_prob(t1)

    # Define name for tree with no branch lengths or supports.
    tree_file2_topo = tree_file2.rsplit('.', 1)[0] + '.topology.tre'

    # Write topology to file path (this is not for use by this script).
    t1.write(format=9, outfile=tree_file2_topo)

    
    # If indicated, append clade names to taxon names in the output trees and
    # conversion table.
    if add_clade_names:
        # Get list of nodes of interest.
        type_seq_dict = get_type_seqs_dict(typeseqsfilename)
        type_seq_list = type_seq_dict.values()
        #print('\n\n')
        #print(typeseqsfilename)
        #print('\n\n')
        #print(type_seq_dict)
        #print('\n\n')
        #print(type_seq_list)
        #print('\n\n')
        orthogroup_nodes = get_nodes_of_interest(t1, type_seq_list)

        # Construct a dictionary storing the clade name corresponding to each
        # sequence name in the input tree.
        seq_clade_name_dict = {}
        for nl in orthogroup_nodes:
            clade_name = get_clade_name_from_model2(nl[0], type_seq_dict)
            for ln in [x.name for x in nl[1].get_leaves()]:
                seq_clade_name_dict[ln] = clade_name

        # Loop over leaf names in tree and append clade names.
        for l in t1.get_leaves():
            if l.name in seq_clade_name_dict.keys():
                l.name = l.name + '__' + seq_clade_name_dict[l.name]
            else:
                l.name = l.name + '__' + 'Unknown'

        # Define path to modified file with clade names appended to taxon
        # names.
        table_file_with_clade_names =\
        tablefilename + '_with_clade_names'

        # Loop over lines in original table.
        with open(tablefilename) as tablefh,\
            open(table_file_with_clade_names,'w') as o:
            for i in tablefh:
                if i.startswith('ZZ') and '_' not in i:
                    o.write(i)
                elif i.strip() == '':
                    o.write(i)
                else:
                    if i.strip() in seq_clade_name_dict.keys():
                        o.write(i.strip() + '__' + seq_clade_name_dict[i.strip()].strip() + '\n')
                    else:
                        o.write(i)



    ## Find a K. nitens sequence if present.
    #first_leaf_from_species = get_first_leaf_from_species(t1, 'Klebsormidium')
    
    ## Root tree on K. nitens sequence if possible.
    #t1.set_outgroup(first_leaf_from_species)
    
    ## If no K. nitens sequence is present, then root on the midpoint node.
    #if first_leaf_from_species is None:
    #    t1.set_outgroup(t1.get_midpoint_outgroup())
    
    # Initiate a dictionary of nodes and labels for the nodes on which to root the
    # tree.
    root_node_dict = {}
    
    # Add midpoint to dict of nodes to root on.
    root_node_dict['midpoint'] = t1.get_midpoint_outgroup()

    # Check that the midpoint outgroup node was found.
    assert root_node_dict['midpoint'] is not None, """No midpoint outgroup was
    identified."""
    
    # Add specific species to dict of nodes to root on.
    if len(taxa_to_root_on) > 0:
        print(taxa_to_root_on)
        for taxon in taxa_to_root_on:
            fl = get_first_leaf_from_species(t1, taxon)
            if fl is not None:
                root_node_dict[taxon] = get_first_leaf_from_species(t1, taxon)
    
    # Loop over nodes too root on, and generate tree files for each.
    list_of_output_pdf_file_paths = []
    list_of_rerooted_tree_objects = []
    for root_node_label in root_node_dict.keys():
        # Root on corresponding node from dictionary.
        t1.set_outgroup(root_node_dict[root_node_label])
    
        # Make a copy of the TreeNode object.
        t2 = t1.copy()
    
        # Customize the node styles generally.
        customize_node_styles_for_visualization(t2)

        if highlight_paralogues:
            # Identify nodes that contain paralogoues.
            nodes_with_paralogues = get_nodes_with_paralogues(t2)

            ## Temp.
            #print('\n\nFor tree n:')
            #for n in nodes_with_paralogues:
            #    print(n)

            # Customize the node styles of clades with paralogues.
            customize_node_styles_for_paralogue_clades(t2, nodes_with_paralogues)

        elif len(highlight_for_removal) > 0:
            # Customize the node styles of leaf nodes identified for removal.
            customize_node_styles_for_clades_to_remove(t2, highlight_for_removal)


        # Initiate a tree style.
        ts = TreeStyle()
        
        # Specify not showing leaf names, because need to define font size, so add as a "face" (see above).
        ts.show_leaf_name = False
        
        # Stretch branches.
        #ts.scale =  100 # 120 pixels per branch length unit (Throws of dimensions).
        #ts.branch_vertical_margin = 6 # 10 pixels between adjacent branches # Not useful.
        
        # Specify showing branch support values (not necessary if showing branch
        # supports (node names) as text faces for each node.
        #ts.show_branch_support = True
    
        # Define path to intermediate pdf file.
        tree_file3 = tree_file2.rsplit('.', 1)[0] + '_' + root_node_label + 'Root.pdf'
        
        # Write tree to intermediate pdf file.
        #print('Rendering PDF image of tree.')
        #t1.render(tree_file3, w=183, units="mm", tree_style=ts)
        #t1.render(tree_file3, tree_style=ts, w=8.5, h=11, units='in', dpi=600)
        #t1.render(tree_file3, tree_style=ts, dpi=600) # Results in prepended blank page, and image exceding length of page.
        t2.render(tree_file3, tree_style=ts, w=8.5, h=11, units='in', dpi=600) 
        # Check that tree was rendered as a pdf file.
        assert os.path.isfile(tree_file3)
        
        subs_model = 'N/A'
        if file_with_subs_model_name is not None:
            if method == 'iqtree':
                # Get model used from iqtree file.
                subs_model = get_subs_model_from_iqtree_file(file_with_subs_model_name)
            elif method == 'mb':
                subs_model = get_subs_model_from_mb_file(file_with_subs_model_name)
            assert subs_model is not 'N/A'
        
        # Define output filename.
        tree_file4 = tree_file3.rsplit('.', 1)[0] + '_labeled.pdf'
        
        # Add a label to the intermediate pdf file with the substitution model and
        # timestamp.
        #print('Labeling PDF image of tree.')
        label_pdf(tree_file3, tree_file4, method, subs_model, timestamp)
        
        # Define file path for pdf with borders.
        latex_file_path = tree_file4.rsplit('.', 1)[0] + '_latex.pdf'
        
        # Make a latex file (so that there are borders around the tree image).
        #print('Adding PDF image of tree to a PDF document.')
        add_borders_to_pdf_with_latex(tree_file4, latex_file_path)
    
        # Remove intermediate files.
        os.remove(tree_file3)
        os.remove(tree_file4)
    
        # Append output file path to list.
        list_of_output_pdf_file_paths.append(latex_file_path)

    # Return list of output file paths.
    return list_of_output_pdf_file_paths


def visualize_tree_in_dir(indp, method, timestamp, taxa_to_root_on,
                          highlight_paralogues, add_clade_names_from_file):
    """Given a path to a directory with output files from a phylogenetic
    analysis, write PDF files with images of the results, and return the paths
    to those files.

    indp is the path to the directory with results.
    method is either 'iqtree', 'mrbayes', or 'raxml'.
    taxa_to_root_on is a list of names such as 'Klebsormidium' or
    'Klebsormidium nitens'
    """
    # Find input tree path based on phylogenetic method used.
    tree_file = None
    file_with_subs_model_name = None
    if method == 'iqtree':
        # Names for input files.
        tree_file = None
        #if os.path.isfile(os.path.join(indp, 'output.treefile')):
        #    tree_file = os.path.join(indp, 'output.treefile')
        #else:
        #    print("""\nWarning: Could not find .treefile file, using .contree
        #    instead.\n""")
        #    tree_file = os.path.join(indp, 'output.contree')
        tree_file = os.path.join(indp, 'output.treefile')

        # Identify file listing which substitution model was used.
        file_with_subs_model_name = '?'
        try:
            file_with_subs_model_name = glob.glob(os.path.join(indp, '*.iqtree'))[0]
        except:
            pass

    elif method == 'mb':
        # Names for input files.
        tree_file = glob.glob(os.path.join(indp, '*.nex'))[0]
        # Identify file listing which substitution model was used.
        file_with_subs_model_name = glob.glob(os.path.join(indp, '*.nex'))[0]
    
    # Check that the files are present in the input directory.
    assert os.path.isfile(tree_file), 'Missing file.'
    
    ## Get tablefilename list.
    #tablefilenames = glob.glob(indp + '/*.table')

    ## Make sure one and only one table file path was found.
    #assert len(tablefilenames) > 0, 'Error: No conversion tables identified in\
    # directory ' + indp
    #
    #assert len(tablefilenames) < 2, 'Error: More than one potential conversion\
    # table files identified in directory '+ indp + ' \nFiles: ' +\
    #str(tablefilenames)
    
    # Get tablefilename.
    #tablefilename = tablefilenames[0]
    tablefilename = find_input_file_in_parent_directory(indp, 'table', [])

    # Get type seqs filename.
    typeseqsfilename = find_input_file_in_parent_directory(indp, 'type_seqs', [])

    # Do not highlight leaves for removal.
    highlight_for_removal = []

    # Call a function to visualize tree represented by the identified files.
    list_of_output_pdf_file_paths = visualize_tree(method,
                                                   timestamp,
                                                   taxa_to_root_on,
                                                   highlight_paralogues,
                                                   highlight_for_removal,
                                                   tree_file,
                                                   tablefilename,
                                                   typeseqsfilename,
                                                   add_clade_names_from_file,
                                                   file_with_subs_model_name
                                                   )

    # Return final output pdf file paths as a list.
    return list_of_output_pdf_file_paths 

