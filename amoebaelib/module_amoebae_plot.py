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
"""Module for functions for plotting results of analyses performed using
AMOEBAE.
"""
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import argparse


def get_text_label_matrix(odf):
    """Takes a dataframe with information to be plotted, and returns a matrix
    of text labels for each data point.
    """
    odf1 = odf.copy()
    row_labels = list(odf1.index)
    column_labels = list(odf1)
    for i in row_labels:
        for j in column_labels:
            if not odf1.at[i, j] == '-' and not odf1.at[i, j] == 'negative':
                odf1.at[i, j] = len(odf1.at[i, j][0])
            elif odf1.at[i, j] == 'negative':
                odf1.at[i, j] = 0
    data_labels = odf1.as_matrix() 
    return data_labels
    

def get_hit_count_matrix(odf):
    """Takes a dataframe with information to be plotted, and returns a matrix
    of hit counts for each data point.
    """
    odf2 = odf.copy()
    row_labels = list(odf2.index)
    column_labels = list(odf2)
    for i in row_labels:
        for j in column_labels:
            if not odf2.at[i, j] == '-' and not odf2.at[i, j] == 'negative':
                odf2.at[i, j] = len(odf2.at[i, j][0])
            elif odf2.at[i, j] == 'negative':
                odf2.at[i, j] = 0
            elif odf2.at[i, j] == '-':
                odf2.at[i, j] = 0
    data_count = odf2.as_matrix()
    return data_count


def get_heat_matrix_from_count_matrix(data_count):
    """Translate matrix of positive hit counts to matrix of heat values for
    heatmap.
    """
    data_heat = data_count.astype(float)
    for y in range(data_count.shape[0]):
        for x in range(data_count.shape[1]):
            if data_count[y,x] >= 1:
                data_heat[y,x] = 0
            elif data_count[y,x] <= -1:
                data_heat[y,x] = 0.5
            else:
                data_heat[y,x] = 1
    return data_heat


def generate_heatmap_figure(column_labels, row_labels, data_labels, data_count,
        data_heat, outpdfpath):
    """Take relevant data, and write heatmap plot to given output path.
    """
    # Initiate plot.
    fig, ax = plt.subplots()
    
    # Set colour scheme for heatmap.
    #heatmap = ax.pcolor(data, cmap=plt.cm.Blues)
    heatmap = ax.pcolor(data_heat, cmap=plt.cm.gray) #, alpha=0.8)
    
    # Format.
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)
    fig.set_size_inches(0.5*len(column_labels), 0.30*len(row_labels))
    
    # Turn off the frame.
    ax.set_frame_on(False)
    
    # Put the major ticks at the middle of each cell.
    #ax.set_xticks(np.arange(data_heat.shape[1]) + 1, minor=False)
    ax.set_xticks(np.arange(data_heat.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(data_heat.shape[0]) + 0.5, minor=False)
    
    # put the minor ticks at the edges of each cell.
    ax.set_xticks(np.arange(data_heat.shape[1]), minor=True)
    ax.set_yticks(np.arange(data_heat.shape[0]), minor=True)
    
    # Add gridlines based on minor ticks (tick marks).
    ax.grid(which='minor', color='gray', linestyle='-', linewidth=2)
    
    # A more natural, table-like display.
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    
    #ax.set_xticklabels(column_labels, minor=False, fontsize=10, fontweight='regular', rotation='45')
    ax.set_xticklabels(column_labels, minor=False, fontsize=10,
            fontweight='regular', rotation='vertical')
    ax.set_yticklabels(row_labels, fontsize=10, fontweight='regular', minor=False)
    
    #if not args.confidence_heat_map:
    # Add text labels to heatmap (number of identified orthologues listed in input spreadsheet).
    for y in range(data_count.shape[0]):
        for x in range(data_count.shape[1]):
            #if data_count[y,x] != 1:
            #plt.text(x + 0.5, y + 0.5, str(abs(data_count[y,x])), #'1', #'%.d' % (data_count[y, x]),
            #         horizontalalignment='center',
            #         verticalalignment='center',
            #         color='black'
            #         )
            font_color = 'black'
            if data_count[y,x] >= 1:
                font_color = 'white'
            plt.text(x + 0.5, y + 0.5, str(data_labels[y,x]), #'1', #'%.d' % (data_count[y, x]),
                     horizontalalignment='center',
                     verticalalignment='center',
                     color=font_color,
                     weight = 'bold'
                     )
    
    # Turn off all the major ticks (tick marks).
    ax = plt.gca()
    
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
        
    # Turn off all the minor ticks (tick marks).
    ax = plt.gca()
    
    for t in ax.xaxis.get_minor_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_minor_ticks():
        t.tick1On = False
        t.tick2On = False
    
    
    # Show plot.
    #plt.show()
    plt.savefig(outpdfpath, bbox_inches='tight')


def generate_stats_file(column_labels, row_labels, data_count,
        out_summary_stats_path):
    """Take relevant data, and write a file with average positive hit counts
    (averaged over genomes and over protein queries).
    """
    # Start writing output file.
    with open(out_summary_stats_path, 'w') as o:
        # Get average positive hits among all genomes and all queries.

        # Sum all positive hits.
        total_sum = 0
        for i in data_count:
            for j in i:
                total_sum += j

        # Divide by number of columns times the number of rows.
        overall_average = total_sum / (len(column_labels) * len(row_labels)) 

        # Write overall average.
        o.write('Overall average number of positive hits among all species and all queries:\n')
        o.write(str(overall_average) + '\n\n')

        # Get average positive hits over all queries for each genome (row).
        o.write('Average number of positive hits among all queries for each species:\n')
        glist = []
        index = -1 
        for genome_name in row_labels:
            index += 1
            glist.append((genome_name, sum(data_count[index]) / len(column_labels)))
        glist.sort(key=lambda x: x[1], reverse=True)
        for i in glist:
            o.write(i[0] + ':' + ' '*(max([len(x) for x in row_labels]) + 1 - len(i[0])))
            o.write(str(i[1]) + '\n')
        o.write('\n')

        # Get average positive hits over all genomes for each query (query title) (column).
        o.write('Average number of positive hits among all species for each query:\n')
        qlist = []
        col = -1 
        for query_name in column_labels:
            col += 1
            qlist.append((query_name, sum(data_count.T[col]) / len(row_labels)))
        qlist.sort(key=lambda x: x[1], reverse=True)
        for i in qlist:
            o.write(i[0] + ':' + ' '*(max([len(x) for x in column_labels]) + 1 - len(i[0])))
            o.write(str(i[1]) + '\n')
        o.write('\n')


def get_sector_label(pct, iterable):
    """Get sector label from iterable list. For constructing coulson plots.
    """
    return next(iterable)


def get_complex_info_list(complex_info_file, column_labels_simple):
    """Take complex info file path and column labels (query title) list, and
    return a list of lists of complex names as the first item followed by all
    the components (query titles) that go in that complex.
    """
    complex_info_list = []
    query_titles_accounted_for = []
    if complex_info_file is not None:
        with open(complex_info_file) as infh:
            for i in infh:
                if not i.startswith('#') and not i.startswith('\n') and not i.startswith(','):
                    spliti = i.strip().replace(' ', '').split(',')
                    complex_name = spliti[0]
                    query_title = spliti[1]

                    # Check that the query title is a valid option.
                    assert query_title in column_labels_simple, """Query title
                    %s listed in input file is not represented in input
                    data.""" % query_title

                    if len([x[0] for x in complex_info_list if x[0] == complex_name]) < 1:
                        # Initiate new list for complex name.
                        complex_info_list.append([complex_name, query_title])
                        query_titles_accounted_for.append(query_title)
                    else:
                        # Add to existing list.
                        for x in complex_info_list:
                            if x[0] == complex_name:
                                x.append(query_title)
                                query_titles_accounted_for.append(query_title)
    else:
        # Assume that each query title is its own separate "complex" (i.e.,
        # none of the proteins form a relevant complex).
        for i in column_labels_simple:
            complex_info_list.append([i, i])
            query_titles_accounted_for.append(i)

    # Check that a list of complexes was compiled.
    assert len(complex_info_list) > 0, """No input complexes."""

    # Add in any query titles that may not have been included in any complexes.
    for query_title in column_labels_simple:
        if query_title not in query_titles_accounted_for:
            complex_info_list.append([query_title, query_title])

    # Check that there are no redundant complex sublists.
    all_complex_titles = [x[0].lstrip().rstrip() for x in complex_info_list]
    assert len(all_complex_titles) == len(list(set(all_complex_titles)))

    # Return complex info list.
    return complex_info_list


def modify_legend_lines(wedges):
    """Take a list of wedge objects, and modify line properties (for pie charts
    used as legends.
    """
    for w in wedges:
        w.set_linewidth(1)
        w.set_edgecolor('black')


def modify_lines(wedges):
    """Take a list of wedge objects, and modify line properties (for pie charts
    with data).
    """
    for w in wedges:
        w.set_linewidth(2)
        w.set_edgecolor('black')


def make_coulson_plot(column_labels_simple,
                      row_labels_simple,
                      data_labels_simple,
                      data_count_simple,
                      data_heat_simple,
                      outpdfpath_coulson,
                      complex_info_file):
    """Make a coulson plot.
    """
    # Re-organize species order based on input file list order.
    #????

    # Define basic input parameters.
    num_species = len(row_labels_simple)

    # Add legend row labels
    row_labels_simple = ['Legend'] + row_labels_simple

    # Extract info from complex info file, and assemble into a list of lists.
    complex_info_list = get_complex_info_list(complex_info_file,
                                              column_labels_simple)

    # Initialize figure.
    num_columns = len(complex_info_list)
    num_rows = len(row_labels_simple)
    fig, axs = plt.subplots(num_rows, num_columns, subplot_kw=dict(aspect="equal"))

    # Scale output figure size based on number of subplots.
    fig.set_size_inches((len(complex_info_list))*2,\
                        (len(row_labels_simple))*2)

    # Decide whether to compare more than one set of paralogue counts.
    compare_counts = False
    if len(data_count_simple) > 1:
        compare_counts = True

    # Iterate over subplots, encoding information:
    for i, species in enumerate(row_labels_simple):
        jnum = -1

        for j, compl in enumerate(complex_info_list):
            jnum += 1

            # Determine number of wedges ("fracs".
            num_fracs = len(compl[1:])
            # Define list of numbers to be used to make wedges/fracs.
            fracs = []
            for x in range(0, num_fracs):
                fracs.append(1)

            # Make legends.
            if species == 'Legend':
                # Define colors.
                cols = []
                for x in range(0, num_fracs):
                    cols.append('gray')

                legend_labels = compl[1:]
                legend_radius = 0.4
                legend_title = compl[0]
                legend_title_size = 15

                # Initiate subplot pie chart.
                if num_columns < 2:
                    wedges, x = axs[i].pie(fracs, labels=legend_labels,\
                            shadow=False, colors=cols, radius=legend_radius)

                    # Add a title for the complex.
                    axs[i].set_title(legend_title, fontsize=legend_title_size)

                    # Modify lines around sectors.
                    modify_legend_lines(wedges)

                else:
                    wedges, x = axs[i,j].pie(fracs, labels=legend_labels,\
                            shadow=False, colors=cols, radius=legend_radius)

                    # Add a title for the complex.
                    axs[i,j].set_title(legend_title, fontsize=legend_title_size)

                    # Modify lines around sectors.
                    modify_legend_lines(wedges)

            # Make subplots for data.
            else:
                if not compare_counts:
                    # Define paralogue counts for each complex component.
                    paralogue_counts_list = []
                    # Determine portion that is relevant to current complex.
                    #for a, x in enumerate(column_labels_simple):
                    #    if x in compl[1:]:
                    #        paralogue_counts_list.append(data_count_simple[i-1][a])
                    for x in compl[1:]:
                        for a, y in enumerate(column_labels_simple):
                            if y == x:
                                paralogue_counts_list.append(data_count_simple[0][i-1][a])
                                break

                    # Make an iterable out of the list.
                    paralogue_counts = iter(paralogue_counts_list)


                    # Define colors for sectors based on paralogue counts.
                    cols = []
                    for x in paralogue_counts_list:
                        if x > 0:
                            # Gray if at least one paralogue found.
                            cols.append('xkcd:green')
                        else:
                            # White if no paralogues found.
                            cols.append('white')

                    # Define which function to use to add appropriate paralogue counts as
                    # labels for each sector.
                    autopct_funct = lambda pct: get_sector_label(pct, paralogue_counts)

                    # Set variables for species/row titles.
                    species_font_size = 15
                    offset = np.array([-1.25, 0.40])

                    sector_text_size = 15

                    # Make data subplots.
                    if num_columns < 2:
                         # Initiate subplot pie chart.
                        wedges, texts, autotexts = axs[i].pie(fracs,\
                                autopct=autopct_funct, textprops=dict(color="w"), shadow=False, colors=cols)

                        # Modify sector label text properties.
                        plt.setp(autotexts, size=sector_text_size, weight="bold")

                        # Modify lines around sectors.
                        modify_lines(wedges)
                        
                        if jnum == 0:
                            if len(complex_info_list) == 1:
                                # Add a title for the complex.
                                title = axs[i].set_title(species,\
                                        fontsize=species_font_size)
                            else:
                                # Add a title for the complex.
                                title = axs[i,j].set_title(species,\
                                        fontsize=species_font_size)
                            #offset = np.array([-0.75, -0.6])
                            #title.set_position(title.get_position() + offset)
                            title.set_position(np.array([0.0, 0.0]) + offset)

                    else:
                        # Initiate subplot pie chart.
                        wedges, texts, autotexts = axs[i,j].pie(fracs,\
                                autopct=autopct_funct, textprops=dict(color="w"), shadow=False, colors=cols)

                        # Modify sector label text properties.
                        plt.setp(autotexts, size=sector_text_size, weight="bold")

                        # Modify lines around sectors.
                        modify_lines(wedges)

                        if jnum == 0:
                            # Define string to make species name italicized.
                            print(species)
                            italic_species = '$\it{' + species.split(' ')[0] + '}$' + ' ' + '$\it{' + species.split(' ')[1] + '}$'
                            # Add a title for the complex.
                            title = axs[i,j].set_title(italic_species,
                                    fontsize=species_font_size)
                            #offset = np.array([-0.75, -0.6])
                            #title.set_position(title.get_position() + offset)
                            title.set_position(np.array([0.0, 0.0]) + offset)

                # Make the plot differently if comparing two sets of paralogue
                # counts.
                elif compare_counts:
                    # Define paralogue counts for each complex component.
                    paralogue_counts_list = []
                    # Determine portion that is relevant to current complex.
                    #for a, x in enumerate(column_labels_simple):
                    #    if x in compl[1:]:
                    #        paralogue_counts_list.append(data_count_simple[i-1][a])
                    for x in compl[1:]:
                        for a, y in enumerate(column_labels_simple):
                            if y == x:
                                first_count = data_count_simple[0][i-1][a]
                                second_count = data_count_simple[1][i-1][a]
                                paralogue_counts_list.append([first_count, second_count])
                                break

                    # Define colors for sectors based on paralogue counts.
                    cols = []
                    for x in paralogue_counts_list:
                        y = max(x)
                        if y > 0:
                            # Gray if at least one paralogue found.
                            if x[0] != x[1]:
                                if x[0] == x[1] -1:
                                    cols.append('xkcd:dark yellow')
                                elif x[0] == x[1] -2:
                                    cols.append('xkcd:dark orange')
                                elif x[0] > x[1]:
                                    cols.append('xkcd:cornflower')
                                else:
                                    cols.append('xkcd:dark red')
                            else:
                                cols.append('xkcd:green')

                        else:
                            # White if no paralogues found.
                            cols.append('white')

                    mod_paralogue_counts_list = []
                    for x in paralogue_counts_list:
                        if x[0] == x[1]:
                            mod_paralogue_counts_list.append(str(x[0]))
                        else:
                            mod_paralogue_counts_list.append(str(x[0]) + '-' + str(x[1]))
                    assert len(mod_paralogue_counts_list) == len(paralogue_counts_list)

                    # Make an iterable out of the modified list.
                    paralogue_counts = iter(mod_paralogue_counts_list)

                    # Define which function to use to add appropriate paralogue counts as
                    # labels for each sector.
                    autopct_funct = lambda pct: get_sector_label(pct, paralogue_counts)

                    # Set variables for species/row titles.
                    species_font_size = 15
                    offset = np.array([-1.25, 0.40])

                    sector_text_size = 15

                    # Make data subplots.
                    if num_columns < 2:
                         # Initiate subplot pie chart.
                        wedges, texts, autotexts = axs[i].pie(fracs,\
                                autopct=autopct_funct, textprops=dict(color="w"), shadow=False, colors=cols)

                        # Modify sector label text properties.
                        plt.setp(autotexts, size=sector_text_size, weight="bold")

                        # Modify lines around sectors.
                        modify_lines(wedges)
                        
                        if jnum == 0:
                            # Add a title for the complex.
                            title = axs[i,j].set_title(species, fontsize=species_font_size)
                            #offset = np.array([-0.75, -0.6])
                            #title.set_position(title.get_position() + offset)
                            title.set_position(np.array([0.0, 0.0]) + offset)

                    else:
                        # Initiate subplot pie chart.
                        wedges, texts, autotexts = axs[i,j].pie(fracs,\
                                autopct=autopct_funct, textprops=dict(color="w"), shadow=False, colors=cols)

                        # Modify sector label text properties.
                        plt.setp(autotexts, size=sector_text_size, weight="bold")

                        # Modify lines around sectors.
                        modify_lines(wedges)

                        if jnum == 0:
                            # Define string to make species name italicized.
                            italic_species = '$\it{' + species.split(' ')[0] + '}$' + ' ' + '$\it{' + species.split(' ')[1] + '}$'
                            # Add a title for the complex.
                            title = axs[i,j].set_title(italic_species,
                                    fontsize=species_font_size)
                            #offset = np.array([-0.75, -0.6])
                            #title.set_position(title.get_position() + offset)
                            title.set_position(np.array([0.0, 0.0]) + offset)

    # Specify text output type so that text can be edited in adobe illustrator.
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # Adjust gap sizes between subplots.
    fig.subplots_adjust(wspace=0, hspace=0)

    # Save figure to pdf file. (transparent=True gets rid of white
    # backgrounds).
    fig.savefig(outpdfpath_coulson, bbox_inches='tight', transparent=True)


def plot_amoebae_res(csv_file, complex_info, outpdfpath, csv_file2=None,
        row_order_file=None):
    """Main function for parsing csv file, translating data into matrices of
    relevant information, and generating plots.
    """
    # Define an output file path if none is input.
    if outpdfpath is None:
        outpdfpath = csv_file + '_plot.pdf'
    outpdfpath_simple = outpdfpath.rsplit('.', 1)[0] + '_simple.pdf'  
    out_summary_stats_path = outpdfpath.rsplit('.', 1)[0] + '_stats.txt' 

    final_output_path = None

    # Define list of csv file paths for input.
    csv_paths = []
    if csv_file2 is None:
        csv_paths = [csv_file]
    else:
        csv_paths = [csv_file, csv_file2]

    # Extract information from each input csv file.
    #outpdfpath_coulson = outpdfpath_simple.rsplit('_', 1)[0] + '_coulson.pdf'
    #make_coulson_plot(column_labels_simple, row_labels_simple, data_labels_simple, data_count_simple,
    #    data_heat_simple, outpdfpath_coulson, complex_info)

    # Define a dictionary for identifying the right columns.
    column_header_dict = {'query title': 'Query title',
                          'query file': 'Query file',
                          'alignment name': 'Alignment for sequence comparison',
                          'taxon name': 'Subject database species (if applicable)',
                          'database filename': 'Subject database file',
                          'positive or not': 'Collective interpretation of reverse search results',
                          'accession': 'Forward hit accession',
                          'program': 'Forward search method',
                          'seq': 'Forward hit sequence',
                          'subseq': 'Forward hit subsequence(s) that align(s) to query',
                          'seq descr': 'Forward hit description',
                          'subseq descr': 'Forward hit description of subsequence(s) that align(s) to query',
                          'evalue': 'Forward hit E-value (top HSP)',
                          'score': 'Forward hit score'
                          }

    # Initiate variables for comparing more than one csv file.
    column_labels_simple_list = []
    row_labels_simple_list = []
    data_labels_simple_list = []
    data_count_simple_list = []
    data_heat_simple_list = []


    file_num = 0
    for csv_file_x in csv_paths:
        file_num += 1

        ######
        # Step 1: Extract relevant information from input csv file and load into
        # a pandas dataframe. 

        # Get data from input file.
        indf = pd.read_csv(csv_file_x)

        # construct an additional dataframe with each datapoint being a list of
        # relevant information.

        # Define a list of index values for new dataframe.
        new_row_index_list = []
        for index, row in indf.iterrows():
            taxon_name = row[column_header_dict['taxon name']]
            db_file = row[column_header_dict['database filename']]
            index_tuple = (taxon_name, db_file)
            if index_tuple not in new_row_index_list:
                new_row_index_list.append(index_tuple)


        # If a database filename list file is provided, then re-order the
        # new_row_index_list list to reflect that order.
        if row_order_file is not None:
            # Get taxon list from input file.
            db_file_list_from_file = []
            with open(row_order_file) as infh:
                for i in infh:
                    if not i.startswith('\n') and not i.startswith('#'):
                        db_file_list_from_file.append(i.strip())
            # Check that the list is the same length as the list from the
            # dataframe.
            #assert len(new_row_index_list) == len(taxon_list_from_file)

            # Make a new list based on order in input file.
            new_db_file_list = []
            for i in db_file_list_from_file:
                for j in new_row_index_list:
                    if j[1] == i:
                        new_db_file_list.append(j)
            # Add remaining db filenames, that were not represented in the
            # input file list.
            for j in new_row_index_list:
                if j[1] not in db_file_list_from_file:
                    new_db_file_list.append(j)
            # Replace list with reordered list.
            new_row_index_list = new_db_file_list

        else:
            # Sort new row indexes.
            new_row_index_list.sort()

        # Define a simplified list of index values for another new dataframe.
        seen = set()
        seen_add = seen.add
        new_row_index_list_simple = [x[0] for x in new_row_index_list if not (x[0] in seen or seen_add(x[0]))]

        # ***Determine column header for column that contains final paralogue count
        # information.
        unique_positive_hit_decis_header = None
        decision_type = None
        for header in indf.columns:
            if header.startswith('Represents an identifiably unique paralogue'):
                unique_positive_hit_decis_header = header
                decision_type = 'nonredundant_positive_hit'
                break
        for header in indf.columns:
            # ***This may need to change if sum_phylo_class is modified.
            if header.startswith('ELW for most likely topology above threshold? '):
                unique_positive_hit_decis_header = header
                decision_type = 'topology_test_result'
                break
        if unique_positive_hit_decis_header is None:
            print('\nWarning: reporting all hits that meet the reverse search criteria.')
            for header in indf.columns:
                if header.startswith('Collective interpretation of reverse search results'):
                    unique_positive_hit_decis_header = header
                    decision_type = 'positive_hit'
                    break
        # Check that the header was found.
        assert unique_positive_hit_decis_header is not None, """Could not identify
        which column contains information regarding unique paralogue counts."""
        assert decision_type is not None


        # Define a list of column header values for new dataframe.
        new_column_header_list = []
        classification_header = None
        if decision_type == 'topology_test_result':
            # Get column with classification (clade names).
            for header in indf.columns:
                if header.startswith('Classification '):
                    classification_header = header
                    break
            assert classification_header is not None

            for index, row in indf.iterrows():
                # Use clade names instead of query titles.
                query_title = row[classification_header]
                query_file = row[column_header_dict['query file']]
                header_tuple = (query_title, query_file)
                if query_title != '-':
                    if header_tuple not in new_column_header_list:
                        new_column_header_list.append(header_tuple)

        # Define a list of column header values for new dataframe.
        else:
            for index, row in indf.iterrows():
                query_title = row[column_header_dict['query title']]
                query_file = row[column_header_dict['query file']]
                header_tuple = (query_title, query_file)
                if header_tuple not in new_column_header_list:
                    new_column_header_list.append(header_tuple)

        # Sort new column headers.
        new_column_header_list.sort()

        # Define a simplified list of header values for another new dataframe.
        seen = set()
        seen_add = seen.add
        new_column_header_list_simple = None
        #if decision_type == 'topology_test_result':
        #    # Make the column headers clade names from the input spreadsheet.
        #    new_column_header_list_simple = clade_names_new_column_header_list
        #else:
        #    new_column_header_list_simple = [x[0] for x in new_column_header_list if not (x[0] in seen or seen_add(x[0]))]
        new_column_header_list_simple = [x[0] for x in new_column_header_list if not (x[0] in seen or seen_add(x[0]))]


        # Initiate new dataframes.
        odf = pd.DataFrame(columns=new_column_header_list, index=new_row_index_list)
        odf_simple = pd.DataFrame(columns=new_column_header_list_simple,
                index=new_row_index_list_simple)

        # Set all values in dataframes to '-' by default.
        for i in new_row_index_list:
            for j in new_column_header_list:
                odf.at[i, j] = '-'
        for i in new_row_index_list_simple:
            for j in new_column_header_list_simple:
                odf_simple.at[i, j] = '-'

        # Temp.
        print(odf)
        print(odf_simple)

        # Iterate over rows in input dataframe and add relevant information to
        # additional dataframes.
        for index, row in indf.iterrows():
            # Determine where in the new dataframes the info should be put.
            index_tuple = (row[column_header_dict['taxon name']], row[column_header_dict['database filename']])
            if decision_type is 'topology_test_result':
                # Skip rows that don't have relevant information.
                if row[classification_header] == '-':
                    continue
                header_tuple = (row[classification_header], row[column_header_dict['query file']])
            else:
                header_tuple = (row[column_header_dict['query title']], row[column_header_dict['query file']])

            #if index_tuple[0] == 'Capsella grandiflora':
            #    print(str(index_tuple) + ' ' + str(header_tuple))

            index_simple = index_tuple[0]
            header_simple = header_tuple[0]

            # Extract relevant information from row.
            #unique_positive_hit = row[unique_positive_hit_decis_header]
            #if unique_positive_hit == 'Yes':
            # ***How to display accurate paralogue counts? Perhaps only if query title is
            # the only row index?
            positive_or_not = row[column_header_dict['positive or not']]
            acc = row[column_header_dict['accession']]
            fwd_evalue = row[column_header_dict['evalue']]
            fwd_score = row[column_header_dict['score']]
            unique_positive_hit_decis = row[unique_positive_hit_decis_header]

            # Add relevant information to dataframe for first plot.
            if positive_or_not == '-':
                # Do not add info if the result is negative anyway.
                # This way if a search was performed it will at least say
                # 'negative', and if not search was performed it will stay as '-'.
                # But, only assign 'negative' if there wasn't already a positive
                # hit recorded.
                if odf.at[index_tuple, header_tuple] == '-':
                    odf.at[index_tuple, header_tuple] = 'negative'
            else:
                if odf.at[index_tuple, header_tuple] == '-' or odf.at[index_tuple,\
                        header_tuple] == 'negative':
                    # First hit to add info for.
                    # Format: [[positive hit acc list], top fwd hit E-value, top fwd hit
                    # score] where the top hit is assumed to be the one that
                    # appears first in the input csv file. 
                    initial_value = [[acc], fwd_evalue, fwd_score] 
                    odf.at[index_tuple, header_tuple] = initial_value
                else:
                    # There are one or more previous hits relevant for this cell in the
                    # dataframe, so update existing value (list).
                    odf.at[index_tuple, header_tuple] = [odf.at[index_tuple, header_tuple][0] +\
                            [acc], odf.at[index_tuple, header_tuple][1], odf.at[index_tuple, header_tuple][2]]

            # Add relevant information to dataframe for first plot.
            if positive_or_not == '-' or unique_positive_hit_decis == 'No':
                if odf_simple.at[index_simple, header_simple] == '-':
                    odf_simple.at[index_simple, header_simple] = 'negative'
            elif unique_positive_hit_decis == 'Yes' or unique_positive_hit_decis == '+':
                if odf_simple.at[index_simple, header_simple] == '-' or odf_simple.at[index_simple, header_simple] == 'negative':
                    # First hit to add info for.
                    # Format: [[positive hit acc list], top fwd hit E-value, top fwd hit
                    # score] where the top hit is assumed to be the one that
                    # appears first in the input csv file. 
                    initial_value = [[acc], fwd_evalue, fwd_score] 
                    odf_simple.at[index_simple, header_simple] = initial_value
                else:
                    # There are one or more previous hits relevant for this cell in the
                    # dataframe, so update existing value (list).
                    odf_simple.at[index_simple, header_simple] = [odf_simple.at[index_simple, header_simple][0] +\
                            [acc], odf_simple.at[index_simple, header_simple][1], odf_simple.at[index_simple, header_simple][2]]

        
        ######
        # Step 2: Translate the dataframe into one or more matrices that can be
        # plotted in various ways. 

        ## Get matrix of positive hit counts (in input spreadsheet).
        #data_count = d.as_matrix()

        # Get matrix for text labels for output plot.
        data_labels = get_text_label_matrix(odf)
        data_labels_simple = get_text_label_matrix(odf_simple)


        # Get matrix of positive hit counts (in input spreadsheet).
        data_count = get_hit_count_matrix(odf)
        data_count_simple = get_hit_count_matrix(odf_simple)

        
        # Get matrix of heat values translated from matrix of counts.
        #if not args.confidence_heat_map:
        data_heat = get_heat_matrix_from_count_matrix(data_count)
        data_heat_simple = get_heat_matrix_from_count_matrix(data_count_simple)

        #elif args.confidence_heat_map:
        #    pass
        #    for y in range(data_count.shape[0]):
        #        for x in range(data_count.shape[1]):
        #            if data_count[y,x] == 0:
        #                data_heat[y,x] = 1
        #            elif 100 > data_count[y,x] > 0:
        #                data_heat[y,x] = (data_heat[y,x] / 100)/1.5
        #            elif data_count[y,x] == 100:
        #                data_heat[y,x] = 0
        
        #print(data_heat)

        # Make rectangular (heatmaps) of only the first csv input file.
        #if file_num == 1:
        #######
        # Step 3: Plot results.
                    
        # Plot.
        column_labels = list(odf)
        row_labels = list(odf.index)
        generate_heatmap_figure(column_labels, row_labels, data_labels, data_count,
            data_heat, outpdfpath)

        # Simple plot.
        column_labels_simple = list(odf_simple)
        row_labels_simple = list(odf_simple.index)
        generate_heatmap_figure(column_labels_simple, row_labels_simple, data_labels_simple, data_count_simple,
            data_heat_simple, outpdfpath_simple)


        if file_num == 1:
            # Generate text file with stats.
            generate_stats_file(column_labels_simple, row_labels_simple, data_count_simple,
                out_summary_stats_path)

            # Only make a coulson plot at this point if there was only one
            # input file.
            if len(csv_paths) == 1:
                # Coulson plot.
                outpdfpath_coulson = outpdfpath_simple.rsplit('_', 1)[0] + '_coulson.pdf'
                make_coulson_plot(column_labels_simple, row_labels_simple,
                        data_labels_simple, [data_count_simple],
                    data_heat_simple, outpdfpath_coulson, complex_info)
                final_output_path = outpdfpath_coulson

        # Define data to use for comparing more than one csv file.
        column_labels_simple_list.append(column_labels_simple)
        row_labels_simple_list.append(row_labels_simple)
        data_labels_simple_list.append(data_labels_simple)
        data_count_simple_list.append(data_count_simple)
        data_heat_simple_list.append(data_heat_simple)


    # If comparing two csv files, then call the make_coulson_plot function
    # with different input.
    if len(csv_paths) > 1:
        # Check that the two input coulson plots are comparable.
        assert len(csv_paths) == 2, """Can only handle up to two input csv files."""
        assert column_labels_simple_list[0] == column_labels_simple_list[1]
        
        species_set_x =  set(row_labels_simple_list[0])
        species_set_y =  set(row_labels_simple_list[1])
        for missing_species in species_set_y - species_set_x:
            print('Species missing from first spreadsheet: %s' % missing_species)
        for missing_species in species_set_x - species_set_y:
            print('Species missing from second spreadsheet: %s' % missing_species)

        assert row_labels_simple_list[0] == row_labels_simple_list[1]
        assert len(data_count_simple_list[0]) == len(data_count_simple_list[1])

        # Call function to make coulson plot.
        outpdfpath_coulson = outpdfpath_simple.rsplit('_', 1)[0] +\
        '_coulson(comparison).pdf'
        make_coulson_plot(column_labels_simple_list[0],
                          row_labels_simple_list[0],
                          data_labels_simple_list[0],
                          data_count_simple_list,
                          data_heat_simple_list[0],
                          outpdfpath_coulson,
                          complex_info)
        final_output_path = outpdfpath_coulson


    # Return main output path.
    return final_output_path



