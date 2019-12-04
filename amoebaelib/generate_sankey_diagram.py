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
"""Module containing functions for generating Sankey diagrams using matplotlib.
"""
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
 

def generate_sankey(title,
                    starting_inflow_unit_count,
                    sankey_outflow_labels,
                    sankey_outflow_proportions,
                    output_pdf_file_path
                    ):
    """Use Matplotlib to generate a Sankey diagram of to illustrate results of
    applying criteria to in similarity searching analyses, given a starting
    value (the input unit count), labels for outflows at each step, and the
    proportion of the total number of input units that flow out at each stage
    (in the same order as for the list of labels). Note that the order of the
    outflow labels and corresponding proportions input to this function can be
    re-arranged to reflect the strategy that was applied, without changing this
    function.
    """
    # Check that the outflow proportions all add to no more than 1.
    sum_of_proportions = sum(sankey_outflow_proportions)
    assert sum_of_proportions >= -1, """Specified outflow proportions for steps
    in the Sankey diagram exceed the total inflow."""

    # Define units.
    #units = '\nsequences'
    units = ''

    # Initiate matplotlib figure.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[], title=title)
    
    # Initiate sankey object.
    sankey = Sankey(ax=ax, unit=None)
    
    # Loop over outflow lists and generate sankey diagram.
    #cur_in = 1.0
    cur_in = 1.0
    count = 0
    for outflow, outflow_label in zip(sankey_outflow_proportions, sankey_outflow_labels):
        #print('\n')
        #print('outflow_label')
        #print(outflow_label)
        #print('outflow')
        #print(outflow)
        count += 1
        #print('count')
        #print(count)
        if count == 1:
            keep = -(cur_in + outflow)
            #print('keep')
            #print(keep)
            sankey.add(flows=[cur_in, outflow, keep],
                       labels=['', outflow_label + '\n' +\
                           str(round(outflow*starting_inflow_unit_count)) + units, ''],
                       orientations=[0, -1, 0]
                       )
            cur_in = -(keep)
        else:
            keep = -(cur_in + outflow)
            #print('keep')
            #print(keep)
            kl = ''
            if count == len(sankey_outflow_proportions):
                kl = 'Positive\n' + str(round((cur_in + outflow)*starting_inflow_unit_count)) + units
            sankey.add(flows=[cur_in, outflow, keep],
                       labels=['', outflow_label + '\n' +\
                           str(round(outflow*starting_inflow_unit_count)) + units, kl],
                       orientations=[0, -1, 0],
                       prior=count-2, 
                       connect=(2, 0)
                       )
            cur_in = -(keep)
    
    # Finish the plot.
    diagrams = sankey.finish()
    #plt.legend() 
    #plt.show() # This has the advantage of stretching the figure so that the outflow labels are more visible.

    # Write plot to output file path.
    fig.savefig(output_pdf_file_path, bbox_inches='tight')


if __name__ == '__main__':
    # This is just example code to generate an example sankey diagram...

    # Define title for output figure.
    title = 'Test generate_sankey function'

    # Define starting inflow unit count.
    starting_inflow_unit_count = 1000

    # Define list of outflow proportions at each stage.
    sankey_outflow_proportions = [-0.15, -0.10, -0.08, -0.05, -0.10]

    # Define list of labels for outflows.
    sankey_outflow_labels = ['Location', 'Internal stops', 'Length', 'Overlap', '%Identity']

    # Define output file path.
    output_pdf_file_path = 'test_sankey_diagram.pdf'
    assert not os.path.isfile(output_pdf_file_path)

    # Call function to generate sankey diagram.
    generate_sankey(title,
                    starting_inflow_unit_count,
                    sankey_outflow_labels,
                    sankey_outflow_proportions,
                    output_pdf_file_path
                    )

    # Open output file.
    subprocess.call(['open', output_pdf_file_path])

    # Delete output file.
    os.remove(output_pdf_file_path)

