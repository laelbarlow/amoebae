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
"""Contains functions for generating histograms.
"""
import sys
import os
#sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#from Bio import SeqIO
import re
import pylab
import subprocess
import numpy as np


#command_line_list = sys.argv
#infilep = str(command_line_list[1])
#infasta = str(command_line_list[2])

## Get values from column as list.
#column_index = 9
#values = []
#with open(infilep) as infh:
#    for i in infh:
#        if not i.startswith('ID') and not i.startswith('\n'):
#            spliti = i.split(',')
#            value = spliti[column_index]
#            values.append(value)

def generate_histogram(title,
                       values,
                       num_bins,
                       output_filename
                       ):
    """Take a list of values, a number of bins, and an output file name, and
    generate a histogram using pylab and write it to the file path.
    """
    # Make histogram of scores.
    pylab.hist(values, bins=num_bins) #specify the number of bins for the histogram
    pylab.title(title)
    pylab.xlabel("Value")
    pylab.ylabel("Number of values in bin")
    #pylab.show() #can do this instead of the savefig method if just want to view
    pylab.savefig(output_filename) #works for pdf or png
    pylab.close()


def generate_double_histogram(title,
                              values1,
                              label1,
                              values2,
                              label2,
                              num_bins,
                              output_filename
                              ):
    """Take two lists of values, a number of bins, and an output file name, and
    generate a histogram using pylab and write it to the file path.
    """
    pylab.style.use('seaborn-deep')
    # Make histogram of scores.
    #pylab.hist(values1, bins=num_bins, label=label1) #specify the number of bins for the histogram
    #pylab.hist(values2, bins=num_bins, label=label2) #specify the number of bins for the histogram
    pylab.hist([values1, values2], bins=num_bins, label=[label1, label2]) #specify the number of bins for the histogram
    pylab.title(title)
    pylab.xlabel("Value")
    pylab.ylabel("Number of values in bin")
    pylab.legend(loc='upper right')
    #pylab.show() #can do this instead of the savefig method if just want to view
    pylab.savefig(output_filename) #works for pdf or png
    pylab.close()
            

if __name__ == '__main__':
    # Generate example plot.

    # Define title for plot.
    title = "Histogram of random values in 30 bins"

    # Define input data for example plot.
    mu, sigma = 0, 0.1 # mean and standard deviation
    s = np.random.normal(mu, sigma, 1000)

    # Define output filepath.
    output_filepath = 'test_histogram_plot.pdf'

    # Call function to generate plot.
    generate_histogram(title,
                       s,
                       30,
                       output_filepath
                       )

    # Open output file.
    #subprocess.call(['open', output_filepath])

    # Delete output file.
    #os.remove(output_filepath)




