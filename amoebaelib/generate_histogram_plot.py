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
import re
import pylab
import subprocess
import numpy as np



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

def autolabel_bars(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def generate_bar_chart(title,
                       categories,
                       labels,
                       num_hits,
                       output_filename
                       ):
    """Take data and use matplotlib to generate a bar chart and write to a
    specified file path.
    """
    ## Simple bar chart.
    #fig, ax = pylab.subplots()
    #pylab.style.use('seaborn-deep')
    #pylab.rcdefaults()
    #fig, ax = pylab.subplots()
    # Example data
    #x_pos = np.arange(len(labels))
    #ax.barh(y_pos, performance, xerr=error, align='center')
    #ax.bar(x_pos, values, align='center')
    #ax.set_xticks(x_pos)
    #ax.set_xticklabels(labels)
    #ax.set_ylabel('Positive hit count')
    #ax.set_title(title)
    #pylab.show()
    #pylab.close()

    pylab.style.use('seaborn-deep')

    #categories = ['Prot', 'Nucl']
    #labels = ['Non-redundant', 'Final positive']
    #num_hits = [[35, 30],
    #            [12,  6]]
    
    assert len(labels) == len(num_hits)
    for sublist in num_hits:
        assert len(sublist) == len(categories)
    
    
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    
    fig, ax = pylab.subplots()
    rects_list = []
    num = 0
    for category, sublist in zip(categories, num_hits):
        num += 1
        if not (num % 2) == 0:
            rects = ax.bar(x - width/2, sublist, width, label=category)
            rects_list.append(rects)
        else:
            rects = ax.bar(x + width/2, sublist, width, label=category)
            rects_list.append(rects)

    # Add numbers to label individual bars.
    for r in rects_list:
        autolabel_bars(r, ax)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of sequences')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    
    fig.tight_layout()
    
    #pylab.show()
    pylab.savefig(output_filename) #works for pdf or png
    pylab.close()

            

if __name__ == '__main__':
    # Generate example plot.

    ## Define title for plot.
    #title = "Histogram of random values in 30 bins"

    ## Define input data for example plot.
    #mu, sigma = 0, 0.1 # mean and standard deviation
    #s = np.random.normal(mu, sigma, 1000)

    ## Define output filepath.
    #output_filepath = 'test_histogram_plot.pdf'

    ## Call function to generate plot.
    #generate_histogram(title,
    #                   s,
    #                   30,
    #                   output_filepath
    #                   )

    ## Open output file.
    #subprocess.call(['open', output_filepath])

    ## Delete output file.
    #os.remove(output_filepath)

    #
    #title = 'test bar chart'
    #values = [20, 10]
    #labels = ['prot', 'nucl']
    #output_filename = 'test_bar_chart.pdf'

    #generate_bar_chart(title,
    #                   values,
    #                   labels,
    #                   output_filename
    #                   )

    ## Open output file.
    #subprocess.call(['open', output_filepath])

    ## Delete output file.
    #os.remove(output_filepath)

    # Test bar chart.

    title = 'test bar chart'
    categories = ['Prot', 'Nucl']
    labels = ['Non-redundant', 'Final positive']
    num_hits = [[35, 30],
                [12,  6]]
    output_filename = 'test_bar_chart.pdf'

    generate_bar_chart(title,
                       categories,
                       labels,
                       num_hits,
                       output_filename
                       )

    # Open output file.
    subprocess.call(['open', output_filename])

    # Delete output file.
    os.remove(output_filename)

