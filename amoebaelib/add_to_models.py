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
"""Module for script: amoebae.
"""
# Import built-in modules.
#import argparse
import sys
import os
import subprocess
#import re
#import settings
import shutil
#import glob
import time
#import pandas as pd

# Import modules from installed libraries/packages.
from Bio import SeqIO
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped

from get_datatype import get_dbtype


# Define functions to be used in amoebae.

def update_models_csv(model_name, csv_file, alignmentfp, topologyfp, subs_model, type_seqsfp,
        taxon):
    """Appends a line to the given spreadsheet with info about the given fasta
    file added to a directory.
    """
    # Get current date.
    cur_date = time.strftime("%Y/%m/%d")

    # Check that no model with the same name already exists.
    existing_model_names = []
    with open(csv_file) as infh:
        for i in infh:
            if not i.startswith('\n'):
                n = i.split(',')[0]
                existing_model_names.append(n)

    ## Get list of type sequences and corresponding clade names.
    #type_seq_clade_name_list = []
    #with open(type_seqsfp) as infh:
    #    for i in infh:
    #        spliti = i.strip().split(',')
    #        type_seq_name = spliti[0]
    #        clade_name = spliti[1]
    #        type_seq_clade_name_list.append('[' + '\"' + type_seq_name +\
    #                '\"' + ',' + '\"' + clade_name + '\"' + ']')

    #type_seq_clade_name_list_string = '[' + ','.join(type_seq_clade_name_list) + ']'

    # Check whether file has a newline character at the end so a new line can
    # start correctly.
    newline = False
    if open(csv_file, 'r').read()[-1] == '\n':
        newline = True

    # Write to file.
    with open(csv_file, 'a') as csvh:
        x = '\n'
        if newline:
            x = ''
        csvh.write(x + model_name + ',' + cur_date + ',' + os.path.basename(alignmentfp) + ',' +\
                os.path.basename(topologyfp) + ',' + subs_model + ','\
                + taxon + ','\
                + os.path.basename(type_seqsfp) + '\n')

