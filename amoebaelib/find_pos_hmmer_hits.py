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
"""Module for find_pos_hmmer_hits.py.
"""
import sys
import os
import argparse
import time
from datapaths import DataPaths
#from get_fas_from_db_dir import get_fas_from_db_dir # This is outdated!


def get_red_acc_list(infilepath):
    """Takes a file that contains redundant accessions and returns a list of
    the accessions it contains.
    """
    # Get redundant accessions from file.
    infilehandle = open(infilepath)
    l = []
    line_num = 0
    for line in infilehandle:
        line_num += 1
        if line_num == 1:
            l = line.rstrip('\n').rstrip(',').split(',')
        else:
            print('More than one line in redundant accessions file.')
    infilehandle.close()

    # Check whether redundant accessions identified.
    if l == []: print('No redundant accessions identified.')
    return l

def write_pos_seqs(infilepath, dbdirpath, outfilepath, prot_name=None):
    """Uses the get_fas_from_db_dir module to get all the sequences
    corresponding to hits in a HMMer_pipeline.py output .csv file.
    """
    with open(infilepath) as i, open(outfilepath, 'w') as o:
        acc_dict = {}
        last_db_name = None
        for line in i:
            split_line = line.split(',')
            db_name = split_line[2]
            acc = split_line[4]
            acc_list = [split_line[4]]
            if db_name != last_db_name:
                acc_dict[db_name] = [acc]
            elif db_name == last_db_name:
                acc_dict[db_name].append(acc)
            last_db_name = db_name
         
        for key in acc_dict.keys():
             o.write(get_fas_from_db_dir(key, acc_dict[key], dbdirpath,
                 prot_name))

def get_csv_name(infilepath, cur_time):
    return infilepath.replace('.csv', '_positive_' + cur_time + '.csv')

def find_pos_hmmer_hitsx(infilepath1, infilepath2, fwdeval, reveval,
        outfilepath, just_evalue):
    """Takes a HMMer_pipeline.py output spreadsheet and finds positive hits
    based on redundant accessions listed in another file.  Writes to another
    spreadsheet.
    """
    #print('\t' + os.path.basename(infilepath1))
    #print('\t' + os.path.basename(infilepath2))
    #print('\t' + str(fwdeval))
    #print('\t' + str(reveval))
    #print('\t' + os.path.basename(outfilepath))
    #print('\n')
    # Get list of redundant accessions from the second infile.

    red_acc_list = get_red_acc_list(infilepath2)

    # Set in and out file variables.
    infilehandle = open(infilepath1)
    outfilehandle = open(outfilepath, 'w')

    # Loop through lines in spreadsheet and write ones with one of the
    # redundant accessions in the top hit position to a new sheet.
    line_num = 0
    for line in infilehandle:
        # Identify the top hit accession from the rBLAST for each HMMer hit.
        line_num += 1
        if not line.startswith('Forward'):
            if not line.startswith('\n'):
                line_list = line.split(',')
                if len(line_list) > 6:
                    top_hit_acc = line_list[6]
                    positive = False

                    # If the top hit accession matches one of the redundant
                    # accessions, then write to the output spreadsheet.
                    for red_acc in red_acc_list:
                        #print('Does ' + red_acc + ' = ' + top_hit_acc + ' ?')
                        if red_acc.strip() == top_hit_acc.strip():
                            #print('Yes\n')
                            positive = True
                        else:
                            #print('No\n')
                            pass

                    # If the just_evalue option is set to True, then ignore
                    # what the top reverse blast hit is.
                    if just_evalue:
                        positive = True

                    if positive:
                        # Only write line if evalue criteria met.
                        fhmmer_e = line_list[5]
                        rblast_e = line_list[7]
                        if (float(fhmmer_e) <= float(fwdeval)) and \
                                (float(rblast_e) <= float(reveval)):
                            outfilehandle.write(line)

    # Close files.
    infilehandle.close()
    outfilehandle.close()


