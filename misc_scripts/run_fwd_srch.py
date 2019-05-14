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
"""Runs forward searches without importing the pandas library or other
dependencies that may not be available on remote servers.

See the run_fwd_srch function in ameobae.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
import argparse
import settings
import module_amoebae_run_searches

parser = argparse.ArgumentParser(
    description='''Perform searches with original queries into subject
    databases.''')
#parser.add_argument('csv_file', help='''Path to summary spreadsheet
#        (CSV) file, which may already contain search summaries.''')
parser.add_argument('fwd_srch_dir', help='''Path to directory that will
        contain forward search output files.''')
args = parser.parse_args() # Different than in 'amoebae'.

# Get current time for timestamp.
#timestamp = time.strftime("%Y_%m_%d_%H_%M_%S")
#timestamp = time.strftime("%Y%m%d%H%M%S")

# Get paths to query and database directories from settings.py file.
query_dir = settings.querydirpath
db_dir = settings.dbdirpath

## Check that csv file is in the right place.
#assert args.srch_dir == os.path.dirname(args.csv_file), """Error: Input csv
#file is not in the specified search directory (srch_dir)."""

## Define output directory path, and make it.
#outdir = None
#if not args.outdir == None:
#    outdir = args.outdir
#else:
#    #outputname = os.path.basename(args.csv_file).rsplit('.', 1)[0] +\
#    #'_fwd_srch_' + timestamp 
#    outputname = 'fwd_srch_' + timestamp 

#    outdir = os.path.join(args.srch_dir, outputname) 
#os.mkdir(outdir)
#assert os.path.isdir(outdir), """Could not create directory: %s""" %\
#outdir

# Get query list from file.
query_file_list =\
module_amoebae_run_searches.get_query_list_from_file(module_amoebae_run_searches.get_out_query_list_path(args.fwd_srch_dir))

# Get database list from file.
db_file_list =\
module_amoebae_run_searches.get_db_list_from_file(module_amoebae_run_searches.get_out_db_list_path(args.fwd_srch_dir))

# Run searches.
module_amoebae_run_searches.run_all_searches(query_file_list, db_file_list,
        args.fwd_srch_dir)
