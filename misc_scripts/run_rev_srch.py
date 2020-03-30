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
"""Runs reverse searches without importing the pandas library or other
dependencies that may not be available on remote servers.

See the run_rev_srch function in ameobae.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
import argparse
import settings
import run_searches


parser = argparse.ArgumentParser(
    description='''Perform searches with forward search hit sequences
    as queries into the original query databases.''')
parser.add_argument('rev_srch_dir', help='''Path to directory that will
        contain output of searches.''')
args = parser.parse_args()


# Get query file path list from file.
query_file_list =\
run_searches.get_query_list_from_file(run_searches.get_out_query_list_path(args.rev_srch_dir))

# Get database file path list from file.
db_file_list =\
run_searches.get_db_list_from_file(run_searches.get_out_db_list_path(args.rev_srch_dir))

# For reverse searches, the queries are in a subdirectory of the
# reverse search directory, not in the query file directory specified
# in the settings.py file.
query_dir = run_searches.get_query_subdir(args.rev_srch_dir)

# Run searches (query_file_list contains full filepaths, so can use
# query_subdir contents).
run_searches.run_all_searches(query_file_list, db_file_list,
        args.rev_srch_dir, query_dir)

# Report output.
print('\n\nReverse search results written to directory:\n\t' +
        args.rev_srch_dir\
        + '\n\n')
