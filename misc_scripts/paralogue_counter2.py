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
"""This script distinguishes predicted protein sequences that are likely
encoded on different loci (i.e., duplicated/paralogous) from those that are
encoded on the same locus (or different alleles of the same gene).

The challenge is to not propogate errors that originated at the genome assembly
and gene prediction stages of genomic analysis, and get an accurate estimate of
the true number of orthologues.

This is particularly relevant for classifying sequences that are likely encoded
on genes that duplicated very recently (within a specific species or genus).
For genes that duplicated relatively early, this task is much easier, because
several species will have both duplicates, as determined via phylogenetic
analysis, so there wouldn't be a question as to whether they were duplicated or
not.

******This script is now redundant with the amoebae file.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import argparse
import settings
from paralogue_counter import count_paralogues2


# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Identifies redundant gene models (predicted peptide
    sequences).""",
    epilog = r"""epilog""")

# Positional arguments.
#parser.add_argument('alignmentdir', help="""Directory path to trimmed alignment
#used for generating backbone tree using phylogenetic analysis, or HMM.""")
parser.add_argument('csv_file', help="""Spreadsheet listing hits to be
classified.""")
parser.add_argument('alignmentdir', help="""Path to directory with alignment
files (nexus format?).""")
parser.add_argument('fastadir', help="""Path to directory with fasta files with
sequences to be classified.""")

# Optional arguments.
parser.add_argument('-f', '--fwdeval', help="""E-value cutoff for similarity
        search hits""", default=0.05)
parser.add_argument('-m', '--metric_name', help="""Name of metric to be used
        for evaluating whether two sequences are redundant or not. Can be
        either 'percent_identity' or 'gene_model_redundancy_index'.""",\
                default='percent_identity')
parser.add_argument('-v', '--metric_value_minimum', help="""Minimum value of
        specified metric for evaluating whether two sequences are redundant or
        not.""", default=98.0)

args = parser.parse_args()


## Define path to alignment directory (Query directory). # (Eventually)
#alignmentdir = settings.querydirpath


# Call the function from the module to count paralogues and output to a
# spreadsheet.
#count_paralogues2(args.alignment, args.fasta, args.fwdeval, args.identity)
count_paralogues2(args.csv_file, args.alignmentdir, args.fastadir,\
        args.fwdeval, args.metric_name, args.metric_value_minimum)


