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

"""Classifies amino acid sequences using a constrained topology (supported by
previous phylogenetic analysis).

This is perhaps a useful alternative to trying to classify sequences using
similarity search results or running many trees manually, it is perhaps even an
alternative to reverse searching as a means of evaluating the specificity of
sequences similar to a query sequence or HMM.

The purpose initially is to simply provide clear evidence for classifying SNARE
sequences that do not make it into my final tree figures to show in
publications.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
#from Bio import AlignIO
#from Bio import SeqIO
#from Bio.Alphabet import IUPAC, Gapped
#from module_afa_to_nex import delete_extra_mesquite_lines, afa_to_nex, nex_to_afa, nex_to_phylip
#import numpy as np
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#import collections
import argparse
#import subprocess
#import time
#import glob
#from ete3 import Tree
#import re
#from trim_nex import trim_nex
from module_amoebae_phylo_clas import classify_seq_with_constrained_tree


# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Classifies sequences using a backbone tree.""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('alignment', help="""Trimmed alignment used for generating
backbone tree using phylogenetic analysis.""")
parser.add_argument('tree', help="""Backbone tree used as basis for constrained
topologies.""")
parser.add_argument('type_seqs', help="""Sequence names for sequences in the
input tree topology that define clades of interest (.txt file with seq names on
lines).""")
parser.add_argument('subs_model', help="""Phylogenetic model determined to be
the best fit for the original alignment.""")
parser.add_argument('fasta', help="""Fasta file with sequences to be
classified.""")
#parser.add_argument('x', help='infilepath2')

# Optional arguments.
# Not done yet:
#parser.add_argument('--boot', help="""Instead of putting each sequence into
#        each clade, runs one ML search with original tree as constraint tree to #        see where the input sequence gets placed.""")

args = parser.parse_args()


classify_seq_with_constrained_tree(args.alignment, args.tree, args.subs_model,
        args.type_seqs, args.fasta)

        



# Ideas:

# split up script so long part can be run on bioinfor without having to install
# too much?

# Add clade names to type_sequences input file (make comma delimited?).

# Need to compare directly against similarity-searching-based approaches to classification.

# What happens when you take out a clade from the input tree? False-positives right?

# Try taking out sequences from input tree one at a time and seeing how well
# the sequences you take out fit back in when running the pipeline. A type of
# cross-validation or jackknifing?

# Do I get significantly different results when I run trees multiple times? Should I do best of several ML searches?
# Add the option of running multiple ML searches and choosing the best one (is that available in IQtree?) so that people can do that if they want to spend more time.

# Would be useful to be able to somehow rank the sequences in the fasta file
# according to how well they might contribute to the tree.

# ***Instead of forcing the sequence into specific clades, try just letting it go anywhere in the tree, by not putting it into the constraint tree. And then parsing the output tree to see where it went. could do bootstrapping in this case too.
