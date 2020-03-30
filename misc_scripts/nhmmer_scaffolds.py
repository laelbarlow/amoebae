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
"""This module is for using nhmmer to find the location of nucleotide sequences
matching an HMM, and outputs to a fasta file.

"""


import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import re
import subprocess
import itertools
import argparse
from get_fas_from_db_dir import get_seq_obj_from_db_fasta
import settings
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import nhmmer_scaffolds

#command_line_list = sys.argv
#hmm_path = str(command_line_list[1])
#prothmmpath = str(command_line_list[2])
#db_name = str(command_line_list[3])

# Set up argument parsing.
parser = argparse.ArgumentParser(
    description = """Searches in scaffolds with nhmmer.""",
    epilog = r"""epilog""")

# Positional arguments.
parser.add_argument('hmm_path', help='infilepath1')
#parser.add_argument('prothmmpath', help='infilepath2')
parser.add_argument('db_name', help='infilepath3')

# Optional arguments.
parser.add_argument('-n', '--no_search', help='Do not run hmmer if already run', action='store_true')

args = parser.parse_args()
            

if __name__ == '__main__':
    # Find query file path from db_name.
    db_path = nhmmer_scaffolds.get_db_fasta_filepath(settings.dbdirpath,
            args.db_name)


    # Run nhmmer on given files.
    nhmmer_outfilepath = args.hmm_path + '_' + args.db_name + '_out.txt'
    if args.no_search: 
        pass
    else:
        nhmmer_scaffolds.run_nhmmer_search(args.hmm_path, db_path, nhmmer_outfilepath)

    # Parse the nhmmer output file.
    hit_align_objs = nhmmer_scaffolds.parse_nhmmer_out(nhmmer_outfilepath)

    # Define path for output sequences.
    seq_filepath = nhmmer_outfilepath.rsplit('.', 1)[0] + '_seqs.fa'

    # Get concatenated sequence for aligning hits that may correspond to the
    # same gene.  nhmmer_scaffolds.get_nhmmer_hit_gene_region_seqs(hit_align_objs,
            seq_filepath)


    # Get temporary three-frame translations to search with the protein HMM.
    # or, add functionality to the get_nhmmer_hit_gene_region_seqs function...
    #three_frame_translation_dirpath = seq_filepath + '_temp_transl'
    #get_three_frame_transl(seq_filepath, three_frame_translation_dirpath)

    # Remove unnecessary translations after searching.





# Notes:

# Problem: Identifying the coding sequence for a gene, and getting an accurate
# translation of that sequence is problematic, because portions of the hmm
# alignment may align to regions of the genomic sequence that do not contain
# homologous codons.

# Maybe try iteratively trimming on either side of whatever aligns (in a local
# alignment) to the HMM so as to extend the local alignment until the whole HMM
# is aligned to the subject sequence. Maybe a bad idea.

# Basically need a way to get useful translatable coding sequences from an
# identified gene based on the nucleotide HMM and not on a huge gene
# re-prediction analysis...

# It may be that re-running nhmmer on a nucleotide region surrounding a gene
# may provide more detailed hit alignments that can be concatenated to
# construct a better sequence (and assuming that the whole HMM is in frame, a
# translation).

    # Assuming that the sequences in the alignment for the HMM were mRNA
    # sequences that started with a start codon, and that you trimmed the
    # alignment of mRNA sequences three nucleotides at a time, you should be
    # able to translate the identified subsequences based on how it aligns to
    # the HMM! Also, have to only trim three bases at a time to preserve
    # reading frame... (make an option for mask_nex.py and trim_nex.py?).

            

# Maybe the simplest way to get protein sequences out of genomic regions is to
# do a six-frame translation of the sequence, and concatenate regions that
# align the best to a protein HMM corresponding to the initial nucleotide HMM.
# This would circumvent the need to generate and search huge amounts of
# translated genomic sequence. The 6-frame translations could be deleted after
# searching, and also, only three frames would be needed if the strand was
# already known. Biopython's translate function could be used to generate the
# translations for searching with hmmer. Perhaps this would be best as a
# separate script that takes the output sequences from nhmmer_scaffolds.py, or
# would info be necessary? One challenge is propogating accurate info about
# what subsequences of the database scaffolds the resulting sequence
# corresponds to (gene model).

# To limit complexity, need to take a stepwise approach whereby functions
# process files generated by functions for the previous step, and output files
# used by the following steps. Instead of deeply nested loops. Maybe output a
# spreadsheet that tracks all the info in a clear way.

