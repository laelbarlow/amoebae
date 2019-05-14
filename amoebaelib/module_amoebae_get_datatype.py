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
"""Module defining function(s) for inferring the data type from a give file.
"""
import re
import sys
from Bio import AlignIO


def get_datatype_for_sequence_string(concat_seq):
    """Take a sequence string and return 'prot' or 'nucl' depending on the
    content of letters from the alphabet.
    """
    # Define regular expressions to identify different data types from sequence
    # strings.
    #alphabets = {'dna': re.compile('^[acgtn]*$', re.I),
    #             'protein': re.compile('^[arndceqghoilkmfpustwyvx*]*$',
    #                          re.I)}
    alphabets = {'dna': re.compile('^[acgtny]*$', re.I),
                 'protein': re.compile('^[arndceqghoilkmfpustwyvx*]*$',
                              re.I)}
    
    # Determine whether the concatenated sequence represents protein or DNA
    # data.
    dbtype = 'Undetermined'
    if alphabets['dna'].search(concat_seq) is not None:
        dbtype = 'nucl'
    elif alphabets['protein'].search(concat_seq) is not None:
        dbtype = 'prot'

    # Check that the data type was determined.
    assert not dbtype == 'Undetermined', "Error: could not identify data type"

    # Return the datatype.
    return dbtype


def get_dbtype(f):
    """Opens input fasta file, examines the first 10 records, and returns
        either 'prot' or 'nucl' depending on what type the sequences are.
        """
    # Get sequence to evaluate.
    concat_seq = ''
    
    # Determine file type.
    filetype = None
    with open(f) as infh:
        for i in infh:
            if i.startswith('>'):
                filetype = 'fasta'
                break
        else:
            # Assume nexus if not fasta.
            filetype = 'nexus'

    # Open fasta file and get first 100 lines as a string.
    if filetype == 'fasta':
        with open(f) as infh:
            line_num = 0
            for i in infh:
                line_num += 1
                if line_num <= 200:
                    if not i.startswith('>'):
                        concat_seq = concat_seq + i.strip().replace('-', '')
                else:
                    break

    # Open nexus file and get concatenation of first 10 sequences (without '-'
    # characters) as a string. Note: This was changed because the previous
    # version was unable to handle larger alignments with interleaved data
    # (which meant that it would only retrieve gaps as sample sequence).
    elif filetype == 'nexus':
        with open(f) as infh:
            alignment = AlignIO.read(infh, "nexus")
            seq_num = 0
            for i in alignment:
                seq_num += 1
                sequence_string = str(i.seq).replace('-', '')
                concat_seq = concat_seq + sequence_string
                if seq_num > 10:
                    break

    # Remove additional characters to be ignored from the concatenated
    # sequence.
    concat_seq = concat_seq.replace('?', '')

    # Determine whether the concatenated sequence represents protein or DNA
    # data.
    dbtype = get_datatype_for_sequence_string(concat_seq)
    
    # Return data type.
    return dbtype



