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
"""Script for converting fastq file to fasta file, taking only a subsample of
the sequences.
"""
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys
import random


# Set values for sampling frequency and max number of samples.
just_convert_all_sequences = True
sampling_frequency = 0.04208 # Fraction of total sequences to sample.
max_samples = 1000000 # Maximum number of sequences to sample.


cmd = sys.argv
in1 = cmd[1]


# Code to check whether two sequence files have the same sequence. This
# verified that just converting text from parsed fastq object attributes is the
# same as converting to fasta with SeqIO.
#with open("first_100_lines.fq") as infh, open("new_file.fa", 'w') as o:
#    seqs = SeqIO.parse(infh, 'fastq')
#    for s in seqs:
#        SeqIO.write(s, o, 'fasta')
#cmd = sys.argv
#in1 = cmd[1]
#in2 = cmd[2]
#with open(in1) as i1, open(in2) as i2:
#    in1_seqs = SeqIO.parse(i1, 'fasta')
#    in2_seqs = SeqIO.parse(i2, 'fasta')
#    for i, j in zip(in1_seqs, in2_seqs):
#        print(i.seq)
#        print(j.seq)
#        assert i.seq == j.seq

# Define some functions.

def decide_whether_to_sample(sampling_frequency):
    """Returns a boolean value given a sampling frequency/probability,
    returning true with that frequency.
    """
    return random.random() < sampling_frequency


def get_random_sample_of_fasta_from_fastq(in1, outfilepath, sampling_frequency,
        max_samples):
    """Open a given FASTQ file and write a random sampling of sequences to a
    FASTA file. Not perfect, because won't always get the right number of
    samples if specifying the exact sampling frequency.
    """
    # Initiate tally of sampled sequences.
    sample_tally = 0
    
    # Initiate tally of sequences observed.
    fastq_seq_tally = 0
    
    # Initiate number of last sequence sampled.
    last_seq_sampled = None

    # Open files.
    with open(in1) as infh, open(outfilepath, 'w') as o:
        # Use lower-level FASTQ parser to increase speed.
        for title, seq, qual in FastqGeneralIterator(infh):
            fastq_seq_tally += 1
            # Determine whether to sample each sequence using a given sampling
            # frequency (increase the sampling frequency slightly so that
            # you get enough).
            sample = decide_whether_to_sample(sampling_frequency * 1.025)
            # Only sample if the function returns True.
            if sample and sample_tally < max_samples:
                o.write('>' + title + '\n' + seq + '\n')
                sample_tally += 1
                last_seq_sampled = fastq_seq_tally 
            #if sample_tally == max_samples:
            #    break
            #assert sample_tally < max_samples
        # Check that the right number of samples were retrieved
        assert sample_tally == max_samples
    
    # Print info.
    print('Total number of sequences sampled: ' + str(sample_tally))
    print('Total number of sequences: ' + str(fastq_seq_tally))
    print('Last sequence sampled: #' + str(last_seq_sampled))


def convert_fastq_to_fasta(in1, outfilepath):
    """Open a given FASTQ file and write all sequences to a fasta file with a
    given path.
    """
    # Initiate tally of sequences observed.
    fastq_seq_tally = 0
    # Open files.
    with open(in1) as infh, open(outfilepath, 'w') as o:
        # Use lower-level FASTQ parser to increase speed.
        for title, seq, qual in FastqGeneralIterator(infh):
            fastq_seq_tally += 1
            o.write('>' + title + '\n' + seq + '\n')

    # Report progress.
    print('\nWrote ' + str(fastq_seq_tally) + ' sequences to ' + outfilepath)


if not just_convert_all_sequences:
    # Randomly select sequences from a fastq file and write them as fasta sequences
    # to a new file.
    
    # Define output file path.
    outfilepath = in1.rsplit('.', 1)[0] + '_sample.fa'
    
    # Randomly sample sequences from the FASTQ file and write to a fasta file.
    sample_info = get_random_sample_of_fasta_from_fastq(in1, outfilepath,
            sampling_frequency, max_samples)

elif just_convert_all_sequences:
    # Define output file path.
    outfilepath = in1.rsplit('.', 1)[0] + '_converted.fa'
    
    # Call function to convert file.
    convert_fastq_to_fasta(in1, outfilepath)










