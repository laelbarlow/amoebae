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
"""This is a script for finding the translation of each sequence in a FASTA
file of nucleotide mRNA sequences that minimizes the number of stop codons.
This may be useful when good quality predicted peptide sequences are not
available yet and you simply wish to do some preliminary searches in the
transcripts with peptide HMMs, for example.

The standard genetic code is used.

Biopython and Python3 must be installed.

If you want to have all three translations for each sequence, then set the
only_best_transl variable to False.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import copy

# Settings:
only_best_transl = True

cmdln = sys.argv
infile = cmdln[1]

# Define output filepath.
outfile = infile.rsplit('.', 1)[0] + '_roughly_translated.faa'
assert not os.path.isfile(outfile), """Make sure that you do not already have
an output file with the path %s.""" % outfile

def make_nucl_divisible_by_three(nucl_seq):
    """Append 'N' to given nucleotide sequence as necessary to make the length
    divisible by 3, and return the modified sequence object.
    """
    length = len(nucl_seq)
    while length % 3 != 0:
        nucl_seq = nucl_seq + 'N'
        length = len(nucl_seq)

    # Return modified (or unmodified sequence object).
    return nucl_seq


def remove_trailing_stops(inseq):
    """Take a peptide sequence object (translation of mRNA) and return a
    sequence object with '*' characters within the last 5 residues removed.
    This is so that the number of internal stop codons can be better compared
    between translations of different frames by ignoring trailing stops.
    """
    # Set number of trailing residues to check.
    num_trailing_residues = 10

    # Define original length.
    original_length = len(inseq)
    #print('\noriginal seq')
    #print(inseq[-10:].seq)

    # Check that the input sequence is at least as long as the number of
    # trailing residues of interest.
    if len(inseq.seq) < num_trailing_residues:
        return inseq

    else:
        # Copy the sequence object.
        inseq2 = copy.deepcopy(inseq)

        # Get the trailing residues.
        last_residues = inseq2[- num_trailing_residues:].seq

        # Remove '*' characters if present.
        other_residues = inseq2[:- num_trailing_residues].seq
        last_residues_without_stops = str(last_residues).replace('*', '') 
        num_stops_removed = len(last_residues) - len(last_residues_without_stops)
        inseq2.seq = other_residues + last_residues_without_stops

        #print('\noriginal length')
        #print(original_length)
        #print('number of stops removed')
        #print(num_stops_removed)
        #print('final length')
        #print(len(inseq))

        #print('non-trailing residues')
        #print(other_residues[-5:])
        #print('trailing residues')
        #print(' '*5 + last_residues)

        # Check that the sequence is the right length.
        assert len(inseq2) == original_length - num_stops_removed

        # Return the modified sequence object.
        return inseq2

def determine_longest_subseq_without_stops(inseq):
    """Take a peptide sequence object(translation of mRNA) and return the
    length of the longest subsequence with no stops.
    """
    seq_string = str(inseq.seq)
    split_seq = seq_string.split('*')
    length_of_longest_subseq = max([len(x) for x in split_seq])
    return length_of_longest_subseq


# Initiate variables to store the lists of internal stops and lengths of
# longest subsequences without stops.
internal_stops = []
longest_subseq_lengths = []

# Parse input FASTA file, translate each sequence, and write to output file.
with open(infile) as infh, open(outfile, 'w') as o:
    for mrna_seq in SeqIO.parse(infh, 'fasta'):
        # Identify relevant info from original sequence.
        original_id = mrna_seq.id
        original_name = mrna_seq.name
        original_description = mrna_seq.description

        # Make sequences that represent each frame for translation, and select
        # the translation with the least stop codons (*).

        # Evaluate translation of frame 1.
        mrna_seq_f1 = make_nucl_divisible_by_three(mrna_seq)
        transl1 = mrna_seq_f1.translate(to_stop=False)
        transl1_without_trailing_stops = remove_trailing_stops(transl1)
        transl1.id = original_id
        transl1.name = original_name
        transl1.description = original_description + ' (translation of frame 1)'
        transl1_stops = transl1_without_trailing_stops.seq.count('*')
        transl1_longest_subseq = determine_longest_subseq_without_stops(transl1)

        # Evaluate translation of frame 2.
        mrna_seq_f2 = make_nucl_divisible_by_three(mrna_seq[1:])
        transl2 = mrna_seq_f2.translate(to_stop=False)
        transl2_without_trailing_stops = remove_trailing_stops(transl2)
        transl2.id = original_id
        transl2.name = original_name
        transl2.description = original_description + ' (translation of frame 2)'
        transl2_stops = transl2_without_trailing_stops.seq.count('*')
        transl2_longest_subseq = determine_longest_subseq_without_stops(transl2)

        # Evaluate translation of frame 3.
        mrna_seq_f3 = make_nucl_divisible_by_three(mrna_seq[2:])
        transl3 = mrna_seq_f3.translate(to_stop=False)
        transl3_without_trailing_stops = remove_trailing_stops(transl3)
        transl3.id = original_id
        transl3.name = original_name
        transl3.description = original_description + ' (translation of frame 3)'
        transl3_stops = transl3_without_trailing_stops.seq.count('*')
        transl3_longest_subseq = determine_longest_subseq_without_stops(transl3)

        ## Set initial values for the translation with the least internal stops
        ## and longest subsequence without stops.
        #transl_with_least_stops = transl1
        #least_stops = transl1_stops
        #transl_with_longest_subseq = transl1
        #longest_subseq = transl1_longest_subseq

        ## Update the translation with the least internal stops and longest
        ## subsequence without stops.
        #if transl2_stops < transl1_stops:
        #    transl_with_least_stops = transl2
        #    least_stops = transl2_stops
        #if transl2_longest_subseq > transl1_longest_subseq:
        #    transl_with_longest_subseq = transl2
        #    longest_subseq = transl2_longest_subseq

        ## Update the translation with the least internal stops and longest
        ## subsequence without stops.
        #if transl3_stops < transl1_stops and transl3_stops < transl2_stops:
        #    transl_with_least_stops = transl3
        #    least_stops = transl3_stops
        #if transl3_longest_subseq > transl1_longest_subseq and\
        #transl3_longest_subseq > transl2_longest_subseq:
        #    transl_with_longest_subseq = transl3
        #    longest_subseq = transl3_longest_subseq

        # Find the translation with the fewest internal stops.
        transl_with_least_stops = sorted([transl1, transl2, transl3], key=lambda x:\
                remove_trailing_stops(x).seq.count('*'),
                reverse=True)[0]

        # Find the translation with the longest stop-free subsequence.
        transl_with_longest_subseq = sorted([transl1, transl2, transl3], key=lambda x:\
                determine_longest_subseq_without_stops(remove_trailing_stops(x)),
                reverse=True)[0]

        # Compare the translation with the fewest stops to the translation with
        # the longest stop-free subsequence to choose the best translation.
        best_translation = None
        if transl_with_least_stops.description == transl_with_longest_subseq.description:
            best_translation = transl_with_least_stops
        else:
            # Determine by what percentage the translation with the
            # longest stop-free subsequence has more internal stops than the
            # one with the fewest
            # stops has.
            percent_more_stops = None
            #print(remove_trailing_stops(transl_with_longest_subseq).seq.count('*') /\
            #            remove_trailing_stops(transl_with_least_stops).seq.count('*')) 

            if remove_trailing_stops(transl_with_least_stops).seq.count('*') == 0:
                percent_more_stops = 10000000 # Just larger than any percentage that you would get from the percent length difference.
            else:
                x = remove_trailing_stops(transl_with_least_stops).seq.count('*')
                percent_more_stops = -100 + 100 *\
                (remove_trailing_stops(transl_with_longest_subseq).seq.count('*') / x) 

            # Determine by what percentage the translation with the longest
            # subsequence has a longer subsequence than the translation with
            # the fewest stops.
            percent_longer_subseq = None
            if determine_longest_subseq_without_stops(remove_trailing_stops(transl_with_least_stops)) == 0:
                percent_longer_subseq = 10000000 # Just larger than any percentage that you would get from the percent length difference.
            else:
                percent_longer_subseq = -100 + 100 *\
                (determine_longest_subseq_without_stops(remove_trailing_stops(transl_with_longest_subseq)) /\
                    determine_longest_subseq_without_stops(remove_trailing_stops(transl_with_least_stops)))

            # Choose the translation that has the greater percent difference.
            if percent_more_stops > percent_longer_subseq:
                best_translation = transl_with_least_stops
            elif percent_more_stops < percent_longer_subseq:
                best_translation = transl_with_longest_subseq
            else:
                # Choose the translation with the fewest stops if the
                # percentages are equal.
                best_translation = transl_with_least_stops


        assert best_translation is not None

        # Add values to the lists of internal stops and lengths of longest
        # subsequences without stops.
        internal_stops.append(remove_trailing_stops(transl_with_least_stops).seq.count('*'))
        longest_subseq_lengths.append(determine_longest_subseq_without_stops(remove_trailing_stops(transl_with_longest_subseq)))

        # Temp.
        #assert transl_with_longest_subseq.description == transl_with_least_stops.description

        # Add text to description of translation with the fewest stops.
        transl_with_least_stops.description =\
        transl_with_least_stops.description + ' (fewest internal stops)'

        # Add text to description of translation with the longest
        # subsequence without stops.
        transl_with_longest_subseq.description =\
        transl_with_longest_subseq.description + ' (longest subsequence without stops)'

        # Either write just the translation with the lowest number of internal
        # stop codons, or write all three translations, depending on settings.
        if only_best_transl:

            # ***Could use different criteria here!

            # Write translation with the least stop codons to the output FASTA
            # file.
            #SeqIO.write([transl_with_least_stops], o, "fasta")
            SeqIO.write([best_translation], o, "fasta")

        else:
            # Write all three translations.
            SeqIO.write([transl1, transl2, transl3], o, "fasta")
        

# Print summary stats.
print('\nAverage number of internal stops in translations:')
print(sum(internal_stops)/len(internal_stops))
print('\nAverage length of longest subsequence without stops in translations:')
print(sum(longest_subseq_lengths)/len(longest_subseq_lengths))




