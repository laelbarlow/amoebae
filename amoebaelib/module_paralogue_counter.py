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
"""This module contains functions to be used in the paralogue_counter.py
script (or other scripts). 
"""
import argparse
import re
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from module_amoebae_nex_to_hmm import nex_to_hmm
import os
import subprocess
from module_afa_to_nex import afa_to_nex, delete_extra_mesquite_lines
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.SeqRecord import SeqRecord
import time
import pandas as pd
import settings
import glob
from module_similarity_score import get_similarity_score, get_score_dataframe_from_file
import module_amoebae
import gffutils
from module_amoebae import mask_nex2
from module_amoebae import get_seq_obj_from_srch_res_csv_info,\
get_hit_range_from_hsp_ranges 
from module_amoebae_trim_nex import trim_nex


def remove_irrelevant_data_from_alignment(ali_plus_top_seq_plus_1):
    """Remove all but the last two sequences, and remove all the positions that
    only contain gaps.

    ***This function is not in use, because as currently written it will crash
    the computer. I don't know why.
    """
    # Parse the input alignment.
    alignment = AlignIO.read(ali_plus_top_seq_plus_1, 'nexus')

    # Remove the original file.
    os.remove(ali_plus_top_seq_plus_1)

    # Remove all but the last two sequences from the alignment.
    alignment = alignment[-2:]


    # Remove positions in the alignment that have only gaps for both
    # sequences.

    # get a list of columns as strings in the original alignment.
    seq_len = alignment.get_alignment_length()
    columns = [alignment[:, col] for col in range(seq_len)] 

    # Get a list of column indexes for columns that need to be remove.
    col_to_remove = []
    col_index = -1
    for col in columns:
        col_index += 1
        if col[0] == '-' and col[1] == '-':
            col_to_remove.append(col_index)

    
    #seq_len = alignment.get_alignment_length()
    #print(seq_len)
    #print(alignment[:,col_to_remove[0]])

    # Remove columns (positions) from alignment.
    num = 0
    for col_index in col_to_remove:
        num += 1
        adj_col_index = col_index - num
        alignment = alignment[:,:adj_col_index] + alignment[:,adj_col_index + 1:]
        #if col_index != (seq_len -1):
        #    alignment = alignment[:,:adj_col_index] + alignment[:,adj_col_index + 1:]
        #else:
        #    alignment = alignment[:,:adj_col_index]

    #seq_len = alignment.get_alignment_length()
    #print(seq_len)
    #print(alignment[:,col_to_remove[0]])

    # Write modified alignment to the same filepath as the original.
    AlignIO.write(alignment, ali_plus_top_seq_plus_1, 'nexus')



def check_for_redundant_seq_names(alignment):
    """Take an alignment, and assert that there are no sequences with identical
    names in the alignment.
    """
    pass


def rank_seqs_by_hmm(inali, inseqs):
    """Takes an HMM and fasta sequences, and ranks the sequences according to
    similarity to the HMM, outputting a list of Bio.SeqIO sequence objects.
    """
    # Convert nexus alignment to HMM.
    hmmname = inali.rsplit('.', 1)[0] + '_TEMP.hmm'
    hmmstdoutname = hmmname + '.stdout.txt'
    nex_to_hmm(inali, hmmname)

    # Search input fasta file with HMM.
    hmmoutput = hmmname + '_out.txt'
    with open(hmmstdoutname, 'w') as o:
        subprocess.call(['hmmsearch', '--tblout', hmmoutput, hmmname, inseqs],\
                stdout=o, stderr=subprocess.STDOUT)

    # Parse output, and get ranked list of accessions.
    space_char = re.compile(r' +')
    ranked_acc_dict = {}
    inum = 0
    with open(hmmoutput) as infh:
        for i in infh:
            if not i.startswith('#'):
                inum += 1
                parsed_i = space_char.split(i.strip())
                acc = parsed_i[0]
                ranked_acc_dict[acc] = inum
    #assert ranked_acc_dict != {}, """Error: could not parse accessions from
    #sequence objects."""

    # Get list of sequence objects from input fasta file.
    num_seqs = 0
    seq_objs = []
    for seq in SeqIO.parse(inseqs, 'fasta'):
        num_seqs += 1
        seq_objs.append(seq)

    # Check that all the sequences were identified as hits with HMMer.
    assert inum == num_seqs, """Error: Not all sequences in input fasta file
    are identified as hits for the HMM."""

    # Rank the list of sequence objects according to similarity to the HMM.
    seq_objs.sort(key=lambda x: ranked_acc_dict[x.id], reverse=False)

    # Delete intermediate files.
    os.remove(hmmname)
    os.remove(hmmoutput)
    os.remove(hmmstdoutname)

    # return ranked list of sequence objects.
    return seq_objs


def get_sig_overlap(inali, score_df):
    """Return true if there is significant overlap between two sequences in an
    alignment.
    """
    minimum_aligning_residues = 20
    minimum_identical_residues = 10
    minimum_similar_residues = 15
    sig_overlap = False

    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(inali)

    #percent_identity = None

    with open(inali) as infh:
        # Check that the input file has the filename extension ".nex".
        assert inali.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')

        # get a list of columns as strings in the original alignment.
        seq_len = alignment.get_alignment_length()
        columns = [alignment[:, col] for col in range(seq_len)] 

        # Count columns where the last two sequences have aligned residues, and
        # where they have identical residues.
        aligned_residues_count = 0
        identical_residues_count = 0
        similar_not_ident_residues_count = 0
        for col in columns:
            last_two = col[-2:]
            if not '-' in last_two:
                aligned_residues_count += 1
                if last_two[0] == last_two[1]:
                    identical_residues_count += 1
                elif get_similarity_score(last_two[0], last_two[1], score_df) > 0:
                    similar_not_ident_residues_count += 1

        # Calculate percent identity.
        #percent_identity = (identical_residues_count/aligned_residues_count)*100

        # Temporary print statements.
        #print('aligned_residues_count: ' + str(aligned_residues_count))
        #print('identical_residues_count: ' + str(identical_residues_count))
        #print('similar_not_ident_residues_count: ' + str(similar_not_ident_residues_count))

        # Determine whether there is significant overlap.
        if aligned_residues_count >= 20 and identical_residues_count >= 10 and\
            (identical_residues_count + similar_not_ident_residues_count) >= 15:
            sig_overlap = True

    # Return percent identity.    
    return sig_overlap


# Unused??
def get_percent_identity(inali):
    """Take an alignment, and return the percent identity between the last two
    sequences in an alignment.
    """
    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(inali)

    percent_identity = None

    with open(inali) as infh:
        # Check that the input file has the filename extension ".nex".
        assert inali.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')

        # get a list of columns as strings in the original alignment.
        seq_len = alignment.get_alignment_length()
        columns = [alignment[:, col] for col in range(seq_len)] 

        # Count columns where the last two sequences have aligned residues, and
        # where they have identical residues.
        aligned_residues_count = 0
        identical_residues_count = 0
        for col in columns:
            last_two = col[-2:]
            if not '-' in last_two:
                aligned_residues_count += 1
                if last_two[0] == last_two[1]:
                    identical_residues_count += 1

        # Calculate percent identity.
        percent_identity = (identical_residues_count/aligned_residues_count)*100

    # Return percent identity.    
    return percent_identity

    # Could calculate similarity as well?


def add_seq_to_alignment2(inseqobj, innexpath, outnexpath):
    """Takes a sequence and adds it to an alignment.
    """
    # Write input sequence object to a file temporarily for input to MUSCLE.
    fa_temp_1 = os.path.join(os.path.dirname(outnexpath), os.path.basename(innexpath) + '_temp1.fa')
    SeqIO.write(inseqobj, fa_temp_1, 'fasta')

    # Temporarily save a copy of the input alignment (nexus format) to fasta
    # file.
    # Check input format.
    informat = None
    with open(innexpath) as infh:
        for i in infh:
            if i.startswith('>'):
                informat = 'fasta'
                break
            else:
                informat = 'nexus'
    # Convert if nexus.
    inali = None
    if informat == 'nexus': 
        # Delete extra lines in input nexus file, if present, because biopython cannot
        # read nexus alignments with these extra lines.
        delete_extra_mesquite_lines(innexpath)

        # Read the alignment file.
        inali = AlignIO.read(innexpath, 'nexus')
    else:
        inali = AlignIO.read(innexpath, 'fasta')
    # Define temporary fasta alignment file path.
    fa_temp_2 = fa_temp_1.replace('_temp1.fa', '_temp2.afa')
    # Write alignment to path in fasta format. 
    AlignIO.write(inali, fa_temp_2, 'fasta')

    # Check that the input sequence does not have the same name as any of the
    # sequences in the input alignment.
    # Get list of sequence names from alignment.
    existing_seq_names = []
    with open(fa_temp_2) as infh:
        for record in AlignIO.read(infh, 'fasta'):
            extra_text = re.compile(r'\.copy$')
            rec_id_no_extra = extra_text.sub('', record.description) 
            existing_seq_names.append(rec_id_no_extra)
            #existing_seq_names.append(record.id)
    # Get additional sequence name.
    additional_seq_name = inseqobj.description

    # Align sequence to alignment (make sure that gap penalties are optimized).
    fa_temp_3 = fa_temp_2.replace('_temp2.afa', '_temp3.afa')
    musclestdoutname = fa_temp_3 + '_musclestdout.txt'
    with open(musclestdoutname, 'w') as o:
        subprocess.call(["muscle", "-profile", "-in1", fa_temp_2, "-in2",
            fa_temp_1, "-out", fa_temp_3], stdout=o, stderr=subprocess.STDOUT)

    assert os.path.isfile(fa_temp_3), """The alignment file %s was not
    produced. Check that the alignment does not contain '?' characters as
    residues for any sequences.""" % os.path.basename(fa_temp_3)

    # Convert alignment to nexus format.
    nex_out_1 = outnexpath
    afa_to_nex(fa_temp_3, nex_out_1) 

    ## Mask alignment (mask out any positions in which the last sequence, the
    ## additional sequence, has a residue but none of the other sequences have
    ## residues at that position).
    #nex_out_2_mask = nex_out_1.rsplit('.', 1)[0] + '_mask.nex'
    #mask_nex2(nex_out_1, nex_out_2_mask) 

    ## Trim alignment.
    #nex_out_3_trim = nex_out_2_mask.rsplit('_', 1)[0] + '_trim.nex'
    #trim_nex(nex_out_2_mask, nex_out_3_trim)

    # Delete temporary files.
    os.remove(fa_temp_1)
    os.remove(fa_temp_2)
    os.remove(fa_temp_3)
    os.remove(musclestdoutname)
    #os.remove(nex_out_1)


def modify_seq_descr_for_tree(inseqobj):
    """Take a SeqIO sequence object, and modify the header for easy parsing in
    alignment and tree files.
    """
    # Set the ID to the a modification of the original description (assuming
    # that the ID is included in the original description).
    if ' ' in inseqobj.description:
        inseqobj.id = inseqobj.description.split(' ', 1)[0] + '__'\
            + inseqobj.description.split(' ', 1)[1].replace(' ', '_')
    else:
        inseqobj.id = inseqobj.description.split(' ', 1)[0] + '__'\
            + inseqobj.description.split('_', 1)[1].replace(' ', '_')
    # Set the description to an empty string.
    inseqobj.description = ''


def add_seq_to_alignment3(inseqobj, innexpath, outnexpath):
    """Takes a sequence and adds it to an alignment. Make sure that the
    sequence header is formatted correctly (e.g., use modify_seq_descr_for_tree
    function).
    """
    # Write input sequence object to a file temporarily for input to MUSCLE.
    fa_temp_1 = os.path.join(os.path.dirname(outnexpath), os.path.basename(innexpath) + '_temp1.fa')
    SeqIO.write(inseqobj, fa_temp_1, 'fasta')

    # Temporarily save a copy of the input alignment (nexus format) to fasta
    # file.
    # Check input format.
    informat = None
    with open(innexpath) as infh:
        for i in infh:
            if i.startswith('>'):
                informat = 'fasta'
                break
            else:
                informat = 'nexus'
    # Convert if nexus.
    inali = None
    if informat == 'nexus': 
        # Delete extra lines in input nexus file, if present, because biopython cannot
        # read nexus alignments with these extra lines.
        delete_extra_mesquite_lines(innexpath)

        # Read the alignment file.
        inali = AlignIO.read(innexpath, 'nexus')
    else:
        inali = AlignIO.read(innexpath, 'fasta')
    # Define temporary fasta alignment file path.
    fa_temp_2 = fa_temp_1.replace('_temp1.fa', '_temp2.afa')
    # Write alignment to path in fasta format. 
    AlignIO.write(inali, fa_temp_2, 'fasta')

    # Check that the input sequence does not have the same name as any of the
    # sequences in the input alignment.
    # Get list of sequence names from alignment.
    existing_seq_names = []
    with open(fa_temp_2) as infh:
        for record in AlignIO.read(infh, 'fasta'):
            extra_text = re.compile(r'\.copy$')
            rec_id_no_extra = extra_text.sub('', record.description) 
            existing_seq_names.append(rec_id_no_extra)
            #existing_seq_names.append(record.id)
    # Get additional sequence name.
    additional_seq_name = inseqobj.description

    # Align sequence to alignment (make sure that gap penalties are optimized).
    fa_temp_3 = fa_temp_2.replace('_temp2.afa', '_temp3.afa')
    musclestdoutname = fa_temp_3 + '_musclestdout.txt'
    with open(musclestdoutname, 'w') as o:
        subprocess.call(["muscle", "-profile", "-in1", fa_temp_2, "-in2",
            fa_temp_1, "-out", fa_temp_3], stdout=o, stderr=subprocess.STDOUT)

    assert os.path.isfile(fa_temp_3), """The alignment file %s was not
    produced. Check that the alignment does not contain '?' characters as
    residues for any sequences.""" % os.path.basename(fa_temp_3)

    # Convert alignment to nexus format.
    nex_out_1 = outnexpath.rsplit('.', 1)[0] + '_untrimmed.nex'
    afa_to_nex(fa_temp_3, nex_out_1) 

    # Mask alignment (mask out any positions in which the last sequence, the
    # additional sequence, has a residue but none of the other sequences have
    # residues at that position).
    nex_out_2_mask = nex_out_1.rsplit('.', 1)[0] + '_mask.nex'
    mask_nex2(nex_out_1, nex_out_2_mask) 

    # Trim alignment.
    #nex_out_3_trim = nex_out_2_mask.rsplit('_', 1)[0] + '_trim.nex'
    #trim_nex(nex_out_2_mask, nex_out_3_trim)
    trim_nex(nex_out_2_mask, outnexpath)

    # Delete temporary files.
    os.remove(fa_temp_1)
    os.remove(fa_temp_2)
    os.remove(fa_temp_3)
    os.remove(musclestdoutname)
    os.remove(nex_out_1)


def get_sequence_comparison_dict(inali, score_df):
    """Takes an alignment in nexus format and determines if the last two
    sequences overlap, etc.

    Overlap could be defined in various ways, depending on what is relevant...
    """
    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(inali)

    # Initiate definitions of relevant variables.
    aligned_residues_count = 0
    identical_residues_count = 0
    similar_not_ident_residues_count = 0
    num_continuous_identical_spans = 0
    length_of_longest_identical_span = 0

    # Extract the relevant information from the alignment file and use to
    # assign values to variables.
    with open(inali) as infh:
        # Check that the input file has the filename extension ".nex".
        assert inali.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')

        # get a list of columns as strings in the original alignment.
        seq_len = alignment.get_alignment_length()
        columns = [alignment[:, col] for col in range(seq_len)] 

        # Count columns where the last two sequences have aligned residues, and
        # where they have identical residues, and where they have similar but
        # not identical residues.
        for col in columns:
            last_two = col[-2:]
            if not '-' in last_two:
                aligned_residues_count += 1
                if last_two[0] == last_two[1]:
                    identical_residues_count += 1
                elif get_similarity_score(last_two[0], last_two[1], score_df) > 0:
                    similar_not_ident_residues_count += 1

        # Get info about continuous spans of identical residues.
        seq2_spans_of_identical_res= []
        prev_identical = False
        cur_span = ''
        col_num = 0
        last_col_num = len(columns)
        for col in columns:
            col_num += 1
            last_two = col[-2:]
            aa1 = last_two[0]
            aa2 = last_two[1] 
            # Ignore positions where there are gaps in both sequences.
            if not aa1 == '-' and not aa2 == '-':

                # Determine whether the position is identical.
                identical = False
                if aa1.upper() == aa2.upper():
                    identical = True

                # Add to existing span or restart a span as necessary.

                if identical:
                    if prev_identical:
                        # Add residue to existing span of identical residues.
                        cur_span = cur_span + aa2

                    elif not prev_identical:
                        # Restart a span of identical residues.
                        cur_span = aa2

                    if col_num == last_col_num:
                        # If on the last residue, then it's the end of the
                        # span, so add to the list.
                        seq2_spans_of_identical_res.append(cur_span)

                elif not identical:
                    if prev_identical:
                        # Then reached end of a span.
                        # Add previous span of identical residues to list.
                        if col_num != 1:
                            seq2_spans_of_identical_res.append(cur_span)

                # Assign current value to "previous" variable.
                prev_identical = identical

            # Account for possibility that sequences are identical and the
            # alignment does not end in a position where the sequences have
            # aligned residues (both are gaps).
            else:
                if col_num == last_col_num:
                    if cur_span != '':
                        seq2_spans_of_identical_res.append(cur_span)

        # Determine relevant information from list of identical spans.
        num_continuous_identical_spans = len(seq2_spans_of_identical_res)
        length_of_longest_identical_span = len(sorted(seq2_spans_of_identical_res,\
                key=lambda x: len(x), reverse=True)[0])

    # Calculate percent identity.
    percent_identity = (identical_residues_count/aligned_residues_count)*100

    # Calculate percent similarity (not counting identical positions).
    percent_similarity = (similar_not_ident_residues_count/aligned_residues_count)*100

    # Calculate percent similarity + identity.
    percent_similarity_plus_identity = ((identical_residues_count + similar_not_ident_residues_count)/aligned_residues_count)*100

    # Calculate gene model redundancy index value.

    #gene_model_redundancy_index =\
    #(identical_residues_count * length_of_longest_identical_span)\
    #    / (1 + (similar_not_ident_residues_count * num_continuous_identical_spans))

    #gene_model_redundancy_index =\
    #(identical_residues_count * length_of_longest_identical_span)\
    #    / ((1 + similar_not_ident_residues_count) * num_continuous_identical_spans)

    #gene_model_redundancy_index =\
    #(identical_residues_count + length_of_longest_identical_span)\
    #    / ((1 + similar_not_ident_residues_count) + num_continuous_identical_spans)

    #gene_model_redundancy_index =\
    #((aligned_residues_count - similar_not_ident_residues_count)*100)\
    #    / aligned_residues_count 

    #gene_model_redundancy_index =\
    #((aligned_residues_count - similar_not_ident_residues_count)\
    #    / aligned_residues_count) + length_of_longest_identical_span 

    #gene_model_redundancy_index =\
    #(identical_residues_count +\
    #        (length_of_longest_identical_span/(identical_residues_count)))\
    #    / ((1 + similar_not_ident_residues_count) +\
    #            (num_continuous_identical_spans\
    #        /identical_residues_count))

    w1 = 10
    w2 = 1
    w3 = 10
    w4 = 1

    gene_model_redundancy_index =\
        (identical_residues_count/aligned_residues_count)*w1 \
        - (similar_not_ident_residues_count/aligned_residues_count)*w2 \
        + (length_of_longest_identical_span/aligned_residues_count)*w3 \
        - (num_continuous_identical_spans/aligned_residues_count)*w4 



    # Construct dictionary with info to return.
    info_dict = {}
    info_dict['aligned_residues_count'] = aligned_residues_count
    info_dict['identical_residues_count'] = identical_residues_count
    info_dict['percent_identity'] = percent_identity 
    info_dict['similar_not_ident_residues_count'] = similar_not_ident_residues_count
    info_dict['percent_similarity'] = percent_similarity
    info_dict['percent_simlarity_plus_identity'] = percent_similarity_plus_identity
    info_dict['num_continuous_identical_spans'] = num_continuous_identical_spans
    info_dict['length_of_longest_identical_span'] = length_of_longest_identical_span
    info_dict['gene_model_redundancy_index'] = gene_model_redundancy_index

    # Return count tuple.    
    return info_dict


def get_all_comparison_output_filepath(outdir):
    """Define path for output file (csv format) that lists all sequence
    comparisons made, and their details.
    """
    return os.path.join(outdir, '0_all_sequence_comparisons.csv')


def find_redun_model_recursively(outdir, cur_ali_num, inali, ranked_seq_objs,
        redundant_gene_model_dict, metric_name, metric_value_minimum,
        minimum_percent_similarity, extra_info_dict, overlap_required):
    """Recursively process sequence objects, classifying them into groups of
    redundant predictions for single loci/alleles.

            info_tuple = (index,
                          acc,
                          seq_obj,
                          evalue, 
                          program,
                          percent_len,
                          hit_range,
                          actual_acc,
                          sequence_filename
                          )
    """
    # Take the top hit as a representative of a true unique locus.
    top_hit_seq_obj = ranked_seq_objs[0]


    redundant_gene_model_dict[top_hit_seq_obj.id][2] =\
        redundant_gene_model_dict[top_hit_seq_obj.id][2] +\
        'likely represents a unique (paralogous) gene locus (see alignment ' + str(cur_ali_num) + '); '

    # Remove top hit sequence object from list.
    ranked_seq_objs.remove(top_hit_seq_obj)

    # Align it to the HMM alignment (using profile-profile alignment).
    ali_plus_top_seq = os.path.join(outdir, str(cur_ali_num) + '.nex')
    cur_ali_num += 1
    add_seq_to_alignment2(top_hit_seq_obj, inali, ali_plus_top_seq)

    # Print input sequence object ids.
    ranked_seq_accs = []
    for i in ranked_seq_objs:
        ranked_seq_accs.append(i.id)
    
    # For each subsequent hit: 
    acc_for_objs_to_remove = []
    for seq_obj in ranked_seq_objs:
        # Get sequence object accessions.
        acc = seq_obj.id

        # Align to the HMM alignment plus top hit.
        ali_plus_top_seq_plus_1 = os.path.join(outdir, str(cur_ali_num) + '.nex')
        cur_ali_num += 1
        add_seq_to_alignment2(seq_obj, ali_plus_top_seq, ali_plus_top_seq_plus_1)

        # Get a dataframe for interpreting amino acid similarity in alignment.
        score_df = get_score_dataframe_from_file(\
                os.path.join(os.path.dirname(os.path.realpath(__file__)), 'blosum62.csv'))

        # Determine whether sequences overlap significantly.
        sig_overlap = get_sig_overlap(ali_plus_top_seq_plus_1, score_df)

        # If no significant overlap, then note as not necessarily represenative of a
        # different locus from the current "top hit".
        if not sig_overlap:
            if overlap_required:
                # Record lack of overlap.
                redundant_gene_model_dict[acc][2] =\
                redundant_gene_model_dict[acc][2] + 'No significant overlap ' +\
                'in alignment with sequence ' + top_hit_seq_obj.id +\
                ' therefore cannot confirm that this represents a distinct' +\
                ' (paralogous) gene locus (see alignment ' + str(cur_ali_num -1) + ');'

                # Remove sequence object from list (because not necessarily
                # representative of a paralogous locus distinct from higher-ranking
                # hits).
                acc_for_objs_to_remove.append(seq_obj.id)

        # If significant overlap, get percent identity to decide whether the
        # sequence is representative of a different locus from the current "top
        # hit".
        else:
            # Determine information about similarity between the sequences in the
            # alignment (would be more efficient if this was done only after it is
            # determined that they do overlap significantly).
            info_dict = get_sequence_comparison_dict(ali_plus_top_seq_plus_1,
                    score_df)

            # Contents of dictionary:
            #info_dict['aligned_residues_count'] = aligned_residues_count
            #info_dict['identical_residues_count'] = identical_residues_count
            #info_dict['percent_identity'] = percent_identity 
            #info_dict['similar_not_ident_residues_count'] = similar_not_ident_residues_count
            #info_dict['percent_similarity'] = percent_similarity
            #info_dict['percent_simlarity_plus_identity'] = percent_similarity_plus_identity
            #info_dict['num_continuous_identical_spans'] = num_continuous_identical_spans
            #info_dict['length_of_longest_identical_span'] = length_of_longest_identical_span
            #info_dict['gene_model_redundancy_index'] = gene_model_redundancy_index

            # Record relevant info to output file with all comparisons
            # (separate from main analysis).
            outputfile = get_all_comparison_output_filepath(outdir)
            with open(outputfile, 'a') as o:
                o.write(','.join([top_hit_seq_obj.id,
                                  acc,
                                  str(info_dict['aligned_residues_count']),
                                  str(info_dict['identical_residues_count']),
                                  str(info_dict['percent_identity']),
                                  str(info_dict['similar_not_ident_residues_count']),
                                  str(info_dict['percent_similarity']),
                                  str(info_dict['num_continuous_identical_spans']),
                                  str(info_dict['length_of_longest_identical_span']),
                                  str(info_dict['gene_model_redundancy_index']),
                                  os.path.basename(ali_plus_top_seq_plus_1),
                                  ali_plus_top_seq_plus_1,
                                  '\n']
                                  ))

            # Decide whether criterion for classifying gene models/predictions
            # as redundant is met.
            criterion_met = False
            metric_value = None
            # Determine value of relevant metric.
            assert metric_name in ['percent_identity',\
                    'gene_model_redundancy_index'], """Invalid metric name
                    provided."""
            if metric_name == 'percent_identity':
                metric_value = info_dict['percent_identity']
            elif metric_name == 'gene_model_redundancy_index':
                metric_value = info_dict['gene_model_redundancy_index']
            # Determine whether criterion met.
            if metric_value >= float(metric_value_minimum):
                criterion_met = True

            # Record percent comparison metric to current "top hit".
            comparison_string1 =\
            metric_name.replace('_', ' ') + ' of ' + str(round(metric_value, 1)) + ' in alignment with sequence ' + top_hit_seq_obj.id + ' '\
            '(longest identical span ' + str(info_dict['length_of_longest_identical_span']) + ' residues)' + ' '
            # Add comparison string to full string.
            redundant_gene_model_dict[acc][2] = redundant_gene_model_dict[acc][2] + comparison_string1
    
            # If criterion met, then classify the sequence as redundant with the top
            # hit.
            if criterion_met:
                # classify as a redundant to top hit.
                redundant_gene_model_dict[acc][2] =\
                redundant_gene_model_dict[acc][2] + 'therefore redundant ' +\
                'with it '

            else:
                # Keep sequence object in list for next recursion.
                # classify as a redundant to top hit.
                redundant_gene_model_dict[acc][2] =\
                redundant_gene_model_dict[acc][2] + 'therefore not redundant ' +\
                'with it '

            # Note alignment file name.
            redundant_gene_model_dict[acc][2] =\
            redundant_gene_model_dict[acc][2] + '(see alignment %s); ' % str(cur_ali_num -1)

            ## Check whether there is overlap with existing gene models in the
            ## corresponding GFF file (if tblastn hit).
            ## Get overlapping entries (genes) from annotation file if available.
            #top_hit_search_program = extra_info_dict[top_hit_seq_obj.id][4]
            #hit_search_program = extra_info_dict[acc][4]
            #overlapping_genes = []
            #if hit_search_program.startswith('tblastn'): 
            #    # Determine whether there is a relevant annotation file to
            #    # check in.
            #    annotation_file = None
            #    sequence_filename = extra_info_dict[acc][8]
            #    with open(settings.db_info_csv) as infh:
            #        df = pd.read_csv(infh)
            #        for index, row in df.iterrows():
            #            if row['Filename'] == sequence_filename:
            #                annotation_file = row['Annotations file']

            #    # If there is an annotation file available, determine whether
            #    # the hit is at the same locus in the nucleotide sequence as
            #    # the peptide sequence hit.
            #    if annotation_file is not None:
            #        annotation_file_path = os.path.join(settings.dbdirpath,
            #                annotation_file)

            #        # Call function for getting IDs of overlapping genes.
            #        overlapping_genes =\
            #        get_overlapping_genes_from_gff3(annotation_file_path,
            #                extra_info_dict[acc][7], extra_info_dict[acc][6])

            ## Add info to info text if the top hit that the current hit is
            ## being compared to is a peptide sequence.
            #at_same_locus_as_top_hit = False
            #if not top_hit_search_program.startswith('tblastn'):
            #    for i in overlapping_genes:
            #        top_hit_acc = top_hit_seq_obj.id
            #        # ***Note: This may not work for all gene ID/accession
            #        # nomenclature schemes.
            #        if i.startswith(top_hit_acc.rsplit('.', 1)[0]):
            #            # Add note to output text string.
            #            redundant_gene_model_dict[acc][2] =\
            #            redundant_gene_model_dict[acc][2] +\
            #            'Overlaps with protein %s at gene locus %s, therefore redundant with it; ' % (top_hit_acc, i)

            #            at_same_locus_as_top_hit = True

            #            # Break the loop.
            #            break

            # Make the decision whether the current hit is redundant with the
            # top hit or not.
            #if at_same_locus_as_top_hit or criterion_met:
            if criterion_met:
                # Add accession to list for current "top hit" sequence.
                if acc not in redundant_gene_model_dict[top_hit_seq_obj.id][3]:
                    redundant_gene_model_dict[top_hit_seq_obj.id][3] =\
                            redundant_gene_model_dict[top_hit_seq_obj.id][3] + [acc]

                # Remove sequence object from list.
                acc_for_objs_to_remove.append(seq_obj.id)

                # Remove sequence object from list (this means that it
                # won't be considered a unique paralogous hit).
                acc_for_objs_to_remove.append(seq_obj.id)

    # Remove objects from the list.
    for acc in acc_for_objs_to_remove:
        for i, o in enumerate(ranked_seq_objs):
            if o.id == acc:
                #print('\tremoving: ' + ranked_seq_objs[i].id)
                del ranked_seq_objs[i]
                break


    # Check that a not was made in the spreadsheet for each sequence removed. 
    #...?

    ranked_seq_accs = []
    for i in ranked_seq_objs:
        ranked_seq_accs.append(i.id)
    
    # Iterate/recurse over the remaining sequences that were not redundant and
    # treat each in turn as a top hit, and consider remaining sequences as either
    # redundant with it, or potential paralogues.
    if len(ranked_seq_objs) > 0:
        return find_redun_model_recursively(outdir, cur_ali_num, inali, ranked_seq_objs,\
                redundant_gene_model_dict, metric_name, metric_value_minimum,\
                minimum_percent_similarity, extra_info_dict, overlap_required)
    else:
        return redundant_gene_model_dict


def count_paralogues(alignment, fasta, fwdeval, identity):
    """Manage input and output details.
    """

    # Check that there are no redundant sequence names or identical sequences in
    # the fasta file.
    extra_text = re.compile(r'\.copy$')
    name_list = []
    seq_list = []
    for record in AlignIO.read(alignment, 'nexus'):
        #print(record.id)
        name = extra_text.sub('', record.id)
        name_list.append(name)
        seq = str(record.seq)
        seq_list.append(seq)

    assert len(name_list) == len(set(name_list)), """Error: Redundant sequence
    names in input alignment."""

    if len(seq_list) != len(set(seq_list)):
        print('\nWarning: Identical sequences in input fasta records.\n')

    # Rank sequences in the fasta file according to similarity to the HMM.
    # use rank_seqs_by_hmm()
    ranked_seq_objs = rank_seqs_by_hmm(alignment, fasta)


    # Construct a dict of ids for true paralogues and ids for sequences that are
    # redundant with each.
    redundant_gene_model_dict = {}
    # Define an output directory for the results (alignments and summary csv file).
    timestamp = time.strftime("%Y_%m_%d_%H_%M_%S")
    outdir = os.path.join(fasta + '__' + os.path.basename(alignment) + '_paralogue_count_' + timestamp)
    os.mkdir(outdir)
    # Add rank info to dict.
    num = 0
    for seq in ranked_seq_objs:
        num += 1
        redundant_gene_model_dict[seq.id] = [num, seq, '', []]
    # Add sequence comparison info to dict.
    #print('\ninitial redundant_gene_model_dict: ' + str(redundant_gene_model_dict))
    full_redundant_gene_model_dict = find_redun_model_recursively(outdir, 1, alignment,\
            ranked_seq_objs, redundant_gene_model_dict, metric_name, metric_value_minimum)

    #print('full_redundant_gene_model_dict: ' + str(full_redundant_gene_model_dict))


    # TEMPORARY:
    # Write assignments in readable format.
    #outfp = args.fasta + '_paralogue_count.csv'
    outfp = os.path.join(outdir, '0_summary.csv')
    with open(outfp, 'w') as o:
        # Write output csv header line.
        o.write(','.join(['Sequence rank against HMM', 
                          'Sequence accession', 
                          'Comparison with other positive hits in the same genome', 
                          'Accessions for sequences that are redundant with the sequence',
                          '\n'])) 
        # Write content of redundancy dict.
        key_num = 0
        for key in sorted(full_redundant_gene_model_dict.keys(),\
                key=lambda x: full_redundant_gene_model_dict[x][0]):
            key_num += 1
            o.write(','.join([str(full_redundant_gene_model_dict[key][0]),
                              key,
                              '\"' + str(full_redundant_gene_model_dict[key][2]) + '\"',
                              '\"' + str(full_redundant_gene_model_dict[key][3]) + '\"',
                              '\n'])
                              )


def get_seq_obj_x(acc_list, fastadir):
    """Take an accession and find the corresponding sequence in fasta file(s)
    in a given directory.
    """
    # Compile a list of sequence object.
    seq_objs = []
    seq_obj_accs = []
    acc_list2 = acc_list
    for f in glob.glob(os.path.join(fastadir, '*.fa')):
        #print('\t\t' + f)
        for seq in SeqIO.parse(f, 'fasta'):
            #print('\t\t\t' + seq.id)
            for acc in acc_list2:
                if seq.id == acc:
                    #print('\t\t\tmatch')
                    seq_objs.append(seq)
                    seq_obj_accs.append(seq.id)
                    break

    #print(acc_list)
    #print(seq_obj_accs)

    # Check that there are no redundant sequences found.
    assert len(seq_obj_accs) <= len(acc_list), """Redundant sequences in input
    sequence file(s)."""

    # Check that there are no missing sequences.
    assert len(seq_obj_accs) == len(acc_list), """Missing sequences in input
    sequence file(s)."""

    # Re-order the list of sequence objects to match the order of corresponding
    # accessions in the input list.
    ranked_seq_obj_list = []
    for acc in acc_list:
        for seq in seq_objs:
            if seq.id == acc:
                ranked_seq_obj_list.append(seq)

    # Return list.
    return ranked_seq_obj_list


def get_seq_obj_x2(acc_list):
    """Take a list of accessions and find the corresponding sequence in fasta file(s)
    in a given directory.
    """
    # Compile a list of sequence object.
    #seq_objs = []
    seq_obj_accs = []
    acc_list2 = acc_list
    ranked_seq_obj_list = []
    for x in acc_list2:
        seq = module_amoebae.get_seqs_from_fasta_db(x[0], [x[1]])[0]
        ranked_seq_obj_list.append(seq)
        acc = seq.id
        seq_obj_accs.append(acc)
    #for f in glob.glob(os.path.join(fastadir, '*.fa')):
    #    #print('\t\t' + f)
    #    for seq in SeqIO.parse(f, 'fasta'):
    #        #print('\t\t\t' + seq.id)
    #        for acc in acc_list2:
    #            if seq.id == acc:
    #                #print('\t\t\tmatch')
    #                seq_objs.append(seq)
    #                seq_obj_accs.append(seq.id)
    #                break

    #print(acc_list)
    #print(seq_obj_accs)

    # Check that there are no redundant sequences found.
    assert len(seq_obj_accs) <= len(acc_list), """Redundant sequences in input
    sequence file(s)."""

    # Check that there are no missing sequences.
    assert len(seq_obj_accs) == len(acc_list), """Missing sequences in input
    sequence file(s)."""

    ## Re-order the list of sequence objects to match the order of corresponding
    ## accessions in the input list.
    #ranked_seq_obj_list = []
    #for acc in acc_list:
    #    for seq in seq_objs:
    #        if seq.id == acc:
    #            ranked_seq_obj_list.append(seq)

    # Return list.
    return ranked_seq_obj_list


def get_cur_ali_num(alignmentdir):
    """Gets number of alignment with largest number in its filename.
    """
    file_paths = glob.glob(os.path.join(alignmentdir, '*.nex'))
    if len(file_paths) == 0:
        return 1
    else:
        return max([int(os.path.basename(x).split('.')[0]) + 1 for x in\
            file_paths])


def count_paralogues2(csv_file, alignmentdir, fastadir, fwdeval, metric_name,
        metric_value_minimum):
    """Manage input and output details.
    fwdeval not used?
    """
    # Define timestamp.
    timestamp = time.strftime("%Y_%m_%d_%H_%M_%S")

    # Define an output directory for the results (alignments and summary csv file).
    outdir = os.path.join(csv_file  + '_paralogue_count_' + timestamp)
    os.mkdir(outdir)


    ## Check that there are no redundant sequence names or identical sequences in
    ## the alignment file.
    #extra_text = re.compile(r'\.copy$')
    #name_list = []
    #seq_list = []
    #for record in AlignIO.read(alignment, 'nexus'):
    #    #print(record.id)
    #    name = extra_text.sub('', record.id)
    #    name_list.append(name)
    #    seq = str(record.seq)
    #    seq_list.append(seq)

    #assert len(name_list) == len(set(name_list)), """Error: Redundant sequence
    #names in input alignment."""

    #if len(seq_list) != len(set(seq_list)):
    #    print('\nWarning: Identical sequences in input fasta records.\n')


    # Use pandas to parse the input csv (spreadsheet) file.
    df = pd.read_csv(csv_file)

    # Define a dictionary of relevant column labels (so I can change the title
    # that it looks for without changing the name that that tile is referred to
    # with). The input spreadsheet needs to have columns with these labels
    # (values).
    column_header_dict = {'query title': 'Protein name',
                          'alignment name': 'Alignment name',
                          'taxon name': 'Subject organism',
                          'accession': 'Positive hit accession'
                          }
    #'rank': 'Hit rank',

    # New column labels (for new columns to be appended).
    new_column_label_list = ['Sequence rank against HMM', 
                             'Sequence accession', 
                             'Comparison with other positive hits in the same genome', 
                             'Accessions for sequences that are redundant with the sequence',
                             'Represents an identifiably unique paralogue']

    # Initiate new dataframe with columns to be appended/joined to existing
    # dataframe.
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list)

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)

    # Parse relevant info from spreadsheet, and construct additional columns to
    # append.
    rank = 0
    alignment_path = None
    acc_list = []
    prev_acc_list = []
    prev_query_title = None
    prev_taxon_name = None
    num_rows = df.shape[0]
    for index, row in df.iterrows():
        rank += 1
        
        # Determine query title, taxon name, and hit accession for result in row.
        query_title = row[column_header_dict['query title']]
        taxon_name = row[column_header_dict['taxon name']]
        acc = row[column_header_dict['accession']]

        # Restart rank counting as appropriate.
        if query_title != prev_query_title:
            rank = 1
        if taxon_name != prev_taxon_name:
            rank = 1

        # Define alignment to use, and add accession to list or restart
        # accession list as appropriate.
        if rank == 1:
            alignment_path = os.path.join(alignmentdir,\
                    row[column_header_dict['alignment name']])
            prev_acc_list = acc_list
            acc_list = [acc]

        else:
            acc_list.append(acc)

        print(' '.join(['Row', str(index), query_title, taxon_name]))


        # Process the accession (sequence) listed in the row, and add info to
        # output dataframe.

        # When appropriate, take a list of accessions, get corresponding
        # sequence objects, process these, and add information to appropriate
        # rows/fields in the dataframe.
        if index != 0 and rank == 1 or index == num_rows -1:
            # Process preceding set of hits for (previous) query title.

            if index == num_rows -1:
                # Reached the end of the dataframe/spreadsheet, so have to set
                # current query title, etc, to prev.
                prev_query_title = query_title
                prev_alignment_path = alignment_path
                prev_taxon_name = taxon_name
                prev_acc_list = acc_list

            # Get list of sequence objects corresponding to list of accessions.
            ranked_seq_objs = get_seq_obj_x(prev_acc_list, fastadir)

            # Construct a dict of ids for true paralogues and ids for sequences that are
            # redundant with each.
            redundant_gene_model_dict = {}

            # Add rank info to dict.
            num = 0
            for seq in ranked_seq_objs:
                num += 1
                redundant_gene_model_dict[seq.id] = [num, seq, '', []]

            # get current alignment number so that it can pick up where it
            # started.
            cur_ali_num = get_cur_ali_num(outdir) 

            # Add sequence comparison info to dict.
            full_redundant_gene_model_dict =\
            find_redun_model_recursively(outdir, cur_ali_num, alignment_path,\
                    ranked_seq_objs, redundant_gene_model_dict, metric_name, metric_value_minimum)

            #print('full_redundant_gene_model_dict: ' + str(full_redundant_gene_model_dict))


            # Add info from dict to dataframe.
            key_num = 0
            for key in sorted(full_redundant_gene_model_dict.keys(),\
                    key=lambda x: full_redundant_gene_model_dict[x][0]):
                key_num += 1
                # Get index for row that info should be added to.
                index_in_dataframe = index - len(full_redundant_gene_model_dict.keys()) + key_num -1
                if index == num_rows -1:
                    index_in_dataframe = index_in_dataframe + 1

                # Determine from dict whether the hit represents the first hit
                # for an identifiably unique paralogous gene locus.
                unique = 'No'
                if 'represents a unique (paralogous) gene locus' in str(full_redundant_gene_model_dict[key][2]):
                    unique = 'Yes'

                # Define dataframe values using info from dict.
                df.at[index_in_dataframe, 'Sequence rank against HMM'] = str(full_redundant_gene_model_dict[key][0])
                df.at[index_in_dataframe, 'Sequence accession'] = key
                df.at[index_in_dataframe, 'Comparison with other positive hits in the same genome'] = str(full_redundant_gene_model_dict[key][2])
                df.at[index_in_dataframe, 'Accessions for sequences that are redundant with the sequence'] = str(full_redundant_gene_model_dict[key][3])
                df.at[index_in_dataframe, 'Represents an identifiably unique paralogue'] = unique

        # Define previous query title and alignment path for next run through loop.
        prev_query_title = query_title
        prev_alignment_path = alignment_path
        prev_taxon_name = taxon_name

    # Write modified dataframe to output file path.
    output_fp = os.path.join(outdir, os.path.basename(csv_file) + '_out.csv')
    df.to_csv(output_fp, index=False)
    print('\nResults written to file:\n' + output_fp)


    # Visualize all comparisons to evaluate whether metric adequately
    # distinguished between sequences.

    # Get filepath for relevant output file. 
    outputfile = get_all_comparison_output_filepath(outdir)

    # Make histogram.
    # ...


def get_overlapping_genes_from_gff3(sql_database, seq_id, cluster_range):
    """Take a path to a gffutils SQL database (based on a gff3 file), the
    sequence id for the relevant sequnce referred to in the gff file, and a
    range (coordinates in the form of a list object with two items). And,
    return IDs for any genes that fall within that range on the sequence. This
    is for determining whether tblastn hits correspond to already annotated
    protein sequences.
    """
    # Parse database.
    db = gffutils.FeatureDB(sql_database, keep_order=True)

    # Find any genes that fall within the range (or overlap at
    # all).
    overlapping_entries = list(db.region(region=(seq_id, cluster_range[0], cluster_range[1]),
                                    completely_within=False, featuretype='gene'
                                    ))

    # Get gene IDs.
    overlap_ids = []
    for gene in overlapping_entries:
        overlap_ids.append(gene.id)

    # Return list of IDs.
    return overlap_ids


def count_paralogues3(csv_file, metric_name, metric_value_minimum,
        minimum_percent_similarity, timestamp, overlap_required,
        minimum_length_of_query_to_be_distinct_paralogue,
        minimum_percent_length_of_query_to_be_distinct_paralogue,
        remove_tblastn_hits_at_annotated_loci, outfp=None):
    """Manage input and output details.
    """
    # Define unique identifier string for this analysis.
    paralogue_count_id = 'paralogue_count_' + timestamp

    # Define directory to acces with alignments.
    alignmentdir = settings.querydirpath

    # Define an output directory for the results (alignments and summary csv file).
    outdir = os.path.join(csv_file.rsplit('.', 1)[0]  + '_' + paralogue_count_id)
    os.mkdir(outdir)

    # Make a file for recording all sequence comparisons made.
    outputfile = get_all_comparison_output_filepath(outdir)
    with open(outputfile, 'a') as o:
        # Write header line.
        o.write(','.join(['ID for sequence 1',
                          'ID for sequence 2',
                          'aligning residues position count',
                          'identical position count',
                          'percent_identity',
                          'similar (not identical) position count',
                          'percent_similarity',
                          'num_continuous_identical_spans',
                          'length_of_longest_identical_span', 'gene_model_redundancy_index',
                          'alignment file name'
                          '\n']
                          ))


    ## Check that there are no redundant sequence names or identical sequences in
    ## the alignment file.
    #extra_text = re.compile(r'\.copy$')
    #name_list = []
    #seq_list = []
    #for record in AlignIO.read(alignment, 'nexus'):
    #    #print(record.id)
    #    name = extra_text.sub('', record.id)
    #    name_list.append(name)
    #    seq = str(record.seq)
    #    seq_list.append(seq)

    #assert len(name_list) == len(set(name_list)), """Error: Redundant sequence
    #names in input alignment."""

    #if len(seq_list) != len(set(seq_list)):
    #    print('\nWarning: Identical sequences in input fasta records.\n')


    # Use pandas to parse the input csv (spreadsheet) file.
    print('Reading input csv file into a pandas dataframe.')
    df = pd.read_csv(csv_file, low_memory=False)

    # Define a dictionary of relevant column labels (so I can change the title
    # that it looks for without changing the name that that tile is referred to
    # with). The input spreadsheet needs to have columns with these labels
    # (values).
    #column_header_dict = {'query title': 'Protein name',
    #                      'alignment name': 'Alignment name',
    #                      'taxon name': 'Subject organism',
    #                      'accession': 'Positive hit accession'
    #                      }

    # Updated for AMOEBAE.
    column_header_dict = {'query title': 'Query title',
                          'alignment name': 'Alignment for sequence comparison',
                          'taxon name': 'Subject database species (if applicable)',
                          'database filename': 'Subject database file',
                          'positive or not': 'Collective interpretation of reverse search results',
                          'accession': 'Forward hit accession',
                          'hit rank': 'Forward hit rank',
                          'subseq coordinates': 'Forward hit coordinates of subsequence(s) that align(s) to query',
                          'program': 'Forward search method',
                          'evalue': 'Forward hit E-value (top HSP)',
                          'hit length': 'Forward hit length',
                          'percent length': 'Forward hit length as a percentage of query length',
                          'seq': 'Forward hit sequence',
                          'subseq': 'Forward hit subsequence(s) that align(s) to query',
                          'seq descr': 'Forward hit description',
                          'subseq descr': 'Forward hit description of subsequence(s) that align(s) to query'
                          }
    #'rank': 'Hit rank',

    # New column labels (for new columns to be appended).
    #new_column_label_list = ['Sequence rank against HMM', 
    #                         'Sequence accession', 
    #                         'Comparison with other positive hits in the same genome', 
    #                         'Accessions for sequences that are redundant with the sequence',
    #                         'Represents an identifiably unique paralogue']
    new_column_label_list = ['Does not have same ID or locus as another hit',
                             'Meets length criteria',
                             'Comparison with other positive hits in the same genome', 
                             'Accessions for sequences that are redundant with the sequence',
                             'Represents an identifiably unique paralogue',
                             'Paralogue name']

    # Check that a column containing the name of the alignment to use is
    # present in the input file.
    assert column_header_dict['alignment name'] in list(df.columns), """Could
    not identify column containing name(s) of alignment(s) to use in input csv
    file %s""" % csv_file

    # Initiate new dataframe with columns to be appended/joined to existing
    # dataframe.
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list)

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)

    # Get sets of indexes for which summaries should be made, organized by
    # query title and taxon combinations as well as corresponding relevant information.
    index_dict = {}
    for index, row in df.iterrows():
        if row[column_header_dict['positive or not']] == '-':
            # Ignore negative hits.
            pass
        else:
            # Get info from row.
            query_title = row[column_header_dict['query title']]
            taxon_name = row[column_header_dict['taxon name']]
            #acc = row[column_header_dict['accession']] + '(fwdhit' +\
            #        str(row[column_header_dict['hit rank']]) + ')'
            alignment_path = os.path.join(alignmentdir,\
                            row[column_header_dict['alignment name']])
            evalue = row[column_header_dict['evalue']]
            percent_len = row[column_header_dict['percent length']]
            length = row[column_header_dict['hit length']]

            # Get hit range.
            hit_range =\
            get_hit_range_from_hsp_ranges(row[column_header_dict['subseq coordinates']])

            # Get sequence using different sources for prot vs. nucl data.
            program = row[column_header_dict['program']]
            description = None
            seq = None
            acc = None
            actual_acc = row[column_header_dict['accession']]
            if program.startswith('tblastn'):
                # Append hit range to accession for tblastn, to distinguish
                # between different clusters of HSPs in the same
                # chromosome/scaffold.
                acc = row[column_header_dict['accession']] +\
                    '[' + str(hit_range[0]) + '..' + str(hit_range[1]) + ']'
                description = row[column_header_dict['subseq descr']]
                seq = row[column_header_dict['subseq']]
            else:
                acc = row[column_header_dict['accession']]
                description = row[column_header_dict['seq descr']]
                seq = row[column_header_dict['seq']]
            seq_obj = get_seq_obj_from_srch_res_csv_info(acc, description, seq)

            sequence_filename = row[column_header_dict['database filename']]

            # Integrate information into dict.
            combo = (query_title, taxon_name)
            # This info tuple is used later as well...
            info_tuple = (index,
                          acc,
                          seq_obj,
                          evalue, 
                          program,
                          percent_len,
                          hit_range,
                          actual_acc,
                          sequence_filename,
                          length
                          )
            if combo not in index_dict.keys():
                # Initiate a new key in dict.
                index_dict[combo] = [alignment_path, [info_tuple]]
            else:
                ## Ignore sequences with redundant accessions.
                #if acc in [x[1] for x in index_dict[combo][1]]: 
                #    df.at[index,\
                #        'Comparison with other positive hits in the same genome']\
                #        = '(Accession (or coordinates) redundant with previous hit(s))'

                ## Ignore hits that are below the minimum percent length.
                ## (Should this be done at a later stage so that it can be
                ## determined whether it is redundant with a protein model due
                ## to location in the genome?)
                #elif int(percent_len) < settings.minimum_percent_length_of_query_to_be_distinct_paralogue:
                #    df.at[index, 'Comparison with other positive hits in the same genome'] =\
                #    '(Hit length is less than %s percent of query length, therefore too short to be a potential distinct paralogue)'\
                #            % str(settings.minimum_percent_length_of_query_to_be_distinct_paralogue)
                ## Ignore hits that are below the absolute minimum length.
                #...

                #else:
                index_dict[combo] = [alignment_path, index_dict[combo][1]\
                        + [info_tuple]]


    # Make the values in the 'Does not have same ID or locus as another
    # hit' and 'Meets length criteria' columns the same as for the
    # 'positive or not' column to start with.
    for index, row in df.iterrows():
        if row[column_header_dict['positive or not']] == '-':
            # Ignore negative hits.
            pass
        else:
            # Change values in the appropriate columns to '+'.
            df.at[index, 'Does not have same ID or locus as another hit'] = '+'
            df.at[index, 'Meets length criteria'] = '+'


    # Iterate over relevant information and feed it into paralogue_counter
    # module for processing, and then put the relevant results in the right place
    # in the dataframe.
    total_combo_count = len(index_dict.keys())
    combo_num = 0
    for combo in index_dict.keys():
        combo_num += 1
        # Print progress.
        print('\n')
        print("""\nFinding redundant sequences for query title %s and taxon %s
            (%s of %s sets, %s percent complete)"""\
            % (combo[0], combo[1], str(combo_num), str(total_combo_count),
                str(round(combo_num*100/total_combo_count))))

        # Sort hit info tuples into those for protein hits and those for
        # nucleotide hits.
        all_tuples = index_dict[combo][1]
        prot_tuples = []
        nucl_tuples = []
        for i in all_tuples:
            if i[4].startswith('tblastn'):
                nucl_tuples.append(i)
            else:
                prot_tuples.append(i)

        # Sort the tuple lists by ascending E-value.
        prot_tuples = sorted(prot_tuples, key=lambda x: x[3])
        nucl_tuples = sorted(nucl_tuples, key=lambda x: x[3])

        # Remove protein hits from the tuples if they have identical accessions
        # to higher-ranking hits.
        reduced_prot_tuples = []
        hit_indexes_to_remove = []
        for num, i in enumerate(prot_tuples):
            if num + 1 < len(prot_tuples):
                for j in prot_tuples[num + 1:]:
                    # Check whether they have the same filename and sequence ID.
                    if j[8] == i[8] and j[7] == i[7]:
                        # Check that the E-values are right.
                        if i[3] <= j[3]:
                            # Mark the one with the lower or equivalent E-value for removal.
                            hit_indexes_to_remove.append(j[0])

                            # Change value in column for ID or locus redundancy
                            # to '-'.
                            df.at[j[0],\
                                'Does not have same ID or locus as another hit']\
                                = '-'

                            # Make a note in the spreadsheet, explaining why
                            # this hit was excluded.
                            df.at[j[0],\
                                'Comparison with other positive hits in the same genome']\
                                = '(This hit has the same ID as protein hit with ID: %s)' % i[1]
        # Remove the identified redundant tuples. 
        if len(hit_indexes_to_remove) > 0:
            for i in prot_tuples:
                if i[0] not in hit_indexes_to_remove:
                    reduced_prot_tuples.append(i)
        else:
            reduced_prot_tuples = prot_tuples

        # Remove nucleotide hits from the tuples, if they represent a gene
        # locus that encodes one of the protein hits, or optionally any
        # annotated gene.
        reduced_nucl_tuples = []
        nucl_tuples_to_remove = []
        if len(prot_tuples) == 0:
            # No protein hits to compare nucleotide hits to, so keep all
            # nucleotide hits.
            reduced_nucl_tuples = nucl_tuples
        else:
            # Look for nucleotide hits that are redundant with protein hits.
            for x in nucl_tuples:
                # Check whether there is overlap with existing gene models in the
                # corresponding GFF file (if tblastn hit).
                # Get overlapping entries (genes) from annotation file if available.
                overlapping_genes = []
                # Determine whether there is a relevant annotation file to
                # check in.
                annotation_file = None
                sequence_filename = x[8]
                with open(settings.db_info_csv) as infh:
                    dfx = pd.read_csv(infh)
                    for index, row in dfx.iterrows():
                        if row['Filename'] == sequence_filename:
                            annotation_file = row['Annotations file']

                            # Change the value back to None if it is NaN in the
                            # DataFrame.
                            if type(annotation_file).__name__ == 'float':
                                annotation_file = None


                # If there is an annotation file available, determine whether
                # the hit is at the same locus in the nucleotide sequence as
                # the peptide sequence hit.
                if annotation_file is not None:
                    # Check that the type is string.
                    assert type(annotation_file).__name__ == 'str'

                    # Join the paths.
                    annotation_file_path = os.path.join(settings.dbdirpath,
                            annotation_file)

                    # Call function for getting IDs of overlapping genes.
                    overlapping_genes =\
                    get_overlapping_genes_from_gff3(annotation_file_path,
                            x[7], x[6])

                # If there are overlapping gene loci with the tblastn hit...
                if len(overlapping_genes) > 0:
                    if remove_tblastn_hits_at_annotated_loci:
                        # If there are any overlapping annotated genes, then remove the
                        # nucleotide tuple and make a note. This only makes sense
                        # if the relevant protein sequences are also being searched
                        # with a similar query.

                        # Change value in column for ID or locus redundancy
                        # to '-'.
                        df.at[x[0],\
                            'Does not have same ID or locus as another hit']\
                            = '-'

                        # Make a note in the dataframe/spreadsheet.
                        df.at[x[0],\
                            'Comparison with other positive hits in the same genome']\
                            = '(This hit overlaps with the following previously annotated loci: %s)'\
                                    % ', '.join(overlapping_genes)

                        # Mark tuple for removal.
                        nucl_tuples_to_remove.append(x)

                    else:
                        # Loop over protein hits and see if any are redundant with the
                        # nucleotide hit.
                        for y in reduced_prot_tuples:
                            # Add info to info text to the spreadsheet if the tblastn hit
                            # is redundant with a protein hit.
                            for i in overlapping_genes:
                                prot_hit_acc = y[7]
                                # ***Note: This may not work for all gene ID/accession
                                # nomenclature schemes.
                                if i.startswith(prot_hit_acc.rsplit('.', 1)[0]):

                                    # Change value in column for ID or locus redundancy
                                    # to '-'.
                                    df.at[x[0],\
                                        'Does not have same ID or locus as another hit']\
                                        = '-'

                                    # Make note in spreadsheet.
                                    df.at[x[0],\
                                        'Comparison with other positive hits in the same genome']\
                                        = '(This hit is at the same locus as protein hit with ID: %s)' % prot_hit_acc

                                    # Mark tuple for removal, because it is redundant
                                    # with a peptide hit.
                                    nucl_tuples_to_remove.append(x)

                                    # Break the loop.
                                    break

        # Generate a reduced list of nucleotide info tuples, if there are any
        # to remove from the original list.
        if len(nucl_tuples_to_remove) > 0:
            nucl_hit_accs_to_remove = [x[1] for x in nucl_tuples_to_remove]
            for i in nucl_tuples:
                if i[1] not in nucl_hit_accs_to_remove:
                    reduced_nucl_tuples.append(i)
        else:
            reduced_nucl_tuples = nucl_tuples

        # Check for redundancy among the tblastn hits in their locations, and
        # remove all but the best hit for each locus.
        # Since the tuples are sorted by E-value, for each tuple in the list,
        # any tuples further down the list that have the same filename and seq
        # ID and overlapping coordinates can be marked for removal.
        further_reduced_nucl_tuples = []
        more_nucl_tuples_to_remove = []
        for num, i in enumerate(reduced_nucl_tuples):
            if num + 1 < len(reduced_nucl_tuples):
                for j in reduced_nucl_tuples[num + 1:]:
                    # Check whether they have the same filename and sequence ID.
                    if j[8] == i[8] and j[7] == i[7]:
                        # Check for overlapping ranges.
                        ir = i[6]
                        jr = j[6]
                        iset = set(range(int(ir[0]), int(ir[1])))
                        jset = set(range(int(jr[0]), int(jr[1])))
                        if len(iset.intersection(jset)) > 0:
                            # Then the hits have overlapping ranges/coordinates on
                            # the same subject sequence.
                            # Check that the E-values are right.
                            if i[3] <= j[3]:
                                # Mark the one with the lower or equivalent E-value for removal.
                                more_nucl_tuples_to_remove.append(j)

                                # Change value in column for ID or locus redundancy
                                # to '-'.
                                df.at[j[0],\
                                    'Does not have same ID or locus as another hit']\
                                    = '-'

                                # Make a note in the spreadsheet, explaining why
                                # this hit was excluded.
                                df.at[j[0],\
                                    'Comparison with other positive hits in the same genome']\
                                    = '(This hit is at the same locus as nucleotide hit with ID: %s)' % i[1]

                                assert j[1] != 'Chr5[17020882..17021032]' # REMOVE THIS???
        # Remove the identified redundant tuples. 
        if len(more_nucl_tuples_to_remove) > 0:
            nucl_hit_indexes_to_remove = [x[0] for x in more_nucl_tuples_to_remove]
            for i in reduced_nucl_tuples:
                if i[0] not in nucl_hit_indexes_to_remove:
                    further_reduced_nucl_tuples.append(i)
        else:
            further_reduced_nucl_tuples = reduced_nucl_tuples


        # Ignore protein hits that are below the minimum length or percent length.
        reduced_prot_tuples2 = []
        prm = []
        for i in reduced_prot_tuples:
            length = i[9]
            percent_length = i[5]
            remove = False
            if int(length) <\
            int(minimum_length_of_query_to_be_distinct_paralogue):
                remove = True
                # Change value in column for length criteria.
                df.at[i[0],\
                    'Meets length criteria']\
                    = '-'
                # Make a note in the spreadsheet.
                df.at[i[0],'Comparison with other positive hits in the same genome'] =\
                '(Hit length is less than %s, therefore too short to be a potential distinct paralogue)'\
                    % str(minimum_length_of_query_to_be_distinct_paralogue)
            if int(percent_length) <\
            int(minimum_percent_length_of_query_to_be_distinct_paralogue):
                remove = True
                # Change value in column for length criteria.
                df.at[i[0],\
                    'Meets length criteria']\
                    = '-'
                # Make a note in the spreadsheet.
                df.at[i[0],'Comparison with other positive hits in the same genome'] =\
                '(Hit length is less than %s percent of query length, therefore too short to be a potential distinct paralogue)'\
                    % str(minimum_percent_length_of_query_to_be_distinct_paralogue)
            if remove:
                prm.append(i[0])
        # Remove the identified redundant tuples. 
        if len(prm) > 0:
            for i in reduced_prot_tuples:
                if i[0] not in prm:
                    reduced_prot_tuples2.append(i)
        else:
            reduced_prot_tuples2 = reduced_prot_tuples

        # Ignore nucleotide hits that are below the minimum length or percent length.
        further_reduced_nucl_tuples2 = []
        nrm = []
        for i in further_reduced_nucl_tuples:
            length = i[9]
            percent_length = i[5]
            remove = False
            if int(length) <\
            int(minimum_length_of_query_to_be_distinct_paralogue):
                remove = True
                # Change value in column for length criteria.
                df.at[i[0],\
                    'Meets length criteria']\
                    = '-'
                # Make a note in the spreadsheet.
                df.at[i[0],'Comparison with other positive hits in the same genome'] =\
                '(Hit length is less than %s, therefore too short to be a potential distinct paralogue)'\
                    % str(minimum_length_of_query_to_be_distinct_paralogue)
            if int(percent_length) <\
            int(minimum_percent_length_of_query_to_be_distinct_paralogue):
                remove = True
                # Change value in column for length criteria.
                df.at[i[0],\
                    'Meets length criteria']\
                    = '-'
                # Make a note in the spreadsheet.
                df.at[i[0],'Comparison with other positive hits in the same genome'] =\
                '(Hit length is less than %s percent of query length, therefore too short to be a potential distinct paralogue)'\
                    % str(minimum_percent_length_of_query_to_be_distinct_paralogue)
            if remove:
                nrm.append(i[0])
        # Remove the identified redundant tuples. 
        if len(nrm) > 0:
            for i in further_reduced_nucl_tuples:
                if i[0] not in nrm:
                    further_reduced_nucl_tuples2.append(i)
        else:
            further_reduced_nucl_tuples2 = further_reduced_nucl_tuples


        # Remove nucleotide hits from the tuples if they have identical accessions
        # to higher-ranking nucleotide hits (This should be redundant, because
        # it is already checked whether they are for the same locus).
        further_reduced_nucl_tuples3 = []
        hit_indexes_to_remove = []
        for num, i in enumerate(further_reduced_nucl_tuples2):
            if num + 1 < len(further_reduced_nucl_tuples2):
                for j in further_reduced_nucl_tuples2[num + 1:]:
                    # Check whether they have the same filename and sequence ID
                    # (seq.id plus coordinates).
                    if j[8] == i[8] and j[1] == i[1]:
                        # Check that the E-values are right.
                        if i[3] <= j[3]:
                            # Mark the one with the lower or equivalent E-value for removal.
                            hit_indexes_to_remove.append(j[0])

                            # Change value in column for length criteria.
                            df.at[j[0],\
                                'Meets length criteria']\
                                = '-'

                            # Make a note in the spreadsheet, explaining why
                            # this hit was excluded.
                            df.at[j[0],\
                                'Comparison with other positive hits in the same genome']\
                                = '(This hit has the same ID as nucleotide hit with ID: %s)' % i[1]
        # Remove the identified redundant tuples. 
        if len(hit_indexes_to_remove) > 0:
            for i in further_reduced_nucl_tuples2:
                if i[0] not in hit_indexes_to_remove:
                    further_reduced_nucl_tuples3.append(i)
        else:
            further_reduced_nucl_tuples3 = further_reduced_nucl_tuples2


        # Check that no accessions listed in the info tuples are identical to
        # each other.
        all_accs = [x[1] for x in reduced_prot_tuples2 + further_reduced_nucl_tuples3]
        assert len(all_accs) == len(list(set(all_accs))), """Info for hits with
        identical accessions is present in the reduced set."""


        # Check that explanations were noted in the spreadsheet (dataframe) for
        # all the tuples that were removed from the original list.
        for i in all_tuples:
            if i[0] not in [x[0] for x in reduced_prot_tuples2\
            + further_reduced_nucl_tuples3]:
                note = df.loc[i[0],'Comparison with other positive hits in the same genome'] 
                assert note != '-', """No justification noted for not
                considering a hit with ID %s and index %s""" % (i[1], str(i[0]))

        # Compile a list of all reduced tuples.
        all_reduced_tuples = reduced_prot_tuples2 + further_reduced_nucl_tuples3

        # Compile lists of just the sequence objects.
        ranked_prot_seq_objs = [x[2] for x in reduced_prot_tuples2]
        ranked_nucl_seq_objs = [x[2] for x in further_reduced_nucl_tuples3]
        # Concatenate the two lists, putting proteins first.
        ranked_seq_objs = ranked_prot_seq_objs + ranked_nucl_seq_objs

        ## Check that there is at least one sequence object in the list.
        #assert len(ranked_seq_objs) >= 1, """There must be at least one
        #sequence object in the list."""

        if len(ranked_seq_objs) >= 1:
            # Make a dict with info for each hit to be used by the
            # find_redun_model_recursively function (info can be accessed by
            # accession).
            extra_info_dict = {}
            for i in all_tuples:
                acc = i[1]
                extra_info_dict[acc] = i

            # Construct a dict of ids for true paralogues and ids for sequences that are
            # redundant with each.
            redundant_gene_model_dict = {}

            # Add rank info to dict.
            num = 0
            for seq in ranked_seq_objs:
                num += 1
                redundant_gene_model_dict[seq.id] = [num, seq, '', []]

            # get current alignment number so that it can pick up where it
            # started.
            cur_ali_num = get_cur_ali_num(outdir) 

            # Define path to alignment to use.
            alignment_path = index_dict[combo][0]

            # Add sequence comparison info to dict.
            full_redundant_gene_model_dict =\
            find_redun_model_recursively(outdir, cur_ali_num, alignment_path,\
                    ranked_seq_objs, redundant_gene_model_dict, metric_name,\
                    metric_value_minimum, minimum_percent_similarity,
                    extra_info_dict, overlap_required)

            # Extract info from full redundant gene model dict, and add to
            # corresponding cells in dataframe.
            key_num = 0
            unique_paralogue_num = 0
            for key in sorted(full_redundant_gene_model_dict.keys(),\
                    key=lambda x: full_redundant_gene_model_dict[x][0]):
                key_num += 1

                # Get index for row that info should be added to.
                # It should be the index that has results for hit with the same
                # accession.
                index_in_dataframe = None

                # Iterate over reduced info tuples for current query title and
                # taxon name (combo).
                for t in all_reduced_tuples:
                    # If the accession in the tuple is the same as for the
                    # current dict key, then use the index from that tuple.
                    if t[1] == key:
                        index_in_dataframe = t[0]

                # Determine from dict whether the hit represents the first hit
                # for an identifiably unique paralogous gene locus.
                unique = 'No'
                if 'represents a unique (paralogous) gene locus' in str(full_redundant_gene_model_dict[key][2]):
                    unique = 'Yes'
                elif 'uncertain whether this represents a unique (paralogous) gene locus' in str(full_redundant_gene_model_dict[key][2]):
                    unique = 'Uncertain'

                # Define dataframe values using info from dict.
                if df.loc[index_in_dataframe, 'Forward search method'].startswith('tblastn'):
                    df.at[index_in_dataframe, 'Comparison with other positive hits in the same genome'] =\
                    str(full_redundant_gene_model_dict[key][1].id) + ' ' + str(full_redundant_gene_model_dict[key][2])
                else:
                    df.at[index_in_dataframe, 'Comparison with other positive hits in the same genome'] = str(full_redundant_gene_model_dict[key][2])
                #print(df.loc[index_in_dataframe, 'Comparison with other positive hits in the same genome'])
                df.at[index_in_dataframe, 'Accessions for sequences that are redundant with the sequence'] = str(full_redundant_gene_model_dict[key][3])
                df.at[index_in_dataframe, 'Represents an identifiably unique paralogue'] = unique

                # Assign paralogue name to "Paralogue name" column, if it is
                # unique.
                if unique == 'Yes':
                    unique_paralogue_num += 1
                    df.at[index_in_dataframe, 'Paralogue name'] =\
                    df.at[index_in_dataframe, column_header_dict['taxon name']].replace(' ', '_')\
                    + '__' + combo[0] + '_' + str(unique_paralogue_num)

                # Check that a note was made for the current index (row) in the
                # dataframe.
                note = df.loc[index_in_dataframe, 'Comparison with other positive hits in the same genome']
                assert note != '-', """No note made for hit with ID %s and index %s"""\
                        % (str(full_redundant_gene_model_dict[key][1].id), str(index_in_dataframe))                


        # Check that explanations were noted in the spreadsheet (dataframe) for
        # all the tuples list processed with the recursive function to check
        # for redundancy.
        for i in all_tuples:
            if i[0] in [x[0] for x in reduced_prot_tuples2\
            + further_reduced_nucl_tuples2]:
                note = df.loc[i[0],'Comparison with other positive hits in the same genome'] 
                assert note != '-', """No note made for hit with ID %s and index %s""" % (i[1], str(i[0]))


    # Change column headers in dataframe so that this function can be called
    # multiple times with the same input spreadsheet to compare results when
    # using different parameters.
    df.columns = list(df.columns[0:-len(new_column_label_list)]) + [x + ' (' + paralogue_count_id + ')' for x in new_column_label_list]

    # Write modified dataframe to output file path.
    output_fp = None
    if outfp == None:
        #output_fp = os.path.join(outdir, os.path.basename(csv_file).rsplit('.',
        #    1)[0] + '_' + paralogue_count_id + '.csv')
        output_fp = os.path.join(os.path.dirname(csv_file), os.path.basename(csv_file).rsplit('.',
            1)[0] + '_' + paralogue_count_id + '.csv')
        print(output_fp)
    else:
        output_fp = outfp
    print('Writing dataframe to csv file.')
    df.to_csv(output_fp, index=False)
    #print('\nResults written to file:\n' + output_fp)


    # Visualize all comparisons to evaluate whether metric adequately
    # distinguished between sequences.

    # Get filepath for relevant output file. 
    outputfile = get_all_comparison_output_filepath(outdir)

    # Return final output filepath for printing.
    return output_fp


def add_alignment_column(incsv, outcsv):
    """Add a column to the input csv spreadsheet which lists the alignment file
    name to use.
    """
    # Define output filepath.
    if outcsv is None:
        outcsv = incsv.rsplit('.', 1)[0] + '_with_ali_col.csv'


    # Loop over rows, and make a dictionary of alignment file names for each
    # query title.

    # Define dict for interpreting headers.
    column_header_dict = {'query title': 'Query title',
                          'alignment name': 'Alignment for sequence comparison',
                          'query filename': 'Query file',
                          'taxon name': 'Subject database species (if applicable)',
                          'database filename': 'Subject database file',
                          'positive or not': 'Collective interpretation of reverse search results',
                          'accession': 'Forward hit accession',
                          'program': 'Forward search method',
                          'seq': 'Forward hit sequence',
                          'subseq': 'Forward hit subsequence(s) that align(s) to query',
                          'seq descr': 'Forward hit description',
                          'subseq descr': 'Forward hit description of subsequence(s) that align(s) to query'
                          }
    # Use pandas to parse the input csv (spreadsheet) file.
    df = pd.read_csv(incsv)

    # Loop over rows and construct dict with query titles as keys and
    # corresponding alignment to use for sequence comparison as values.
    complete_query_title_list = []
    query_title_alignment_dict = {}
    for index, row in df.iterrows():
        #if row[column_header_dict['positive or not']] == '+':
        #    query_title = row[column_header_dict['query title']]
        #    query_file = row[column_header_dict['query filename']]
        #    complete_query_title_list.append(query_title)
        #    if query_file.endswith('.afaa') or query_file.endswith('.afna'):
        #        if query_title not in query_title_alignment_dict.keys():
        #            query_title_alignment_dict[query_title] = [query_file]
        #        else:
        #            if not query_file in query_title_alignment_dict[query_title]:
        #                query_title_alignment_dict[query_title] = query_title_alignment_dict[query_title] + [query_file]

        #    # TEMPORARY!!:
        #    #elif query_title in ['AP5_mu', 'TTRAY1', 'AP3_sigma', 'COPI_beta','ClathrinHC']:
        #    #        query_title_alignment_dict[query_title] = [query_file]
        #    #if query_title == 'COPI_alpha':
        #    #        query_title_alignment_dict[query_title] = ['COPI_Alpha_3.afaa']
        #    elif query_title in ['AP5_sigma']:
        #            query_title_alignment_dict[query_title] = [query_file]

        query_title = row[column_header_dict['query title']]
        query_file = row[column_header_dict['query filename']]
        complete_query_title_list.append(query_title)
        if query_file.endswith('.afaa') or query_file.endswith('.afna'):
            if query_title not in query_title_alignment_dict.keys():
                query_title_alignment_dict[query_title] = [query_file]
            else:
                if not query_file in query_title_alignment_dict[query_title]:
                    query_title_alignment_dict[query_title] = query_title_alignment_dict[query_title] + [query_file]

    
    # Check that an alignment file name was assigned to each query title.
    titles_without_alignments = []
    for query_title in list(set(complete_query_title_list)):
        if query_title not in query_title_alignment_dict.keys():
            titles_without_alignments.append(query_title)
    assert len(titles_without_alignments) == 0, """No alignment file names
    could be found to use for the following query titles:
        %s""" % (str(titles_without_alignments))

    # Check that one and only one alignment file name has been assigned to each
    # query title.
    for query_title in query_title_alignment_dict.keys():
        assert len(query_title_alignment_dict[query_title]) > 0, """Could not
        identify an alignment file to use for query title: %s""" % query_title

        assert len(query_title_alignment_dict[query_title]) < 2, """Identified
        more than one potential alignment file to use for query title %s:
            %s"""\
        % (query_title, str(query_title_alignment_dict[query_title]))
    

    # Use the dict to add a column to the dataframe.

    # New column labels (for new columns to be appended).
    new_column_label_list = [column_header_dict['alignment name']]

    # Initiate new dataframe with columns to be appended/joined to existing
    # dataframe.
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list)

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)


    # Update the new column with alignment names.
    for index, row in df.iterrows():
        if row[column_header_dict['positive or not']] == '+':
            query_title = row[column_header_dict['query title']]
            alignment_to_use = query_title_alignment_dict[query_title][0]
            df.at[index, column_header_dict['alignment name']] = alignment_to_use

    # Write updated dataframe to output csv file.
    df.to_csv(outcsv, index=False)

    return outcsv


def add_model_column(incsv, outcsv):
    """Add a column to the input csv spreadsheet which lists the alignment file
    name to use.

    ***Very similar to the add_alignment_to_column function defined above.
    ***Not yet customized for purpose!
    """
    # Define output filepath.
    if outcsv is None:
        outcsv = incsv.rsplit('.', 1)[0] + '_with_ali_col.csv'


    # Loop over rows, and make a dictionary of alignment file names for each
    # query title.

    # Define dict for interpreting headers.
    column_header_dict = {'query title': 'Query title',
                          'alignment name': 'Alignment for sequence comparison',
                          'query filename': 'Query file',
                          'taxon name': 'Subject database species (if applicable)',
                          'database filename': 'Subject database file',
                          'positive or not': 'Collective interpretation of reverse search results',
                          'accession': 'Forward hit accession',
                          'program': 'Forward search method',
                          'seq': 'Forward hit sequence',
                          'subseq': 'Forward hit subsequence(s) that align(s) to query',
                          'seq descr': 'Forward hit description',
                          'subseq descr': 'Forward hit description of subsequence(s) that align(s) to query'
                          }
    # Use pandas to parse the input csv (spreadsheet) file.
    df = pd.read_csv(incsv)

    # Loop over rows and construct dict.
    complete_query_title_list = []
    query_title_alignment_dict = {}
    for index, row in df.iterrows():
        if row[column_header_dict['positive or not']] == '+':
            query_title = row[column_header_dict['query title']]
            query_file = row[column_header_dict['query filename']]
            complete_query_title_list.append(query_title)
            if query_file.endswith('.afaa') or query_file.endswith('.afna'):
                if query_title not in query_title_alignment_dict.keys():
                    query_title_alignment_dict[query_title] = [query_file]
                else:
                    if not query_file in query_title_alignment_dict[query_title]:
                        query_title_alignment_dict[query_title] = query_title_alignment_dict[query_title] + [query_file]

            ## TEMPORARY!!:
            #elif query_title in ['AP5_mu', 'TTRAY1', 'AP3_sigma', 'COPI_beta','ClathrinHC']:
            #        query_title_alignment_dict[query_title] = [query_file]

    
    # Check that an alignment file name was assigned to each query title.
    titles_without_alignments = []
    for query_title in list(set(complete_query_title_list)):
        if query_title not in query_title_alignment_dict.keys():
            titles_without_alignments.append(query_title)
    assert len(titles_without_alignments) == 0, """No alignment file names
    could be found to use for the following query titles:
        %s""" % (str(titles_without_alignments))

    # Check that one and only one alignment file name has been assigned to each
    # query title.
    for query_title in query_title_alignment_dict.keys():
        assert len(query_title_alignment_dict[query_title]) > 0, """Could not
        identify an alignment file to use for query title: %s""" % query_title

        assert len(query_title_alignment_dict[query_title]) < 2, """Identified
        more than one potential alignment file to use for query title %s:
            %s"""\
        % (query_title, str(query_title_alignment_dict[query_title]))
    

    # Use the dict to add a column to the dataframe.

    # New column labels (for new columns to be appended).
    new_column_label_list = [column_header_dict['alignment name']]

    # Initiate new dataframe with columns to be appended/joined to existing
    # dataframe.
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list)

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)


    # Update the new column with alignment names.
    for index, row in df.iterrows():
        if row[column_header_dict['positive or not']] == '+':
            query_title = row[column_header_dict['query title']]
            alignment_to_use = query_title_alignment_dict[query_title][0]
            df.at[index, column_header_dict['alignment name']] = alignment_to_use

    # Write updated dataframe to output csv file.
    df.to_csv(outcsv, index=False)

    return outcsv
