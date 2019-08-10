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
"""Module for amoebae script with functions for the replace_seqs command.
"""
# Import built-in modules.
import sys
import os
import subprocess
import re
import shutil
import glob
import copy

# Import major package modules.
import pandas as pd

# Import amoebae modules.
import settings
from module_afa_to_nex import nex_to_afa, afa_to_nex, determine_alphabet
from module_add_to_db import make_blast_db
from misc_functions import get_fa_record_text_from_obj
from module_paralogue_counter import add_seq_to_alignment3
from module_align_to_profile_iter import do_align_iteratively
from module_mask_nex import mask_nex
from module_amoebae_trim_nex import trim_nex
from module_afa_to_nex import delete_extra_mesquite_lines
from module_amoebae import get_seqs_from_fasta_db

# Import functions for working with sequences, alignments, and trees.
#from ete3 import Tree, TreeStyle, Tree, TextFace, add_face_to_node
from Bio import SearchIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment


# Define functions to be used in amoebae.

def afa_to_fa(infilepath, outfilepath):
    """Takes an aligned fasta file (protein) and writes a .fa file.
    """
    inhandle = open(infilepath)
    outhandle = open(outfilepath, 'w')
    afa_file = SeqIO.parse(inhandle, "fasta")
    for record in afa_file:
        x = str(record.seq)
        record.seq = x.replace('-', '')
        outhandle.write(get_fa_record_text_from_obj(record))
    inhandle.close()
    outhandle.close()


def replace_seqs_in_alignment_with_seqs_from_fasta(alignment, fasta=None):
    """Takes sequences from an alignment, finds the top blastp hit in a given
    fasta file for each one, then aligns the top hits to the alignment, then
    removes the original sequences from the alignment, writes the alignemnt,
    and returns the path to the output alignment file.

    Note: Might be better to explicitly mask the new sequences to only the
    columns that align with the original trimmed sequences instead of applying
    general masking criteria, because some input alignments might be masked
    more stringently.
    """
    # Make a folder to contain temporary files.
    tempdirpath = alignment + '_temp'
    if os.path.isdir(tempdirpath):
        shutil.rmtree(tempdirpath)
    os.mkdir(tempdirpath)

    # Copy nexus alignment into temporary directory.
    #nex_ali_copy = os.path.join(tempdirpath, os.path.basename(alignment))
    #shutil.copyfile(alignment, nex_ali_copy) 

    # Remove extra mesquite lines if present.
    delete_extra_mesquite_lines(alignment)

    # Get number of sequences and accessions in the input alignment.
    alignment_seq_num = 0
    original_seq_accs = []
    for i in AlignIO.read(alignment, 'nexus'):
        alignment_seq_num += 1
        original_seq_accs.append(i.id)

    # Convert alignment to fasta.
    afa_alignment = os.path.join(tempdirpath,
            os.path.basename(alignment).rsplit('.', 1)[0] + '.afaa')
    nex_to_afa(alignment, afa_alignment) 

    # Remove '?' characters from alignment, and replace space characters in
    # fasta headers with underscores..
    modified_original_seq_accs = []
    afa_alignment2 = afa_alignment.rsplit('.', 1)[0]\
        + '_withoutquestionmarks.afaa'
    with open(afa_alignment) as infh, open(afa_alignment2, 'w') as o:
        for i in infh:
            if i.startswith('>'):
                #new_header_line = i.replace(' ', '_')
                new_header_line = i.rstrip().replace(' ', '_') + '__ORIGNALSEQ\n'
                o.write(new_header_line)
                # Assumes that there are no multi-line fasta headers.
                modified_original_seq_accs.append(new_header_line.strip()[1:])
            else:
                o.write(i.replace('?', '-'))

    # Check that the correct number of headers was found.

    # Make a copy of the alignment in nexus format.
    nex_ali_copy = os.path.join(tempdirpath, os.path.basename(alignment))
    afa_to_nex(afa_alignment2, nex_ali_copy)

    # Initiate a variable to store the path to a fasta file with top hits or
    # full-length sequences with the same IDs from the same genomes.
    top_hit_fasta = None

    # Initiate variable to store a list of all headers for sequences to be
    # added to the alignment.
    nonredun_top_hit_headers = None

    # If there is a fasta file provided, then find top hits for each sequence
    # in the input alignment.
    if fasta is not None:
        # Remove hyphens from fasta file (make it an unaligned fasta).
        fa_query_file = afa_alignment2.rsplit('.', 1)[0] + '.faa'
        afa_to_fa(afa_alignment2, fa_query_file)

        # Copy the input fasta file into the temporary directory.
        fasta_copy = os.path.join(tempdirpath, os.path.basename(fasta.rsplit('.',
            1)[0] + '.faa'))
        shutil.copyfile(fasta, fasta_copy) 

        # Make input fasta file a blastable database.
        make_blast_db(fasta_copy)

        # BLASTP search with the sequences from the alignment as the query into the
        # input fasta file.
        print('\nRunning blastp search')
        blast_outpath = fasta_copy.rsplit('.', 1)[0] + '_blastp_out.txt'
        subprocess.call(['blastp', '-query', fa_query_file, '-db', fasta_copy,
        '-out', blast_outpath, '-outfmt', '5', '-evalue', str(0.000000001),
        '-num_threads', str(4)]) 

        # Parse the output, and make a list of all the top hits for each
        # sequence/query/searchrecord.
        # Iterate through each search record (one for each query sequence) and
        # compile a list of headers for top hits.
        print('\nParsing blastp output file\n')
        top_hit_headers = []
        queries_retrieving_no_hits = copy.copy(modified_original_seq_accs)
        log_file_path = alignment.rsplit('.', 1)[0] + '_with_newseqs.csv' 
        with open(log_file_path, 'w') as o:
            o.write('Header of original sequence,Header of replacement sequence,E-value of top HSP\n')
            idx = SearchIO.index(blast_outpath, 'blast-xml')
            qnum = 0
            for key in sorted(idx.keys()):
                qnum += 1
                query_id = idx[key].id
                o.write(query_id + ',')
                print('\tQuery ' + str(qnum) + ':' + query_id)
                found_a_hit = False
                for hit in idx[key]:
                    found_a_hit = True
                    print('\t\tTop hit:' + hit.id)
                    print('\t\tTop HSP E-value:' + str(hit[0].evalue))
                    top_hit_headers.append(hit.id)
                    #if query_id in queries_retrieving_no_hits:
                    #    queries_retrieving_no_hits.remove(query_id)
                    queries_retrieving_no_hits.remove(query_id)
                    o.write(hit.id + ', ' + str(hit[0].evalue) + '\n')
                    break
                if not found_a_hit:
                    o.write('N/A,N/A\n')
            idx.close()

        # Remove redundant headers from the list.
        nonredun_top_hit_headers = list(set(top_hit_headers))

        # Calculate number of redundant sequences among top hits.
        num_redun_top_hits = len(top_hit_headers) - len(nonredun_top_hit_headers)

        # Report stats.
        print('\nNumber of sequences in input alignment:')
        print(alignment_seq_num)
        print('\nNumber of top hits:')
        print(len(top_hit_headers))
        print('\nNumber of nonredundant top hits:')
        print(len(nonredun_top_hit_headers))
        print('\nNumber of redundant sequences among top hits:')
        print(num_redun_top_hits)

        # Check that the number of top hits is not below 90 percent of the number
        # of query sequences.
        percent_of_queries_without_hits = round((len(queries_retrieving_no_hits) /\
            alignment_seq_num) * 100, 2)
        print("""\nTotal queries that did not retrieve any hits in the input fasta
        file: %s""" % str(len(queries_retrieving_no_hits)))
        print("""\nPercent of queries that did not retrieve any hits in the input
        fasta file: %s""" % str(percent_of_queries_without_hits))
        assert percent_of_queries_without_hits <= 10, """More than 10% of queries
        retrieved no hits in the input fasta file.""" 

        ## Check that the number of top hits is the same as the number of query
        ## sequences.
        #assert len(top_hit_headers) == alignment_seq_num, """Different number
        #of top hits than query sequences. There may be no hits for some sequences
        #in the input alignment."""

        # Retrieve the top hit sequences from the input fasta file, and align them
        # to the input alignment file.
        #print('\nAligning top hits to input alignment')
        #with open(fasta_copy) as infp:
        #    for seq in SeqIO.parse(infp, 'fasta'):
        #        if seq.id in nonredun_top_hit_headers:
        #            # Remove '?' characters from sequence?
        #            #seq.seq = str(seq.seq).replace('?', '')
        #            assert not '?' in str(seq.seq)
        #            # Align to input alignment.
        #            new_nex_path = nex_ali_copy + '_temp.nex' 
        #            add_seq_to_alignment3(seq, nex_ali_copy, new_nex_path)
        #            # Delete old nex.
        #            os.remove(nex_ali_copy)
        #            # Rename new nex as old nex.
        #            os.rename(new_nex_path, nex_ali_copy)
        print('\nGetting top hits from input fasta file')
        top_hit_fasta = fasta_copy + '_top_hits.faa'
        with open(fasta_copy) as infp, open(top_hit_fasta, 'w') as o:
            for seq in SeqIO.parse(infp, 'fasta'):
                if seq.id in nonredun_top_hit_headers:
                    # Write to file.
                    SeqIO.write(seq, o, 'fasta')

    # Otherwise look in relevant files in the Genomes directory.
    else:
        # Iterate through sequences in the input alignment and compile a list
        # of sequence objects corresponding to each sequence.
        seq_objects = []
        nonredun_top_hit_headers = []
        with open(afa_alignment) as infh:
            for seq in SeqIO.parse(infh, 'fasta'):
                # Extract the species name and sequence ID from the sequence header.
                header = seq.id
                #print('\nFinding full-length sequence for sequence: %s' % header)
                species_name = header.split('__')[0]
                #print('species name: ' + species_name)
                seq_id = header.split('__')[1]
                #print('sequence identifier: ' + seq_id)

                # Determine which databases could correspond to the species name,
                # given information in the genome info spreadsheet.
                db_names_with_species_name = []
                with open(settings.db_info_csv) as infh2:
                    dfx = pd.read_csv(infh2)
                    for index, row in dfx.iterrows():
                        if species_name.replace('_', ' ') in row['Species (if applicable)']:
                            filename = row['Filename']
                            if filename.endswith('.faa'):
                                #print(filename)
                                db_names_with_species_name.append(filename)

                # Check that at least one database might have the full-length sequence.
                assert len(db_names_with_species_name) > 0, """No databases
                identified for species name %s.""" % species_name

                # For each possible database/file, try to retrieve a sequence with
                # the ID.
                full_length_seq_obj = None
                for db_name in db_names_with_species_name:
                    # Call a function to retrieve the sequence from a database
                    # given an ID.
                    try:
                        # Try fast retrieval, and if that doesn't work then try
                        # the slow method.
                        full_length_seq_obj =\
                        get_seqs_from_fasta_db(db_name, [seq_id], False)[0]

                    except:
                        print('Looking in another file.')
                        continue

                # Check that a sequence was found.
                assert full_length_seq_obj is not None

                # Change sequence header to be identical to that for the
                # sequence that it is replacing.
                full_length_seq_obj.description = ''
                full_length_seq_obj.id = header

                # Add sequence object to list. 
                seq_objects.append(full_length_seq_obj)

            # Append each sequence object in the list to a fasta file.
            top_hit_fasta = alignment.rsplit('.', 1) [0] + '_full_length_seqs.faa'
            with open(top_hit_fasta, 'w') as o:
                for seq_object in seq_objects:
                    # Add header to list of headers for sequences to add.
                    nonredun_top_hit_headers.append(seq_object.id)

                    # Write to file.
                    SeqIO.write(seq_object, o, 'fasta')

            # Make list of sequence headers for sequences to add non-redundant.
            nonredun_top_hit_headers = list(set(nonredun_top_hit_headers))

    # Align top hits to copy of input alignment.
    print('\nAligning replacement sequences to input alignment')
    new_nex_path = alignment.rsplit('.', 1)[0] + '_with_newseqs.nex' 
    do_align_iteratively(nex_ali_copy, top_hit_fasta, new_nex_path)
    
    # Mask the alignment and write to the directory that the input alignment is
    # in.
    masked_nex_path = new_nex_path.rsplit('.', 1)[0] + '.mask.nex'
    mask_nex(new_nex_path, masked_nex_path)

    # Trim the alignment.
    trimmed_nex_path = masked_nex_path.rsplit('.', 1)[0] + '.trim.nex'
    trim_nex(masked_nex_path, trimmed_nex_path)
    
    # Remove all the original sequences from the trimmed alignment.
    new_ali_obj = AlignIO.read(trimmed_nex_path, 'nexus')
    #new_ali_obj_copy = copy.deepcopy(new_ali_obj)
    alphabet_to_use = determine_alphabet(alignment)
    new_ali_obj_reduced = MultipleSeqAlignment([], alphabet=alphabet_to_use)
    inum = -1
    for i in new_ali_obj:
        inum += 1
        if i.id in nonredun_top_hit_headers:
            new_ali_obj_reduced.append(new_ali_obj[inum])

    # Write reduced alignment to directory that input alignment is in.
    final_output_nex = trimmed_nex_path.rsplit('.', 1)[0] + '.newseqs.nex'
    AlignIO.write(new_ali_obj_reduced, final_output_nex, 'nexus')

    # Remove temporary directory and contents.
    shutil.rmtree(tempdirpath)

    if fasta is not None:
        # Report number of sequences for which top hits were not identified.
        if len(queries_retrieving_no_hits) > 0:
            print("""\nThe following sequences (%s percent)from the input alignment retrieved no
            blastp hits from the input fasta file:""" % percent_of_queries_without_hits)
            for i in queries_retrieving_no_hits:
                print('\t' + i)

        # Report stats again.
        print('\nNumber of sequences in input alignment:')
        print(alignment_seq_num)
        print('\nNumber of top hits:')
        print(len(top_hit_headers))
        print('\nNumber of nonredundant top hits:')
        print(len(nonredun_top_hit_headers))
        print('\nNumber of redundant sequences among top hits:')
        print(num_redun_top_hits)

    # Print paths to main output files.
    print('\nFasta file with identified replacement sequences:')
    print(top_hit_fasta)
    print('\nTrimmed nexus alignment with only replacement sequences:')
    print(final_output_nex)


    # Return the output alignment file path.
    return final_output_nex

