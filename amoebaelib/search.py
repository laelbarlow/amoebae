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
"""Module for functions for amoebae script that are for performing
similarity searches and parsing the output.
"""
# Import built-in modules.
#import argparse
#import sys
import os
import subprocess
import re
import settings
#import shutil
import glob
#import time
import pandas as pd
import math
#import copy

# Import modules from installed libraries/packages.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped

#from Bio import SearchIO

# Import SearchIO and suppress experimental warning
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

import amoebae_m
from srchresfile import SrchResFile
import column_header_lists 
from run_searches import get_query_list_from_file,\
get_db_list_from_file, get_out_query_list_path, get_out_db_list_path,\
determine_search_method, search_result_filepath, run_any_search,\
run_all_searches, get_out_hmm_path, get_query_subdir
from search_scaffolds import split_tblastn_hits_into_separate_genes,\
get_hit_seq_record_and_coord, get_hit_seq_record_and_coord2, get_cluster_range


def intersect(a, b):
    """ return the intersection of two lists """
    return list(set(a) & set(b))

#from get_datatype import get_dbtype



# Define functions to be used in amoebae.


def get_search_software_version(method):
    """Parses output from running the -h option to find what version a given
    program is.
    """
    # Find out version.
    version = None
    if 'blast' in method:
        blast_vers_output = subprocess.check_output([method, "-version"])
        version = str(blast_vers_output).replace('b\'blast',\
                'blast').replace('\\n\'', '').replace('+\\nP', '+\nP').split('\
 ', 2)[2]

    #else:
    #    version = 'SearchMethodHere'
    
    # Check that it worked.
    assert version is not None, """Version of software used (%s) could not be
    determined.""" % method

    # Return string specifying software version.
    return version


def get_search_output_file_type(filepath):
    """Take a search output file and return the name of the program that was
    used to generate it (and perhaps the format style).
    """
    with open(filepath) as infh:
        ...


def fwd_search_results_to_csv(search_outdir, csv_file):
    """Parse forward search output files and add a summary of the results to
    a csv file.
    """
    # Get list of search filenams (start with numbers in order).
    outfiles = glob.glob(os.path.join(search_outdir, '*out.txt'))

    # Loop through files and extract the relevant info for writing to the csv
    # file. 
    for f in outfiles:
        with open(f) as fh, open(csv_file, 'a') as o:
            # Get file type (name of program that was used to generate it).
            outfiletype = get_search_output_file_type(f)

            # Parse file according to its format.
            if outfiletype == 'blast' or outfiletype == 'hmmer':
                # Use the SearchIO module to parse.
                search_obj = SearchIO.read(fh)

                # Loop through forward search hits, and write info to csv file.
                # ***Use the pandas library....

            # Only blast and hmmer working for now... (later hhpred).
            else:
                pass


def get_corr_evalue(evalue):
    if evalue == 0:
        evalue = 1e-300
    return evalue


def get_redun_hits_in_dbs(query_title,
                          query_file_list,
                          db_file_list,
                          csv_file,
                          outdir,
                          timestamp,
                          blast_report_evalue_cutoff,
                          blast_max_target_seqs,
                          hmmer_report_evalue_cutoff,
                          hmmer_report_score_cutoff,
                          num_threads_similarity_searching
                          ):
    """Finds all items (sequences or profiles) in each given database file
    (fasta or hmm databases) that may be positive hits for query title (e.g.,
    'VAMP7' or 'R-SNAREs') for each given query file (sequence or profile) and
    appends a record of these hits to a given spreadsheet file path for later
    manual annotation to distinguish positive(redundant) hits from negative
    hits. 
    """
    # Make an output spreadsheet file if it does not already exist.
    if not os.path.isfile(csv_file):
        # Define string for header line in output spreadseheet.
        header_line = ','.join(['Query Title',
                                'Query File',
                                'Database File',
                                'Search Method',
                                'Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)',
                                'Hit Number',
                                'E-value',
                                'Order of magnitude of E-value difference compared to top hit',
                                'Length as percent of top hit sequence length',
                                'Accession',
                                'Description',
                                '\n'
                                ])

        # Write header line to spreadsheet.
        with open(csv_file, 'w') as csvh:
            csvh.write(header_line)

    # Use another function to run searches and write results to outputdir.
    run_all_searches(query_file_list,
                     db_file_list,
                     outdir,
                     blast_report_evalue_cutoff,
                     blast_max_target_seqs,
                     hmmer_report_evalue_cutoff,
                     hmmer_report_score_cutoff,
                     num_threads_similarity_searching)

    # Manage assignment of query title(s).
    find_query_titles = True
    if query_title is not None:
        find_query_titles = False

    # Parse search results and write summary to output.
    # ***The first part of this code block could perhaps be abstracted to
    # another function, because it will be needed as the first step in
    # interpreting results of forward and reverse searches as well.
    with open(csv_file, 'a') as o:
        # Loop over query file names.
        for q in query_file_list:
            # Define query title to use.
            if find_query_titles:
                query_title = amoebae_m.get_query_title_from_csv(q)
            # Loop over database file names.
            for d in db_file_list:
                # Find relevant result file in output directory.
                search_result_path = search_result_filepath(q, d, outdir)

                # Determine file format of search result file.
                # For certain formats the search method is required, so may
                # need to figure out how to get that info later (see
                # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html)
                #srch_file_prog_vers_fmt = srchresfile.get_srch_file_info(search_result_path)
                #srch_file_prog = srch_file_prog_vers_fmt[0]
                #srch_file_prog_vers = srch_file_prog_vers_fmt[1]
                #srch_file_format = srch_file_prog_vers_fmt[2]

                parsed_file_obj = SrchResFile(search_result_path)
                srch_file_prog = parsed_file_obj.program
                srch_file_prog_vers = parsed_file_obj.version
                srch_file_format = parsed_file_obj.format

                # Loop over hits (objects) in search result file.
                hit_num = -1
                top_hit_evalue = parsed_file_obj.hit_evalue(0)
                top_hit_len = len(parsed_file_obj.hit_sequence(0))
                for hit in SearchIO.read(search_result_path, srch_file_format): 
                    hit_num += 1

                    # Get E-value difference between current hit and top hit
                    # top HSPs.
                    e_top = top_hit_evalue
                    e_cur = parsed_file_obj.hit_evalue(hit_num)
                    evaldiff = get_evaldiff(get_corr_evalue(e_top),\
                            get_corr_evalue(e_cur))

                    # Get hit length as a percentage of the top hit sequence
                    # length.
                    cur_hit_len = len(parsed_file_obj.hit_sequence(hit_num))
                    percent_len = round((cur_hit_len / top_hit_len) * 100)

                    # Get description.
                    description = parsed_file_obj.hit_description(hit_num)
                    if not description.startswith('\"'):
                        description = '\"' + description
                    if not description.endswith('\"'):
                        description = description + '\"'

                    # Append line with relevant info to spreadsheet.
                    line_to_write = ','.join([query_title, 
                                      q,
                                      d, 
                                      srch_file_prog + ' ' + srch_file_prog_vers,
                                      '-', 
                                      str(hit_num + 1), 
                                      str(e_cur),
                                      str(evaldiff),
                                      str(percent_len),
                                      '\"' + parsed_file_obj.hit_id(hit_num) + '\"',
                                      description
                                      ])
                    o.write(line_to_write)
                    o.write('\n')
    # Remove output directory?
    #os.remove(...)

    print('DONE!')

    # Return path to main output file.
    return csv_file 


def get_query_len(query_filename):
    """Takes a query filename, gets the full path, get the query type, and
    returns the query length (assumes fasta, because even for HMMer searches
    the query listed will be .afaa).
    """
    query_length = None

    # Get full path.
    query_path = os.path.join(settings.querydirpath, query_filename)

    # Get query length.
    query_length = None
    for seq in SeqIO.parse(query_path, 'fasta'):
        query_length = len(seq)
        break

    # Check that it worked.
    assert not query_length is None, """Could not determine query length for
    query file: %s""" % query_filename

    # Return length.
    return query_length


def get_subsequences_matching_query(hit_obj):
    """Takes a SearchIO hit object, and returns a list with the first element
    being a string representing the coordinates of the query to which the
    subject sequence aligns, and the second being the corresponding subject
    subsequence(s).
    """
    # ***Eventually add more sophisticated approach similar to
    # search_scaffolds.py, so that multiple HSPs can be represented.
    return [str(list(hit_obj[0].query_range)).replace(' ', ''),\
            str(hit_obj[0].hit.seq).replace('-', '')]


def get_hit_obj_for_hsp_cluster(query_res_obj, cluster):
    """Take a SearchIO QueryResult object and a list of HSPs, and return the
    Hit object from the QueryResult object that contains the HSPs.
    """
    # Search for hits that contain the HSPs.
    hit_obj = None
    found = False
    for hit in query_res_obj:
        for hsp in hit:
            if hsp.hit_id == cluster[0].hit_id:
                if hsp.hit_description == cluster[0].hit_description:
                    hit_obj = hit
                    found = True
                    break
        if found:
            break

    # Check that a hit was found.
    assert hit_obj is not None, """Could not identify a Hit object that
    contains the given list of HSPs."""

    # Return the Hit object.
    return hit_obj


def get_rows_for_fwd_srch_df(df, 
                             q, 
                             d, 
                             search_result_path, 
                             column_label_list,
                             max_evalue, 
                             max_gap_setting, 
                             exonerate_score_threshold,
                             do_not_use_exonerate,
                             max_hits_to_sum,
                             max_length_diff
                             ):
    """Parses search result file and adds a corresponding row to an output
    pandas dataframe.
    """
    # Parse given search output file.
    parsed_file_obj = SrchResFile(search_result_path)
    srch_file_prog = parsed_file_obj.program
    srch_file_prog_vers = parsed_file_obj.version
    srch_file_format = parsed_file_obj.format

    # Extract some initial info.
    query_len = get_query_len(q)
    query_res_obj = SearchIO.read(search_result_path,\
            srch_file_format)

    # Determine the number of hits in the search result file.
    num_hits = parsed_file_obj.num_hits

    ## Modify the Bio.SearchIO query result object if the search program was
    ## tblastn, such that HSPs are grouped as separate hits if they likely
    ## represent separate genes. ???
    ##hit_list = query_res_obj.hits
    #if query_res_obj.program == 'tblastn' and num_hits > 0:
    #    query_res_obj = split_tblastn_hits_into_separate_genes(query_res_obj, settings.max_gap)
    #    #hit_list = split_tblastn_hits_into_separate_genes(query_res_obj, settings.max_gap)
        
    # Define a sub-dataframe to populate and return.
    subdf = pd.DataFrame(columns=column_label_list) # Define query title to use.
    query_title = amoebae_m.get_query_title_from_csv(q)

    # Initiate dataframe for a new row to append to dataframe
    # for output.
    new_row_df = pd.DataFrame(columns=column_label_list)
    new_row_df.loc[0] = ['-'] * len(column_label_list)
    #new_row_df.loc[0] = ['not???'] * len(column_label_list) # Temp

    # Fill in info to row.
    new_row_df.loc[0]['Query title'] = query_title
    new_row_df.loc[0]['Query file'] = q
    new_row_df.loc[0]['Query species (if applicable)'] = amoebae_m.get_query_taxon_from_csv(q)

    # Get "taxon" from query filename.
    #taxon = module_amoebae.get_query_taxon_from_filename(q)
    taxon = amoebae_m.get_query_taxon_from_csv(q)

    # Use query "taxon" to look up the database, if any, that
    # the query came from (in the query info csv file).
    new_row_df.loc[0]['Query database name'] = amoebae_m.get_db_filename_for_query_from_db_csv(taxon)

    new_row_df.loc[0]['Query accession (if applicable)'] = query_res_obj.id
    new_row_df.loc[0]['Query description'] = query_res_obj.description
    new_row_df.loc[0]['Query length'] = query_len
    new_row_df.loc[0]['Subject database species (if applicable)'] =\
        amoebae_m.get_species_for_db_filename(d) ### PROBLEM?
    new_row_df.loc[0]['Subject database file'] = d
    new_row_df.loc[0]['Forward search method'] = srch_file_prog + ' ' + srch_file_prog_vers


    if num_hits == 0:
        # No hits to record.
        new_row_df.loc[0]['Forward hit rank'] = '0'

        # Append row dataframe to dataframe for output.
        subdf = subdf.append(new_row_df, ignore_index=True)

    elif num_hits > 0:
        # ***Note: TBLASTN results are parsed very differently from results of
        # searches into protein databases, because a single subject sequence
        # (chromosome, scaffold, etc) may contain HSPs that represent different
        # genes. To account for this, the code below allows these to be
        # recorded as separate hits in the forward search summary.
        if parsed_file_obj.program == 'tblastn':
            # Get list of lists of HSPs corresponding to each potential gene in
            # the subject sequences.
            max_gap = max_gap_setting
            hsp_clusters = split_tblastn_hits_into_separate_genes(query_res_obj, max_gap)

            hit_num = -1
            # These lines are problematic...:
            top_hit_evalue = parsed_file_obj.hit_evalue(0) # CHANGE!!!
            top_hit_score = parsed_file_obj.hit_score(0) # CHANGE!!!
            top_hit_len = len(parsed_file_obj.hit_sequence(0)) # CHANGE!!!
            # Iterate over clusters of HSPs.
            #for hit in hit_list: 
            #print('\n\tProcessing clusters of HSPs within hits')
            for clusterplus in hsp_clusters: 
                cluster = clusterplus[0]
                #row_num += 1
                hit_num += 1

                # Only process hits up to the maximum number specified.
                if max_hits_to_sum != 0:
                    if hit_num + 1 > int(max_hits_to_sum):
                        break # *******

                # Get SearchIO Hit object from QueryResult object that
                # corresponds to the HSPs in the cluster.
                hit = get_hit_obj_for_hsp_cluster(query_res_obj, cluster)

                # Print message about the cluster.
                #print('\t\tCluster ' + str(hit_num + 1) + ': ' + str(len(cluster)) +\
                #' HSPs in subject sequence ' + hit.id + ' ' + hit.description)

                # Get score difference between current hit and top hit.
                score_cur = max([x.bitscore for x in cluster])
                scorediff = round(abs(score_cur - top_hit_score))

                # Get E-value difference between current hit and top hit
                # top HSPs.
                e_top = top_hit_evalue
                e_cur = min([x.evalue for x in cluster])
                evaldiff = get_evaldiff(get_corr_evalue(e_top), get_corr_evalue(e_cur))

                # Get a sequence record.
                hit_seq_record_and_coord = None
                if do_not_use_exonerate:
                    # Get a sequence record, without exonerate (just based on
                    # TBLASTN HSPs).
                    hit_seq_record_and_coord = get_hit_seq_record_and_coord(hit, cluster)
                else:
                    # Get a sequence record using exonerate to define the
                    # coding sequence and translation.
                    genetic_code_number = '1'
                    hit_seq_record_and_coord =\
                    get_hit_seq_record_and_coord2(search_result_path,
                                                  hit,
                                                  cluster,
                                                  d,
                                                  q,
                                                  genetic_code_number,
                                                  exonerate_score_threshold
                                                  )
                    if hit_seq_record_and_coord == None:
                        # Do not include the cluster in the results, because
                        # exonerate could not identify any sequence to
                        # translate (given the specified parameters).
                        continue

                # Define hit sequence object.
                hit_seq = hit_seq_record_and_coord[0] 

                # Get hit length as a percentage of the top hit sequence
                # length.
                cur_hit_len = len(hit_seq)
                percent_len = round((cur_hit_len / top_hit_len) * 100)

                new_row_df.loc[0]['Forward hit rank'] = hit_num + 1
                new_row_df.loc[0]['Forward hit score'] = score_cur
                new_row_df.loc[0]['Forward hit score difference from top hit score'] = scorediff
                new_row_df.loc[0]['Forward hit E-value (top HSP)'] = e_cur
                new_row_df.loc[0]['Forward hit E-value (top HSP) order of magnitude difference compared to top hit'] = evaldiff
                new_row_df.loc[0]['Forward hit length'] = cur_hit_len
                new_row_df.loc[0]['Forward hit length as a percentage of query length'] =\
                        round((cur_hit_len/query_len) * 100)
                new_row_df.loc[0]['Forward hit percent query cover'] =\
                        round((cur_hit_len/query_len) * 100) # Same as percent length, because all residues identified align to query sequence.
                #new_row_df.loc[0]['Forward hit accession'] = hit.id + '_(fwdhit' + str(hit_num + 1) + ')' 
                new_row_df.loc[0]['Forward hit accession'] = hit.id
                new_row_df.loc[0]['Forward hit description'] = hit.description


                # Need to handle nucleotide and protein sequences
                # differently!
                new_row_df.loc[0]['Forward hit sequence'] = '-'

                # Get subsequence(s) to use for reverse search.
                #range_and_seq = get_subsequences_matching_query(hit)
                #new_row_df.loc[0]['Forward hit coordinates of subsequence(s) that align(s) to query'] =\
                #range_and_seq[0] 
                #new_row_df.loc[0]['Forward hit subsequence(s) that align(s) to query'] = range_and_seq[1]

                subseq_and_coord = hit_seq_record_and_coord
                new_row_df.loc[0]['Forward hit coordinates of subsequence(s) that align(s) to query'] =\
                str(subseq_and_coord[1]).replace(' ', '')
                new_row_df.loc[0]['Forward hit description of subsequence(s) that align(s) to query']\
                    = subseq_and_coord[0].description
                new_row_df.loc[0]['Forward hit subsequence(s) that align(s) to query']\
                = str(subseq_and_coord[0].seq)

                # Fill in info about how close the hit is to the ends of the
                # subject sequence (searching in nucleotide data, do nothing
                # otherwise).
                cluster_range = get_cluster_range(cluster)
                subject_seq_len = hit.seq_len
                proximity = min([cluster_range[0], subject_seq_len - cluster_range[1]])
                new_row_df.loc[0]['Proximity (bp) to end of subject sequence (if searching in nucleotide sequences)']\
                = str(proximity)

                # Calculate difference in length between query and hit
                # sequences.
                length_diff = int(abs(query_len - cur_hit_len))

                # Determine whether criteria have been met.
                decis = '-'
                if e_cur <= float(max_evalue) and length_diff <= int(max_length_diff):
                    decis = '+'
                new_row_df.loc[0]['Positive/redundant (+) or negative (-) hit based on E-value criterion'] = decis

                # Make sure that all the info was identified and added to
                # the row.
                #assert not '???' in list(new_row_df.loc[0]), """Could not find
                #all information."""


                # Append row dataframe to dataframe for output.
                subdf = subdf.append(new_row_df, ignore_index=True)

        else:
            # Loop over hits in query_res_obj...
            hit_num = -1
            top_hit_evalue = parsed_file_obj.hit_evalue(0)
            top_hit_score = parsed_file_obj.hit_score(0)
            top_hit_len = len(parsed_file_obj.hit_sequence(0))
            #for hit in hit_list: 
            for hit in query_res_obj: 
                #row_num += 1
                hit_num += 1

                # Only process hits up to the maximum number specified.
                if max_hits_to_sum != 0:
                    if hit_num + 1 > int(max_hits_to_sum):
                        break # *******


                # Get score difference between current hit and top hit.
                score_cur = parsed_file_obj.hit_score(hit_num)
                scorediff = round(abs(score_cur - top_hit_score))

                # Get E-value difference between current hit and top hit
                # top HSPs.
                e_top = top_hit_evalue
                e_cur = parsed_file_obj.hit_evalue(hit_num)
                evaldiff = get_evaldiff(get_corr_evalue(e_top), get_corr_evalue(e_cur))

                # Get hit length as a percentage of the top hit sequence
                # length.
                hit_seq = parsed_file_obj.hit_sequence(hit_num) 
                cur_hit_len = len(hit_seq)
                percent_len = round((cur_hit_len / top_hit_len) * 100)

                new_row_df.loc[0]['Forward hit rank'] = hit_num + 1
                new_row_df.loc[0]['Forward hit score'] = score_cur
                new_row_df.loc[0]['Forward hit score difference from top hit score'] = scorediff
                new_row_df.loc[0]['Forward hit E-value (top HSP)'] = e_cur
                new_row_df.loc[0]['Forward hit E-value (top HSP) order of magnitude difference compared to top hit'] = evaldiff
                new_row_df.loc[0]['Forward hit length'] = cur_hit_len
                new_row_df.loc[0]['Forward hit length as a percentage of query length'] =\
                        round((cur_hit_len/query_len) * 100)

                # Get subsequence that aligns to query, and its coordinates
                # within the target sequence.
                subseq_and_coord = parsed_file_obj.hit_subsequence_and_coord(hit_num)

                new_row_df.loc[0]['Forward hit percent query cover'] =\
                        round((len(subseq_and_coord[0])*100)/query_len)
                new_row_df.loc[0]['Forward hit accession'] = parsed_file_obj.hit_id(hit_num)
                new_row_df.loc[0]['Forward hit description'] = parsed_file_obj.hit_description(hit_num)

                # Need to handle nucleotide and protein sequences
                # differently!
                if parsed_file_obj.program == 'tblastn':
                    new_row_df.loc[0]['Forward hit sequence'] = '-'
                else:
                    new_row_df.loc[0]['Forward hit sequence'] = str(hit_seq.seq)

                # Get subsequence(s) to use for reverse search.
                #range_and_seq = get_subsequences_matching_query(hit)
                #new_row_df.loc[0]['Forward hit coordinates of subsequence(s) that align(s) to query'] =\
                #range_and_seq[0] 
                #new_row_df.loc[0]['Forward hit subsequence(s) that align(s) to query'] = range_and_seq[1]

                new_row_df.loc[0]['Forward hit coordinates of subsequence(s) that align(s) to query'] =\
                str(subseq_and_coord[1]).replace(' ', '')
                new_row_df.loc[0]['Forward hit description of subsequence(s) that align(s) to query']\
                    = subseq_and_coord[0].description
                new_row_df.loc[0]['Forward hit subsequence(s) that align(s) to query']\
                = str(subseq_and_coord[0].seq)

                # Calculate difference in length between query and hit
                # sequences.
                length_diff = int(abs(query_len - cur_hit_len))

                # Determine whether search criteria have been met.
                decis = '-'
                if e_cur <= float(max_evalue) and length_diff <= int(max_length_diff):
                    decis = '+'
                new_row_df.loc[0]['Positive/redundant (+) or negative (-) hit based on E-value criterion'] = decis

                # Make sure that all the info was identified and added to
                # the row.
                #assert not '???' in list(new_row_df.loc[0]), """Could not find
                #all information."""


                # Append row dataframe to dataframe for output.
                subdf = subdf.append(new_row_df, ignore_index=True)

    # Return populated sub-dataframe to append to full dataframe.
    return subdf


def write_fwd_srch_res_to_csv(outdir,
                              query_file_list,
                              db_file_list,
                              csv_file, 
                              timestamp, 
                              max_evalue, 
                              max_gap_setting, 
                              exonerate_score_threshold,
                              do_not_use_exonerate,
                              max_hits_to_sum,
                              max_length_diff):
    """Parse output of a forward search (from running the fwd_srch command of
    amoebae) and append rows to input csv with information for interpreting
    the forward results. 
    """
    # Get database info spreadsheet path from settings module.
    db_info_csv = settings.db_info_csv

    # Define string for header line in output spreadsheet.
    # ***If this is changed, then make sure to change the corresponding code
    # for filling in data into the rows (below).
    column_label_list = column_header_lists.fwd_column_label_list

    # Make an output spreadsheet file if it does not already exist.
    if not os.path.isfile(csv_file):
        ## Write header line to spreadsheet.
        with open(csv_file, 'w') as csvh:
            csvh.write(','.join(column_label_list))


    # Get data from input file.
    df = pd.read_csv(csv_file) #, index_col=0)

    # Check that the headers match up the main spreadsheet and
    # the row to be appended.
    assert len(list(df.columns)) ==\
    len(column_label_list), """Unexpected number of column
    headers in input spreadsheet."""
    assert set(df.columns) == set(column_label_list),\
    """Input spreadsheet has different column header names than
    expected."""

    num_searches = len(query_file_list) * len(db_file_list)

    # Parse search results and write summary to output.
    #row_num = -1
    with open(csv_file, 'a') as o: 
        num = 0
        # Loop over query file names.
        for q in query_file_list:
            # Define query title to use.
            #query_title = module_amoebae.get_query_title_from_csv(q)
            # Loop over database file names.
            for d in db_file_list:
                num += 1
                print('Result ' + str(num) + ' of ' + str(num_searches))
                # Temporary 'if' statement???
                if q.rsplit('.', 1)[1] == 'afaa' and d.rsplit('.', 1)[1] == 'fna':
                    warning_text = """\nWARNING: Assuming that did not search
                    with profile query %s in nucleotide data %s""" % (q, d)
                    print(warning_text)
                    #o.write(warning_text)

                ## TEMPORARY elif statement:
                #elif d.startswith('Zmays'):
                #    print('Skipping Zea mays result to save time')

                else:
                    # Find relevant result file in output directory.
                    search_result_path = search_result_filepath(q, d, outdir)

                    # Print activity.
                    print('Extracting information from search result file ' +\
                            os.path.basename(search_result_path))

                    # Determine file format of search result file.
                    # For certain formats the search method is required, so may
                    # need to figure out how to get that info later (see
                    # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html)
                    parsed_file_obj = SrchResFile(search_result_path)
                    srch_file_prog = parsed_file_obj.program
                    srch_file_prog_vers = parsed_file_obj.version
                    srch_file_format = parsed_file_obj.format

                    # Get rows to add to forward search dataframe from info from
                    # the search output file.
                    #subdf = get_rows_for_fwd_srch_df(df, q, d, search_result_path, srch_file_prog,
                    #        srch_file_prog_vers, srch_file_format,
                    #        column_label_list)
                    subdf = get_rows_for_fwd_srch_df(df, 
                                                     q, 
                                                     d, 
                                                     search_result_path,
                                                     column_label_list, 
                                                     max_evalue, 
                                                     max_gap_setting,
                                                     exonerate_score_threshold,
                                                     do_not_use_exonerate,
                                                     max_hits_to_sum,
                                                     max_length_diff)

                    # Append sub-dataframe to full dataframe.
                    df = df.append(subdf, ignore_index=True)

        # Write updated dataframe to output spreadsheet.
        print('Writing dataframe to csv file')
        df.to_csv(csv_file + '_out.csv', index=False)

    # Return main output path.
    return csv_file 


#def get_query_subdir(outdir):
#    """Takes a reverse-search output directory and returns a path for a reverse
#    search query subdirectory.
#    """
#    return os.path.join(outdir, '0_rev_srch_queries')


def get_rev_query_path_from_fwd_srch_row(query_subdir, df_row, aasubseq=False):
    """Take a path for a directory containing reverse search queries (or to
    write such files to), and a row from a forward search summary csv parsed
    with pandas, and returns a path for the corresponding reverse search query
    file to be written or read.
    """
    # Only use the first 50 characters of the coordinates string.
    description = df_row['Forward hit coordinates of subsequence(s) that align(s) to query'][0:50]
    # Unless using the full amino acid sequence.
    if df_row['Subject database file'].endswith('faa'):
        if not aasubseq:
            description = 'full'

    acc = df_row['Forward hit accession']
    acc = acc.replace('__', '_')
    acc = acc.replace('/', '_')
    #acc = acc.replace(':', '_')
    #acc = acc.replace('+', '_')

    # Define the reverse search query filepath.
    rev_query_path = os.path.join(query_subdir,
            acc + '__' +\
            description\
            + '__' + df_row['Subject database file'])

    # Change filename extension to .faa if .fna, because reverse searches only
    # performed with amino acid sequences. TEMPORARY, because need to
    # accomodate hhsearch, etc in the future.
    rev_query_path = rev_query_path.rsplit('.', 1)[0] + '.faa'

    return rev_query_path


def get_rev_queries(csv_file, query_subdir, aasubseq, nafullseq):
    """Take a summary of forward search results, extract sequences to use as
    reverse search queries, write them to output files in a given directory.
    """
    # Read spreadsheet with forward search results.
    df = pd.read_csv(csv_file)

    # Loop over forward search hits, and select appropriate ones to write to
    # query files.
    reverse_search_query_filepaths = []
    already_written_query_list = []
    for index, row in df.iterrows():
        # Only consider forward search hits that meet the E-value cutoff.
        if row['Positive/redundant (+) or negative (-) hit based on E-value criterion'] == '+':
            # Treat hits that are individual fasta sequences different than
            # hits that are profiles.
            if row['Subject database file'].endswith('.faa')\
                    or row['Subject database file'].endswith('.fna'): 
                # Then hit is a single fasta sequence.

                # Determine whether prot or nucl.
                dtype = amoebae_m.get_dbtype_from_file_exten(row['Subject database file'])

                # Determine alphabet to use.
                #?

                # Assign attributes from csv file differently depending on
                # whether the full sequence or a subsequence is to be used for
                # reverse search.
                sequence = None
                description = None
                outexten = None

                if dtype == 'prot':
                    #outexten = '.faa'
                    if aasubseq:
                        description = row['Forward hit coordinates of subsequence(s) that align(s) to query']
                        sequence = row['Forward hit subsequence(s) that align(s) to query']
                    # By default, use full sequence.
                    else:
                        description = row['Forward hit description']
                        sequence = row['Forward hit sequence']

                elif dtype == 'nucl':
                    #outexten = '.fna'
                    if nafullseq:
                        description = row['Forward hit description']
                        sequence = row['Forward hit sequence']
                    # By default, use identified subsequence.
                    else:
                        description = row['Forward hit coordinates of subsequence(s) that align(s) to query']
                        sequence = row['Forward hit subsequence(s) that align(s) to query']

                # Check that it worked.
                assert sequence is not None, """Could not find which sequence
                to use for reverse search query."""
                assert description is not None, """Could not find which
                description to use for reverse search query."""

                # Long sequence descriptions may be problematic, because used
                # in output filenames...
                assert len(description) <= 1000, """Sequence description is
                really long."""
                
                # Check whether sequence already added to reverse search
                # queries.
                #seq_ident_string = row['Subject database file'] + row['Forward hit accession'] + description
                seq_ident_string =\
                get_rev_query_path_from_fwd_srch_row(query_subdir, row,
                        aasubseq)

                if seq_ident_string not in already_written_query_list:
                    # Initiate sequence object to write to reverse search query
                    # file.
                    row_seq_obj = SeqRecord(Seq(sequence))
                    row_seq_obj.id = row['Forward hit accession']
                    row_seq_obj.description = description

                    # Define reverse search query file path (use description as
                    # name).
                    # (Use the 'Forward hit description' and 'Forward hit
                    # description of sequence portion that aligns to query' as
                    # applicable to determine which reverse search result needs
                    # to be summarized in which row of the spreadsheet for the
                    # sum_rev_srch command).

                    # Make sure there are no slash characters in the sequence
                    # description before using it to name a file (no longer
                    # necessary, because the description is no longer used to
                    # name the file).
                    #print(description)
                    #assert '/' not in description, """Can't use sequence
                    #description to name query file, because it has a slash
                    #character in it."""

                    #description = description.replace('__', '_')

                    # Make sure there are no instances of '__' in the sequence
                    # description before using it to name a file (because that
                    # will be used to parse the description out).
                    #assert '__' not in description, """Can't use sequence
                    #description to name query file, because it has a '__'
                    #(double underscore) in it."""

                    # Define output filepath.
                    #rev_query_path = os.path.join(query_subdir,
                    #        row['Forward hit accession'] + '__' +\
                    #        row['Forward hit coordinates of subsequence(s) that align(s) to query']\
                    #        + '__' + row['Subject database file'])
                    rev_query_path =\
                    get_rev_query_path_from_fwd_srch_row(query_subdir, row,
                            aasubseq)

                    # Add query filepath to list for searching with.
                    reverse_search_query_filepaths.append(rev_query_path)

                    # Write sequence to filepath.
                    SeqIO.write(row_seq_obj, rev_query_path, 'fasta')

                    # Add query identity to list of those queries already used.
                    already_written_query_list.append(seq_ident_string)


            else:
                # Database consisted of profiles of some sort...
                x = False
                assert x, """Can't handle forward search hits in profile
                databases yet..."""

    # Check that the list of queries written is nonredundant.
    assert len(already_written_query_list) ==\
    len(set(already_written_query_list)), """Redundant queries written to
    reverse search query dir."""

    # Return reverse search query filepath list (absolute file paths).
    return reverse_search_query_filepaths


def get_csv_with_rev_path(csv_file):
    """Takes a csv file path and returns the name for a new version of the
    file, after addition of columns with info about reverse search results.
    """
    return csv_file.rsplit('.', 1)[0] + '_1.csv'


def get_redun_hit_dict(redun_hit_csv):
    """Takes a csv file listing redundant hits, and returns a nested dictionary
    with the structure:
    Query title > query filename > database filename > [redun hit accs]

    Only works with information from searches into sequence data, not profile
    databases.
    """
    # Load csv as a pandas dataframe.
    df = pd.read_csv(redun_hit_csv)

    # Construct dictionary.
    redun_hit_dict = {}
    for index, row in df.iterrows():
        # Define information to use.
        query_title = row['Query Title']
        #query_file = row['Query File']
        db_file = row['Database File']
        decis = row['Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)']
        acc = row['Accession']

        # Make an empty list to begin with.
        if query_title not in redun_hit_dict.keys():
            redun_hit_dict[query_title] = {}
        #if query_file not in redun_hit_dict[query_title].keys():
        #    redun_hit_dict[query_title][query_file] = {}
        #if db_file not in redun_hit_dict[query_title][query_file].keys():
        #    redun_hit_dict[query_title][query_file][db_file] = []
        if db_file not in redun_hit_dict[query_title].keys():
            redun_hit_dict[query_title][db_file] = []

        # Update dictionary if the row corresponds to a redundant hit.
        if decis == '+':
            # Otherwise just add accession to list.
            #redun_hit_dict[query_title][query_file][db_file] =\
            #    redun_hit_dict[query_title][query_file][db_file] + [acc]
            redun_hit_dict[query_title][db_file] =\
                redun_hit_dict[query_title][db_file] + [acc]

    # Return the constructed dict.
    return redun_hit_dict


def get_evaldiff(evalue1, evalue2):
    """Return absolute value of the order of magnitude evalue difference for
    two given E-values.  """
    return round(abs(math.log(float(evalue1),10) - math.log(float(evalue2),10)))


def write_rev_srch_res_to_csv(rev_srch_id,
                              outdir,
                              query_file_list,
                              db_file,
                              csv_file,
                              redun_hit_csv,
                              min_evaldiff,
                              timestamp,
                              aasubseq,
                              max_num_reverse_searches_per_database
                              ):
    """Parse output of a forward search (from running the fwd_srch command of
    amoebae) and append rows to input csv with information for interpreting
    the forward results. 
    """
    redun_hit_dict = None
    if redun_hit_csv is not None:
        # Get redundant hit dict from from input redundant hit csv file.
        redun_hit_dict = get_redun_hit_dict(redun_hit_csv)

    assert rev_srch_id is not None
    assert os.path.isdir(outdir)
    assert len(query_file_list) > 0
    assert os.path.isfile(csv_file)
    assert redun_hit_csv is None or os.path.isfile(redun_hit_csv)

    # Get database info spreadsheet path from settings module.
    db_info_csv = settings.db_info_csv
    assert os.path.isfile(db_info_csv)

    # Get same list of header titles that the sum_fwd_srch command makes (for
    # parsing).
    fwd_column_label_list = column_header_lists.fwd_column_label_list 

    ## Make an output spreadsheet file if it does not already exist.
    #if not os.path.isfile(csv_file):
    #    ## Write header line to spreadsheet.
    #    with open(csv_file, 'w') as csvh:
    #        csvh.write(','.join(column_label_list))


    # Get data from input file.
    print('\tReading input csv file into a pandas dataframe.')
    df = pd.read_csv(csv_file) #, index_col=0)

    ## Check that the headers match up the main spreadsheet and
    ## the row to be appended.
    #assert len(list(new_row_df.columns)) ==\
    #len(column_label_list), """Unexpected number of column
    #headers in input spreadsheet."""
    #assert set(new_row_df.columns) == set(column_label_list),\
    #"""Input spreadsheet has different column header names than
    #expected."""

    # Define a list of new column headers to be appended.
    new_column_label_list = column_header_lists.rev_column_label_list
    #assert len(new_column_labels_list) > 0

    # Initiate new dataframe with columns to be appended/joined to spreadsheet.
    #num_rows = len(list(df.index))
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list) 

    # Get a copy of all the original column headers, which may include some for
    # previous reverse search information.
    #original_csv_col_headers = df.columns.copy()

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)

    num_rev_srchs = df['Query title'].count()
    assert num_rev_srchs > 0

    # Iterate over rows in input spreadsheet with forward search results, and
    # determine which reverse search data needs to be found, and update the row
    # by adding this information.
    row_num = 0
    for index, row in df.iterrows():
        row_num += 1

        # Check that a reverse search needed to be done for this row.
        if row['Positive/redundant (+) or negative (-) hit based on E-value criterion'] == '-':
            pass

        # Prevent analyzing reverse searches when the forward hit rank is over
        # a certain number.
        elif max_num_reverse_searches_per_database != 0 and row['Forward hit rank'] > max_num_reverse_searches_per_database:
            pass

        else:
            print('\n\tParsing reverse search ' + str(row_num) + ' of ' +\
                    str(num_rev_srchs) + '(' +\
                    str(round(row_num*100/num_rev_srchs)) +\
                    '% complete for this reverse search db)')
        
            row['Reverse search database name'] = db_file

            # Get query info.
            query_title = row['Query title']
            query_file = row['Query file']

            # Get species for db_file.
            db_file_sp = amoebae_m.get_species_for_db_filename(db_file)
            row['Reverse search species (if applicable)'] = db_file_sp


            ## ***Assume that if no redun_hit_csv is provided, then hmmscan was
            ## used.
            #hmmscan_used = False
            #if redun_hit_csv is None:
            #    hmmscan_used = True

            ## Check that a redundant hit csv file was provided.
            #assert redun_hit_csv is not None, """Must provide a path to a
            #csv file that lists redundant hits in this case."""

            # Get redundant hit list.
            print('\n')
            print(query_title)
            print(query_file)
            print(db_file)

            redun_hit_list = None
            if redun_hit_csv is not None:
                #redun_hit_list =\
                #redun_hit_dict[query_title][query_file][db_file]
                redun_hit_list =\
                list(set(redun_hit_dict[query_title][db_file]))
                row['Redundant hit list applied'] = redun_hit_list

            # Get corresponding query file path.
            query_subdir = get_query_subdir(outdir)
            corresp_query_path =\
            get_rev_query_path_from_fwd_srch_row(query_subdir, row, aasubseq)
            assert os.path.isfile(corresp_query_path), """Could not figure out
            where the corresponding reverse search query is for a row. Expected
            file does not exist: %s""" % corresp_query_path

            # Get corresponding reverse search result file path.
            search_result_path = search_result_filepath(corresp_query_path,
                    db_file, outdir)
            print('\t\t' + os.path.basename(search_result_path))

            # Determine file format of search result file.
            # For certain formats the search method is required, so may
            # need to figure out how to get that info later (see
            # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html)
            parsed_file_obj = SrchResFile(search_result_path)
            srch_file_prog = parsed_file_obj.program
            srch_file_prog_vers = parsed_file_obj.version
            srch_file_format = parsed_file_obj.format


            ## Loop over hits (objects) in reverse search result file, extract
            ## the necessary information, add it to the dataframe.

            # ***This needs to be customized for different search programs...

            #hit_num = 0
            #top_hit_evalue = None
            #top_hit_score = None
            #top_hit_len = None

            # Get top hit and first negative/nonredun hit (if present) by parsing
            # search result file.
            #top_hit_obj = None
            top_hit_acc = None
            first_neg_hit_rank = None
            query_res_obj = SearchIO.read(search_result_path,\
                    srch_file_format)
            if len(query_res_obj) >= 1:
                top_hit_acc = parsed_file_obj.hit_id(0)
                first_neg_hit_rank = None
                if redun_hit_list is not None: 
                    first_neg_hit_rank = parsed_file_obj.rank_of_first_nonredun_hit(redun_hit_list)
                else:
                    # The redun_hit_list is None when there is no redundant hit
                    # list CSV file input via the --redun_hit_csv option.

                    #######################################
                    # ***This old code doesn't make sense and results in
                    # misleading summary of results:

                    # Assume any hits after the top reverse search hit are not
                    # equivalent/orthologous to the original query.
                    #if len(query_res_obj) > 1:
                    #    first_neg_hit_rank = 1 # (the first rank is 0)

                    ## Assume that the top reverse search hit matches the
                    ## original query.
                    #redun_hit_list = [top_hit_acc]
                    #######################################

                    # Indicate that there are no hits to be considered
                    # nonredundant with the original query in the reverse
                    # search results.
                    first_neg_hit_rank = None

                    # All reverse search hits should be considered redundant,
                    # so all their accessions should be added to the list.
                    redun_hit_list = parsed_file_obj.all_hit_ids()


            #if len(query_res_obj) >= 1:
            #    hit_num = 0
            #    for hit in query_res_obj: 
            #        hit_num += 1
            #        if hit_num == 1:
            #            top_hit_obj = hit
            #        if hit.id not in redun_hit_list:
            #            first_neg_hit_obj = hit
            #            break

            # Extract info from hits, if any, and interpret. 

            if top_hit_acc is not None:
                # Then there is at least one reverse search hit. 
                # Fill in info for top hit.
                row['Top reverse hit description'] = parsed_file_obj.hit_description(0)
                row['Top reverse search hit E-value'] = parsed_file_obj.hit_evalue(0)
                row['Top reverse search hit score'] = parsed_file_obj.hit_score(0)

                if top_hit_acc in redun_hit_list:
                    # Then potential positive.

                    if first_neg_hit_rank is not None:
                        # Then there is at least one non-redundant hit.
                        row['First non-redundant hit'] = parsed_file_obj.hit_description(first_neg_hit_rank)
                        evaldiff =\
                            get_evaldiff(get_corr_evalue(parsed_file_obj.hit_evalue(0)),\
                                    get_corr_evalue(parsed_file_obj.hit_evalue(first_neg_hit_rank)))
                        row['Order of magnitude E-value difference between top hit and first non-redundant hit'] =\
                                evaldiff
                        scorediff =\
                            parsed_file_obj.hit_score(0) - parsed_file_obj.hit_score(first_neg_hit_rank)
                        row['Difference in score between top hit and first non-redundant hit'] =\
                                scorediff
                        # Determine whether the criteria are met.
                        if float(evaldiff) >= float(min_evaldiff):
                            # Result is positive.
                            row['Note'] = 'The order of magnitude difference is greater than or equal to %s'\
                                    % str(min_evaldiff)
                            row['Positive (+) or negative (-) hit based on reverse search criteria'] = '+'

                    else:
                        # Then there are no non-redundant hits. 
                        # Result is positive.
                        row['Note'] = 'All reverse search hits are redundant with the original query'
                        row['Positive (+) or negative (-) hit based on reverse search criteria'] = '+'
                else:
                    # Then definitely not positive.
                    row['Note'] = 'Top reverse search hit is not the original query'
                    row['Positive (+) or negative (-) hit based on reverse search criteria'] = '-'

            else:
                # Then there are no reverse search hits.
                # Result is negative.
                row['Note'] = 'No hits in reverse search'
                row['Positive (+) or negative (-) hit based on reverse search criteria'] = '-'



            # ...
            # Update row in dataframe with new information.
            df.loc[index] = row


    # ...
    # For end of function:
    assert row_num > 0

    # Give the newly added columns unique names so that they do not get
    # confused with other similar sets of columns for other reverse search
    # results.
    #df.columns = fwd_column_label_list + [x + ' (' + rev_srch_id + ')' for x in new_column_label_list]
    #print('df.columns slice:')
    #print(str(list(df.columns[0:-len(new_column_label_list)])))
    #print('altered column headers:')
    #print(str([x + ' (' + rev_srch_id + ')' for x in new_column_label_list]))
    df.columns = list(df.columns[0:-len(new_column_label_list)]) + [x + ' (' + rev_srch_id + ')' for x in new_column_label_list]

    # Write joined to output path.
    final_outfp = get_csv_with_rev_path(csv_file)
    df.to_csv(final_outfp, index=False)

    # Return main output path.
    return final_outfp 


def write_interp_csv(csv_file, outfp, fwd_evalue_cutoff, rev_evalue_cutoff):
    """Take a csv file and write a new one with an additional column with
    interpretation of which forward search results are positive based on all
    the reverse search results performed.
    """
    # Read the csv file into a pandas dataframe.
    df = pd.read_csv(csv_file)

    # Compile a list of headers for columns that contain relevant +/-
    # information, and E-values.
    relev_col_headers = []
    fwd_evalue_header = None 
    rev_evalue_headers = []
    for i in df.columns:
        # Record if the column has +/- info.
        if i.startswith('Positive (+) or negative (-) hit based on reverse search criteria'):
            relev_col_headers.append(i)
        # Record if the column contains forward search E-values.
        if i == 'Forward hit E-value (top HSP)':
            fwd_evalue_header = i
        # Record if the column contains reverse search E-values.
        if i.startswith('Top reverse search hit E-value'):
            rev_evalue_headers.append(i)
    assert len(relev_col_headers) >= 1, """Could not identify columns
    containing relevant information in input spreadsheet: %s""" % csv_file
    assert fwd_evalue_header is not None, """Could not identify column
    containing forward search E-values in input spreadsheet: %s""" % csv_file
    assert len(rev_evalue_headers) >= 1, """Could not identify columns
    containing reverse search E-values in input spreadsheet: %s""" % csv_file

    # Define a list of new column headers to be appended.
    new_column_label_list = column_header_lists.interp_column_label_list

    # Initiate new dataframe with columns to be appended/joined to spreadsheet.
    #num_rows = len(list(df.index))
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list) 

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)

    # Iterate over rows in input spreadsheet with forward search results, and
    # determine which reverse search data needs to be found, and update the row
    # by adding this information.
    for index, row in df.iterrows():
        decis = '+'
        for header in relev_col_headers:
            if row[header] == '-':
                decis = '-'
                break
        if fwd_evalue_cutoff is not None:
            if decis == '+':
                if float(row[fwd_evalue_header]) > float(fwd_evalue_cutoff):
                    decis = '-'
        if rev_evalue_cutoff is not None:
            if decis == '+':
                for header in rev_evalue_headers:
                    if float(row[header]) > float(rev_evalue_cutoff):
                        decis = '-'
                        break

        # Check that a decision could be made.
        assert decis == '+' or decis == '-', """Could not determine whether the
        results are positive or negative for row number %s.""" % str(index + 1)

        # Record decision in dataframe.
        df.at[index, 'Collective interpretation of reverse search results'] =\
        decis

    # Write new dataframe to output csv file.
    df.to_csv(outfp, index=False)


def write_fwd_srch_interp_csv(csv_file, outfp, score_cutoff):
    """Take a csv file and write a new one with an additional column with
    interpretation of which forward search results are positive based on all
    the reverse search results performed.
    """
    # Read the csv file into a pandas dataframe.
    df = pd.read_csv(csv_file)

    ## Compile a list of headers for columns that contain relevant +/-
    ## information.
    #relev_col_headers = []
    #for i in df.columns:
    #    if i.startswith('Positive (+) or negative (-) hit based on reverse search criteria'):
    #        relev_col_headers.append(i)
    #assert len(relev_col_headers) >= 1, """Could not identify columns
    #containing relevant information in input spreadsheet: %s""" % csv_file

    # Define a list of new column headers to be appended.
    new_column_label_list = column_header_lists.interp_column_label_list

    # Initiate new dataframe with columns to be appended/joined to spreadsheet.
    #num_rows = len(list(df.index))
    df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    # Set default value of all fields to '-'.
    for index, row in df_new_cols.iterrows():
        df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list) 

    # Join constructed dataframe to input dataframe.
    df = df.join(df_new_cols)

    # Iterate over rows in input spreadsheet with forward search results, and
    # determine which reverse search data needs to be found, and update the row
    # by adding this information.
    for index, row in df.iterrows():
        decis = None
        if row['Forward hit score'] is not '-' and float(row['Forward hit score']) > float(score_cutoff):
            decis = '+'
        else:
            decis = '-'

        # Check that a decision could be made.
        assert decis == '+' or decis == '-', """Could not determine whether the
        results are positive or negative for row number %s.""" % str(index + 1)

        # Record decision in dataframe.
        df.at[index, 'Collective interpretation of reverse search results'] =\
        decis

    # Write new dataframe to output csv file.
    df.to_csv(outfp, index=False)




def write_redun_hit_interp_csv(csv_file, outfp):
    """Take a csv file and write a new one with an additional column with
    interpretation of which forward search results redundant for the query.

    ***This currently very inefficient.

    ***Currently assumes that different query titles refer to
    non-overlapping/non-nested sets of related proteins.
    """
    # Read the csv file into a pandas dataframe.
    df = pd.read_csv(csv_file)

    ## Compile a list of headers for columns that contain relevant +/-
    ## information.
    #relev_col_headers = []
    #for i in df.columns:
    #    if i.startswith('Positive (+) or negative (-) hit based on reverse search criteria'):
    #        relev_col_headers.append(i)
    #assert len(relev_col_headers) >= 1, """Could not identify columns
    #containing relevant information in input spreadsheet: %s""" % csv_file

    ## Define a list of new column headers to be appended.
    #new_column_label_list = ['interp']

    ## Initiate new dataframe with columns to be appended/joined to spreadsheet.
    ##num_rows = len(list(df.index))
    #df_new_cols = pd.DataFrame(columns=new_column_label_list, index=df.index)

    ## Set default value of all fields to '-'.
    #for index, row in df_new_cols.iterrows():
    #    df_new_cols.loc[index] = [] + ['-'] * len(new_column_label_list) 

    ## Join constructed dataframe to input dataframe.
    #df = df.join(df_new_cols)

    # Iterate over rows and decide whether it is redun or not.
    combos_for_which_neg_found = []
    for index, row in df.iterrows():
        interp = '+'

        query_title = row['Query Title']
        query_file = row['Query File']
        acc = row['Accession']
        hit_rank = row['Hit Number']
        combo = (query_title, query_file)
        #print(' '.join([query_title, query_file, str(hit_rank), acc]))

        if combo in combos_for_which_neg_found:
            # A higher-ranking hit in the same search was non-redundant, so
            # this one is non-redundant as well.
            interp = '-'

        else:
            # Look for evidence that the hit is non-redundant.
            for index2, row2 in df.iterrows():
                query_title2 = row2['Query Title']
                query_file2 = row2['Query File']
                acc2 = row2['Accession']
                hit_rank2 = row2['Hit Number']
                #combo2 = (query_title2, query_file2)
                if query_title2 != query_title and acc2 == acc and hit_rank2 == 1:
                    #print('\t' + ' '.join([query_title2, query_file2, str(hit_rank2), acc2]))
                    interp = '-'
                    combos_for_which_neg_found.append(combo)
                    break
     
        # Record whether redundant or not.
        #df.at[index, 'interp'] = interp
        df.at[index, 'Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)'] = interp

    # Deal with hits that are now listed as redundant for more than one query
    # title.
    for index, row in df.iterrows():
        query_title = row['Query Title']
        query_file = row['Query File']
        acc = row['Accession']
        hit_rank = row['Hit Number']
        evalue = row['E-value']
        redun = row['Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)']

        if redun == '+': 
            # Loop over all rows, and determine whether the current hit is
            # retrieved by another query title with a better e-value.
            for index2, row2 in df.iterrows():
                query_title2 = row2['Query Title']
                query_file2 = row2['Query File']
                acc2 = row2['Accession']
                hit_rank2 = row2['Hit Number']
                evalue2 = row2['E-value']
                redun2 = row2['Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)']

                if acc2 == acc and query_title2 != query_title and redun2 == '+':
                    if float(evalue2) < float(evalue):
                        redun = '-'
                        print(acc + ' is redundant for query title ' +\
                                query_title2 + ' not for query title ' +\
                                query_title + ' (E-values of ' + str(evalue2)\
                                + ' and ' + str(evalue) + ' respectively).')
                        break
                    else:
                        pass

            # Change the decision if necessary.
            df.at[index, 'Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)'] = redun


    # Check that the lists of redundant accessions are mutually exclusive sets.
    # Generate dict of accession lists for each query_title.
    query_title_acc_list_dict = {}
    for index, row in df.iterrows():
        query_title = row['Query Title']
        acc = row['Accession']
        redun = row['Positive/redundant (+) or negative (-) hit for queries with query title (edit this column)']
        if redun == '+':
            if query_title not in query_title_acc_list_dict.keys():
                query_title_acc_list_dict[query_title] = [acc]
            else:
                query_title_acc_list_dict[query_title] = query_title_acc_list_dict[query_title] + [acc]
    # Check that the lists of redundant accessions in the dict are non-overlapping.
    for key in query_title_acc_list_dict.keys():
        for key2 in query_title_acc_list_dict.keys():
            if key != key2:
                assert len(intersect(query_title_acc_list_dict[key], query_title_acc_list_dict[key2])) == 0, """Redundant
                    accession lists overlap."""


    # Write new dataframe to output csv file.
    df.to_csv(outfp, index=False)




