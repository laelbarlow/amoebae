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
"""Module for amoebae script that are general purpose (used in modules that
define more specialized functions).
"""
# Import built-in modules.
#import argparse
import sys
import os
sys.path.append(os.path.dirname(sys.path[0]))
import subprocess
import re
from datapaths import DataPaths
#import shutil
import glob
import time
import datetime
import pandas as pd
import platform
import collections

# Import modules from installed libraries/packages.
from Bio import SeqIO
from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from get_datatype import get_dbtype
from afa_to_nex import delete_extra_mesquite_lines 
#from paralogue_counter import get_seq_obj_from_srch_res_csv_info, get_hit_range_from_hsp_ranges
#from module_amoebae import get_seq_obj_from_srch_res_csv_info,\
#get_hit_range_from_hsp_ranges 

# Define functions to be used in amoebae.

def get_corr_fasta_exten(infp):
    """Determines the correct filename extension that a given fasta file should
    have.
    """
    # Automatically determine what the dbtype is.
    dbtype = get_dbtype(infp)
    assert dbtype == 'prot' or dbtype == 'nucl', """Error: Could not determine
    data type in file " + f + " using get_datatype module."""

    # Determine appropriate filename extension to use.
    exten = None
    if dbtype == 'prot':
        exten = 'faa'
    elif dbtype == 'nucl':
        exten = 'fna'

    return exten


def get_dbtype_from_file_exten(infp):
    """Determine data type from fasta filename extension.
    """
    dbtype = None
    exten = infp.rsplit('.', 1)[1]
    if exten == 'faa':
        dbtype = 'prot'
    elif exten == 'fna':
        dbtype = 'nucl'
    assert dbtype is not None, """Could not determine data type for file:\n\t
        %s""" % infp
    return dbtype


def get_seqs_from_fasta_db(db_name, accs, main_data_dir, slow=False):
    """Returns a list of SeqRecord objects corresponding to the given accessions
    in the given database file. 
    """
    # Get database directory from DataPaths(main_data_dir) module.
    db_dir = DataPaths(main_data_dir).dbdirpath

    # Get database filepath.
    #db_path = os.path.join(db_dir, db_name)
    db_path = os.path.join(db_dir, os.path.basename(db_name))
    assert os.path.isfile(db_path), """Path is not a file: %s""" % db_path

    # Old way:
    ## Parse database file and compile a list of sequence objects to return. 
    #seq_objs = []
    #with open(db_path) as dbh:
    #    for seq in SeqIO.parse(dbh, 'fasta'):
    #        acc = seq.id.strip()
    #        if acc in accs:
    #            seq_objs.append(seq)
    #        if len(seq_objs) == len(accs):
    #            break


    # Retrieve sequences from fasta file and write to a temporary file.
    temp_fa_path = db_path + '_TEMP_FASTA.fa'
    if not slow:
        # Use esl-sfetch to retrieve the sequences and write to a temporary file.
        if os.path.isfile(temp_fa_path):
            os.remove(temp_fa_path)
        with open(temp_fa_path, 'a') as o:
            for acc in accs:
                # Get sequence as text.
                subprocess.call(['esl-sfetch', db_path, acc], stdout=o)
    elif slow:
        # Parse sequence using a slower method that does not make use of
        # esl-sfetch.
        if os.path.isfile(temp_fa_path):
            os.remove(temp_fa_path)
        with open(temp_fa_path, 'a') as o, open(db_path) as db_handle:
            for acc in accs:
                # Get sequence as text.
                all_seq_ids = None
                with open(db_path) as db_handle:
                    all_seq_ids = [x.id for x in SeqIO.parse(db_handle, 'fasta')]  
                if acc in all_seq_ids:
                    with open(db_path) as db_handle:
                        for x in SeqIO.parse(db_handle, 'fasta'):
                            if x.id == acc:
                                # Write to temp fasta file.
                                SeqIO.write([x], o, 'fasta')
                                break
                else:
                    accs_that_start_with_acc = []
                    with open(db_path) as db_handle:
                        for x in SeqIO.parse(db_handle, 'fasta'):
                            if x.id.startswith(acc):
                                accs_that_start_with_acc.append(x.id)
                    # Check that only one accession starts with.
                    #assert len(accs_that_start_with_acc) == 1, """More than one
                    #accession in file starts with %s""" % acc
                    if len(accs_that_start_with_acc) < 1:
                        print("No accessions start with %s" % acc)
                    elif len(accs_that_start_with_acc) > 1:
                        print("More than one accession starts with %s" % acc)

                    if len(accs_that_start_with_acc) >= 1:
                        with open(db_path) as db_handle:
                            for x in SeqIO.parse(db_handle, 'fasta'):
                                if x.id.startswith(acc):
                                    # Write to temp fasta file.
                                    SeqIO.write([x], o, 'fasta')
                                    break

    # Parse the fasta file to get Seq objects.
    seq_objs = []
    for s in SeqIO.parse(temp_fa_path, 'fasta'):
        seq_objs.append(s) 

    # Remove the temporary fasta file.
    os.remove(temp_fa_path)

    # Return the list of sequence objects.
    return seq_objs


def get_subseq_from_fasta_db(db_name, acc, subseq_coord, main_data_dir):
    """Returns a SeqRecord object corresponding to the subsequence with the
    given coordinates in the sequence with the given accessions in the given
    database file. 

    Note: the input subsequence coordinates ('subseq_coord') are the start and
    end residue numbers for the subsequence, not python-style slices.
    """
    # Get database directory from DataPaths(main_data_dir) module.
    db_dir = DataPaths(main_data_dir).dbdirpath

    # Get database filepath.
    #db_path = os.path.join(db_dir, db_name)
    db_path = os.path.join(db_dir, os.path.basename(db_name))
    assert os.path.isfile(db_path), """Path is not a file: %s""" % db_path

    # Old way:
    ## Parse database file and compile a list of sequence objects to return. 
    #seq_objs = []
    #with open(db_path) as dbh:
    #    for seq in SeqIO.parse(dbh, 'fasta'):
    #        acc = seq.id.strip()
    #        if acc in accs:
    #            seq_objs.append(seq)
    #        if len(seq_objs) == len(accs):
    #            break

    # Use esl-sfetch to retrieve the sequences and write to a temporary file.
    temp_fa_path = db_path + '_TEMP_FASTA.fa'
    if os.path.isfile(temp_fa_path):
        os.remove(temp_fa_path)
    with open(temp_fa_path, 'a') as o:
        # Get sequence as text.
        subprocess.call(['esl-sfetch', db_path, acc], stdout=o)

    # Parse the fasta file to get Seq objects.
    seq_obj = None
    seq_obj = SeqIO.read(temp_fa_path, 'fasta')
    seq_obj.description = seq_obj.description.rstrip('\"') + ' ' + str(subseq_coord) + '\"'

    # ***Re-use code from search scaffolds to verify validity of input
    # subseq_coord...?
    #...

    # Construct new sequence.
    new_seq = ''
    for subseq in subseq_coord:
        start = subseq[0]
        end = subseq[1]
        new_seq = new_seq + seq_obj.seq[start: end + 1] # Double-check this!

    seq_obj.seq = new_seq

    # Remove the temporary fasta file.
    os.remove(temp_fa_path)

    # Return the list of sequence objects.
    return seq_obj
                

def get_query_taxon_from_filename(query_filename):
    """Takes a query file name and extracts the database name (abbreviation of
    species name usually) or taxon name/indicator from the file name.

    Assumes that the file is named as:
        query_title + '_' + dbname/taxon + '_' + accession (if appl.) + extension

        where neither query_title nor dbname have underscores in them.
    """
    try:
        return query_filename.split('_')[1]
    except:
        return query_filename.rsplit('.', 1)[0]


def get_query_title_from_filename(query_filename):
    """Takes a query file name and extracts the query title/name from the file
    name.

    Assumes that the file is named as:
        query_title + '_' + dbname/taxon + '_' + accession (if appl.) + extension

        where neither query_title nor dbname have underscores in them.
    """
    return query_filename.split('_')[0]


def get_query_title_from_csv(query_filename, main_data_dir):
    """Take a query filename, look up corresponding query title in the query
    directory csv specified in the DataPaths(main_data_dir) module, and return that.
    """
    # Parse query info csv file.
    query_info_csv_path = DataPaths(main_data_dir).query_info_csv
    df = pd.read_csv(query_info_csv_path, encoding='utf-8')

    # Get query title.
    df.set_index('Filename', inplace=True)
    query_title = df.loc[query_filename]['Query title']

    # Check that query title is not a dataframe.
    assert not isinstance(query_title, pd.DataFrame), """Could not properly
    parse query info file (%s), there may be duplicate entries in this file.""" \
        % query_info_csv_path 

    # Check that query title is a string.
    assert isinstance(query_title, str), """Could not properly identify query
    title for query file %s in CSV file %s. Please ensure that the query title
    contains some alphabetic characters. Query title identified: %s""" \
            % (query_filename, query_info_csv_path, query_title) 

    # Return query title.
    return query_title


def get_query_taxon_from_csv(query_filename, main_data_dir):
    """Take a query filename, look up corresponding query taxon in the query
    directory csv specified in the DataPaths(main_data_dir) module, and return that.
    """
    # Parse query info csv file.
    df = pd.read_csv(DataPaths(main_data_dir).query_info_csv, encoding='utf-8')

    # Return query title.
    df.set_index('Filename', inplace=True)
    return df.loc[query_filename]['Query taxon (species if applicable)']


def get_species_from_db_csv(taxon, main_data_dir):
    """Take a database name/species abbreviation from a taxon name extracted
    from a query filename. If there is a corresponding species name in the
    database directory information csv file specified in the DataPaths(main_data_dir) module
    return that.
    """
    df = pd.read_csv(DataPaths(main_data_dir).db_info_csv, encoding='utf-8')
    # Species name to return is not applicable by default.
    sp = '-'
    for f in list(df['Filename']):
        if f.rsplit('.', 1)[0] == taxon:
            df.set_index('Filename', inplace=True)
            sp = df.loc[f]['Species (if applicable)']
            break
    return sp


def get_db_filename_for_query_from_db_csv(taxon, main_data_dir):
    """Take a database name/species abbreviation or taxon name extracted from a
    query filename, and if there is a corresponding database file name in the
    database directory information csv file specified in the DataPaths(main_data_dir) module,
    then return that. Otherwise, just return 'N/A'.
    """
    # Define database file name as not applicable, by default.
    dbfn = '-'

    # Check whether the given "taxon" name exists in the database info csv.
    try:
        # Try loading the dataframe.
        df = pd.read_csv(DataPaths(main_data_dir).db_info_csv, encoding='utf-8')
    except:
        # Print an error message.
        print("""Error: Could not load contents of csv file as pandas dataframe:\n\n
        \t%s\n\nCheck that the file was saved properly in comma separated value
        format (UTF-8 encoding).""" % DataPaths(main_data_dir).db_info_csv)

        # Exit the script.
        print('Quitting script.')
        sys.exit()

    col_list = df['Taxon'].tolist()
    # If it does, then get filename for corresponding database.
    if taxon in col_list:
        df.set_index('Taxon', inplace=True)
        dbfn = df.loc[taxon]['Filename']

    # Return database file name or 'N/A'.
    return dbfn


def get_species_for_db_filename(db_filename, main_data_dir):
    """Takes a database filename, and returns the species name that appears in
    the database info csv file (may be '-' if not applicable).
    """
    df = pd.read_csv(DataPaths(main_data_dir).db_info_csv, encoding='utf-8')
    df.set_index('Filename', inplace=True)
    sp = df.loc[db_filename]['Species (if applicable)']

    #print('\nTrying to get species name from genome info csv file.')
    #print('genome info file path: ' + DataPaths(main_data_dir).db_info_csv)
    #print('db_filename: ' + db_filename)
    #print('value in species column: ' + sp)

    # Check that the value retrieved makes sense.
    assert type(sp) is str, """There is more than one entry (row) for the
    filename %s in the file %s.""" % (db_filename, DataPaths(main_data_dir).db_info_csv)

    # Return the species name from the spreadsheet.
    return sp


def write_seqs_to_fasta(csv_file, 
                        output_dir, 
                        abbrev=False,
                        paralogue_names=False, 
                        only_descr=False, 
                        subseq=False,
                        all_hits=False, 
                        split_by_top_rev_srch_hit=None,
                        split_by_query_title=None,
                        split_to_query_fastas=False):
    """Takes a csv file (output by amoebae) and writes listed sequences to a
    fasta file.
    """
    # Define output subdirectories.
    #subdirpath = output_fasta.rsplit('.', 1)[0] + '_fasta_files'
    #subdirpath = output_dir + '_fasta_files'
    subdirpath = output_dir
    assert not os.path.isdir(subdirpath), """Output directory already
    exists: %s""" % subdirpath

    # Read csv file into a pandas dataframe.
    #df = pd.read_csv(csv_file, low_memory=False)
    df = pd.read_csv(csv_file, encoding='utf-8')

    # Do it differently depending on if all hits are to be included or not.
    seq_objs = []
    #if not all_hits: # ******
    # Get column header for column with Yes/No data.
    yes_no_col = None
    if all_hits:
        # Use all forward hits.
        print('\nWriting all forward search hits to FASTA files.')
        pass
    else:
        # Look a column to use to indicate whether a hit meets the criteria or
        # not.
        for header in df.columns:
            if header.startswith('Represents a potential paralogue'):
                yes_no_col = header
                print('\nWriting non-redundant sequences to FASTA files.')
                break
        if yes_no_col is None:
            for header in df.columns:
                if header == 'Collective interpretation of reverse search results':
                    yes_no_col = header
                    print("""\nWriting all sequences meeting forward and reverse
                    search criteria to FASTA files.""")
                    break
        assert yes_no_col is not None, """Could not determine which column to
        use to decide whether each sequence should be written."""

    # Get column header for column listing top reverse search hit id or
    # description if necessary.
    top_rev_hit_id_col = None
    if split_by_top_rev_srch_hit is not None:
        for header in df.columns:
            if header.startswith('Top reverse hit description'):
                if split_by_top_rev_srch_hit in header:
                    top_rev_hit_id_col = header
                    break
        assert top_rev_hit_id_col is not None, """Could not determine which column
        contains top reverse search hit sequence IDs."""

    # Iterate through rows and construct a list of relevant sequence
    # objects for each query title or for each top reverse search hit.
    query_title_seq_obj_dict = {}
    rev_srch_top_hit_seq_obj_dict = {}
    for index, row in df.iterrows():
        # Determine whether the forward hit sequence recorded in this row
        # should be included in the output file(s).
        use_this_fwd_hit = False
        if row['Forward hit rank'] == 0:
            pass
        elif yes_no_col is not None:
            result = row[yes_no_col]
            if result == 'Yes' or result == '+': 
                use_this_fwd_hit = True
        else:
            use_this_fwd_hit = True

        # Find the information for constructing a sequence object for this
        # forward hit.
        if use_this_fwd_hit:
            program = row['Forward search method']
            db_file = row['Subject database species (if applicable)']
            acc = row['Forward hit accession']
            description = None
            seq = None
            if not program.startswith('tblastn'):
                if subseq:
                    description = row['Forward hit description of subsequence(s) that align(s) to query']
                    seq = row['Forward hit subsequence(s) that align(s) to query'] 
                else:
                    description = row['Forward hit description']
                    seq = row['Forward hit sequence'] 
            else:
                # If the hit is a tblastn (nucleotide) hit, then always
                # use the subsequence.
                description = row['Forward hit description of subsequence(s) that align(s) to query']
                seq = row['Forward hit subsequence(s) that align(s) to query'] 


            # Use the species name if writing the file for input to
            # phylogenetics software.
            if abbrev:
                description = db_file.rsplit('.', 1)[0]
                if program.startswith('tblastn'):
                    hit_range = get_hit_range_from_hsp_ranges(row['Forward hit coordinates of subsequence(s) that align(s) to query'])
                    description = str(hit_range[0]) + '...' + str(hit_range[1]) + '__'\
                            + description

                # The output sequences may be aligned and visualized using
                # Mesquite, and if an alignment is saved in Mesquite, then
                # space characters will be randomly (not in all cases) replaced
                # with underscores, because Mesquite does not like nexus
                # alignments to have space characters in taxon names. This
                # prevents downstream parsing of the sequence names. So, if the
                # abbrev option is used, and the sequences are destined for an
                # alignment for phylogenetic analysis, etc. then it is best to
                # avoid spaces, and use a double-underscore to separate the
                # accession number from the species name (and convert that back
                # to a space in the final output tree file).
                description = description.replace(' ', '_')

            elif paralogue_names:
                col_name = None
                for x in df.columns:
                    if x.startswith('Paralogue name '):
                        col_name = x
                assert col_name is not None, """Could not find column with
                paralogue names."""
                description = df.at[index, col_name]

            elif only_descr:
                description = ''
                orig_descr = row['Forward hit description'].strip('\"')
                if len(orig_descr.split(' ', 1)) < 2:
                    acc = orig_descr.split('_', 1)[0] + '__' +\
                    orig_descr.split('_', 1)[1].replace(' ', '_')
                else:
                    acc = orig_descr.split(' ', 1)[0] + '__' +\
                    orig_descr.split(' ', 1)[1].replace(' ', '_')
            
            # Instantiate a sequence object for the sequence of interest to
            # be written to output.
            seq_obj = get_seq_obj_from_srch_res_csv_info(acc, description,
                    seq, abbrev, paralogue_names)

            if split_by_query_title or split_to_query_fastas:
                # Determine what the query title is for this hit.
                query_title = row['Query title']

                # Check that a sequence with the same description is not already present in
                # the dictionary value (list).
                found_seq_with_same_descr = False
                if query_title in query_title_seq_obj_dict.keys():
                    for seq_descr in [(x.id, x.description) for x in query_title_seq_obj_dict[query_title]]:
                        if seq_obj.id == seq_descr[0] and seq_obj.description == seq_descr[1]:
                            found_seq_with_same_descr = True 
                            break

                # Determine if the top hit id is already a key in the dict,
                # and add the sequence object to the list that is the value
                # for the appropriate key.
                if not found_seq_with_same_descr:
                    if query_title not in query_title_seq_obj_dict.keys():
                        query_title_seq_obj_dict[query_title] = [seq_obj]
                    else:
                        # Don't want to write more than one sequence to each
                        # file if splitting the hits into files to use as
                        # queries.
                        if not split_to_query_fastas:
                            query_title_seq_obj_dict[query_title] =\
                            query_title_seq_obj_dict[query_title] + [seq_obj]

            elif split_by_top_rev_srch_hit:
                # Define sequence ID for top reverse search hit.
                top_rev_hit_id = row[top_rev_hit_id_col]

                # Check that a sequence with the same description is not already present in
                # the dictionary value (list).
                found_seq_with_same_descr = False
                if top_rev_hit_id in rev_srch_top_hit_seq_obj_dict.keys():
                    for seq_descr in [(x.id, x.description) for x in
                            rev_srch_top_hit_seq_obj_dict[top_rev_hit_id]]:
                        if seq_obj.id == seq_descr[0] and seq_obj.description == seq_descr[1]:
                            found_seq_with_same_descr = True 
                            break

                # Determine if the top hit id is already a key in the dict,
                # and add the sequence object to the list that is the value
                # for the appropriate key.
                if not found_seq_with_same_descr:
                    if top_rev_hit_id not in rev_srch_top_hit_seq_obj_dict.keys():
                        rev_srch_top_hit_seq_obj_dict[top_rev_hit_id] = [seq_obj]
                    else:
                        rev_srch_top_hit_seq_obj_dict[top_rev_hit_id] =\
                        rev_srch_top_hit_seq_obj_dict[top_rev_hit_id] + [seq_obj]

    # Write sequences to different output files according to the
    # reverse search hit that they retrieved.
    os.mkdir(subdirpath)
    print('\nOutput fasta files:')
    if split_by_query_title:
        for seq_obj_list_key in query_title_seq_obj_dict.keys():
            specific_output_fasta_path = os.path.join(subdirpath,\
                    seq_obj_list_key.strip() + '_matches.fa')
            print(specific_output_fasta_path)
            with open(specific_output_fasta_path, 'w') as o:
                SeqIO.write(query_title_seq_obj_dict[seq_obj_list_key], o, 'fasta')

    elif split_by_top_rev_srch_hit:
        for seq_obj_list_key in rev_srch_top_hit_seq_obj_dict.keys():
            # Only write sequences that retrieved a hit in the reverse search
            # results.
            if seq_obj_list_key != '-':
                # Define output filename.
                #specific_output_fasta_name = output_fasta.rsplit('.', 1)[0] + '_' +\
                #seq_obj_list_key.strip().replace(' ', '_').replace('\"',\
                #    '').replace('[', '').replace(']', '') + '_matches.fa' 
                #specific_output_fasta_name =\
                #seq_obj_list_key.strip().replace(' ', '_').replace('\"',\
                #    '').replace('[', '').replace(']', '') + '_matches.fa' 

                specific_output_fasta_name = re.sub('[^0-9a-zA-Z]+', '_',
                        seq_obj_list_key.strip().strip('\"')) + '_matches.fa'

                #print(seq_obj_list_key)
                #print(seq_obj_list_key.strip().strip('\"'))
                #print(specific_output_fasta_name)

                # Define output filepath.
                specific_output_fasta_path = os.path.join(subdirpath, specific_output_fasta_name)
                print(specific_output_fasta_path)

                # Write sequences to file.
                with open(specific_output_fasta_path, 'w') as o:
                    SeqIO.write(rev_srch_top_hit_seq_obj_dict[seq_obj_list_key], o, 'fasta')

    elif split_to_query_fastas:
        # Note: Currently only writes one file for each query title...
        for seq_obj_list_key in query_title_seq_obj_dict.keys():
            species = 'SPECIES'
            #species = 'Ddiscoideum'
            seq_id = query_title_seq_obj_dict[seq_obj_list_key][0].id
            specific_output_fasta_path = os.path.join(subdirpath,\
                    seq_obj_list_key.strip() + '_' + species + '_' + seq_id + '.faa')
            print(specific_output_fasta_path)
            with open(specific_output_fasta_path, 'w') as o:
                SeqIO.write(query_title_seq_obj_dict[seq_obj_list_key], o, 'fasta')


def record_amoebae_info_in_log_file(commandline, outdir, start_time, end_time,
        output_timestamp, main_output_path):
    """Record info including full commandline command, git commit ID, time to
    run, etc in a log file.
    """
    # Identify log file.
    log_basename = '0_amoebae_log.txt'
    log_path = os.path.join(outdir, log_basename)
    if not os.path.isfile(log_path):
        with open(log_path, 'w') as o:
            o.write('Recording of AMOEBAE commands run with output to directory:\n%s\n\n\n'\
                    % os.path.basename(outdir))

    # Open log file and append information.
    with open(log_path, 'a') as o:
        # Record time.
        readable_time = time.ctime()
        o.write('\n' + readable_time + '\n\n')

        # Record command run.
        o.write('AMOEBAE run as:\n\n\t')
        o.write(commandline + '\n\n')

        # Record main output file.
        o.write('Main output:\n' + main_output_path + '\n\n')

        # Record output timestamp.
        o.write('Output timestamp: ' + output_timestamp + '\n\n')

        # Record the time it took for the process to complete.
        elapsed = end_time - start_time
        o.write('Total run time: ' + str(datetime.timedelta(seconds=elapsed)) + '\n\n')

        # Record git repository version information.
        script_dir = os.path.dirname(os.path.realpath(__file__)) 
        git_hash = str(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=script_dir).strip())
        git_branch = str(subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=script_dir).strip())  
        o.write('Git repository (code) version: ' + git_hash + ' (branch name: ' + git_branch + ')\n\n')

        # Record system information.
        o.write('System info: ' + str(platform.uname()) + '\n')

        # Write extra line break.
        o.write('\n')


def mask_alignment2(alignment):
    """Takes an alignment object and adds a 'MASK' sequence using certain
    criteria for inclusion of positions.
    """
    # Get length of sequences and number of sequences in alignment.
    seq_len = alignment.get_alignment_length()
    num_seqs = len(alignment)

    # get a list of columns as strings in the original alignment.
    columns = [alignment[:, col] for col in range(seq_len)] 

    # Iterate over columns and make a mask sequence to append.
    mask_seq = ''
    for col in columns:
        mask_char = apply_mask_criteria2(col) 
        mask_seq = mask_seq + mask_char

    # Generate a SeqRecord object with the mask_seq.
    #empty_mask_seq = Seq(mask_seq, IUPAC.protein)
    empty_mask_seq = Seq(mask_seq)
    empty_mask_rec = SeqRecord(empty_mask_seq, id='MASK', name='MASK')

    # Add the mask_seq sequence to the alignment.
    masked_alignment = alignment
    masked_alignment.append(empty_mask_rec)

    return masked_alignment


def mask_nex2(infilepath, outfilepath):
    """Takes a filepath and adds a MASK sequence.
    ***Includes all positions!!!
    """
    # Delete extra lines in input nexus file, if present, because biopython cannot
    # read nexus alignments with these extra lines.
    delete_extra_mesquite_lines(infilepath)

    # Define the name of the output file.
    outfilename = outfilepath

    with open(infilepath) as infh, open(outfilename, 'w') as o:
        # Check that the input file has the filename extension ".nex".
        assert infilepath.endswith('.nex'), "Error: Input file name must have the\
 extension '.nex'."

        # Read the alignment file.
        alignment = AlignIO.read(infh, 'nexus')

        masked_alignment = mask_alignment2(alignment)
        AlignIO.write(masked_alignment, o, 'nexus')


def get_seq_obj_from_srch_res_csv_info(acc, description, seq, abbrev=False,
        paralogue_names=False):
    """Takes components and assembles a Bio.Seq sequence object.
    If abbreviating, now puts the description before the assession, but used to be the other way
    around.
    """
    seq_obj = Seq(seq)
    seq_rec_obj = SeqRecord(seq_obj)

    if abbrev:
        #seq_rec_obj.id = acc + '__' + description
        seq_rec_obj.id = description + '__' + acc
        seq_rec_obj.description = ''

    elif paralogue_names:
        seq_rec_obj.id = description
        seq_rec_obj.description = ''

    else: 
        seq_rec_obj.id = acc
        seq_rec_obj.description = description

    return seq_rec_obj


def get_hit_range_from_hsp_ranges(subseq_coord):
    """Take subsequence coordinates as a string, and return a list containing
    the lowest and highest numbers.
    """
    #print('\nget_hit_range_from_hsp_ranges:')
    #print('subseq_coord:')
    #print(subseq_coord)
    # Compile regular expressions to identify relevant substrings.
    low = re.compile(r'\[\d+,')
    high = re.compile(r',\d+\]')

    # Compile lists of numbers.
    low_list = [x[1:-1] for x in low.findall(subseq_coord)]
    #print('low_list:')
    #print(low_list)
    high_list = [x[1:-1] for x in high.findall(subseq_coord)]
    #print('high_list:')
    #print(high_list)

    assert '-' not in subseq_coord, """There appears to be unexpected
    characters in the subsequence coordinates."""

    # Return a list containing the lowest and highest numbers.
    return [min(low_list), max(high_list)]


def apply_mask_criteria2(column):
    """Apply simple masking criteria to a single column, and return '-' if the
    column does not meet the criteria, and 'I' if it does.
    """
    # Return '-' by default.
    mask_char = '-'

    # Get column features.
    num_seqs = len(column)
    half_num_seqs = num_seqs / 2
    num_gaps_in_col = column.count('-')
    column_no_gaps = column.replace('-', '')

    # Check that the column is not entirely composed of gaps.
    #assert not column_no_gaps == '', "Error: Empty positions in input alignment."

    if column_no_gaps == '':
        return mask_char

    elif not column_no_gaps == '':

        most_common_residue = collections.Counter(column_no_gaps).most_common(1)[0]
        most_common_residue_count = most_common_residue[1]
        percent_identity = most_common_residue_count * 100 / num_seqs 

        # If less than half the sequences have a gap at this position of the
        # alignment, then include the position. 
        #if num_gaps_in_col < half_num_seqs:
        #    mask_char = 'I'

        #if num_gaps_in_col < (num_seqs * 0.30):
        #    mask_char = 'I'

        # If percent identity is at least 50, then include position.
        #if percent_identity >= 50:
        #    mask_char = 'I'

        # If the last row in the column has a residue and at least one of the
        # other rows in the column has a residue instead of a gap, then include
        # the column/position in the trimmed alignment.
        last_row_value = column[-1]
        original_row_values = column[0:-1]
        if len(original_row_values.replace('-', '')) >= 1:
            mask_char = 'I'
        #if last_row_value != '-' and len(column_no_gaps) > 1:
        #    mask_char = 'I'

        # If at least one original sequence has a residue at the position, then
        # it should be included.

        return mask_char


def disqualifying_string_in_file_basename(f, disqualifying_strings):
    """Returns True if one or more of the strings in the input list is in the
    basename of the input file path.
    """
    in_basename = False
    for s in disqualifying_strings:
        # Only look in the last 15 characters for the disqualifying strings.
        if s in os.path.basename(f)[-15:]:
            in_basename = True
            break
    return in_basename


def find_input_file_in_parent_directory(indp, extension, disqualifying_strings):
    """Given a directory path (e.g., path to a directory that was downloaded
    after running a phylogenetic analysis on CIPRES) and the filename extension
    that your target file will have (e.g., 'table' or 'nex'), return the path
    to a file in the parent directory of the input directory path that is the
    most similar to the input directory basename and has the specified
    extension.
    """
    # Identify parent directory.
    parent_dir = os.path.dirname(indp)

    # Identify all files in parent directory that have the given extension.
    all_rel_files_in_parent_dir = []
    for f in glob.glob(os.path.join(parent_dir, '*.' + extension)):
        if os.path.isfile(f):
            # Only look at files that don't have disqualifying strings (e.g.,
            # masked and trimmed alignments instead of original alignment).
            if not disqualifying_string_in_file_basename(f, disqualifying_strings):
                all_rel_files_in_parent_dir.append(f)

    # Check that there are any relevant files.
    assert len(all_rel_files_in_parent_dir) > 0, """No files in parent
    directory %s with filename extension %s.""" % (parent_dir, extension)

    # If there is only one file in the list, then it's the only option so you
    # could just return that...
    
    # Find the file with the longest identical prefix to the input directory
    # basename.
    indp_bn = os.path.basename(indp)
    # Loop over prefixes from longest to shortest.
    file_found = False
    for i in range(0, len(indp_bn) + 1)[::-1]:
        prefix = indp_bn[:i] 
        files_with_prefix = []
        for f in all_rel_files_in_parent_dir:
            if os.path.basename(f).startswith(prefix):
                files_with_prefix.append(f)
        if len(files_with_prefix) > 0:
            file_found = True
        # Return path if there is only one.
        if len(files_with_prefix) == 1:
            #print(files_with_prefix[0])
            return files_with_prefix[0]
        elif len(files_with_prefix) > 1:
            # Return the path to the file with the shortest basename.
            file_with_shortest_basename = sorted(files_with_prefix, key=lambda x: len(x))[0]
            return file_with_shortest_basename
        ## Report an error if there is more than one.
        #assert len(files_with_prefix) > 1, """More than one relevant file
        #found with extension %s in directory %s.""" % (extension, parent_dir)

    # Check that a file was found.
    assert file_found, """Could not identify an appropriate file with extension
    %s in directory %s.""" % (extension, parent_dir)


