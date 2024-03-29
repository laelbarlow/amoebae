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
similarity searches, but do not require installation of pandas library.
"""
from datapaths import DataPaths
import os
import subprocess
import datetime 
import time
import genetic_code_config


def get_query_subdir(outdir):
    """Takes a reverse-search output directory and returns a path for a reverse
    search query subdirectory.
    """
    return os.path.join(outdir, '0_rev_srch_queries')


def get_query_list_from_file(infp):
    """Get query list from file.
    """
    query_list = []
    with open(infp) as infh:
        for i in infh:
            if not i.startswith('#') and not i.startswith('\n'):
                query_list.append(i.strip())
    return query_list


def get_db_list_from_file(infp):
    """Get database list from file.
    """
    db_list = []
    with open(infp) as infh:
        for i in infh:
            if not i.startswith('#') and not i.startswith('\n'):
                db_list.append(i.strip())
    #print('db_list')
    #print(db_list)
    return db_list


def get_out_bash_path(outdir):
    """Takes an output directory that search results are written to, and
    returns a path for a file to write query file list to.
    """
    return os.path.join(outdir, '0_run_searches.sh')


def get_out_query_list_path(outdir):
    """Takes an output directory that search results are written to, and
    returns a path for a file to write query file list to.
    """
    return os.path.join(outdir, '0_queries.txt')


def get_out_db_list_path(outdir):
    """Takes an output directory that search results are written to, and
    returns a path for a file to write database file list to.
    """
    return os.path.join(outdir, '0_databases.txt')


def determine_search_method(query_exten, db_exten):
    """Determine which similarity search method to use based on the filename
    extensions of the query and database files.

    Eventually integrate hhsuite for HMM-HMM comparison.

    Assuming blastp rather than phmmer for now.
    """
    # Import default methods, for cases where several could be used, from the
    # DataPaths(main_data_dir) module.
    #faa_in_faa = DataPaths(main_data_dir).faa_in_faa # phmmer, blastp.

    # Decide on method.
    method = None

    if query_exten == 'faa': 
        if db_exten == 'faa':
            method = 'blastp'
        elif db_exten == 'fna':
            method = 'tblastn'
        elif db_exten == 'hmmdb':
            method = 'hmmscan'

    elif query_exten == 'fna': 
        if db_exten == 'faa':
            method = 'blastx'
        elif db_exten == 'fna':
            method = 'blastn'
        elif db_exten == 'hmmdb':
            pass

    elif query_exten == 'afaa': 
        if db_exten == 'faa':
            method = 'hmmsearch'
        elif db_exten == 'fna':
            pass
        elif db_exten == 'hmmdb':
            pass

    elif query_exten == 'afna': 
        if db_exten == 'faa':
            pass
        elif db_exten == 'fna':
            method = 'nhmmer'
        elif db_exten == 'hmmdb':
            pass

    assert method is not None, """Error: Could not determine which similarity
    searching method to apply for query with extension '%s' and database with
    extension '%s'""" % (query_exten, db_exten)

    # Return chosen method.
    return method


def search_result_filepath(query_filename, db_filename, dirpath):
    """Returns a filepath for a search result file given a directory path, query
    filename, and database filename. This is for both naming and identifying
    filepaths after searches have been performed. 
    """
    #outfile = os.path.join(dirpath, os.path.basename(query_filename).rsplit('.', 1)[0] + '__' +\
    #        db_filename.rsplit('.', 1)[0] + '_srch_out.txt')
    outfile = os.path.join(dirpath, os.path.basename(query_filename).rsplit('.', 1)[0] + '__' +\
            db_filename.rsplit('.', 1)[0] + '_' + db_filename.rsplit('.', 1)[1] + '_srch_out.txt')
    return outfile


def get_out_hmm_path(new_afa_path):
    """Define an hmm file path for a given aligned fasta file path.
    """
    new_exten = None
    old_exten = new_afa_path.rsplit('.', 1)[1]
    if old_exten == 'afaa':
        new_exten = 'hmm'
    elif old_exten == 'afna':
        new_exten = 'nhmm'

    # Check that it worked.
    assert new_exten is not None, """Error: HMM extension could not be
    determined for input file: %s""" % query_file

    # Get new path.
    hmm_file_path = new_afa_path.rsplit('.', 1)[0] + '.' + new_exten

    # Return the new path.
    return hmm_file_path


def run_any_search(queryfile,
                   dbfile,
                   outfile,
                   blast_report_evalue_cutoff,
                   blast_max_target_seqs,
                   hmmer_report_evalue_cutoff,
                   hmmer_report_score_cutoff,
                   num_threads_similarity_searching
                   ):
    """Run similarity search.

    Import info from DataPaths(main_data_dir).py to specify options for running external
    software? For example, number of threads.

    ***Need to refactor functions in the module_nhmmer_search module?
    """
    # Determine method to use based on the input file types.
    query_exten = queryfile.rsplit('.', 1)[1]
    dbfile_exten = dbfile.rsplit('.', 1)[1]
    method = determine_search_method(query_exten, dbfile_exten)

    # Get version number for software.
    #version = 'version'
    #version = get_search_software_version(method)

    # Get relevant DataPaths(main_data_dir).
    # Get cutoffs for recording hits.
    blast_evalcut = str(blast_report_evalue_cutoff)
    blast_max_target_seqs = str(blast_max_target_seqs)
    hmmer_evalcut = str(hmmer_report_evalue_cutoff)
    hmmer_scorecut = str(hmmer_report_score_cutoff)
    # Get number of threads to use.
    num_threads = str(num_threads_similarity_searching)

    # Construct search command.
    run_command = []
    if method == 'blastp':
        run_command = [method, '-query', queryfile, '-db', dbfile, '-out',
                outfile, '-num_threads', num_threads, '-outfmt', '5', '-evalue',
                blast_evalcut, '-max_target_seqs',
                blast_max_target_seqs]
    elif method == 'tblastn':
        # Set the genetic code.
        tblastn_ncbi_gen_code = genetic_code_config.ncbi_genetic_code_number

        run_command = [method, '-query', queryfile, '-db', dbfile, '-out',
                outfile, '-num_threads', num_threads, '-outfmt', '5', '-evalue',
                blast_evalcut, '-max_target_seqs',
                blast_max_target_seqs, '-db_gencode',
                tblastn_ncbi_gen_code]
    elif method == 'blastx':
        run_command = [method, '-query', queryfile, '-db', dbfile, '-out',
                outfile, '-num_threads', num_threads, '-outfmt', '5', '-evalue',
                blast_evalcut, '-max_target_seqs',
                blast_max_target_seqs]
    elif method == 'blastn':
        run_command = [method, '-query', queryfile, '-db', dbfile, '-out',
                outfile, '-num_threads', num_threads, '-outfmt', '5', '-evalue',
                blast_evalcut, '-max_target_seqs',
                blast_max_target_seqs]
    elif method == 'hmmsearch':
        # Use HMM file rather than '.afaa' file.
        actual_queryfile = get_out_hmm_path(queryfile)
        #run_command = [method, "-T", hmmer_scorecut, "--cpu", num_threads,
        #        '--tblout', outfile, actual_queryfile, dbfile]
        run_command = [method, "-T", hmmer_scorecut, "--cpu", num_threads,
                '-o', outfile, actual_queryfile, dbfile]
    elif method == 'hmmscan':
        #run_command = [method, "-T", hmmer_scorecut, "--cpu", num_threads,
        #        '--tblout', outfile, queryfile, dbfile]
        run_command = [method, "-T", hmmer_scorecut, "--cpu", num_threads,
                '-o', outfile, dbfile, queryfile]
    elif method == 'nhmmer':
        # Use HMM file rather than '.afna' file.
        actual_queryfile = get_out_hmm_path(queryfile)
        #run_command = [method, "-T", hmmer_scorecut, "--cpu", num_threads,
        #        '--tblout', outfile, actual_queryfile, dbfile]
        run_command = [method, "-T", hmmer_scorecut, "--cpu", num_threads,
                '-o', outfile, actual_queryfile, dbfile]

    # Prepend program name with directory path if necesssary and specified in
    # the DataPaths(main_data_dir).py file (this is a work-around for a particular remote
    # server).
    #server_program_dirpath = DataPaths(main_data_dir).server_program_dirpath
    #run_command = [os.path.join(server_program_dirpath, run_command[0])] +\
    #run_command[1:]

    # Run command.
    subprocess.call(run_command)

    # Check that the output file is not empty.
    assert os.path.getsize(outfile) != 0, """Error: Search output file is
    empty. This may be due to low memory available on this system.\nEmpty file:
    %s""" % outfile

    # Return string with command used to run search.
    #search_descr =  method + ' (' + version + ')' + ' run with command:\n\t' +\
    #        ' '.join(run_command) + '\n'
    search_descr = ' '.join(run_command)
    return search_descr


def run_all_searches(query_file_list,
                     db_file_list,
                     outdir,
                     blast_report_evalue_cutoff,
                     blast_max_target_seqs,
                     hmmer_report_evalue_cutoff,
                     hmmer_report_score_cutoff,
                     num_threads_similarity_searching,
                     main_data_dir,
                     query_dir=None
                     ):
    """Search with every query file in a given list into every database file in
    another given list using appropriate methods.
    """
    # Current time.
    start_time = time.time()

    # Get query and database directories from DataPaths(main_data_dir).
    if query_dir == None:
        query_dir = DataPaths(main_data_dir).querydirpath
    db_dir = DataPaths(main_data_dir).dbdirpath

    # Write a query file list file to output directory.
    out_query_file = get_out_query_list_path(outdir)
    with open(out_query_file, 'w') as o:
        for q in query_file_list:
            o.write(q + '\n')

    # Write a database file list file to output directory.
    out_db_file = get_out_db_list_path(outdir)
    with open(out_db_file, 'w') as o:
        for d in db_file_list:
            o.write(d + '\n')

    # Create a log file.
    logfile = os.path.join(outdir, '0_search_log.txt')

    # Loop over each query-database pair.
    with open(logfile, 'w') as o:
        srch_num = 0
        # Loop over query files.
        for q in query_file_list:
            # Loop over database files.
            for d in db_file_list:

                # Check that database file is a single FASTA file.
                assert '.faa' not in d.rsplit('.', 1)[0], """The database file
                name %s does not appear to be formatted correctly. This may
                have resulted from a file parsing error.""" % d

                if q.rsplit('.', 1)[1] == 'afaa' and d.rsplit('.', 1)[1] == 'fna':
                    warning_text = """\nWARNING: Not searching with profile query %s
                    in nucleotide data %s\n\n""" % (q, d)
                    print(warning_text)
                    o.write(warning_text)
                else:
                    srch_num += 1

                    # Get name of output file.
                    #print(outdir)
                    outfile = search_result_filepath(q, d, outdir)
                    #print(outfile)

                    # Get full filepaths, and verify existence.
                    qfull = None
                    #if os.path.isfile(q):
                    #    qfull = q
                    #else:
                    #    qfull = os.path.join(query_dir, q)
                    qfull = os.path.join(query_dir, q)
                    assert os.path.isfile(qfull), """Specified query file path is
                    not a file: %s""" % qfull

                    dfull = None
                    #if os.path.isfile(d):
                    #    dfull = d
                    #else:
                    #    dfull = os.path.join(db_dir, d)
                    dfull = os.path.join(db_dir, d)
                    assert os.path.isfile(dfull), """Specified database file path
                    is not a file: %s\n Please ensure that a FASTA file with
                    the filename %s exists in the input data.""" % (dfull, d)

                    # Search start time.
                    search_start_time = time.time()

                    # Run the similarity search and get a description of the search
                    # command (write to a log file?).
                    command_descr = run_any_search(qfull, dfull, outfile,
                                                   blast_report_evalue_cutoff,
                                                   blast_max_target_seqs,
                                                   hmmer_report_evalue_cutoff,
                                                   hmmer_report_score_cutoff,
                                                   num_threads_similarity_searching)

                    # Write description of search to log file.
                    o.write(command_descr + '\n')

                    # End time.
                    search_end_time = time.time()
                    # Record time elapsed.
                    search_elapsed = search_end_time - search_start_time
                    o.write('Run time: ' + str(datetime.timedelta(seconds=search_elapsed)) + '\n')

        # End time.
        end_time = time.time()
        # Record time elapsed.
        elapsed = end_time - start_time
        o.write('Total run time: ' + str(datetime.timedelta(seconds=elapsed)) + '\n')
        
