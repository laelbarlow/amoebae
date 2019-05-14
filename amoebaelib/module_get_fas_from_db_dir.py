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
"""This module contains functions for retrieving fasta sequences from the
database directory used by run_BLAST_pipeline.py and HMMer_pipeline.py.

This will need to be updated to work with nucleotide data in the database
directory system.
"""

from Bio import SeqIO


def get_seq_obj_from_db_fasta(acc_list, fa_path):
    """Takes a list of accessions and a fasta file path, and returns seq
    objects for corresponding sequences.
    """
    # Construct a list of sequences (objects) that correspond to accessions.
    # Need a way to check them off the list as it goes, and break the loop when
    # they are all identified.
    record_list = []
    acc_not_found_list = acc_list
    with open(fa_path) as fa:
        rec_num = 0
        for record in SeqIO.parse(fa, 'fasta'):
            rec_num += 1
            for acc in acc_list:
                if acc == record.id:
                    record_list.append(record)
                    acc_not_found_list.remove(acc)

    # Check whether any sequences not identified.
    if len(acc_not_found_list) > 0:
        print('\nThe sequences corresponding to the following accessions were \
not identified:')
        for acc in acc_not_found_list: print(acc)
        print('\n')

    return record_list


def get_fas_from_db_dir(db_name, acc_list, dbdirpath, prot_name=None):
    """Takes a dir path (e.g., /Users/Lael/Documents/Data/Genomes_2016), a
    database name (e.g., Gintestinalis), and one or more accession numbers, and
    returns the corresponding fasta sequence(s) as a string.
    """
    # print('db_name: ' + db_name)
    # print('acc_list: ' + str(acc_list))
    # print('dbdirpath: ' + dbdirpath)
    # print('prot_name: ' + str(prot_name))
    fa_paths = glob.glob(os.path.join(dbdirpath, db_name + '/db_*/*.fa')) #'/db_Prot/*.fa'))
    if len(fa_paths) > 1:
        print('\nProblem: more than one fasta file identified in the database\
path.')
    elif len(fa_paths) < 1:
        print('\nProblem: No fasta file identified in the database path.')

    # Check that fasta file paths were identified.
    assert len(fa_paths) >= 1, "Error: Could not identify any fasta file paths."

    # Define fasta file path.
    fa_path = fa_paths[0]

    # Get sequence objects.
    record_list = get_seq_obj_from_db_fasta(acc_list, fa_path)

    # Get the necessary info from the list of objects.
    fas_string = ''
    for record in record_list:
        fas_string = str(fas_string + get_abbrev_fa_record_text_from_obj(record,
                db_name, prot_name))

    return fas_string
