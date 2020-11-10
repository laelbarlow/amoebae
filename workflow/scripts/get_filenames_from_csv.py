#!/bin/env python3
"""Module defining functions for retrieving FASTA filenames from input CSV
files.
"""
import sys
import os
import requests

def get_query_file_dict_from_csv(input_csv_path):
    """Parse input CSV file to define a dictionary with desination filenames as
    values and NCBI sequence accessions keys.
    """
    query_dict = {}
    with open(input_csv_path) as infh:
        # Iterate over lines in CSV file.
        lnum = 0
        for i in infh:
            lnum += 1
            # Ignore irrelevant lines.
            if not i.startswith('Filename') \
                    and not i.startswith('\n') \
                    and not i.startswith(','):
                # Parse info in line.
                spliti = i.strip().split(',')
                filename = spliti[0]
                seq_id = str(spliti[1])

                # Check that the filename extension is acceptable.
                assert filename.rsplit('.', 1)[1] in ['faa', 'afaa', 'fna'],\
                """Error: Input file %s line %s: Data filename must have extension
                faa, afaa, or fna not %s.""" % (input_csv_path, str(lnum),\
                        filename.rsplit('.', 1)[1])

                # Add info to query dictionary.
                if filename not in query_dict.keys():
                    if seq_id == "":
                        query_dict[filename] = []
                    else:
                        query_dict[filename] = [seq_id]
                else:
                    if seq_id == "":
                        pass
                    else:
                        query_dict[filename] = query_dict[filename] + [seq_id]

    # Return the dict.
    return query_dict


def get_query_filenames_from_csv(input_csv_path, source_type):
    """Parse input CSV file and return a list of input FASTA filenames.
    """
    # Check that source_type is valid.
    assert source_type in ['all', 'ncbi', 'local']

    query_dict = {}
    with open(input_csv_path) as infh:
        # Iterate over lines in CSV file.
        lnum = 0
        for i in infh:
            lnum += 1
            # Ignore irrelevant lines.
            if not i.startswith('Filename') \
                    and not i.startswith('\n') \
                    and not i.startswith(','):
                # Parse info in line.
                spliti = i.strip().split(',')
                filename = spliti[0]
                seq_id = str(spliti[1].strip())

                # Check that the filename extension is acceptable.
                assert filename.rsplit('.', 1)[1] in ['faa', 'afaa', 'fna'],\
                """Error: Input file %s line %s: Data filename must have extension
                faa, afaa, or fna not %s.""" % (input_csv_path, str(lnum), \
                        filename.rsplit('.', 1)[1])

                # Decide whether to include.
                include = False
                if source_type == 'all':
                    include = True
                elif source_type == 'ncbi' and seq_id != '':
                    include = True
                elif source_type == 'local' and seq_id == '':
                    include = True

                if include:
                    # Add info to query dictionary.
                    if filename not in query_dict.keys():
                        query_dict[filename] = [seq_id]
                    else:
                        query_dict[filename] = query_dict[filename] + [seq_id]

    # Generate lists of filenames for files with single or multiple seqs.
    single_fasta_filenames = []
    multi_fasta_filenames = []
    for filename in query_dict.keys():
        num_ids = len(query_dict[filename])
        if source_type == 'local' and filename.endswith('.afaa'):
            multi_fasta_filenames.append(filename)
        elif num_ids == 1:
            single_fasta_filenames.append(filename)
        else:
            multi_fasta_filenames.append(filename)

    # Check that each list is nonredundant.
    assert len(single_fasta_filenames) == len(list(set(single_fasta_filenames)))
    assert len(multi_fasta_filenames) == len(list(set(multi_fasta_filenames)))

    # Check that the two output lists are non-overlapping.
    for i in single_fasta_filenames:
        for j in multi_fasta_filenames:
            assert i != j, """Error: %s identified as both a single FASTA query
            and as a multi FASTA query."""

    # Return the list.
    return (single_fasta_filenames, multi_fasta_filenames)


def get_database_filenames_from_csv(input_csv_path, source_type):
    """Parse input CSV file and return a list of input FASTA filenames.
    """
    # Check that source_type is valid.
    assert source_type in ['all', 'ncbi', 'local']

    filenames = []
    with open(input_csv_path) as infh:
        # Iterate over lines in CSV file.
        lnum = 0
        for i in infh:
            lnum += 1
            # Ignore irrelevant lines.
            if not i.startswith('Filename') \
                    and not i.startswith('\n') \
                    and not i.startswith(','):
                # Parse info in line.
                spliti = i.strip().split(',')
                filename = spliti[0]
                location = spliti[4].strip()

                # Check that the filename extension is acceptable.
                assert filename.rsplit('.', 1)[1] in ['faa', 'fna', 'gff3', 'hmmdb'],\
                """Error: Input file %s line %s: Data filename must have extension
                .faa, .fna, or .gff3, not %s: %s""" % (input_csv_path, str(lnum), filename)

                # Decide whether to include.
                include = False
                if source_type == 'all':
                    include = True
                elif source_type == 'ncbi' and location != '':
                    include = True
                elif source_type == 'local' and location == '':
                    include = True

                if include:
                    # Add filename to list.
                    filenames.append(filename)

    # Return the list.
    return filenames

