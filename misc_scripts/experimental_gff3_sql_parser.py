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
"""This script is for testing the gffutils package to figure out how to extract
information from sql files with information from gff3 files.
"""
import os
import sys
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import gffutils
from module_search_scaffolds import check_if_two_hsp_ranges_overlap


command_line_list = sys.argv
sqldb = str(command_line_list[1]) # Path to sql file
pid1 = str(command_line_list[2]) # Protein accession 1
pid2 = str(command_line_list[3]) # Protein accession 2



def cds_record_matches_prot_id(cds_obj, prot_id):
    """Take a CDS record object (constructed using gffutils) and a protein
    sequence ID, and return true if the CDS record is linked to the protein
    sequence ID (the CDS is part of the sequence encoding that peptide).
    """
    if cds_obj.id[4:] == prot_id or \
       cds_obj.attributes['Name'][0] == prot_id or \
       cds_obj.attributes['protein_id'][0] == prot_id:
        return True
    else:
        return False


def two_proteins_same_gene_in_gff3_2(sql_database, prot_id_1, prot_id_2):
    """Take a path to a gffutils SQL database (based on a gff3 file), the
    sequence id for the relevant sequence referred to in the gff file, and a
    range (coordinates in the form of a list object with two items). And,
    return IDs for any genes that fall within that range on the sequence. This
    is for determining whether two protein hits with different IDs correspond
    to the same gene.
    """
    # Parse the SQL database that was generated by gffutils for the relevant
    # GFF3 file.
    db = gffutils.FeatureDB(sql_database, keep_order=True)

    # Initiate variable to store True if the two proteins are encoded at the
    # same gene locus.
    both_prot_encoded_at_same_locus = False

    # Initiate variables to store nucleotide sequence IDs corresponding to the
    # first coding sequences for each protein sequence.
    nucl_seq_id_for_cds_1 = None
    nucl_seq_id_for_cds_2 = None

    # Initiate a list of nucleotide sequence IDs.
    nucl_seq_ids = []

    # Initiate variables to store first CDS objects for each protein.
    prot_id_1_cds_1 = None
    prot_id_2_cds_1 = None

    # Initiate lists of all the CDS record start and stop positions for CDS
    # records associated with the two protein sequences.
    prot_id_1_cds_starts_stops = []
    prot_id_2_cds_starts_stops = []

    # Loop over all coding sequence (CDS) records in the database.
    stop_looping_through_cds_records = False
    for cds in db.features_of_type('CDS', order_by='start'):
        # Loop over the two protein IDs.
        for prot_id in [prot_id_1, prot_id_2]:
            # Determine whether the record matches one of the protein IDs. 
            if cds_record_matches_prot_id(cds, prot_id):
                # Set trigger to break the loop if a relevant sequence is
                # found, and then a different nucleotide sequence is
                # considered (only coding sequences in the same nucleotide
                # sequence should be considered).
                if len(nucl_seq_ids) == 0:
                    nucl_seq_ids = [cds.seqid]
                else:
                    if cds.seqid != nucl_seq_ids[0]:
                        stop_looping_through_cds_records = True

                # Store CDS objects if they are the first for either of the two
                # proteins.
                if prot_id_1_cds_1 is None and prot_id == prot_id_1:
                    prot_id_1_cds_1 = cds
                elif prot_id_2_cds_1 is None and prot_id == prot_id_2:
                    prot_id_2_cds_1 = cds 

                # Record start and stop positions.
                if prot_id == prot_id_1:
                    prot_id_1_cds_starts_stops.append(cds.start)
                    prot_id_1_cds_starts_stops.append(cds.end)
                elif prot_id == prot_id_2:
                    prot_id_2_cds_starts_stops.append(cds.start)
                    prot_id_2_cds_starts_stops.append(cds.end)

        # See whether CDS for the two proteins correspond to the same gene
        # record (this information may not be available in the database). If
        # so, then there is no need to continue looping over the CDS records.
        if prot_id_1_cds_1 is not None and prot_id_2_cds_1 is not None:
            # Identify gene IDs, if present.
            gene_ids = []
            for cds_obj in [prot_id_1_cds_1, prot_id_2_cds_1]:
                parents = db.parents(cds, featuretype='gene')
                parents_list = list(parents)
                if len(parents_list) != 1:
                    # Break loop over the CDS objects if the first object does
                    # not have a gene record associated with it.
                    break
                else:
                    gene_ids.append(parents_list[0].id)
            # Trigger breaking the loop over CDS records if the two CDS records
            # are associated with the same gene record. 
            if len(gene_ids) == 2 and gene_ids[0] == gene_ids[1]:
                stop_looping_through_cds_records = True
                # In this case, both proteins are encoded at the same gene
                # locus.
                both_prot_encoded_at_same_locus = True

        # Stop looping over CDS records, if triggered.
        if stop_looping_through_cds_records:
            break

    # Check that at least one CDS record linked to at least one of the protein
    # sequence IDs was identified.
    assert prot_id_1_cds_1 is not None or prot_id_2_cds_1 is not None, """No
    CDS records corresponding to protein sequence IDs %s or %s could be
    identified in the gffutils SQL database file %s. This suggests that either
    the database does not contain relevant information, or there was an error
    in retrieval of the information.""" % (prot_id_1, prot_id_2, sql_database)

    # If necessary, determine whether the coding sequences corresponding to the
    # two proteins overlap on the nucleotide sequence.
    if not both_prot_encoded_at_same_locus:
        if len(prot_id_1_cds_starts_stops) > 0 and \
                len(prot_id_2_cds_starts_stops) > 0:
            # Get ranges for both CDS.
            prot_id_1_cds_range = [min(prot_id_1_cds_starts_stops),
                                   max(prot_id_1_cds_starts_stops)]
            prot_id_2_cds_range = [min(prot_id_2_cds_starts_stops),
                                   max(prot_id_2_cds_starts_stops)]
            # Determine if the CDS record ranges overlap.
            overlap = check_if_two_hsp_ranges_overlap([prot_id_1_cds_range,
                                                       prot_id_2_cds_range])
            if overlap:
                both_prot_encoded_at_same_locus = True


    # Return True if the two proteins have been determined to be encoded at the
    # same gene locus.
    if both_prot_encoded_at_same_locus:
        return True
    else:
        return False




if __name__ == '__main__':
    x = two_proteins_same_gene_in_gff3_2(sqldb, pid1, pid2)
    print(x)