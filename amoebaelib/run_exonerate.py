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
"""Code for running exonerate as a subprocess, and extracting relevant
information from the output.

exonerate can be installed via Anaconda:
    conda install -c bioconda exonerate

This module can be run as a script:

    run_exonerate.py QUERY.faa DATABASE.fna ACCESSION SUBSEQSTART SUBSEQENCE GENETICCODE

    For example: 

    run_exonerate.py Athaliana_COPIGamma.faa PabiesSc.fna MA_93046 37000 39000 1
"""
# Import modules.
import sys
import os
sys.path.append(os.path.dirname(sys.path[0]))
from datapaths import DataPaths
#from get_fas_from_db_dir import get_seq_obj_from_db_fasta
from amoebae_m import get_seqs_from_fasta_db
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO
import subprocess


# Define classes.

class ExonerateLocusResult:
    """Parses relevant information as needed from a file containing output of
    running exonerate on a nucleotide sequence containing potential exons of
    one gene locus.
    """
    def __init__(self,
                 exonerate_output_file_path,
                 subject_seq_fasta,
                 position_of_subject_seq_start_in_original,
                 genetic_code_number,
                 rough_start,
                 rough_end
                 ):
        # Check that there is a usable nucleotide sequence in the file, and
        # extract a list of FASTA sequence strings for the HSPs in the order
        # that they appear in the file.
        with open(exonerate_output_file_path) as exonerate_outputfh:
            # Get file contents as a string.
            file_text_string = exonerate_outputfh.read()

            # Check that there is at least one FASTA sequence in the file.
            fasta_string = None
            try:
                fasta_string = file_text_string.split('\n------------')[4].strip()
            except:
                print("Could not identify FASTA sequence in exonerate output file %s" % exonerate_output_file_path)
            if fasta_string is None:
                # If no FASTA sequence can be found in the exonerate output file,
                # then nothing was identified by exonerate, despite TBLASTN
                # identifying an HSP in the sequence. In this case, simply set both
                # attributes (the sequence and location string) to None.
                self.location_string = None
                self.seq_record = None
            else:
                # Collect information regarding the sequence identified by
                # exonerate.

                # Check that the FASTA string is in FASTA format.
                assert fasta_string.startswith('>')

                # Compile a list of FASTA strings from the exonerate output
                # file.
                fasta_strings = []
                for i in file_text_string.split('\n------------')[4:]:
                    if i.strip().startswith('>'):
                        fasta_strings.append(i.strip())

                # Use Biopython to parse the exonerate output file, and identify the
                # HSP to use.
                hsp_to_use = None
                fasta_string_to_use = ''
                #with open(exonerate_output_file_path) as exonerate_outputfh:
                try:
                    #parsed_exonerate = list(SearchIO.parse(exonerate_output_file_path, 'exonerate-text'))
                    parsed_exonerate =\
                    list(SearchIO.parse(exonerate_output_file_path,\
                        'exonerate-vulgar')) # Use the vulgar string to avoid erroneous start and stop positions.
                    #print(parsed_exonerate)
                    #for hsp in sorted(list(parsed_exonerate), key=lambda x: x.score, reverse=True):

                    # Make a list of HSPs.
                    hsp_list = []
                    for x in parsed_exonerate:
                        for y in x:
                            for z in y:
                                hsp_list.append(z)

                    # Note: HSPs must be in the same order as the FASTA
                    # sequences (in order of appearance in the exonerate output
                    # file).

                    assert len(hsp_list) == len(fasta_strings), """Different
                    number of HSPs and FASTA sequences identified in the exonerate
                    output file %s""" % exonerate_output_file_path
                    for hsp, fasta_string in zip(hsp_list, fasta_strings):
                        if len(fasta_string) > len(fasta_string_to_use):
                            fasta_string_to_use = fasta_string
                            hsp_to_use = hsp
                    assert hsp_to_use is not None
                except:
                    print("""\n\nWarning: This exonerate output file could
                    not be parsed using Bio.SearchIO: %s""" % exonerate_output_file_path)
                    # Just use a rough estimate of location in the sequnce
                    # header. Maybe Biopython will be updated at some point... 

                    # Choose the longest FASTA sequence in the exonerate
                    # output.
                    for fasta_string in fasta_strings:
                        if len(fasta_string) > len(fasta_string_to_use):
                            fasta_string_to_use = fasta_string

                assert fasta_string_to_use is not ''

                # Check that the start and end positions are not negative
                # numbers.
                assert rough_start >= 0
                assert rough_end >= 0

                location_string = 'approx[' + str(rough_start) + ',' + str(rough_end) + ']'
                if not hsp_to_use is None:
                    #print('\n')
                    #print(exonerate_output_file_path)
                    #print(fasta_string_to_use)

                    # Define the locations of the fragments within the selected HSP for
                    # recording in the output summary and sequence header.
                    locations = []
                    concatenated_fragment_seqs = ''
                    additional_seq_len = position_of_subject_seq_start_in_original - 1
                    #for hspfragments in hsp_to_use:
                    #    #print('\n')
                    #    #print(hspfragments)
                    for fragment in hsp_to_use:
                        #print(fragment)
                        #print(fragment.hit.seq)

                        #location = '[' + str(fragment.hit_start + 1) + ',' + str(fragment.hit_end) + ']'
                        #location = [min([fragment.hit_start + 1, fragment.hit_end]), max([fragment.hit_start + 1, fragment.hit_end])]

                        # Check that the start and end positions are not
                        # negative numbers.
                        assert fragment.hit_range[0] >= 0
                        assert fragment.hit_range[1] >= 0

                        location = [fragment.hit_range[0] + 1 +\
                                additional_seq_len, fragment.hit_range[1] + additional_seq_len]
                        locations.append(location)
                        
                            #concatenated_fragment_seqs = concatenated_fragment_seqs + fragment.hit

                    #print(concatenated_fragment_seqs.seq)
                    #print(transl_exonerate_seq_obj.seq)
                    #assert str(concatenated_fragment_seqs.seq) == str(transl_exonerate_seq_obj.seq)

                    # Sort locations from 5' to 3'. 
                    locations.sort(key=lambda x: x[0])
                    # Construct a string with locations.
                    location_string = '[' + ','.join([str(x).replace(' ', '') for x in locations]) + ']'
                    #print(location_string)

                # Define location description strings as an attribute of instances of
                # this class.
                self.location_string = location_string

                # Write FASTA sequence to a file.
                temp_fastafp = exonerate_output_file_path.rsplit('.', 1) [0] + '_seq.fna'
                with open(temp_fastafp, 'w') as o:
                    ## Get file contents as a string.
                    #file_text_string = exonerate_outputfh.read()

                    ## Check that there is at least one FASTA sequence in the file.
                    #fasta_string = None
                    #seq_position = XXXX + XX
                    #try:
                    #    fasta_string = file_text_string.split('\n------------')[4].strip()
                    #except:
                    #    print("Could not parse exonerate output file %s" % exonerate_output_file_path)
                    #assert fasta_string is not None
                    #assert fasta_string.startswith('>')

                    #print('\n\n\n')
                    #print('Sequence identified by exonerate in file %s:' % exonerate_output_file_path)
                    #print(fasta_string)
                    #print('\n\n\n')

                    o.write(fasta_string_to_use)

                # Parse the nucleotide sequence file.
                exonerate_seq_obj = SeqIO.read(temp_fastafp, 'fasta')
                #print(exonerate_seq_obj.description)

                # Change alphabet type for the sequence to ambiguous DNA, so that it
                # can be translated (if it is translatable).
                exonerate_seq_obj.seq.alphabet = IUPAC.ambiguous_dna

                ## Code block for debugging issues with translating sequences using
                ## Biopython.
                #print("\nexonerate_seq_obj:")
                #print(exonerate_seq_obj)
                #print("\nexonerate_seq_obj (" + exonerate_seq_obj.id + ' '\
                #        + location_string + ":")
                #print(exonerate_seq_obj)
                #print('exonerate_seq_obj.seq:')
                #print(exonerate_seq_obj.seq)
                #print("\nexonerate_seq_obj.seq.alphabet:")
                #print(exonerate_seq_obj.seq.alphabet)

                # Translate the nucleotide sequence with the appropriate genetic code,
                # in the appropriate strand.
                transl_exonerate_seq_obj = exonerate_seq_obj.translate(table=genetic_code_number)
                transl_exonerate_seq_obj.id = exonerate_seq_obj.id
                transl_exonerate_seq_obj.description = exonerate_seq_obj.description
                #print(transl_exonerate_seq_obj)

                #print('\nCoding sequence identified by exonerate:')
                #print(exonerate_seq_obj.seq)
                #print('\nTranslation of the coding sequence:')
                #print(transl_exonerate_seq_obj.seq)
                #print('\n')
                
                # Add locations to header.
                transl_exonerate_seq_obj.description =\
                transl_exonerate_seq_obj.description + ' ' + location_string
                #print(transl_exonerate_seq_obj.description)

                # Define sequence object as an attribute of instances of this class.
                self.seq_record = transl_exonerate_seq_obj



# Define functions.

def get_subseq_from_nucl(subj_seq_id, 
                         db_filepath,
                         target_subseq_start,
                         target_subseq_end,
                         output_fasta_path,
                         main_data_dir,
                         additional_flanking_basepairs=0
                         ):
    """Take a sequence ID, FASTA filepath, start position, end position,
    and output filepath, and write a FASTA subsequence to the filepath.
    """
    # Retrieve object for full hit sequence.
    full_seq_obj = get_seqs_from_fasta_db(db_filepath,
                                          [subj_seq_id],
                                          main_data_dir,
                                          False)[0]
    #full_seq_obj = get_seq_obj_from_db_fasta([subj_seq_id], db_filepath)[0]
    #full_seq_obj = get_seq_obj_from_db_fasta([subj_seq_id], db_filepath)[0]

    # Redefine more inclusive subsequence of interest, including flanking
    # regions.
    target_subseq_start = max([target_subseq_start - additional_flanking_basepairs, 1])
    target_subseq_end = min([target_subseq_end + additional_flanking_basepairs, len(full_seq_obj)])

    # Get appropriate subsequence.
    #subseq_seq = str(full_seq_obj.seq)[int(target_subseq_start) -1: int(target_subseq_end)]
    subseq_seq = full_seq_obj.seq[int(target_subseq_start) -1: int(target_subseq_end)]

    # Check that sequence is the right length.
    assert len(subseq_seq) == abs(int(target_subseq_start) - int(target_subseq_end)) + 1, "Error 1"

    # Update sequence attribute.
    full_seq_obj.seq = subseq_seq

    ## Append position info to the header.
    #full_seq_obj.description = full_seq_obj.description + ' [' +\
    #target_subseq_start + '..' + target_subseq_end + ']'

    # Write the subsequence to a new FASTA file.
    #assert not os.path.isfile(output_fasta_path), """Specified path for output
    #FASTA file already exists: %s""" % output_fasta_path
    with open(output_fasta_path, 'w') as o:
        SeqIO.write(full_seq_obj, o, 'fasta')

    # Return the start position of the identified subsequence as an integer, to
    # keep track of where the exons start and end in the original sequence.
    return target_subseq_start


# Run exonerate as a subprocess.
def run_exonerate_as_subprocess(query_prot_faa,
                                subj_subseq_fna,
                                exonerate_output_filepath,
                                exonerate_score_threshold,
                                genetic_code=1
                                ):
    """Take a query peptide FASTA sequence file (with a single sequence
    record), a subject sequence subsequence nucleotide FASTA file to search in,
    and an output file path, then run exonerate and write output to the
    specified path.
    """
    # Define exonerate command list.
    exonerate_command_list = ['exonerate',
                              query_prot_faa,
                              subj_subseq_fna,
                              #'--gappedextension', # Do gapped HSP extension (default).
                              '--extensionthreshold',
                              '50', # Set gapped extension threshold (default=50).
                              '--gapextend',
                              '-4', # For Affine models (default=-4).
                              '-i',
                              '-15', # Specify penalty for adding introns default -30.
                              '-m',
                              'protein2genome',
                              '--ryo',
                              '\"\n------------\n>%ti\n%tcs\n------------\"',
                              #'--showalignment',
                              #'False', # Could save space.
                              '--geneticcode',
                              str(genetic_code),
                              '--score',
                              exonerate_score_threshold,
                              '--proteinhspthreshold', 
                              '25', # Default setting is 30.
                              '--proteinwordlimit', 
                              '4', # Default setting is 4.
                              '--frameshift',
                              '-28', # Default penalty for adding a frameshift is -28.
                              '--splice3',
                              'primate', # Default setting.
                              '--splice5',
                              'primate' # Default setting.
                             ]

    # Call exonerate in a subprocess.
    with open(exonerate_output_filepath, 'w') as o:
        subprocess.call(exonerate_command_list, stdout=o, stderr=subprocess.STDOUT)



if __name__ == '__main__':
    # Parse input.
    command_line_list = sys.argv
    query_faa = str(command_line_list[1])
    target_fna_name = str(command_line_list[2])
    target_seq_id = str(command_line_list[3])
    target_subseq_start = str(command_line_list[4])
    target_subseq_end = str(command_line_list[5])
    genetic_code = str(command_line_list[6])

    # Get filepath for specified query FASTA filename.
    query_dir = DataPaths(main_data_dir).querydirpath
    query_faa_path = os.path.join(query_dir, query_faa)
    assert os.path.isfile(query_faa_path), """Specified query file path is
    not a file: %s""" % query_faa_path

    # Get filepath for specified subject FASTA filename.  
    db_dir = DataPaths(main_data_dir).dbdirpath
    target_fna_path = os.path.join(db_dir, target_fna_name)
    assert os.path.isfile(target_fna_path), """Specified database file path is
    not a file: %s""" % target_fna_path

    # Define path to FASTA file with subsequence of interest from target
    # nucleotide sequence.
    subseq_fasta_path = query_faa.rsplit('.', 1)[0] + '_subject_subseq.fna'

    # Extract relevant subsequence from input target sequence (region identified in
    # a previous step using TBLASTN).
    get_subseq_from_nucl(target_seq_id, 
                         target_fna_path,
                         int(target_subseq_start),
                         int(target_subseq_end),
                         subseq_fasta_path
                         )

    # Define path to exonerate output file.
    exonerate_output_filepath = subseq_fasta_path.rsplit('.', 1)[0] + '_exonerate_out.txt'
    
    # Run exonerate as a subprocess.
    run_exonerate_as_subprocess(query_faa_path,
                                subseq_fasta_path,
                                exonerate_output_filepath,
                                '20',
                                genetic_code
                                )
    
    # Parse output of exonerate.
    #parse_exonerate_output(exonerate_output_filepath)
    genetic_code_number = '1'
    position_of_subject_seq_start_in_original = int(target_subseq_start)
    exonerate_locus_result_obj = ExonerateLocusResult(exonerate_output_filepath,
                                                      subseq_fasta_path,
                                                      position_of_subject_seq_start_in_original,
                                                      genetic_code_number
                                                      )
    
    # Write output to file.
    exonerate_output_seq_filepath = subseq_fasta_path.rsplit('.', 1)[0] +\
    '_exonerate_out_seq.faa'
    with open(exonerate_output_seq_filepath, 'w') as o:
        SeqIO.write(exonerate_locus_result_obj.seq_record, o, 'fasta')

    # Delete temporary intermediate files.
    os.remove(subseq_fasta_path)



