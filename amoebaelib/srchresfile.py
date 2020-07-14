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
"""Module for defining a module for parsing similarity search output files with
various formats.
"""
# Import built-in modules.
#import argparse
#import sys
import os
#import subprocess
import re
from datapaths import DataPaths
#import shutil
#import glob
#import time
#import pandas as pd

# Import modules from installed libraries/packages.
#from Bio import SeqIO
#from Bio import SearchIO
# Import SearchIO and suppress experimental warning
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

#from Bio import AlignIO
#from Bio.Alphabet import IUPAC, Gapped
from Bio.Blast import NCBIXML

from amoebae_m import get_seqs_from_fasta_db, get_subseq_from_fasta_db
from search_scaffolds import get_blastp_hit_seq_obj_and_coord, get_tblastn_hit_seq_obj_and_coord

# How to iterate over all hits in a result file?


class SrchResFile:
    """Parses relevant information as needed from similarity search output
    file, which may be one of several types.
    """
    def __init__(self, filepath, main_data_dir):
        # Check that the input file path exists.
        assert os.path.isfile(filepath), """Input filepath does not exist:
        %s""" % filepath

        self.filepath = filepath

        self.main_data_dir = main_data_dir

        # Get basic info from file.
        info = get_srch_file_info(filepath)

        # Determine program used to generate input file.
        self.program = None
        self.program = info[0]
        assert self.program is not None, """Could not determine the name of
        the program that produced the similarity search result file: %s"""\
        % filepath

        # Determine version of program used to generate input file.
        self.version = None
        self.version = info[1]
        assert self.version is not None, """Could not determine the version of
        the program that produced the similarity search result file: %s"""\
        % filepath

        # Determine format type of input file.
        self.format = None
        self.format = info[2]
        assert self.format is not None, """Could not determine the name of
        the format type of the similarity search result file: %s"""\
        % filepath

        # Check that the file contains only results for a search with a single
        # query.
        if not self.format == 'hhsearch':
            assert len(list(SearchIO.parse(filepath, self.format))) == 1, """More than
            one search result contained in input file: %s""" % filepath
        else:
            pass # ...

        # Determine number of hits in input file.
        self.num_hits = None
        if not self.format == 'hhsearch':
            self.num_hits = len(SearchIO.read(filepath, self.format))
            assert self.num_hits is not None, """Could not determine the number of
            hits listed in the similarity search result file: %s"""\
            % filepath
        else:
            pass # ...

        #assert self.num_hits >= 0, """Could not determine the number of
        #hits listed in the similarity search result file: %s"""\
        #% filepath

        # Get the query and database file paths.
        #self.query_file = None
        self.db_file = None

        assert self.format != 'hmmer3-tab', """Does not work with tabular
        format."""
        if not self.format == 'hhsearch':
            p = SearchIO.read(filepath, self.format)
            self.db_file = os.path.basename(p.target)

        else:
            pass # ...

        #assert self.query_file is not None, """Could not determine query file
        #name listed in the similarity search result file: %s"""\
        #% filepath
        assert self.db_file is not None, """Could not determine database file
        name listed in the similarity search result file: %s"""\
        % filepath

        # Define full path to database file.
        self.db_file_path = None
        self.db_file_path = os.path.join(DataPaths(self.main_data_dir).dbdirpath, self.db_file)
        # Check that it is a real file.
        assert os.path.isfile(self.db_file_path), """Path to database is not a file:
        %s""" % self.db_file_path


    def hit_id(self, hit_rank):
        """Return id (accession) (if present) for sequence/profile with a given
        rank in the results.
        """
        hit_id = None
        if not self.format == 'hhsearch':
            hit_id = SearchIO.read(self.filepath, self.format)[hit_rank].id

        else:
            pass # ...

        # Check that it worked.
        assert hit_id is not None, """Could not retrieve id for hit with rank
        %s in search result file: %s""" % (str(hit_rank), filepath)

        # Return the hit id.
        return hit_id


    def all_hit_ids(self):
        """Return a list of ids (accessions) (if any are present) for all the
        hits.
        """
        hit_ids = []

        if not self.format == 'hhsearch':
            hit_ids = [x.id for x in SearchIO.read(self.filepath, self.format)]

        else:
            pass # ...

        return hit_ids


    def hit_description(self, hit_rank):
        """Return description (if present) for sequence/profile with a given
        rank in the results.
        """
        hit_descr = None
        if not self.format == 'hhsearch':
            hit_descr = SearchIO.read(self.filepath,\
                    self.format)[hit_rank].description
            if hit_descr == '':
                hit_descr = SearchIO.read(self.filepath,\
                        self.format)[hit_rank].id

        else:
            pass # ...

        # Check that it worked.
        assert hit_descr is not None, """Could not retrieve description for hit
        with rank %s in search result file: %s""" % (str(hit_rank), filepath)

        # Return the hit description.
        return hit_descr


    def hit_evalue(self, hit_rank):
        """Return E-value for hit sequence/profile with a given rank in the
        results.
        """
        hit_evalue = None
        if self.format == 'blast-xml':
            # Get the E-value for the top hit.
            with open(self.filepath) as infh:
                for blast_record in NCBIXML.parse(infh):
                    hit_num = 0
                    for hit in blast_record.descriptions:
                        if hit_num == hit_rank: 
                            hit_evalue = hit.e
                            break
                        hit_num += 1
                    break

        elif self.format == 'hmmer3-text':
            hit_evalue = SearchIO.read(self.filepath,\
                    self.format)[hit_rank].evalue

        elif self.format == 'hhsearch':
            pass # ...

        else:
            pass # ...?

        # Check that a value for the E-value was found.
        assert hit_evalue is not None, """Could not retrieve E-value for hit
        with rank %s in search result file: %s""" % (str(hit_rank), self.filepath)

        # Check that the E-value looks like an E-value
        assert hit_evalue is 0 or isinstance(hit_evalue, float), """The E-value
        %s identified in the file %s does not look like an E-value.""" %\
        (hit_evalue, self.filepath)

        # Return the hit E-value.
        return hit_evalue


    def hit_score(self, hit_rank):
        """Return bitscore (from BLAST programs) or sequence score (from HMMer
        output) (or something else for hhsearch?) for hit sequence/profile with
        a given rank in the results.
        """
        hit_score = None
        if self.format == 'blast-xml':
            # Get the score for the top hit.
            with open(self.filepath) as infh:
                for blast_record in NCBIXML.parse(infh):
                    for hit in blast_record.descriptions:
                        hit_score = hit.score
                        break
                    break

        elif self.format == 'hmmer3-text':
            hit_score = SearchIO.read(self.filepath,\
                    self.format)[hit_rank].bitscore

        elif self.format == 'hhsearch':
            pass # ...

        else:
            pass # ...?

        # Check that it worked.
        assert hit_score is not None, """Could not retrieve score for hit
        with rank %s in search result file: %s""" % (str(hit_rank), self.filepath)

        # Return the hit score.
        return hit_score


    def hit_sequence(self, hit_rank):
        """Return a Seq object for full sequence of subject sequence.
        """
        # Get path for databases directory.
        dbdir_path = DataPaths(self.main_data_dir).dbdirpath

        # Get sequence object.
        seq_obj = None
        seq_id = self.hit_id(hit_rank) 
        db_path = os.path.join(dbdir_path, self.db_file)
        assert os.path.isfile(db_path), """Given path is not a file."""
        #seq_obj = get_seqs_from_fasta_db(db_path, [seq_id])[0]
        seq_obj = get_seqs_from_fasta_db(self.db_file, [seq_id], self.main_data_dir)[0]

        # Check that it worked.
        assert seq_obj is not None, """Could not retrieve sequence for hit."""

        # Return sequence object.
        return seq_obj


    def hit_subsequence_and_coord(self, hit_rank, max_gap=10000):
        """Return a Seq object for subsequence of subject sequence that
        actually aligns to the query sequence/profile.
        """
        subseq_obj = None
        subseq_coord = None

        # Get hit object with SearchIO parser.
        searchio_hit_obj = None
        with open(self.filepath) as infh:
            hit_num = -1 
            for hit in SearchIO.read(infh, self.format):
                hit_num += 1
                if hit_num == hit_rank:
                    searchio_hit_obj = hit
                    break

        # Check that hit object was retreived.
        assert searchio_hit_obj != None, """Could not retrieve Bio.SearchIO
        hit object from file."""

        # Process SearchIO hit object to get subsequence object and
        # coordinates differently depending on format.
        if self.format == 'blast-xml':
            # Need to concatenate HSPs in a logical manner, and differently for
            # blastp vs. tblastn.
            if self.program == 'blastp':
                # Use the search_scaffolds module.
                xlist =\
                get_blastp_hit_seq_obj_and_coord(searchio_hit_obj, max_gap)
                subseq_obj = xlist[0]
                subseq_coord = xlist[1]

            elif self.program == 'tblastn':
                # Use the search_scaffolds module.
                xlist =\
                get_tblastn_hit_seq_obj_and_coord(searchio_hit_obj, max_gap)
                subseq_obj = xlist[0]
                subseq_coord = xlist[1]

        elif self.format == 'hmmer3-text':
            # More straight-forward, because the sequences can be taken
            # directly from the database file.

            # Get path for databases directory.
            dbdir_path = DataPaths(self.main_data_dir).dbdirpath

            # Get coordinates.
            subseq_coord = get_hmmer_hit_seq_coord(searchio_hit_obj,
                                                   self.db_file,
                                                   self.main_data_dir)

            # Get sequence object.
            seq_id = self.hit_id(hit_rank) 
            db_path = os.path.join(dbdir_path, self.db_file)
            subseq_obj = get_subseq_from_fasta_db(db_path,
                                                  seq_id,
                                                  subseq_coord,
                                                  self.main_data_dir)

        else:
            pass # ...?

        # Check that it worked.
        assert subseq_obj is not None, """Could not retrieve sequence for hit."""
        assert subseq_coord is not None, """Could not retrieve sequence
        coordinates for hit."""

        # Return sequence object and coordinates.
        return [subseq_obj, subseq_coord]


    def hit_query_subsequence(self, hit_rank):
        """Return Seq object for subsequence of query sequence (or consensus of
        profile?) that aligns to the subject sequence.

        Is this even a necessary method?
        """
        pass


    def rank_of_first_nonredun_hit(self, redun_hit_id_list):
        """Efficiently determine the rank of the first hit that does not have
        an ID that is in the given list, or return None if there aren't any.
        """
        first_nonredun_hit_rank = None

        with open(self.filepath) as infh:
            hit_num = -1 
            for hit in SearchIO.read(infh, self.format):
                hit_num += 1
                if hit.id not in redun_hit_id_list:
                    first_nonredun_hit_rank = hit_num
                    break

        # Return the rank.
        return first_nonredun_hit_rank






################################
# Define functions used by the SrchResFile class that are best separated from
# methods of the class.


def get_srch_file_info(search_result_path):
    """Returns the format string necessary for the SearchIO module to parse a
    given file.
    """
    # Define regular expressions that will only match any line in a file if it
    # is of a certain format. See the SearchIO documentation here:
    # http://biopython.org/DIST/docs/api/Bio.SearchIO-module.html
    fmt_expr_dict = {'blast-xml': [re.compile(r'^<!DOCTYPE BlastOutput')],
                     'hmmer3-tab': [re.compile(r'^# Option settings: hmmsearch --tblout')],
                     'hmmer3-text': [re.compile(r'^# HMMER ')]
                     }
    
    # See if any of the regular expressions match the file.
    ident_fmt = None
    with open(search_result_path) as srh:
        for line in srh:
            for fmt in fmt_expr_dict.keys():
                for r in fmt_expr_dict[fmt]:
                    if r.search(line) is not None:
                        ident_fmt = fmt
                        break
                if ident_fmt is not None:
                    break

    # Check that a format was found.
    assert ident_fmt is not None, """Error: Could not identify format for
    search result file: %s""" % search_result_path

    ident_prog = None
    with open(search_result_path) as srh:
        ident_prog = SearchIO.read(srh, ident_fmt).program
    assert ident_prog is not None, """Could not identify program name from file
    format type."""

    # Decide which dict key to use based on identified program version, and
    # then identify version in file.
    ident_vers = None
    #vers_key = None
    with open(search_result_path) as srh:
        ident_vers = SearchIO.read(srh, ident_fmt).version

    # Check that the program version was found.
    assert ident_vers is not None, """Error: Could not identify name of program
    used to generate search result file: %s""" % search_result_path

    # Return a list defining the input file format.
    return [ident_prog, ident_vers, ident_fmt]


def get_hmmer_hit_seq_coord(searchio_hit_obj,
                            db_file,
                            main_data_dir,
                            extra=0):
    """Takes a Bio.SearchIO hmmer3-text hit object and returns specially
    formatted coordinates for the subsequence (top HSP) that match the query.

    Note: The coordinates are the literal start and end residue numbers for the
    relevant portion of the subject sequence.
    """

    # Check that the hit object is usable.
    #print('\n')
    #print(searchio_hit_obj.__dict__)
    #assert searchio_hit_obj.domain_obs_num >= 1, """No domains identified, so
    #no hit sequence coordinates exist."""

    if searchio_hit_obj.domain_obs_num < 1:
        # Just return coordinates that will give you an empty string when
        # slicing the sequence.
        return [[0, 0]]


    else:
        #domainstart = int(linelist[10]) - 1
        domainstart = searchio_hit_obj[0].hit_start - 1

        #domainend = int(linelist[11])
        domainend = searchio_hit_obj[0].hit_end

        #seqlen = len(str(seq_record.seq))
        #seqlen = get_seq_len_from_db(searchio_hit_obj.id,\
        #        db_file)
        # Get sequence length from sequence in database.
        seqlen = len(get_seqs_from_fasta_db(db_file,
                                            [searchio_hit_obj.id],
                                            main_data_dir)[0])
        nterm_margin = domainstart
        cterm_margin = seqlen - domainend
        nterm_extra = int(extra)
        cterm_extra = int(extra)
        if nterm_margin < nterm_extra:
            nterm_extra = nterm_margin
        if cterm_margin < cterm_extra:
            cterm_extra = cterm_margin
        cut_start = domainstart - nterm_extra
        cut_end = domainend + cterm_extra #+ 1
        #domain = seq_record.seq[cut_start:cut_end]
        #domain = seq_record.seq[cut_start:cut_end]
        #seq_record.seq = domain
        #o.write('>' + seq_record.id + '\n' 
        #        + str(seq_record.seq) + '\n')

        # Return coordinates.
        return [[cut_start, cut_end]]
        



