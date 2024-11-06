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
"""This module is for using nhmmer to find the location of nucleotide sequences
matching an HMM, and outputs to a fasta file.

"""


import sys
import os
import re
import subprocess
import itertools
from get_fas_from_db_dir import get_seq_obj_from_db_fasta
from datapaths import DataPaths
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC

# Maybe it would be better to get the whole gene region and align that to the
# nucleotide HMM, and then make a translation based on that alignment...


class NhmmerHitAlig:
    def __init__(self):
        self.query_name = None
        self.db_filepath = None
        self.subj_seq_id = None
        self.strand = None
        self.score = None
        self.evalue = None
        self.hmmfrom = None
        self.hmmto = None
        self.alifrom = None
        self.alito = None
        self.envfrom = None
        self.envto = None

    def sequ(self):
        # Retrieve object for full hit sequence.
        # Maybe there is a more efficient way to do this???
        full_seq_obj = get_seq_obj_from_db_fasta([self.subj_seq_id],
                self.db_filepath)[0]

        # Get appropriate subsequence that aligned to HMM.
        alig_seq = ''
        if self.strand == '+':
            alig_seq = str(full_seq_obj.seq)[int(self.alifrom) -1: int(self.alito)]

        elif self.strand == '-':
            #***********alifrom comes after alito if on neg strand!!!!!!!!!!!!
            #***********Have to reverse them to get the right subsequence.
            subseqstring = str(full_seq_obj.seq)[int(self.alito) -1: int(self.alifrom)]
            #x = Seq(subseqstring, IUPAC.ambiguous_dna)
            x = Seq(subseqstring)
            y = x.reverse_complement()
            alig_seq = str(y)

        # Check that sequence is the right length.
        assert len(alig_seq) == abs(int(self.alito) - int(self.alifrom)) + 1, "Error 1"

        # Return subsequence.
        return alig_seq

    def sequprot(self):
        """Translate the sequence based on how it aligns to the HMM.

        ***Assumes that the HMM starts with the start codon of the aligned
        homologues, and that the rest of the alignment is in frame.

        ***Assumes indels keep the whole hit in frame.

        ***Assumes standard genetic code.

        ***Does not always to a good job of translating in the right frame...

        May not be sufficient. May need to align to whole query sequence to
        determine frame along length.
        """
        # Determine frame by how it aligns to the HMM, and remove bases
        # accordingly before translating.
        seq_in_frame = frameify_seq(self.sequ(), self.hmmfrom)

        # Make the in frame sequence a Seq object, and translate.
        #seq_in_frame_obj = Seq(seq_in_frame, IUPAC.ambiguous_dna)
        seq_in_frame_obj = Seq(seq_in_frame)
        translation = str(seq_in_frame_obj.translate())

        # Return translation of sequence.
        #print('found seq for ' + self.subj_seq_id)
        return translation

    def db_name(self):
        # Return the name of the db of origin.
        return os.path.basename(self.db_filepath).rsplit('_', 1)[0]

    def hit_range(self):
        x = [self.alifrom, self.alito]
        return (min(x), max(x))


#class HmmerHitAlig:
#    def __init__(self):
#        self.query_name = None
#        self.db_filepath = None
#        self.subj_seq_id = None
#        self.original_strand = None
#        self.score = None
#        self.evalue = None
#        self.hmmfrom = None
#        self.hmmto = None
#        self.alifrom = None
#        self.alito = None
#        self.envfrom = None
#        self.envto = None
#        self.subj_seq_position = None
#
#    def position(self):
#        """Returns a string with the position in the subject sequence.
#        """
#        return '[' + self.alifrom + '..' + self.alito + ']'
#
#    def actual_position(self):
#        """Returns a string with the position in the original db_filepath.
#        """
#        subj_seq_start = self.subj_seq_position.strip('[]').split('..')[0]
#        subj_seq_end = self.subj_seq_position.strip('[]').split('..')[1]
#        # This gets really complex when you have to consider the
#        # strand!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#        if self.original_strand == '+':
#            return '[' + self.alifrom + 
#
#        elif self.original_strand == '-':
#
#
#
#    #def sequ(self):
#    #    # Retrieve object for full hit sequence.
#    #    # Maybe there is a more efficient way to do this???
#    #    full_seq_obj = get_seq_obj_from_db_fasta([self.subj_seq_id],
#    #            self.db_filepath)[0]
#
#    #    # Get appropriate subsequence that aligned to HMM.
#    #    alig_seq = ''
#    #    if self.strand == '+':
#    #        alig_seq = str(full_seq_obj.seq)[int(self.alifrom) -1: int(self.alito)]
#
#    #    elif self.strand == '-':
#    #        #***********alifrom comes after alito if on neg strand!!!!!!!!!!!!
#    #        #***********Have to reverse them to get the right subsequence.
#    #        subseqstring = str(full_seq_obj.seq)[int(self.alito) -1: int(self.alifrom)]
#    #        x = Seq(subseqstring, IUPAC.ambiguous_dna)
#    #        y = x.reverse_complement()
#    #        alig_seq = str(y)
#
#    #    # Check that sequence is the right length.
#    #    assert len(alig_seq) == abs(int(self.alito) - int(self.alifrom)) + 1, "Error 1"
#
#    #    # Return subsequence.
#    #    return alig_seq
#
#    #def sequprot(self):
#    #    """Translate the sequence based on how it aligns to the HMM.
#
#    #    ***Assumes that the HMM starts with the start codon of the aligned
#    #    homologues, and that the rest of the alignment is in frame.
#
#    #    ***Assumes indels keep the whole hit in frame.
#
#    #    ***Assumes standard genetic code.
#
#    #    ***Does not always to a good job of translating in the right frame...
#
#    #    May not be sufficient. May need to align to whole query sequence to
#    #    determine frame along length.
#    #    """
#    #    # Determine frame by how it aligns to the HMM, and remove bases
#    #    # accordingly before translating.
#    #    seq_in_frame = frameify_seq(self.sequ(), self.hmmfrom)
#
#    #    # Make the in frame sequence a Seq object, and translate.
#    #    seq_in_frame_obj = Seq(seq_in_frame, IUPAC.ambiguous_dna)
#    #    translation = str(seq_in_frame_obj.translate())
#
#    #    # Return translation of sequence.
#    #    #print('found seq for ' + self.subj_seq_id)
#    #    return translation
#
#    #def db_name(self):
#    #    # Return the name of the db of origin.
#    #    return os.path.basename(self.db_filepath).rsplit('_', 1)[0]
#
#    #def hit_range(self):
#    #    x = [self.alifrom, self.alito]
#    #    return (min(x), max(x))


def frameify_seq(inseq, hmmfrom):
    """Take a string, and return a modified string that is in fram with the
    HMM.
    """
    seq_in_frame = None
    hmmfromposition = int(hmmfrom)
    remainder = hmmfromposition % 3
    if remainder == 1:
        #frame = 1
        seq_in_frame = inseq 
    elif remainder == 2:
        #frame = 2
        seq_in_frame = inseq[2:]
    elif remainder == 0:
        #frame = 3
        seq_in_frame = inseq[1:]

    return seq_in_frame


def check_if_hits_overlapping(hit_objs):
    """Takes a list of nhmmer hit objects, and returns true if they overlap.
    """
    # Compile a list of the ranges.
    ranges = []
    for hit in hit_objs:
        r = hit.hit_range() 
        ranges.append(r)
    
    # Generate list of all possible combinations of two ranges.
    combos = itertools.combinations(ranges, 2)

    # Loop over combinations and check for overlap.
    overlap = False
    for combo in combos:
        x = set(combo[0])
        y = set(combo[1])
        intersect = x.intersection(y)
        if intersect != set():
            overlap = True

    return overlap


def get_db_fasta_filepath(dir_string, dbname):
    """Return the name of the appropriate fasta file."""
    subdir_list = [subdir for subdir in os.listdir(dir_string)]
    fasta_filepath = None
    for subdir in subdir_list:
        if subdir == dbname.strip():
            deep_subdir = os.path.join(dir_string, subdir)
            deep_subdir_list = [deepsubdir for deepsubdir in os.listdir(deep_subdir) if\
                    os.path.isdir(os.path.join(deep_subdir, deepsubdir))]
            for deeper_subdir in deep_subdir_list:
                if re.search(r'_Nucl', deeper_subdir):
                    deeper_subdir = os.path.join(deep_subdir, deeper_subdir)
                    onlyfiles = [f for f in os.listdir(deeper_subdir) if\
                    os.path.isfile(os.path.join(deeper_subdir, f))]
                    for f in onlyfiles:
                        if f.endswith('.fa'):
                            fasta_filepath = os.path.join(deeper_subdir, f)
                            assert not fasta_filepath == None, "Error: Fasta file not identified."
    return fasta_filepath


def run_nhmmer_search(query_path, db_path, nhmmer_outfilepath):
    """Runs nhmmer."""
    subprocess.call(["nhmmer", "-E", "0.05", "--cpu", "2", "--tblout",\
            nhmmer_outfilepath, query_path, db_path])


def run_hmmsearch(query_path, db_path, hmmsearch_outfilepath):
    """Runs hmmsearch (protein queries)."""
    # Could use a bit score cutoff instead using the "-T" option.
    subprocess.call(["hmmsearch", "-E", "0.05", "--cpu", "2", "--domtblout",\
            hmmsearch_outfilepath, query_path, db_path])
    

def parse_nhmmer_out(nhmmer_outfilepath):
    """Opens an nhmmer output file (tabular output), extracts the information for location of
    domains in the sequences, and returns a list of hit alignment objects.
    """
    # Loop through the lines and identify the relevant info applying to all hits.
    target_filepath = None
    hmm_name = None
    with open(nhmmer_outfilepath) as i:
        # Get info from end of file.
        for line in i:
            split_mod_line = re.split(r' +', line)

            if line.startswith('# Target file:'):
                target_filepath = line.replace(' ', '').split(':')[1].strip()

            if line.startswith('# Query file:'):
                hmm_name = line.replace(' ', '').split(':')[1].strip()

    # Loop through the lines again and make objects for each hit alignement.
    hit_alignments = []
    with open(nhmmer_outfilepath) as i:
        for line in i:
            write = None
            
            # Identify lines with relevant info.
            if not line.startswith('#'):
                split_mod_line = re.split(r' +', line)

                # Construct object.
                a = NhmmerHitAlig()
                a.query_name = hmm_name
                a.db_filepath = target_filepath
                a.subj_seq_id = split_mod_line[0]
                a.strand = split_mod_line[11]
                a.score = split_mod_line[13]
                a.evalue = split_mod_line[12]
                a.hmmfrom = split_mod_line[4]
                a.hmmto = split_mod_line[5]
                a.alifrom = split_mod_line[6]
                a.alito = split_mod_line[7]
                a.envfrom = split_mod_line[8]
                a.envto = split_mod_line[9]

                # Add object to list.
                hit_alignments.append(a)

    # Check that at least one hit was found.
    assert len(hit_alignments) >= 1, "Error: No hits in file " +\
 nhmmer_outfilepath

    # Return list of objects.
    return hit_alignments


#def parse_hmmer_out(hmmer_outfilepath):
#    """Opens a hmmer output file (tabular output), extracts the information for location of
#    domains in the sequences, and returns a list of hit alignment objects.
#    """
#    # Loop through the lines and identify the relevant info applying to all hits.
#    target_filepath = None
#    hmm_name = None
#    with open(hmmer_outfilepath) as i:
#        # Get info from end of file.
#        for line in i:
#            split_mod_line = re.split(r' +', line)
#
#            if line.startswith('# Target file:'):
#                target_filepath = line.replace(' ', '').split(':')[1].strip()
#
#            if line.startswith('# Query file:'):
#                hmm_name = line.replace(' ', '').split(':')[1].strip()
#
#    # Loop through the lines again and make objects for each hit alignement.
#    hit_alignments = []
#    with open(hmmer_outfilepath) as i:
#        for line in i:
#            write = None
#            
#            # Identify lines with relevant info.
#            if not line.startswith('#'):
#                split_mod_line = re.split(r' +', line)
#
#                # Construct object.
#                a = HmmerHitAlig()
#                a.query_name = hmm_name a.db_filepath = target_filepath
#                a.subj_seq_id = split_mod_line[0]
#                a.original_strand = split_mod_line[24].strip('()')
#                a.score = split_mod_line[13]
#                a.evalue = split_mod_line[12]
#                a.hmmfrom = split_mod_line[15]
#                a.hmmto = split_mod_line[16]
#                a.alifrom = split_mod_line[17]
#                a.alito = split_mod_line[18]
#                a.envfrom = split_mod_line[19]
#                a.envto = split_mod_line[20]
#                a.subj_seq_position = split_mod_line[23]
#
#                #print(a.db_filepath)
#                #print(a.subj_seq_id)
#                #print(a.strand)
#                #print(a.score)
#                #print(a.evalue)
#                #print(a.hmmfrom)
#                #print(a.hmmto)
#                #print(a.alifrom)
#                #print(a.alito)
#                #print(a.envfrom)
#                #print(a.envto)
#                #print(a.pos_in_orig_db_seq)
#
#                # Add object to list.
#                hit_alignments.append(a)
#
#    # Check that at least one hit was found.
#    #assert len(hit_alignments) >= 1, "Error: No hits in file " +\
#    #hmmer_outfilepath
#
#    # Return list of objects.
#    if len(hit_alignments) < 1:
#        # Return an empty list if no alignments were in input file.
#        return []
#    else:
#        return hit_alignments


# Obsolete?
def write_nhmmer_hit_fasta(nhmmer_hit_objs, fasta_filepath):
    """Takes a list of nhmmer hit alignment objects, and writes corresponding
    fasta records to a fasta file.
    """
    with open(fasta_filepath, 'w') as o:
        for i in nhmmer_hit_objs:
            fasta_header = '>' + i.subj_seq_id + ' [' + i.alifrom + '..' +\
            i.alito + '] ' + '(' + i.strand + ') ' + i.db_name() + '\n'

            o.write(fasta_header + i.sequ() + '\n')


def cluster_ali_recur(sublists, origlist):
    """Recursively cluster each object in the list with those that are within
    10,000 bp. Return a dict.
    """
    # Problem: cannot assume that the whole gene will be no more than 10kb
    # long. The furthest exon may have a start position further than 10kb from
    # the start position of some other exons. This function is a little sloppy,
    # and may yeild inaccurate results in some cases.
    #***The assumption should be that individual introns should be no more than
    #~10 kb.
    
    # Assume that gene coding sequence does not exceed a certain value.
    max_coding_sequence_len = 1000000

    if len(origlist) == 0:
        return sublists

    elif len(origlist) == 1:
        return sublists + [origlist]

    else:
        sublist = []
        new_orig_list = []
        for i in origlist:
            # If start positions are within 10000 bp, add to sublist.
            # (The first one will be included as well).
            if abs(int(origlist[0].alifrom) - int(i.alifrom)) -1 <= max_coding_sequence_len:
                sublist.append(i)
            else:
                new_orig_list.append(i)
        return [sublist] + cluster_ali_recur(sublists, new_orig_list)


def get_concat_seq(nhmmer_hit_objs, dtype='nucl'):
    """Takes a list of nhmmer hit alignment objects that should correspond to
    the same gene, and returns a concatenated sequence.

    May not be sufficient if the hit alignments do overlap.
    """
    # Make sure the hits aren't overlapping.
    assert not check_if_hits_overlapping(nhmmer_hit_objs), "Error: Overlapping\
 nhmmer hits should not be concatenated."

    # Determine strand.
    strand = nhmmer_hit_objs[0].strand
    db = nhmmer_hit_objs[0].db_name()

    if strand == '+':
        # Sort by start position.
        nhmmer_hit_objs.sort(key=lambda x: int(x.alifrom))
    elif strand == '-':
        # Sort by start position in reverse.
        nhmmer_hit_objs.sort(key=lambda x: int(x.alifrom), reverse=True)

    # Assemble concatenated file.
    header = '>' + nhmmer_hit_objs[0].subj_seq_id + ' ' + db + ' '
    concat_seq = ''
    for i in nhmmer_hit_objs:
        fromtolist = [i.alifrom, i.alito]
        header = header + '[' + min(fromtolist) + '..' + max(fromtolist) + ']' + ',' 
        if dtype == 'nucl':
            concat_seq = concat_seq + i.sequ()
        elif dtype == 'prot':
            concat_seq = concat_seq + i.sequprot()
    header = header.rstrip(',') + ' (' + strand + ')' + '\n' 

    return header + concat_seq + '\n' 


def get_gene_region(nhmmer_hit_objs):
    """Takes a list of nhmmer hit alignment objects that should correspond to
    the same gene, and returns a concatenated sequence (of appropriate strand).

    So, if the gene region is on the negative strand, then it returns the
    sequence of the negative strand (5' to 3') by getting the reverse
    complement of the forward strand.

    ***Does not account for paralogous genes possibly being very close to each
    other on the same strand of the same chromosome...
    """
    # Determine strand.
    strand = nhmmer_hit_objs[0].strand
    db = nhmmer_hit_objs[0].db_name()

    # Determine the inclusive range of this hits.
    incl_range = None
    if strand == '+':
        # Left to right.
        incl_range = (int(min(nhmmer_hit_objs, key=lambda x: int(x.alifrom)).alifrom),
                int(max(nhmmer_hit_objs, key=lambda x: int(x.alito)).alito))
    elif strand == '-':
        # Right to left.
        incl_range = (int(min(nhmmer_hit_objs, key=lambda x: int(x.alito)).alito),
                int(max(nhmmer_hit_objs, key=lambda x: int(x.alifrom)).alifrom))

    # Get full sequence object.
    full_seq_obj = get_seq_obj_from_db_fasta([nhmmer_hit_objs[0].subj_seq_id],
            nhmmer_hit_objs[0].db_filepath)[0]

    # Define how much extra to include on either side of the inclusive range.
    extra = 100

    # Define extended range (with extra on either side).
    left_margin = incl_range[0] -1
    left_extra = min([extra, left_margin])
    right_margin = len(full_seq_obj) - incl_range[1]
    right_extra = min([extra, right_margin])
    exten_range = (incl_range[0] - left_extra, incl_range[1] + right_extra)

    # Assemble header string. 
    header = '>' + nhmmer_hit_objs[0].subj_seq_id + ' ' + db + ' [' +\
 str(exten_range[0]) + '..' + str(exten_range[1]) + '] (' + strand + ')'

    # Get sequence.
    sequence = None
    if strand == '+':
        sequence = full_seq_obj.seq[exten_range[0] -1:exten_range[1]]
    elif strand == '-':
        sequence = full_seq_obj.seq[exten_range[0] -1:exten_range[1]].reverse_complement()

    # Return a printable fasta record string.
    #fasta_string = header + '\n' + str(sequence) + '\n' 
    #return fasta_string
    
    # Return a sequence object.
    return [header, sequence] 


def display_dict(nhmmer_hit_objs_dict2):
    """Prints relevant info."""
    print('\n')
    for key1 in nhmmer_hit_objs_dict2.keys():
        print(key1)
        for key2 in nhmmer_hit_objs_dict2[key1].keys():
            print('\t' + key2)
            obj_sublists = nhmmer_hit_objs_dict2[key1][key2]
            num = 0
            for i in obj_sublists:
                num += 1
                print('\t\tgene ' + str(num))
                for j in i:
                    print('\t\t\t' + j.alifrom + '..' + j.alito)


def get_three_frame_transl_string(inputs):
    """Takes a list comprising the header string and sequence string of an
    input sequence.
    """
    # Needs to transfer coordinates from original while taking into account the
    # strand, etc.
    header = inputs[0]
    strand = header.rsplit(' ', 1)[1].strip('()')
    position_in_original = re.compile(r'\[\d+\.\.\d+\]')
    start_end = position_in_original.findall(header)[0].strip('[]').split('..')
    start_end.sort(key=lambda x: int(x))
    #seq = Seq(str(inputs[1]), IUPAC.ambiguous_dna)
    seq = Seq(str(inputs[1]))

    frame1 = header + ' frame 1' + '\n' + str(seq.translate()) + '\n'

    x = 0
    y = 0
    if strand == '+':
        x = 1
    elif strand == '-':
        y = -1

    frame2_start_end = '[' + str(int(start_end[0]) + x) + '..' + str(int(start_end[1]) + y) + ']'
    frame2_header = position_in_original.sub(frame2_start_end, header)
    frame2 = frame2_header + ' frame 2' + '\n' + str(seq[1:].translate()) + '\n' 

    x = 0
    y = 0
    if strand == '+':
        x = 2
    elif strand == '-':
        y = -2

    frame3_start_end = '[' + str(int(start_end[0]) + x) + '..' + str(int(start_end[1]) + y) + ']'
    frame3_header = position_in_original.sub(frame3_start_end, header)
    frame3 = frame3_header + ' frame 3' + '\n' + str(seq[2:].translate()) + '\n' 

    outstring = frame1 + frame2 + frame3

    return outstring


# Obsolete?
def get_three_frame_transl(seq_filepath, three_frame_translation_dirpath):
    """Produce three-frame translations of forward strand DNA sequences in a
    given file, and write to files in a given directory path.
    """
    pass 



def get_nhmmer_hit_gene_region_seqs(nhmmer_hit_objs, seq_filepath):
    """Take a list of nhmmer hit alignment objects and return one or more
    nucleotide subsequences that may correspond to the coding sequence for a
    single gene.
    """
    # Cluster the objects into groups that may be from the same gene.
    # (on same strand, within a certain number of nucleotides of each other)

    # Split the list of objects into separate sub-lists for each fasta record
    # for which there are hits in the subject database.
    nhmmer_hit_objs_ids = []
    nhmmer_hit_objs_dict = {}
    for i in nhmmer_hit_objs:
        if i.subj_seq_id in nhmmer_hit_objs_ids:
            nhmmer_hit_objs_dict[i.subj_seq_id] = nhmmer_hit_objs_dict[i.subj_seq_id] + [i]
        else:
            nhmmer_hit_objs_dict[i.subj_seq_id] = [i]
            nhmmer_hit_objs_ids.append(i.subj_seq_id)

    # Split the list of object into those that are on the positive strand, and
    # those that are on the negative strand.    
    nhmmer_hit_objs_dict2 = {}
    for key in nhmmer_hit_objs_dict:
        # Identify objects with alignments on either positive or negative strand.
        obj_list = nhmmer_hit_objs_dict[key]
        pos = []
        neg = []
        for i in obj_list:
            assert i.strand == '+' or i.strand == '-', "Error: Strand not specified."
            if i.strand == '+':
                pos.append(i)
            elif i.strand == '-':
                neg.append(i)

        # Make a dict with sub-lists.
        obj_dict = {}
        obj_dict['positive strand'] = pos
        obj_dict['negative strand'] = neg

        # Replace original list with dict.
        nhmmer_hit_objs_dict2[key] = obj_dict

    # Further divide the objects based on which are likely to belong to the
    # same gene.
    for key1 in nhmmer_hit_objs_dict2.keys():
        for key2 in nhmmer_hit_objs_dict2[key1].keys():
            obj_list = nhmmer_hit_objs_dict2[key1][key2]
            obj_sublists = cluster_ali_recur([], obj_list)
            nhmmer_hit_objs_dict2[key1][key2] = obj_sublists

    ## Display contents of dictionary structure.
    display_dict(nhmmer_hit_objs_dict2)   

    # Iterate over all the lists of objects corresponding to potential
    # individual genes, and write consensus sequences to an output fasta file.
    # This step seems to take a really long time for some reason, longer than
    # running the HMMer search in the first place...
    # (Would be nice to sort the output sequences according to e-value somehow)
    gene_num = 0
    with open(seq_filepath, 'w') as o:
        for key1 in nhmmer_hit_objs_dict2.keys():
            for key2 in nhmmer_hit_objs_dict2[key1].keys():
                obj_lists = nhmmer_hit_objs_dict2[key1][key2]
                if len(obj_lists) > 0:
                    # Loop over lists of objects corresponding to exons for
                    # single genes.
                    for l in obj_lists:
                        gene_num += 1

                        # Get concatenated sequence for identified gene.
                        c = get_concat_seq(l, 'nucl')
                        o.write(c)

                        # Get gene region.
                        #s = get_gene_region(l)


                        # Sort objects in hmmer alignment list by start
                        # alignment.
                        



                        # Generate temporary 3-frame translation of gene
                        # region.
                        #temp_filepath = seq_filepath + '_temp_transl' + str(gene_num) + '.fa'
                        #with open(temp_filepath, 'w') as tempoutfh:
                        #    tempoutfh.write(get_three_frame_transl_string(s))

                        # Run hmmer search on 3-frame translations.
                        #tab_hmmsearch_outpath = temp_filepath + '_hmmsearchout' + str(gene_num) + '.txt'
                        #run_hmmsearch(prothmmpath, temp_filepath, tab_hmmsearch_outpath)                           

                        # Parse the hmmsearch output, and extract the relevant
                        # gene models.
                        # ...
                        #hmmer_objs = parse_hmmer_out(tab_hmmsearch_outpath)

                        # Remove temporary file with 3-frame translations.
                        #os.remove(temp_filepath)

                        # Remove temporary hmmsearch output file.
                        #os.remove(tab_hmmsearch_outpath)


                        # get three-frame translations of gene region.
                        #o.write(get_three_frame_transl_string(s))
