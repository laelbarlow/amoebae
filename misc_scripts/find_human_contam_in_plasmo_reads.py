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
"""Script to search Plasmodium paired-end RNA-seq reads against human and
plasmodium coding sequences, and determine what percentage retrieve a human
coding sequence as the top hit.

Reads are ~150bp. 6GB file for each of the two ends.
Search in Plasmodium 3D7 coding sequences.
nHMMer might be faster than blastn, and could work because looking for
identical sequences anyway, but would take time to make seqs into HMMs.
Perhaps deleting search output files right away before creating new ones would
be good, and just record in a very simple manner (update a running tally).
Account for reads that do not retrieve any hits.
Write test code to use fragments of Pfal and Hsap coding sequences as test
queries.
Combine Hsap and Pfal coding sequences into a single blastable database (or not
if using nHMMer).

"""
import sys
import os
import time
import datetime
import subprocess
from Bio import SeqIO

# Import SearchIO and suppress experimental warning
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

from Bio.Blast import NCBIXML

# Current time.
start_time = time.time()
timestamp = time.strftime("%Y%m%d%H%M%S")

cmd = sys.argv
in1 = cmd[1] # FASTA file with queries.
in2 = cmd[2] # BLASTable database.

# Define some function.
def hit_is_human(hit):
    """Take an NCBIXML.parse hit object description and return true if the
    sequence is human (the accession has an underscore as the third character).
    """
    # Have to remove the prefixes in the form "gnl|BL_ORD_ID|0", which blastn
    # adds to the description.
    if str(hit).split(' ', 1)[1][2] == '_':
        return True
    else:
        return False

# Using a different script, create reasonably sized files containing reads to
# use as queries. Perhaps a random sample.

## Open input file.
#with open(in1) as infh:
#    # Parse input fasta file.  
#    seqs = SeqIO.parse(in1, "fasta")
#    # Iterate over sequences in input fasta file.
#    for seq in seqs:



# Run blastn search with the first input as query file and the second as
# database. 
method = 'blastn'
queryfile = in1
dbfile = in2
num_threads = str(4)
blast_evalcut = str(0.05)
outfile = in1.rsplit('.', 1)[0] + '_blastn_out.txt'
# Report on progress.
print('\nSearching with query sequences from %s in the database %s.' % (in1, in2))
# Run the command to search with blastn.
subprocess.call([method, '-query', queryfile, '-db', dbfile, '-out',
        outfile, '-num_threads', num_threads, '-outfmt', '5', '-evalue',
        blast_evalcut])

# Check that output file exists.
assert os.path.isfile(outfile), """No search output file."""

# Report on progress.
print('\nSearched with query sequences from %s in the database %s.' % (in1, in2))

        
# Open and parse search output file.
total_queries = 0
human_hit_tally = 0
#with open(outfile) as infh:
#    # Iterate through each search record (one for each query sequence).
#    for blast_record in NCBIXML.parse(infh):
#        total_queries += 1
#        #for hit in blast_record.descriptions:
#        for hit in blast_record.descriptions:
#            human = hit_is_human(hit)
#            if human:
#                human_hit_tally += 1
#            break
with open(outfile) as infh:
    # Iterate through each search record (one for each query sequence).
    idx = SearchIO.index(outfile, 'blast-xml')
    for key in sorted(idx.keys()):
        total_queries += 1
        for hit in idx[key]:
            human = hit_is_human(hit)
            if human:
                human_hit_tally += 1
            break
    idx.close()

# Remove unnecessary files.
os.remove(outfile)

# End time.
end_time = time.time()
# Record time elapsed.
elapsed = end_time - start_time
print('\nRun time: ' + str(datetime.timedelta(seconds=elapsed)))


# Print results.
print('\nResults:')
print('\nTotal number of query sequences was %s' % str(total_queries))
print('\nNumber of queries that retrieved a human transcript as the top hit was '\
    + str(human_hit_tally) + ' (' + str((human_hit_tally/total_queries)*100) +\
    '%)')
