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
"""Module for dacksify_pos_hmmer_hits.py
"""
import sys
import time
from datapaths import DataPaths
import math
import re
import os


def get_rblast_string(relev_blast_res, red_acc_file):
    """Takes relevant blast results, and generates a string for writing to the
    output spreadsheet.
    """
    string = None
    
    # Define features of top rblast hit.
    top_hit = relev_blast_res[0][3]
    top_hit_eval_raw = relev_blast_res[0][2]
    top_hit_eval = relev_blast_res[0][2]
    if float(top_hit_eval) == 0:
        top_hit_eval = 1e-300

    # Determine whether the top hit has a redun acc, and calculate certain
    # values if it does.
    if is_redun_acc(red_acc_file, relev_blast_res[0][1]):

        top_hit = relev_blast_res[0][3]
        top_hit_eval_raw = relev_blast_res[0][2]
        top_hit_eval = relev_blast_res[0][2]
        if float(top_hit_eval) == 0:
            top_hit_eval = 1e-300

        last_redun_hit_eval = relev_blast_res[1][2]
        if float(last_redun_hit_eval) == 0:
            last_redun_hit_eval = 1e-300

        if relev_blast_res[2] != [None]: # i.e. if there are nonredundant hits.
            # Define first nonredundant hit evalue.
            first_nonred_eval = relev_blast_res[2][2]
            if float(first_nonred_eval) == 0:
                first_nonred_eval = 1e-300

            # Determine diff btw top hit and first nonredundant hit.
            eval_dif1 = abs(math.log(float(top_hit_eval), 10)) -\
            abs(math.log(float(first_nonred_eval), 10))

            # Determine diff btw last redundant hit and first nonredundant hit.
            eval_dif2 = abs(math.log(float(last_redun_hit_eval), 10)) -\
            abs(math.log(float(first_nonred_eval), 10))

            # Define the first nonredundant hit sequence description.
            first_nonredun_descr = relev_blast_res[2][3]

            # Define the string to return.
            string = top_hit + ',' + str(top_hit_eval_raw) + ',' + str(eval_dif1)\
            + ',' + str(eval_dif2) + ',' + str(first_nonredun_descr)

        else:
            string = top_hit + ',' + str(top_hit_eval_raw) + ',' + 'N/A'\
            + ',' + 'N/A' + ',' + 'N/A' 
    else:
        string = top_hit + ',' + str(top_hit_eval_raw) + ',' + 'N/A'\
            + ',' + 'N/A' + ',' + '(Top blast hit is not the original query)' 
        
    return string


def find_relev_blast(blast_res, red_acc_file):
    """Finds the hit corresponding to the original query, and the first hit
    that is not redundant with the original query in the rblast searches.

    Returns a list composed of lists for the hit corresponding to the original
    query and the first result that is not redundant with the original query.
    """
    x = None

    if blast_res is not None:
        ## TEMPORARY for determing whether the top hit has a redun acc.
        #if is_redun_acc(red_acc_file, blast_res[0][0]):
        #    print('***Yes')
        #else:
        #    print('***No')

        # First hit (redundant with original query, because HMMer_pipeline).
        orig_query_sublist = [1] + blast_res[0] 

        # Find first blast hit that is not redundant with the original query, and
        # the last redundant hit.
        first_non_red_sublist = [None]
        last_red_sublist = [None]
        len_blast_res = len(blast_res)
        res_num = 0
        for sublist in blast_res:
            res_num += 1
            if not is_redun_acc(red_acc_file, sublist[0].strip()):
                #sublist is the first non-redundant hit
                first_non_red_sublist = [res_num] + sublist
                #print('first non-redundant blast result is ' +\
                #        str(first_non_red_sublist))

                #determine the sublist corresponding to the last redundant hit
                last_red_sublist = [res_num - 1] + blast_res[res_num - 2]
                #print('last redundant blast result is ' + str(last_red_sublist))

                break
                
            elif res_num == len_blast_res:
                #reached last hit without finding any nonredun hits set first non-redun hit as None.
                #print('last blast result is redundant with original query, and is ' + str(sublist))
                first_non_redun_sublist = [None]

                #determine the sublist corresponding to the last redundant hit, the
                # last hit in this case.
                last_red_sublist = [res_num] + blast_res[res_num -1]
                #print('last (redundant) blast result is ' + str(last_red_sublist))

        # Return the relevant lists in a list.
        x = [orig_query_sublist, last_red_sublist, first_non_red_sublist]

    else:
        pass

    return x

def get_acc_for_prot_name(prot_name, key_file):
    """Searches in a table file for the chosen representative (refseq?)
    accession for the original (Hsap) query that corresponds to the given
    protein name.
    """
    with open(key_file) as keyfilehandle:
        for i in keyfilehandle:
            x = ''
            acc = i.split(',')[0].strip()
            if acc == prot_name.strip():
                x = i.split(',')[1]
        return x


def get_rblast_list(linestring):
    """Takes a line from the HMMer_pipeline output and returns a list of lists
    containing the info for a reverse blast hit.
    """
    # Split the string.
    split = linestring.rstrip('\n').split(',', maxsplit=6)

    x = None

    # Only return an rblast string if there are rblast results.
    if len(split) > 6:

        # Define the portion of the line that contains the rblast results.
        blast_res = split[6]
        
        # Identify fields with internal commas, and remove the internal commas.
        matches = re.findall(r'\"(.+?)\"', blast_res)
        for match in matches:
            blast_res = blast_res.replace(match, match.replace(',', ''))

        # Split the modified line on commas.
        blast_res_split = blast_res.split(',')

        len_res = len(blast_res_split)
        x = [blast_res_split[i:i+3] for i in list(range(100))[0:len_res:3]]

    return x 


def is_redun_acc(redun_acc_file, acc):
    """Returns True if a given accession is in a list of accessions appearing
    in a given file.
    """
    with open(redun_acc_file) as infilehandle:
        line = infilehandle.read()
        redun_list = line.rstrip('\n').split(',') 
        #print('\n\n' + str(redun_list) + '\n\n')

        result = False
        for i in redun_list:
            if acc.strip() == i.strip():
                #print(acc.strip() + ' == ' + i.strip())
                result = True
                break
            else:
                #print(acc.strip() + ' != ' + i.strip())
                pass
        return result 


def find_unabbrev_name(abbrev, tablefilepath):
    """Takes an abbreviation, and returns the full name.
    The table file is a csv file based on the table of genome databases
    accessed.
    """
    # print(abbrev) # Temporary
    with open(tablefilepath) as table:
        x = ''
        for i in table:
            isplit = i.split(',')
            if isplit[3].strip() == abbrev.strip():
                x = isplit[2]
                break

        assert x != '', "Error: Unable to identify abbreviation in table."
        return x 


def get_red_acc_listx(infilepath):
    """Takes a file that contains redundant accessions and returns a list of
    the accessions it contains.
    """
    # Get redundant accessions from file.
    infilehandle = open(infilepath)
    l = []
    line_num = 0
    for line in infilehandle:
        line_num += 1
        if line_num == 1:
            l = line.rstrip('\n').rstrip(',').split(',')
        else:
            print('More than one line in redundant accessions file.')
    infilehandle.close()

    # Check whether redundant accessions identified.
    if l == []: print('No redundant accessions identified.')
    return l


def query_retrieved(hmmer_string, red_acc_file):
    """Returns true if a match for the original query was retrieved in the
    rblast results as a top hit.
    """
    red_accs = get_red_acc_listx(red_acc_file)
    
    split = hmmer_string.split(',')

    if len(split) > 6:
        if split[6].strip() in red_accs:
            return True


def dacksify(infp, outfp, tablefilepath, red_acc_file, evaldiffset):
    """Dacksifies
    """
    outfile = open(outfp, 'w')

    with open(infp) as infilehandle, open(outfp, 'w') as outfilehandle:
        # Define a header line for output spreadsheet.
        header_string = """Protein name, Reverse BLAST database, HMM name, Subject organism, Positive hit \
accession, Positive hit HMMer Evalue, Top reverse BLAST hit in query database, Reverse BLAST E-value, Order of \
magnitude e-value difference top hit and first negative hit, Order of \
magnitude e-value difference last positive hit and first negative hit,\
First negative reverse BLAST hit""" 
        outfilehandle.write(header_string)
        outfilehandle.write('\n')

        # Get rblast database name.
        #rblastdbname = os.path.basename(red_acc_file).split('_',
        #1)[0].strip('*')

        l = 0
        for i in infilehandle:
            l += 1

            # Only process lines containing data.
            if not i.startswith('Forward'):
                if not i.startswith('\n'):

                    # Make sure that a match for the query was retrived in the
                    # rBLAST. This is important for dacksifying raw
                    # HMMer_pipeline.py outputs.
                    if query_retrieved(i, red_acc_file): 

                        # Split the line on commas. 
                        spliti = i.split(',')

                        # Determine the protein name from the HMM name.
                        protein_name = spliti[0].split('_', 1)[0]

                        # Get a list of blast result accessions and evalues.
                        blast_res = get_rblast_list(i)

                        # Find the relevant rblast results.
                        relev_blast_res = find_relev_blast(blast_res, red_acc_file)  

                        # Make blast result string for output
                        rblast_string = get_rblast_string(relev_blast_res, red_acc_file)

                        # Replace abbrev names with full names. 
                        spliti[2] = find_unabbrev_name(spliti[2], tablefilepath)
                        rblastdbname = find_unabbrev_name(spliti[1], tablefilepath)

                        # Remove hit number and rblast hits.
                        spliti = [spliti[0], spliti[2], spliti[4], spliti[5]]

                        # Join the list with back into a string with commas.
                        i = ','.join(spliti)

                        # Add back the selected rblast results to the string, and add the
                        # protein name as the first column.
                        i = protein_name + ',' + rblastdbname + ',' + i + ',' + rblast_string + '\n'

                        # Apply e-value difference criterion (only write lines
                        # if they the hit meets the criterion). This assumes
                        # that there are no forward or top reverse hits that
                        # were retrieved with evalues above the threshold (see
                        # find_pos_hmmer_hits.py if not).
                        evaldiff = i.split(',')[8].strip()
                        if evaldiff == 'N/A' or float(evaldiff) >= evaldiffset:

                            # Write the string as a line in output.
                            outfilehandle.write(i)


                    # If the query was not retrieved, then proceed differently.***
                    else: 

                        # Split the line on commas. 
                        spliti = i.split(',')

                        # Determine the protein name from the HMM name.
                        protein_name = spliti[0].split('_', 1)[0]

                        # Get a list of blast result accessions and evalues.
                        blast_res = get_rblast_list(i)

                        # Find the relevant rblast results.
                        relev_blast_res = find_relev_blast(blast_res, red_acc_file)  

                        # If no blast hits, then fill in the fields.
                        rblast_string = '(No hits),N/A,N/A,N/A,N/A' 

                        # If there are blast results, then fill in the fields.
                        if blast_res is not None:
                            # Make blast result string for output
                            rblast_string = get_rblast_string(relev_blast_res, red_acc_file)

                        # Replace abbrev names with full names. 
                        spliti[2] = find_unabbrev_name(spliti[2], tablefilepath)
                        rblastdbname = find_unabbrev_name(spliti[1], tablefilepath)

                        # Remove hit number and rblast hits.
                        spliti = [spliti[0], spliti[2], spliti[4],\
                                spliti[5].rstrip()]

                        # Join the list with back into a string with commas.
                        i = ','.join(spliti)

                        # Add back the selected rblast results to the string, and add the
                        # protein name as the first column.
                        i = protein_name + ',' + rblastdbname + ',' + i + ',' + rblast_string + '\n'

                        # Apply e-value difference criterion (only write lines
                        # if they the hit meets the criterion). This assumes
                        # that there are no forward or top reverse hits that
                        # were retrieved with evalues above the threshold (see
                        # find_pos_hmmer_hits.py if not).

                        # evaldiff = i.split(',')[8].strip()
                        # if evaldiff == 'N/A' or float(evaldiff) >= evaldiffset:

                        # Write the string as a line in output.
                        outfilehandle.write(i)
