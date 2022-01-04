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
"""This module contains functions for automatically predicting redundant hits
in the CSV file output by the get_redun_hits command.

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import itertools
import functools
import random


def compare_two_hit_id_lists(ids1, ids2):
    """Take two lists of sequence IDs, and remove ids so that each list
    contains only unique ids.
    """
    # Copy the input lists.
    ids1a = ids1
    ids2a = ids2

    # Make the lists the same length.
    if len(ids1a) > len(ids2a):
        ids2a = ids2a + ["Filler"]*(len(ids1a) - len(ids2a))
    elif len(ids1a) < len(ids2a):
        ids1a = ids1a + ["Filler"]*(len(ids2a) - len(ids1a))
    else:
        pass
    assert len(ids1a) == len(ids2a)

    # Zip through the lists and generate non-redundant lists.
    ids1a2 = []
    ids2a2 = []
    stop_adding_to_ids1a2 = False
    stop_adding_to_ids2a2 = False
    for x, y in zip(ids1a, ids2a):
        if not stop_adding_to_ids1a2:
            if not x in ids2a2 and not x == "Filler":
                ids1a2.append(x)
            else:
                stop_adding_to_ids1a2 = True
        if not stop_adding_to_ids2a2:
            if not y in ids1a2 and not y == "Filler":
                ids2a2.append(y)
            else:
                stop_adding_to_ids2a2 = True
        if stop_adding_to_ids1a2 and stop_adding_to_ids2a2:
            break

    # Remove all instances of "Filler" from both lists.
    ids1a2 = list(filter(lambda a: a != "Filler", ids1a2))
    ids2a2 = list(filter(lambda a: a != "Filler", ids2a2))

    # Return the new lists.
    return (ids1a2, ids2a2)


def get_nonredun_dict_key_pairs(in_dict):
    """Take a dict and return a list of all unique two-key combinations of keys
    of the dict.
    """
    # Don't modify the dict if there is only one key.
    if len(in_dict.keys()) == 1:
        return in_dict

    # Copy the input dict.
    in_dict2 = in_dict.copy()

    # Get a list of lists of all unique combinations of keys in the input dict.
    two_key_combos1 = list(itertools.permutations(in_dict2.keys(), 2))
    two_key_combos2 = [sorted(list(x)) for x in two_key_combos1]
    two_key_combos3 = []
    for elem in two_key_combos2:
        if elem not in two_key_combos3:
            two_key_combos3.append(elem)

    # Return list of uniqe combos.
    return two_key_combos3


def process_query_title_hit_id_list_dict(in_dict):
    """Take a dict with keys as query titles and values as lists of sequence
    IDs"""
    # Don't modify the dict if there is only one key.
    if len(in_dict.keys()) == 1:
        return in_dict

    # Copy the input dict.
    in_dict2 = in_dict.copy()

    # Get a list of lists of all unique combinations of keys in the input dict.
    two_key_combos = get_nonredun_dict_key_pairs(in_dict2)

    ## Check whether the number of combos is super big.
    #assert len(two_key_combos) <= 500, """Attempting to iterate over more than
    #500 combinations of hit ID lists for different query titles. This is
    #probably impractical."""

    # Iterate over all pairs of keys and generate a new dict.
    for combo in two_key_combos:

        # Retrieve lists from input dict.
        ids1 = in_dict2[combo[0]]
        ids2 = in_dict2[combo[1]]

        # Get non-redundant lists.
        if not set(ids1).isdisjoint(ids2):
            ids1a, ids2a = compare_two_hit_id_lists(ids1, ids2)

            # Assign non-redundant lists to values in output dict.
            in_dict2[combo[0]] = ids1a
            in_dict2[combo[1]] = ids2a

    # Return new dict with nonredundant lists/values.
    nonredun_dict = in_dict2
    return nonredun_dict


def parse_input_redun_hit_csv_line(l):
    sl = l.split(',')
    query_title = sl[0].strip('\"')
    hit_id = sl[9]
    evalue = sl[6]
    return (query_title, hit_id, evalue)


def cmp_paths(lhs, rhs):
    return 0 if any(key in lhs for key in rhs) else -1


def cluster_overlapping_query_title_lists(query_title_lists):
    """Take a list of lists of query titles and return a list of lists of
    query titles in which the sublists are non-overlapping.
    """
    # Make each sub-list nonredundant.
    query_title_lists = [list(set(x)) for x in query_title_lists]

    sublists_overlapping = True
    while sublists_overlapping:
        # Cluster query title lists together in sublists.
        random.shuffle(query_title_lists)
        paths = query_title_lists
        key_paths = functools.cmp_to_key(cmp_paths)
        groups = [list(group) for key, group in itertools.groupby(paths, key_paths)]

        # Chain the groups together.
        chained_groups = []
        for x in groups:
            xlist = []
            for y in x:
                xlist = xlist + y
            chained_groups.append(xlist)

        # Make the chained subgroups nonredundant.
        nonredun_chained_groups = [list(set(x)) for x in chained_groups]

        # Check whether the sublists are overlapping.
        total_list = []
        for l in nonredun_chained_groups:
            total_list = total_list + l
        if len(total_list) == len(list(set(total_list))):
            sublists_overlapping = False

        # Update the list of query title lists.
        query_title_lists = nonredun_chained_groups

    # Return a new list of query title lists.
    return query_title_lists


def generate_dicts_from_redun_hit_csv(in_csv):
    """Take the output CSV file from get_redun_hits command and return a dict
    with query titles as keys and lists of sequence IDs as values.
    """
    out_dict = {}
    id_dict = {}
    with open(in_csv) as infh:
        inum = 0
        for i in infh:
            inum += 1
            if not inum == 1:
                # Add to out_dict.
                query_title, hit_id, evalue = parse_input_redun_hit_csv_line(i)
                hit_plus_evalue = (hit_id, evalue)
                if query_title not in out_dict.keys():
                    out_dict[query_title] = [hit_plus_evalue]
                else:
                    out_dict[query_title] = out_dict[query_title] + [hit_plus_evalue]

                # Add to id_dict.
                if hit_id not in id_dict.keys():
                    id_dict[hit_id] = [query_title]
                else:
                    id_dict[hit_id] = id_dict[hit_id] + [query_title]

    # Sort all the lists (values) by ascending E-value, and remove the
    # E-values.
    for query_title in out_dict.keys():
        #print(out_dict[query_title])
        #print('\n\n')
        out_dict[query_title] = sorted(out_dict[query_title], key=lambda x: float(x[1]))
        #print(out_dict[query_title])
        #print('\n\n')
        out_dict[query_title] = [x[0] for x in out_dict[query_title]]
        #print(out_dict[query_title])
        #print('\n\n')

    # Divide dict into separate dicts with overlapping elements in their lists
    # (values).
    out_dicts = []
    query_title_lists = list(id_dict.values())
    nonredun_chained_groups = cluster_overlapping_query_title_lists(query_title_lists)
    #for x in nonredun_chained_groups:
    #    print('\n\n')
    #    print(x)
    #assert 2!=2
    for query_title_list in nonredun_chained_groups:
        out_dict_x = {}
        for query_title in query_title_list:
            for key in out_dict.keys():
                if key == query_title:
                    out_dict_x[query_title] = out_dict[query_title]
        out_dicts.append(out_dict_x)

    # Return the dict with just ordered lists of sequence IDs.
    return out_dicts


def get_csv_with_redun_hit_predictions(in_csv, out_csv):
    """Take a CSV file and write a new CSV file to the given path with '+' in
    the appropriate rows.
    """
    # Parse the input CSV file, and generate a dict of query titles and
    # associated IDs.
    dicts1 = generate_dicts_from_redun_hit_csv(in_csv)

    # Remove ids from the lists so that they are non-redundant.
    dicts2 = []
    for d in dicts1:
        dicts2.append(process_query_title_hit_id_list_dict(d))
    xnum = 0

    # Put the dicts together in one big dict with all the query titles as keys.
    final_dict = {}
    for d in dicts2:
        for k in d.keys():
            print(k)
            assert k not in final_dict.keys()
            final_dict[k] = d[k]

    # Write new CSV file.
    with open(in_csv) as in_csvh, open(out_csv, 'w') as o:
        inum = 0
        for i in in_csvh:
            inum += 1
            if not inum == 1:
                query_title, hit_id, evalue = parse_input_redun_hit_csv_line(i)
                if hit_id in final_dict[query_title]:
                    si = i.split(',')
                    o.write(','.join(si[:4] + ['+'] + si[5:]))
                else:
                    o.write(i)
            else:
                o.write(i)


