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
"""Contains functions used in search_scaffolds.py.
"""

import sys
import os
sys.path.append(os.path.dirname(sys.path[0]))
import copy
import subprocess
from Bio import SeqIO
#from Bio import SearchIO
# Import SearchIO and suppress experimental warning
import warnings
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO

from Bio.Blast import NCBIXML
from Bio.Blast import Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import itertools
import settings
from run_exonerate import ExonerateLocusResult, get_subseq_from_nucl,\
run_exonerate_as_subprocess


def check_if_two_hsp_ranges_overlap(lists):
    """Takes two lists which contain the start and stop positions of two HSPs,
    and returns True if they overlap.

    Note: This assumes that the SearchIO coordinates are provided using Python
    string slicing conventions (this is suggested to be the case here:
    https://biopython.org/DIST/docs/api/Bio.SearchIO-module.html).

    """
    x = set(range(lists[0][0], lists[0][1]))
    y = set(range(lists[1][0], lists[1][1]))
    intersect = x.intersection(y)

    overlap = False
    if len(list(intersect)) > 0:
        overlap = True

    return overlap


def check_if_two_hsps_overlap(hsp_obj1, hsp_obj2):
    """Takes two HSP objects and returns True if their ranges overlap, and
    False if not.
    """
    x = hsp_obj1.hit_range
    y = hsp_obj2.hit_range

    overlap = check_if_two_hsp_ranges_overlap([x, y])

    return overlap


def check_if_hsps_overlapping(hsp_objs):
    """Takes a list of HSP objects, and returns true if they overlap.
    """
    # Compile a list of the ranges.
    ranges = []
    for hsp in hsp_objs:
        r = hsp.hit_range
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


def check_strands(hsp_objs, fwd_strand):
    """Return True if hsps in the list are not all on the same strand.
    """
    hsps_on_different_strands = False
    for hsp in hsp_objs:
        x = True
        if hsp.hit_frame < 0:
            x = False
        if x != fwd_strand:
            hsps_on_different_strands = True

    return hsps_on_different_strands


def whether_fwd_strand(hsp):
    if hsp.hit_frame < 0:
        return False
    else:
        return True


def check_for_duplicate_ranges(range_list):
    """Return True if the given list of tuples contains duplicates.
    """
    non_redun_range_list = list(set(range_list))
    x = None
    if len(non_redun_range_list) < len(range_list):
        x = True
    else:
        x = False

    return x


def sort_left_ranges(left_ranges):
    return list(reversed(sorted(left_ranges, key=lambda x: x[1])))


def reduce_left_ranges(sorted_left_ranges, top_hsp_range, max_gap):
    # Reduce the list of left ranges to only those that meet the criteria.
    reduced_left_ranges = []
    last_range = top_hsp_range
    for r in sorted_left_ranges:
        # Only include range if it ends before last range starts and
        # within max_gap.
        if (last_range[0] - max_gap - 2) < r[1] < last_range[0]:
            reduced_left_ranges.append(r)
            last_range = r
        else:
            break

    return reduced_left_ranges
    

def sort_right_ranges(right_ranges):
    return sorted(right_ranges, key=lambda x: x[0])


def reduce_right_ranges(sorted_right_ranges, top_hsp_range, max_gap):
    # Reduce the list of right ranges to only those that meet the criteria.
    reduced_right_ranges = []
    last_range = top_hsp_range
    for r in sorted_right_ranges:
        # Only include range if it starts after last range ends, and within max_gap.
        if (last_range[1] + max_gap + 2) > r[0] > last_range[1]:
            reduced_right_ranges.append(r)
            last_range = r

    return reduced_right_ranges


def get_proximate_hsp_ranges(range_list, top_hsp_range, max_gap):
    """Takes a list of ranges and returns another list, but corresponding to
    the top HSP and proximate HSPs that do not overlap.

    Potential problem: Does not take into account the E-value of overlapping
    HPSs, and may exclude some HSPs with better E-values because of their
    position. It would perhaps be better to take E-values into account?
    """
    # Compile a list of proximate HSP ranges upstream (left regardless of which
    # way 5' or 3' is) of the range for the top HSP.
    left_ranges = [x for x in range_list if x[1] < top_hsp_range[0]]

    reduced_left_ranges = []

    if len(left_ranges) > 0:

        # Sort the ranges in descending order of their higher value.
        sorted_left_ranges = sort_left_ranges(left_ranges)

        # Reduce the list of left ranges to only those that meet the criteria.
        reduced_left_ranges = reduce_left_ranges(sorted_left_ranges, top_hsp_range, max_gap)

    # Compile a list of proximate HSPs downstream (right regardless of which
    # way 5' or 3' is) of top HSP.
    right_ranges = [x for x in range_list if x[0] > top_hsp_range[1]]

    reduced_right_ranges = []

    if len(right_ranges) > 0:

        # Sort the ranges in ascending order of their lower value.
        sorted_right_ranges = sort_right_ranges(right_ranges)

        # Reduce the list of right ranges to only those that meet the criteria.
        reduced_right_ranges = reduce_right_ranges(sorted_right_ranges,
                top_hsp_range, max_gap)

    # Concatenate left, top, and right range lists to make the final list.
    final_proximate_ranges = list(reversed(reduced_left_ranges))\
                             + [top_hsp_range]\
                             + reduced_right_ranges

    # Compile list of ranges that fall within the range of HSPs in the final
    # list, but were not included because they overlap with higher ranking hsps
    # in the final list. 
    leftover_proximate_ranges = []
    for r in range_list:
        if r not in final_proximate_ranges:
            for fpr in final_proximate_ranges:
                if check_if_two_hsp_ranges_overlap([r, fpr]):
                    leftover_proximate_ranges.append(r)

    # Return the final list of ranges.
    return [final_proximate_ranges, leftover_proximate_ranges]

    
def remove_hsps_with_redun_ranges(range_list, hsp_objs):
    """Remove HSPs with redundant ranges.
    """
    non_redun_hsps = []

    # Get list of ranges of which there are duplicates.
    duplicate_ranges = []
    for r in range_list:
        if range_list.count(r) > 1 and r not in duplicate_ranges:
            duplicate_ranges.append(r)

    # Get corresponding list of HSP objects.
    hsps_with_duplicate_ranges = []
    for r in duplicate_ranges:
        hsp_sublist = []
        for h in hsp_objs:
            if h.hit_range == r:
                hsp_sublist.append(h)
        hsps_with_duplicate_ranges.append(hsp_sublist)

    # For each list of HSPs with identical ranges, identify top HSP.
    top_hsps_from_sublists = []
    for hsp_sublist in hsps_with_duplicate_ranges:
        top_hsp_in_sublist = sorted(hsp_sublist, key=lambda x: x.evalue)[0]
        top_hsps_from_sublists.append(top_hsp_in_sublist)

    # Remove the duplicate HSPs from the complete list.
    for h in hsp_objs:
        # If range is in the duplicate ranges, only append if it is the top
        # HSP.
        append = False
        if h.hit_range in duplicate_ranges:
            if h in top_hsps_from_sublists:
                append = True
        else:
            append = True

        if append:
            non_redun_hsps.append(h)

    return non_redun_hsps


def reduce_to_best_hsps(sorted_hsps, max_gap, recursion_num, name, leftright):
    """Recursively select HSPs upstream of the top HSP for inclusion in the
    cluster.

    Note: This function may be inefficient.
    """
    ## Check that no HSPs with identical hit_ranges exist in the input list.
    #for hsp1 in sorted_hsps:
    #    for hsp2 in sorted_hsps:
    #        if not hsp1 == hsp2:
    #            assert hsp1.hit_range != hsp2.hit_range, """HSPs with identical
    #            hit ranges exist in input hsp list."""
    
    ## Display.
    #print('\nHSPs input to reduce_to_best_hsps R' + str(recursion_num) + ' ' +\
    #name + ':')
    #print('Complete list:')
    #for hsp in sorted_hsps:
    #    print(hsp.hit_range)
    #print('\n')

    reduced_sorted_hsps = []
    discard_list = []
    # Iterate over the HSPs in the input list, and find up to one HSP to remove.
    hsp_index = -1
    for hsp in sorted_hsps: 
        hsp_index += 1
        current_hsp = sorted_hsps[hsp_index]
        if len(sorted_hsps) > hsp_index + 1:
            # Consider the most proximal hsp to the first hsp.
            next_hsp = sorted_hsps[hsp_index + 1]
            
            # If the next hsp does not represent an appropriate upstream or
            # downstream portion of the query sequence (depending on whether on
            # the plus or minus strand, upstream or downstream of top hsp),
            # then discard it.
            stream = None
            if current_hsp.hit_frame > 0:
                # On plus strand.
                #print('on plus strand')
                if leftright == 'left':
                    #print('going left')
                    # Subsequent HSPs should be upstream of the current HSP.
                    stream = 'up'
                else:
                    #print('going right')
                    # Subsequent HSPs should be downstream of the current HSP.
                    stream = 'down'
            elif current_hsp.hit_frame < 0:
                # On minus strand.
                #print('on minus strand')
                if leftright == 'left':
                    #print('going left')
                    # Subsequent HSPs should be downstream of the current HSP.
                    stream = 'down'
                else:
                    #print('going right')
                    # Subsequent HSPs should be upstream of the current HSP.
                    stream = 'up'
            # Check whether next hsp is appropriately upstream or downstream in
            # terms of its query range.
            if stream == 'up':
                # Compare the end positions of the HSPs in the query sequence
                # as an indication of whether they are in an appropriate order.
                if next_hsp.query_end >= current_hsp.query_end:
                    # Then not upstream as it should be.
                    discard_list.append(next_hsp)
                    # Break the loop so that the HSP can be removed before
                    # proceeding with the next recursion.
                    break
                else:
                    pass

                    # TEMPORARY
                    #assert not next_hsp.hit_start == 8467356
                    #assert not next_hsp.hit_start == 8466585
                    #assert current_hsp.hit_start != 8466585 or\
                    #next_hsp.hit_start != 8467356
                    #assert current_hsp.hit_start != 8467356 or\
                    #next_hsp.hit_start != 8466585

            elif stream == 'down':
                # Compare the start positions of the HSPs in the query sequence
                # as an indication of whether they are in an appropriate order.
                if next_hsp.query_start <= current_hsp.query_start:
                    # Then not downstream as it should be.
                    discard_list.append(next_hsp)
                    # Break the loop so that the HSP can be removed before
                    # proceeding with the next recursion.
                    break

                    # TEMPORARY
                    #print('\n')
                    #print(current_hsp.hit_range)
                    #print('current_hsp.query_start: ' + str(current_hsp.query_start))
                    #print(next_hsp.hit_range)
                    #print('next_hsp.query_start: ' + str(next_hsp.query_start))
                    #assert not next_hsp.hit_start == 8467356
                    #assert not next_hsp.hit_start == 8466585
                    #assert current_hsp.hit_start != 8466585 or\
                    #next_hsp.hit_start != 8467356
                    #assert current_hsp.hit_start != 8467356 or\
                    #next_hsp.hit_start != 8466585



            # If the first hsp overlaps with the next hsp, then keep the hsp that
            # has the better E-value.
            if check_if_two_hsps_overlap(next_hsp, current_hsp):
                ## Display.
                #print(str(current_hsp.hit_range) + ' overlaps with ' +\
                #        str(next_hsp.hit_range))

                if next_hsp.evalue < current_hsp.evalue:
                    discard_list.append(current_hsp)

                    ## Display.
                    #print(str(next_hsp.hit_range) + ' ' + str(next_hsp.evalue) +\
                    #' < ' + str(current_hsp.hit_range) + ' ' +\
                    #str(current_hsp.evalue)) 
                    #print('so discarding ' + str(current_hsp.hit_range))

                elif next_hsp.evalue == current_hsp.evalue:
                    # Use bitscore to break the tie.
                    if next_hsp.bitscore > current_hsp.bitscore:
                        discard_list.append(current_hsp)
                    else:
                        discard_list.append(next_hsp)
                        # TEMPORARY
                        #assert not next_hsp.hit_start == 8467356
                        #assert not next_hsp.hit_start == 8466585
                        #assert current_hsp.hit_start != 8466585 or\
                        #next_hsp.hit_start != 8467356
                        #assert current_hsp.hit_start != 8467356 or\
                        #next_hsp.hit_start != 8466585

                    #assert next_hsp.bitscore != current_hsp.bitscore,\
                    #"""Houston, we have a problem: Which HSP to keep?."""

                else:
                    # The current HSP has the better E-value, so discard the
                    # other one.
                    discard_list.append(next_hsp)

                    # TEMPORARY
                    #assert not next_hsp.hit_start == 8467356
                    #assert not next_hsp.hit_start == 8466585
                    #assert current_hsp.hit_start != 8466585 or\
                    #next_hsp.hit_start != 8467356
                    #assert current_hsp.hit_start != 8467356 or\
                    #next_hsp.hit_start != 8466585

                    ## Display.
                    #if next_hsp.query_range == (274,341):
                    #    print('YYYYYYYYYYY')
                    #print(str(next_hsp.hit_range) + ' ' + str(next_hsp.evalue) +\
                    #' > ' + str(current_hsp.hit_range) + ' ' +\
                    #str(current_hsp.evalue)) 
                    #print('so discarding ' + str(next_hsp.hit_range))

                # Break the loop so that the HSP can be removed before
                # proceeding with the next recursion.
                break
            # If they do not overlap, then keep both, and move onto the next
            # HSP in the input list, unless the next HSP is too far away from
            # the current HSP.
            else:
                # Check that the next HSP is within max_gap of the current
                # HSP. If not, then discard it and all remaining HSPs that are
                # even further.

                ## Display.
                #print(str(current_hsp.hit_range) + ' does not overlap with ' +\
                #        str(next_hsp.hit_range))

                gap_list = [abs(current_hsp.hit_range[0] - next_hsp.hit_range[1]),
                            abs(next_hsp.hit_range[0] - current_hsp.hit_range[1])
                           ]
                if min(gap_list) > max_gap:
                    for i in sorted_hsps[hsp_index + 1:]:
                        discard_list.append(i) 

                    ## Display.
                    #print(str(next_hsp.hit_range) + ' is not within ' +\
                    #str(max_gap) + ' bp of ' + str(current_hsp.hit_range)\
                    #+ ' so removing ' + str(next_hsp.hit_range) +\
                    #' and all further HSPs.')

                    # Break the loop so that the HSP can be removed before
                    # proceeding with the next recursion.
                    break
                else:
                    ## Display.
                    #print(str(next_hsp.hit_range) + ' is within ' +\
                    #str(max_gap) + ' bp of ' + str(current_hsp.hit_range)\
                    #+ ' so keeping ' + str(next_hsp.hit_range))

                    pass

    ## Display.
    #for hsp in discard_list:
    #    print('\tremoved ' + str(hsp.hit_range))


    # If there are no more HSPs to remove, then return the list. 
    if len(discard_list) == 0:
        ## Display.
        #print('No more HSPs to remove.')

        return sorted_hsps 

    # Otherwise, remove the HSPs and do another recursion with the reduced
    # list.
    else:
        # Generate a reduced list of HSPs.
        for i in sorted_hsps:
            append = True
            for j in discard_list:
                # Note: Sometimes different HSPs have identical hit ranges.
                if i.hit_range == j.hit_range and i.query_range == j.query_range:
                    append = False
            if append:
                reduced_sorted_hsps.append(i)
        # Do another recursion on the reduced list.
        return reduce_to_best_hsps(reduced_sorted_hsps, max_gap, recursion_num\
                + 1, name, leftright)


def find_missed_hsp_ranges(complete_ranges,
                     final_proximate_ranges,
                     leftover_proximate_ranges):
    """Testable function to test for problems.
    """
    missed_hsps = False

    for r in complete_ranges:
        total_prox = final_proximate_ranges + leftover_proximate_ranges
        if r not in total_prox:
            rset = set(range(min([x[0] for x in total_prox]), max([x[1] for x in total_prox]))) 
            rset2 = set(range(min([x[0] for x in final_proximate_ranges]), max([x[1] for x in final_proximate_ranges]))) 
            #if r[0] in rset or r[1] in rset:
            in_range = False
            for x in range(r[0], r[1]):
                if x in rset2:
                    in_range = True
            if in_range:
                print('missed hsp')
                print(r)
                missed_hsps = True

    return missed_hsps


# Obsolete function:
def get_best_proximate_hsps(hsp_objs, top_hsp, max_gap):
    """Return a list of HSP objects that are sufficiently proximate to the top
    HSP, and among which there is no overlap in the sequence ranges. And, also
    return a list of HSP objects that are also proximate to these HSPs, but
    were not included in the final list of HSPs because they overlap with the
    ones that were. Instead of simply selecting HSPs based on hit range,
    consider E-values, favoring HSPs with better E-values if HSPs overlap and
    some have to be excluded.
    """
    # Compile a list of proximate HSPs upstream (left regardless of which
    # way 5' or 3' is) of the the top HSP.
    left_hsps = [x for x in hsp_objs if x.hit_range[0] < top_hsp.hit_range[1]]
    # Remove top HSP from list.
    left_hsps.remove(top_hsp)

    reduced_left_hsps = []

    if len(left_hsps) > 0:

        # Sort the hsps in descending order of their 3' position.
        sorted_left_hsps = sorted(left_hsps, key=lambda x: x.hit_range[1], reverse=True)

        # Reduce the list of left ranges to only those that meet the criteria.
        reduced_left_hsps = reduce_to_best_hsps([top_hsp] + sorted_left_hsps,
                max_gap, 0, 'Left', 'left')
        assert len(reduced_left_hsps) <= len([top_hsp] + sorted_left_hsps)

    # Compile a list of proximate HSPs downstream (right regardless of which
    # way 5' or 3' is) of top HSP.
    right_hsps = [x for x in hsp_objs if x.hit_range[1] > top_hsp.hit_range[0]]
    # Remove top HSP from list.
    right_hsps.remove(top_hsp)

    reduced_right_hsps = []

    if len(right_hsps) > 0:

        # Sort the HSPs in ascending order of their 5' position.
        sorted_right_hsps = sorted(right_hsps, key=lambda x: x.hit_range[0], reverse=False)

        # Reduce the list of right HSPs to only those that meet the criteria.
        reduced_right_hsps = reduce_to_best_hsps([top_hsp] + sorted_right_hsps,
                max_gap, 0, 'Right', 'right')
        assert len(reduced_right_hsps) <= len([top_hsp] + sorted_right_hsps)

    # Check that all HSPs in the input list are accounted for in the left
    # and/or right lists.
    assert len(hsp_objs) <= (len(left_hsps) + 1 + len(right_hsps)), """Did not
    account for all HSPs in input list."""

    # Concatenate left, top, and right HSP lists to make the final list.
    final_proximate_hsps = list(reversed(reduced_left_hsps[1:]))\
                             + [top_hsp]\
                             + reduced_right_hsps[1:]

    ## Unnecessary???
    ## Check that top hsp was included.
    #found = False
    #for hsp in final_proximate_hsps:
    #    if hsp.hit_range == top_hsp.hit_range:
    #        found = True
    #assert found

    # Compile list of HSPs that fall within the range of HSPs in the final
    # list, but were not included because they overlap with higher ranking HSPs
    # in the final list. 
    #leftover_proximate_hsps = []
    #for h in hsp_objs:
    #    if h.hit_range not in [x.hit_range for x in final_proximate_hsps]:
    #        for fh in final_proximate_hsps:
    #            if check_if_two_hsps_overlap(h, fh): 
    #                leftover_proximate_hsps.append(h)
    #                break
    #for h in hsp_objs:
    #    for fh in final_proximate_hsps:
    #        if h.hit_range != fh.hit_range or h.query_range != fh.query_range:
    #            if check_if_two_hsps_overlap(h, fh): 
    #                leftover_proximate_hsps.append(h)
    #                break

    leftover_proximate_hsps = []
    cluster_range = get_cluster_range(final_proximate_hsps)
    cluster_range_set = set(range(cluster_range[0] - max_gap, cluster_range[1] + max_gap))
    #print('\ncluster_range')
    #print(cluster_range)
    #print('hsps in complete list')
    for h in hsp_objs:
        #print('\t' + str(h.hit_range))
        # Determine whether the HSP is in the list of proximate HSPs.
        possible = True
        for fh in final_proximate_hsps:
            if h.hit_range == fh.hit_range and h.query_range == fh.query_range:
                possible = False
                #print('\t\talready in proximate hsp list for cluster')
                break
        # If not, then determine whether it is within the range of the
        # proximate HSPs.
        if possible:
            # If it is, then add it to the list of leftover HSPs.
            #if h.hit_start in cluster_range_set or h.hit_end in cluster_range_set: 
            for x in range(h.hit_start,h.hit_end):
                if x in cluster_range_set:
                    #print('\t\tadded to leftover list for cluster')
                    leftover_proximate_hsps.append(h)
                    break

    ## Redundant.
    ## Check that there are no HSPs in the list that are within the range of the
    ## cluster that are not in either the final or leftover lists.
    #cluster_range = get_cluster_range(final_proximate_hsps)
    #for h in hsp_objs:
    #    if h.hit_range not in [x.hit_range for x in (final_proximate_hsps + leftover_proximate_hsps)]:
    #        r = range(cluster_range[0], cluster_range[1])
    #        in_range = False
    #        if h.hit_range[0] in r:
    #            in_range = True
    #        elif h.hit_range[1] in r:
    #            in_range = True
    #        assert not in_range, """A HSP still falls
    #        within the range of identified proximate HSPs."""

    # Check whether any HSPs that should have been identified and added
    # to the leftover_proximate_hsps list were missed.
    complete_hsp_ranges = [x.hit_range for x in hsp_objs]
    final_proximate_hsp_ranges = [x.hit_range for x in final_proximate_hsps]
    leftover_proximate_hsp_ranges = [x.hit_range for x in leftover_proximate_hsps]
    missed = find_missed_hsp_ranges(complete_hsp_ranges, final_proximate_hsp_ranges,leftover_proximate_hsp_ranges)
    assert not missed

    # Return the final and leftover (excluded proximate) hsp lists.
    return [final_proximate_hsps, leftover_proximate_hsps] 

# Obsolete function:
def get_proximate_hsp_objs(hsp_objs, top_hsp, max_gap):
    """Return a list of HSP objects that are sufficiently proximate to the top
    HSP, and among which there is no overlap in the sequence ranges. And, also
    return a list of HSP objects that are also proximate to these HSPs, but
    were not included in the final list of HSPs because they overlap with the
    ones that were.
    """
    # Define maximum gap length (bp) between consecutive HSPs upstream and
    # downstream of the top HSP.
    max_gap = max_gap

    # Define the range for the top HSP.
    top_hsp_range = top_hsp.hit_range

    # Get a list of hit ranges
    range_list = []
    for h in hsp_objs:
        range_list.append(h.hit_range)

    # Check whether any of the ranges are identical (and can therefore be used
    # to identify different HSPs).
    has_duplicate_ranges = check_for_duplicate_ranges(range_list)

    if has_duplicate_ranges:
        # Remove HSPs with redundant ranges.
        non_redun_hsps = remove_hsps_with_redun_ranges(range_list, hsp_objs)

        # Get a list of hit ranges again.
        range_list = []
        for h in non_redun_hsps:
            range_list.append(h.hit_range)

        # Check whether any of the ranges are identical again.
        has_duplicate_ranges = check_for_duplicate_ranges(range_list)
        assert not has_duplicate_ranges, "Error: Some HSPs have identical ranges."

    # Get a list of ranges for HSPs proximate to the top HSP, and not
    # overlapping with each other. 
    phr = get_proximate_hsp_ranges(range_list,
                                   top_hsp_range,
                                   max_gap)
    proximate_hsp_ranges = phr[0]
    leftover_proximate_hsp_ranges = phr[1]

    # Generate a list of proximate HPS objects, and leftover HSP objects (those
    # that don't make it into the final list because they overylap with
    # higher-ranking hsps) from the list of their ranges.
    final_proximate_hsps = []
    for r in proximate_hsp_ranges:
        for h in hsp_objs:
            if h.hit_range == r:
                final_proximate_hsps.append(h)
    leftover_proximate_hsps = []
    for r in leftover_proximate_hsp_ranges:
        for h in hsp_objs:
            if h.hit_range == r:
                leftover_proximate_hsps.append(h)

    # Check whether the HSPs overlap. Unnecessary?
    overlap = check_if_hsps_overlapping(hsp_objs)

    #num = 0
    #for h in final_proximate_hsps:
    #    num += 1
    #    #print('>' + str(h.hit_range[0]) + ' ' + str(h.hit_range[1]) + '\n' + str(h.hit.seq).replace('-', '') + '\n')

    return [final_proximate_hsps, leftover_proximate_hsps]


def get_proximate_hsp_objs2(hsp_objs, top_hsp, max_gap):
    """Return a list of HSP objects that are sufficiently proximate to the top
    HSP, and among which there is no overlap in the sequence ranges. And, also
    return a list of HSP objects that are also proximate to these HSPs, but
    were not included in the final list of HSPs because they overlap with the
    ones that were. Instead of simply selecting HSPs based on hit range,
    consider E-values, favoring HSPs with better E-values if HSPs overlap and
    some have to be excluded.
    """
    # Define the range for the top HSP.
    top_hsp_range = top_hsp.hit_range

    # Get a list of hit ranges
    range_list = []
    for h in hsp_objs:
        range_list.append(h.hit_range)

    # Check whether any of the ranges are identical (and can therefore be used
    # to identify different HSPs).
    has_duplicate_ranges = check_for_duplicate_ranges(range_list)

    if has_duplicate_ranges:
        # Remove HSPs with redundant ranges.
        non_redun_hsps = remove_hsps_with_redun_ranges(range_list, hsp_objs)

        # Get a list of hit ranges again.
        range_list = []
        for h in non_redun_hsps:
            range_list.append(h.hit_range)

        # Check whether any of the ranges are identical again.
        has_duplicate_ranges = check_for_duplicate_ranges(range_list)
        assert not has_duplicate_ranges, "Error: Some HSPs have identical ranges."

    ## Get a list of ranges for HSPs proximate to the top HSP, and not
    ## overlapping with each other. 
    #phr = get_proximate_hsp_ranges(range_list,
    #                               top_hsp_range,
    #                               max_gap)
    #proximate_hsp_ranges = phr[0]
    #leftover_proximate_hsp_ranges = phr[1]

    ## Generate a list of proximate HPS objects, and leftover HSP objects (those
    ## that don't make it into the final list because they overylap with
    ## higher-ranking hsps) from the list of their ranges.
    #final_proximate_hsps = []
    #for r in proximate_hsp_ranges:
    #    for h in hsp_objs:
    #        if h.hit_range == r:
    #            final_proximate_hsps.append(h)
    #leftover_proximate_hsps = []
    #for r in leftover_proximate_hsp_ranges:
    #    for h in hsp_objs:
    #        if h.hit_range == r:
    #            leftover_proximate_hsps.append(h)

    ## Get a list of HSPs proximate to the top HSP, and not
    ## overlapping with each other. 
    #bph = get_best_proximate_hsps(hsp_objs,
    #                              top_hsp,
    #                              max_gap)
    #final_proximate_hsps = bph[0]
    #leftover_proximate_hsps = bph[1]

    ## Check whether the HSPs overlap. Unnecessary?
    #overlap = check_if_hsps_overlapping(hsp_objs)

    #return [final_proximate_hsps, leftover_proximate_hsps]

    # Compile a list of proximate HSPs upstream (left regardless of which
    # way 5' or 3' is) of the the top HSP.
    left_hsps = [x for x in hsp_objs if x.hit_range[0] < top_hsp.hit_range[1]]
    # Remove top HSP from list.
    left_hsps.remove(top_hsp)

    reduced_left_hsps = []

    if len(left_hsps) > 0:

        # Sort the hsps in descending order of their 3' position.
        sorted_left_hsps = sorted(left_hsps, key=lambda x: x.hit_range[1], reverse=True)

        # Reduce the list of left ranges to only those that meet the criteria.
        reduced_left_hsps = reduce_to_best_hsps([top_hsp] + sorted_left_hsps,
                max_gap, 0, 'Left', 'left')
        assert len(reduced_left_hsps) <= len([top_hsp] + sorted_left_hsps)

    # Compile a list of proximate HSPs downstream (right regardless of which
    # way 5' or 3' is) of top HSP.
    right_hsps = [x for x in hsp_objs if x.hit_range[1] > top_hsp.hit_range[0]]
    # Remove top HSP from list.
    right_hsps.remove(top_hsp)

    reduced_right_hsps = []

    if len(right_hsps) > 0:

        # Sort the HSPs in ascending order of their 5' position.
        sorted_right_hsps = sorted(right_hsps, key=lambda x: x.hit_range[0], reverse=False)

        # Reduce the list of right HSPs to only those that meet the criteria.
        reduced_right_hsps = reduce_to_best_hsps([top_hsp] + sorted_right_hsps,
                max_gap, 0, 'Right', 'right')
        assert len(reduced_right_hsps) <= len([top_hsp] + sorted_right_hsps)

    # Check that all HSPs in the input list are accounted for in the left
    # and/or right lists.
    assert len(hsp_objs) <= (len(left_hsps) + 1 + len(right_hsps)), """Did not
    account for all HSPs in input list."""

    # Concatenate left, top, and right HSP lists to make the final list.
    final_proximate_hsps = list(reversed(reduced_left_hsps[1:]))\
                             + [top_hsp]\
                             + reduced_right_hsps[1:]

    ## Unnecessary???
    ## Check that top hsp was included.
    #found = False
    #for hsp in final_proximate_hsps:
    #    if hsp.hit_range == top_hsp.hit_range:
    #        found = True
    #assert found

    # Compile list of HSPs that fall within the range of HSPs in the final
    # list, but were not included because they overlap with higher ranking HSPs
    # in the final list. 
    #leftover_proximate_hsps = []
    #for h in hsp_objs:
    #    if h.hit_range not in [x.hit_range for x in final_proximate_hsps]:
    #        for fh in final_proximate_hsps:
    #            if check_if_two_hsps_overlap(h, fh): 
    #                leftover_proximate_hsps.append(h)
    #                break
    #for h in hsp_objs:
    #    for fh in final_proximate_hsps:
    #        if h.hit_range != fh.hit_range or h.query_range != fh.query_range:
    #            if check_if_two_hsps_overlap(h, fh): 
    #                leftover_proximate_hsps.append(h)
    #                break

    leftover_proximate_hsps = []
    cluster_range = get_cluster_range(final_proximate_hsps)
    cluster_range_set = set(range(cluster_range[0] - max_gap, cluster_range[1] + max_gap))
    #print('\ncluster_range')
    #print(cluster_range)
    #print('hsps in complete list')
    for h in hsp_objs:
        #print('\t' + str(h.hit_range))
        # Determine whether the HSP is in the list of proximate HSPs.
        possible = True
        for fh in final_proximate_hsps:
            if h.hit_range == fh.hit_range and h.query_range == fh.query_range:
                possible = False
                #print('\t\talready in proximate hsp list for cluster')
                break
        # If not, then determine whether it is within the range of the
        # proximate HSPs.
        if possible:
            # If it is, then add it to the list of leftover HSPs.
            #if h.hit_start in cluster_range_set or h.hit_end in cluster_range_set: 
            for x in range(h.hit_start,h.hit_end):
                if x in cluster_range_set:
                    #print('\t\tadded to leftover list for cluster')
                    leftover_proximate_hsps.append(h)
                    break

    ## Redundant.
    ## Check that there are no HSPs in the list that are within the range of the
    ## cluster that are not in either the final or leftover lists.
    #cluster_range = get_cluster_range(final_proximate_hsps)
    #for h in hsp_objs:
    #    if h.hit_range not in [x.hit_range for x in (final_proximate_hsps + leftover_proximate_hsps)]:
    #        r = range(cluster_range[0], cluster_range[1])
    #        in_range = False
    #        if h.hit_range[0] in r:
    #            in_range = True
    #        elif h.hit_range[1] in r:
    #            in_range = True
    #        assert not in_range, """A HSP still falls
    #        within the range of identified proximate HSPs."""

    # Check whether any HSPs that should have been identified and added
    # to the leftover_proximate_hsps list were missed.
    complete_hsp_ranges = [x.hit_range for x in hsp_objs]
    final_proximate_hsp_ranges = [x.hit_range for x in final_proximate_hsps]
    leftover_proximate_hsp_ranges = [x.hit_range for x in leftover_proximate_hsps]
    missed = find_missed_hsp_ranges(complete_hsp_ranges, final_proximate_hsp_ranges,leftover_proximate_hsp_ranges)
    assert not missed

    # Return the final and leftover (excluded proximate) hsp lists.
    return [final_proximate_hsps, leftover_proximate_hsps] 


def get_hit_seq_record_and_coord(hit, proximate_hsp_objs):
    """Given a list of SearchIO HSP objects (arguably representing a single
    gene), return a Bio.Seq record composed of concatenated HSP sequences.
    """
    # Get strand that the HSPs are on.
    strand = get_strand_word(proximate_hsp_objs[0].hit_frame)

    # Concatenate the (translated peptide) hsp sequences in order from 5' to 3'
    # as they appear in the nucleotide subject sequence.
    hit_sequence = None
    if strand == 'plus': 
        hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))
    elif strand == 'minus': 
        # Reverse the order of the HSPs.
        hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs][::-1]))

    # Make a description of where the hit sequence came from.
    seq_origin = ','.join(['[' + str(x.hit_range[0]) + '..' + str(x.hit_range[1])\
        + ']' for x in proximate_hsp_objs])

    # Define sequence coordinates as a list.
    seq_coord = [[x.hit_range[0],x.hit_range[1]]\
        for x in proximate_hsp_objs]

    # Generate a fasta sequence record.
    hit_seq_record = SeqRecord(hit_sequence)
    hit_seq_record.id = hit.id
    descr_without_quotation_marks = hit.description.replace('\"', '')
    hit_seq_record.description = '\"' + descr_without_quotation_marks + ' ' +\
    strand + ' strand ' + seq_origin + '\"'

    return [hit_seq_record, seq_coord]


def get_hit_seq_record_and_coord2(search_result_path,
                                  hit,
                                  proximate_hsp_objs,
                                  db_filename,
                                  query_filename,
                                  genetic_code_number,
                                  exonerate_score_threshold
                                  ):
    """Given a list of SearchIO HSP objects (arguably representing a single
    gene), retrieve the corresponding region of the subject sequence, use
    exonerate to predict exon locations and extract translated sequences,
    return a list with the first element being a Bio.Seq record from the
    exonerate analysis and the second element being a string representing the
    position of predicted exons.
    """
    ## Get strand that the HSPs are on.
    #strand = get_strand_word(proximate_hsp_objs[0].hit_frame)

    ## Concatenate the (translated peptide) hsp sequences in order from 5' to 3'
    ## as they appear in the nucleotide subject sequence.
    #hit_sequence = None
    #if strand == 'plus': 
    #    hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))
    #elif strand == 'minus': 
    #    # Reverse the order of the HSPs.
    #    hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs][::-1]))

    ## Make a description of where the hit sequence came from.
    #seq_origin = ','.join(['[' + str(x.hit_range[0]) + '..' + str(x.hit_range[1])\
    #    + ']' for x in proximate_hsp_objs])

    ## Define sequence coordinates as a list.
    #seq_coord = [[x.hit_range[0],x.hit_range[1]]\
    #    for x in proximate_hsp_objs]

    ## Generate a fasta sequence record.
    #hit_seq_record = SeqRecord(hit_sequence)
    #hit_seq_record.id = hit.id
    #descr_without_quotation_marks = hit.description.replace('\"', '')
    #hit_seq_record.description = '\"' + descr_without_quotation_marks + ' ' +\
    #strand + ' strand ' + seq_origin + '\"'


    # Get full path to query FASTA file.
    query_dir = settings.querydirpath
    query_faa = os.path.join(query_dir, query_filename)
    assert os.path.isfile(query_faa), """Specified query file path is
    not a file: %s""" % query_faa

    # Get full path to database FASTA file.
    db_dir = settings.dbdirpath
    target_fna_path = os.path.join(db_dir, db_filename)
    assert os.path.isfile(target_fna_path), """Specified database file path is
    not a file: %s""" % target_fna_path

    # Retrieve subject sequence of interest.

    # Determine ID of subject sequence of interest from the Biopython hit
    # object.
    target_seq_id = hit.id

    # Get start and end of the subsequence of interest within the subject
    # sequence based on the range that the HSPs are found in.
    start_of_5prime_hsp = min([x.hit_range[0] for x in proximate_hsp_objs])
    end_of_3prime_hsp = max([x.hit_range[1] for x in proximate_hsp_objs])

    # Define the number of additional basepairs to include on either side.
    additional_flanking_basepairs = 500

    # Define path to FASTA file with subsequence of interest from target
    # nucleotide sequence.
    subseq_fasta_path = search_result_path.rsplit('.', 1)[0] +\
    '_subject_subseq_' + target_seq_id + '_' + str(start_of_5prime_hsp) + '-' +\
    str(end_of_3prime_hsp) + '.fna'

    # Extract relevant subsequence from input target sequence (region identified in
    # a previous step using TBLASTN).
    target_subseq_start = get_subseq_from_nucl(target_seq_id, 
                                               target_fna_path,
                                               start_of_5prime_hsp,
                                               end_of_3prime_hsp,
                                               subseq_fasta_path,
                                               additional_flanking_basepairs
                                               )

    # Define path to exonerate output file.
    exonerate_output_filepath = subseq_fasta_path.rsplit('.', 1)[0] + '_exonerate_out.txt'
    
    # Run exonerate as a subprocess.
    genetic_code_number = '1'
    run_exonerate_as_subprocess(query_faa,
                                subseq_fasta_path,
                                exonerate_output_filepath,
                                exonerate_score_threshold,
                                genetic_code_number
                                )
    
    # Parse output of exonerate.
    #parse_exonerate_output(exonerate_output_filepath)
    rough_start = start_of_5prime_hsp
    rough_end = end_of_3prime_hsp
    position_of_subject_seq_start_in_original = int(target_subseq_start)
    exonerate_locus_result_obj = ExonerateLocusResult(exonerate_output_filepath,
                                                      subseq_fasta_path,
                                                      position_of_subject_seq_start_in_original,
                                                      genetic_code_number,
                                                      rough_start,
                                                      rough_end
                                                      )
    
    # Delete temporary intermediate files.
    os.remove(subseq_fasta_path)

    if exonerate_locus_result_obj.seq_record is None:
        # Return None to indicate that exonerate did not identify anything
        # (even though TBLASTN did).
        return None
    else:
        # Write output to file.
        exonerate_output_seq_filepath = subseq_fasta_path.rsplit('.', 1)[0] +\
        '_exonerate_out_seq.faa'
        with open(exonerate_output_seq_filepath, 'w') as o:
            SeqIO.write(exonerate_locus_result_obj.seq_record, o, 'fasta')

        # Return the sequence record and location/coordinate string.
        return [exonerate_locus_result_obj.seq_record, exonerate_locus_result_obj.location_string]


def get_hit_seq_obj(hit, max_gap):
    """Takes a SearchIO hit object (containing one or more HSPs) and returns a
    fasta sequence object likely to correspond to a single gene (the most
    similar gene to the query in the nucleotide subject sequence).
    """
    # Make a list of HSP objects within the given hit object.
    hsp_objs = []
    for hsp in hit:
        hsp_objs.append(hsp)

    # Determine whether the top HSP is on the forward or reverse strand.
    top_hsp = hsp_objs[0]
    fwd_strand = whether_fwd_strand(top_hsp) 

    # Check that all the hsps are on the same strand.
    hsps_on_different_strands = check_strands(hsp_objs, fwd_strand)

    # If HSPs are on different strands, then make a reduced list of HSP objects
    # composed of only those that are on the same strand as the top HSP.
    hsp_objs_red = hsp_objs 
    if hsps_on_different_strands:
        hsp_objs_red = []
        for hsp in hsp_objs:
            x = whether_fwd_strand(hsp)
            # Only include in reduced list if on same strand as top HSP.
            if x == fwd_strand:
                hsp_objs_red.append(hsp)

    # (If more than one hsp)

    # Only include HSPs that are sufficiently proximate to the top HSP, and not
    # overlapping.
    proximate_hsp_objs = get_proximate_hsp_objs2(hsp_objs_red, top_hsp,
            max_gap)[0]

    # Concatenate the hsp sequences in order from 5' to 3' as they appear
    # in the nucleotide subject sequence.
    # Make a description of where the hit sequence came from.
    # Generate a fasta sequence record.
    hit_seq_record = get_hit_seq_record_and_coord(hit, proximate_hsp_objs)[0]

    return hit_seq_record


def get_tblastn_hit_seq_obj_and_coord(hit, max_gap):
    """Takes a SearchIO hit object (containing one or more HSPs) and returns a
    fasta sequence object likely to correspond to a single gene (the most
    similar gene to the query in the nucleotide subject sequence).
    """
    # Make a list of HSP objects within the given hit object.
    hsp_objs = []
    for hsp in hit:
        hsp_objs.append(hsp)

    # Determine whether the top HSP is on the forward or reverse strand.
    top_hsp = hsp_objs[0]
    fwd_strand = whether_fwd_strand(top_hsp) 

    # Check that all the hsps are on the same strand.
    hsps_on_different_strands = check_strands(hsp_objs, fwd_strand)

    # If HSPs are on different strands, then make a reduced list of HSP objects
    # composed of only those that are on the same strand as the top HSP.
    hsp_objs_red = hsp_objs 
    if hsps_on_different_strands:
        hsp_objs_red = []
        for hsp in hsp_objs:
            x = whether_fwd_strand(hsp)
            # Only include in reduced list if on same strand as top HSP.
            if x == fwd_strand:
                hsp_objs_red.append(hsp)

    # (If more than one hsp)

    # Only include HSPs that are sufficiently proximate to the top HSP, and not
    # overlapping.
    proximate_hsp_objs = get_proximate_hsp_objs2(hsp_objs_red, top_hsp,
            max_gap)[0]

    # Concatenate the hsp sequences in order from 5' to 3' as they appear
    # in the nucleotide subject sequence.
    #hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))
    hit_sequence = None
    if strand == 'plus': 
        hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))
    elif strand == 'minus': 
        # Reverse the order of the HSPs.
        hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs][::-1]))


    # Make a description of where the hit sequence came from.
    #seq_origin = ','.join(['[' + str(x.hit_range[0]) + '..' + str(x.hit_range[1])\
    #    + ']' for x in proximate_hsp_objs])

    # Define sequence coordinates as a list.
    seq_coord = [[x.hit_range[0],x.hit_range[1]]\
        for x in proximate_hsp_objs]

    # Generate a fasta sequence record.
    hit_seq_record = SeqRecord(hit_sequence)
    hit_seq_record.id = hit.id
    descr_without_quotation_marks = hit.description.replace('\"', '')
    hit_seq_record.description = '\"' + descr_without_quotation_marks + ' ' +\
    str(seq_coord) + '\"'

    # Return the sequence record and sequence coordinates.
    return [hit_seq_record, seq_coord]


def recursively_add_clusters(hsp_list, hsp_clusters, max_gap):
    """take a list of SearchIO HSP objects and a list of HSP clusters
    (representing distinct genes), and add clusters recursively. 
    """
    # Check that there are hsps in the hsp list.
    assert len(hsp_list) > 0, """The input list is empty."""

    # Visualize input hsps.
    #print_cluster(clusterplus, hit_num, cluster_num, num_dots, startend=None)
    #print_cluster([hsp_list,[]], 'Hit X', 'recursively_add_clusters recursion x', 150, startend=None)

    # Check that the clusters are listed in order of ascending E-value.
    hsp_list = sorted(hsp_list, key=lambda x: x.evalue)
    #print('\t\thsp_list: ' + str([x.hit_range for x in hsp_list]))

    # Define top HSP.
    top_hsp = hsp_list[0]
    #print('\t\ttop_hsp: ' + str(top_hsp.hit_range))

    # Only include HSPs that are sufficiently proximate to the top HSP, and not
    # overlapping.
    #pho = get_proximate_hsp_objs(hsp_list, top_hsp, max_gap)
    pho = get_proximate_hsp_objs2(hsp_list, top_hsp, max_gap)
    proximate_hsp_objs = pho[0]
    assert len(proximate_hsp_objs) > 0, """Did not properly identify HSPs
    within hit for clustering at locus."""
    leftover_proximate_hsp_objs = pho[1]
    #print('\t\tproximate_hsp_objs: ' + str([x.hit_range for x in proximate_hsp_objs]))

    # Visualize proximate HSPs.
    #print_cluster([proximate_hsp_objs, leftover_proximate_hsp_objs], 'Hit X', 'recursively_add_clusters recursion x proximate hsps', 150, startend=None)

    #print('XXXXXX')
    #print(top_hsp)
    #print(sorted(proximate_hsp_objs, key=lambda x: x.evalue)[0])

    ## Check that top hsp included in proximate list. Unecessary?
    #found = False
    #for hsp in proximate_hsp_objs:
    #    if hsp.hit_range == top_hsp.hit_range:
    #        found = True
    #assert found

    # Add proximate HPS cluster to list.
    hsp_clusters.append([proximate_hsp_objs, leftover_proximate_hsp_objs])

    # Note: need to clean up redundancy in code below.

    # Remove proximate hsp objects (and leftover proximate hsp objects) from
    # the initial list.
    new_hsp_list = copy.deepcopy(hsp_list)
    new_new_hsp_list = []
    hsps_to_remove = []
    if len(proximate_hsp_objs) > 0:
        # Compile a list of HSPs to remove.
        for hsp1 in new_hsp_list:
            for hsp2 in proximate_hsp_objs + leftover_proximate_hsp_objs:
                if str(hsp1.hit) == str(hsp2.hit):
                    if hsp1.hit_range == hsp2.hit_range:
                        # Note: Different HSPs may have identical hit_range, but
                        # different query_range if the query contains repeats.
                        if hsp1.query_range == hsp2.query_range:
                            hsps_to_remove.append(hsp1)
                            break
        # Remove from list.
        #for hsp in hsps_to_remove:
        #    new_hsp_list.remove(hsp)
        for hsp1 in new_hsp_list:
            append = True
            for hsp2 in hsps_to_remove:
                if str(hsp1.hit) == str(hsp2.hit):
                    if hsp1.hit_range == hsp2.hit_range:
                        if hsp1.query_range == hsp2.query_range:
                            append = False 
                            break
            if append:
                new_new_hsp_list.append(hsp1)

        # Check whether the removal worked.
        assert len(hsp_list) > len(new_new_hsp_list), """Did not remove any HSPs
        from the list."""
        assert len(hsps_to_remove) == len(set(hsps_to_remove)), """Redundant
        HSPs were added to the hsps_to_remove list."""
        # Display.
        if not len(set(hsps_to_remove)) == (len(set(new_hsp_list)) - len(set(new_new_hsp_list))):
            print('Number of HSPs to remove: ' + str(to_remove_count))
            print('Number of HSPs removed: ' + str(len(new_hsp_list) -\
                len(new_new_hsp_list)))
            print('HSPs not removed that should have been:')
            hsps_that_should_have_been_removed_but_were_not =\
                list(set(hsps_to_remove)\
                - (set(new_hsp_list) - set(new_new_hsp_list)))

            print(len(hsps_that_should_have_been_removed_but_were_not))

            for hsp in hsps_that_should_have_been_removed_but_were_not:
                print(hsp.hit_range)

        assert len(set(hsps_to_remove)) == (len(set(new_hsp_list)) - len(set(new_new_hsp_list)))

    #print('new_hsp_list: ' + str(new_hsp_list))

    # Check that there are no redundant HSPs in the clusters.
    for clusterplus in hsp_clusters:
        cluster = clusterplus[0]
        assert len(cluster) == len(set(cluster)), """Redundant HSPs in list for
        cluster."""

    # Check that there are no overlaps in the sets of HSPs in the clusters
    # (redundancy).
    for cluster1 in hsp_clusters:
        for cluster2 in hsp_clusters:
            if cluster1[0] != cluster2[0]:
                for hsp1 in cluster1[0]:
                    for hsp2 in cluster2[0]:
                        assert hsp1.hit_range != hsp2.hit_range,\
                        """Different clusters of HSPs contain overlapping sets of HSPs.""" 

    ## Redundant:
    ## Check that there are no HSPs in the reduced list that fall within the
    ## range of the cluster just identified.
    #current_cluster_range1 = get_cluster_range(proximate_hsp_objs)
    #for hsp1 in new_new_hsp_list:
    #    r1 = set(range(current_cluster_range1[0], current_cluster_range1[1]))
    #    if not set([hsp1.hit_start, hsp1.hit_end]).intersection(r1) == set():
    #        print('problem leftover HSP:')
    #        print(hsp1)
    #        print('HSPs in cluster:')
    #        for hsp in proximate_hsp_objs:
    #            print(hsp)
    #    assert set([hsp1.hit_start, hsp1.hit_end]).intersection(r1) == set(),\
    #    """HSP %s falls within the range of the current identified cluster %s."""\
    #    % (str(hsp1.hit_range), str(current_cluster_range1))

    #if len(leftover_proximate_hsp_objs) > 0:
    #    current_cluster_range2 = get_cluster_range(leftover_proximate_hsp_objs)
    #    for hsp1 in new_new_hsp_list:
    #        r2 = set(range(current_cluster_range2[0], current_cluster_range2[1]))
    #        assert set([hsp1.hit_start, hsp1.hit_end]).intersection(r2) == set(),\
    #        """HSP %s falls within the range of the HSPs excluded from the current
    #        identified cluster due to overlap."""\
    #        % str(hsp1.hit_range)

    # Check that there are no gaps between HSPs larger than max_gap.
    last_hsp = None 
    for hsp in proximate_hsp_objs:
        if last_hsp is not None:
            assert last_hsp.hit_end - hsp.hit_start <= max_gap, """A gap between
            adjacent HSPs is larger than max_gap."""

    ## Redundant:
    ## Check that top hsp is not included in remaining list.
    #found = False
    #for hsp in new_new_hsp_list:
    #    if hsp.hit_range == top_hsp.hit_range:
    #        found = True
    #assert not found

    # Redundant check?:
    # Check that none of the HSPs in the proximate HSP list/cluster have
    # overlapping ranges with each other.
    for hsp1 in proximate_hsp_objs:
        for hsp2 in proximate_hsp_objs:
            if not hsp1.hit_range == hsp2.hit_range:
                overlap = False
                if check_if_two_hsps_overlap(hsp1, hsp2):
                    overlap = True
                if overlap:
                    print(hsp1)
                    print(hsp2)
                assert not overlap, """HSPs in proximate HSP list have
                overlapping ranges."""

    # Check that the order of HSPs in the proximate HSPs selected makes sense
    # given their query ranges.
    strand = get_strand_word(proximate_hsp_objs[0].hit_frame)
    l = None
    if strand == 'plus':
        l = proximate_hsp_objs
    else:
        l = reversed(proximate_hsp_objs)
    last_hsp_query_start = -1
    last_hsp_query_end = -1
    #print('\nhsps in proximate_hsp_objs:')
    for hsp in l:
        #print(str(hsp.query_start) + ' ' + str(hsp.hit_start))
        assert hsp.query_start > last_hsp_query_start\
               or hsp.query_end > last_hsp_query_end,\
               """Query ranges of the proximate HSPS are not in the appropriate
               order."""
        last_hsp_query_start = hsp.query_start
        last_hsp_query_end = hsp.query_end

    # Stop rule.
    if len(new_new_hsp_list) == 0:
        return hsp_clusters
    # Recursion.
    else:
        return recursively_add_clusters(new_new_hsp_list, hsp_clusters, max_gap)


def get_hsp_clusters(hit, max_gap):
    """Take a SearchIO Hit object and return a list of lists, where the
    sublists are HSPs that cluster together on the subject sequence such that
    it is likely that they represent a distinct gene.
    """
    # Make a list of HSP objects within the given hit object.
    hsp_objs = []
    for hsp in hit:
        hsp_objs.append(hsp)

    # Split HSPs into those on the forward strand and those on the reverse
    # strand.
    hsp_objs_fwd = []
    hsp_objs_rev = []
    for hsp in hsp_objs:
        on_fwd_strand = whether_fwd_strand(hsp)
        if on_fwd_strand:
            hsp_objs_fwd.append(hsp)
        else:
            hsp_objs_rev.append(hsp)

    # Define list of HSP lists
    list_of_hsp_lists = []
    templist = [hsp_objs_fwd, hsp_objs_rev]
    for i in templist:
        if len(i) > 0:
            list_of_hsp_lists.append(i)
    assert len(list_of_hsp_lists) > 0

    #for i in list_of_hsp_lists:
    #    for j in i:
    #        print(j)

    ## Order by lowest E-value of best hit.
    #if len(list_of_hsp_lists) > 1:
    #    print('\nlist_of_hsp_lists: ' + str(list_of_hsp_lists))
    #    list_of_hsp_lists = sorted(list_of_hsp_lists, key=lambda x: min([y.evalue for y in x]))

    # Add to list of clusters. 
    hsp_clusters = []
    for hsp_list in list_of_hsp_lists:
        if len(hsp_list) > 0:
            #hsp_clusters = recursively_add_clusters(hsp_list, hsp_clusters, max_gap)
            #print('\n\nRunning recursively_add_clusters on a strand of ' + hit.id)
            strand_hsp_clusters = recursively_add_clusters(hsp_list, [], max_gap)

            ## Redundant check, because done in recursively_add_clusters
            ## function:
            ## Check that the clusters do not overlap with each other on this
            ## strand of the subject sequence.
            #for cluster1 in strand_hsp_clusters:
            #    for cluster2 in strand_hsp_clusters:
            #        if cluster1[0] != cluster2[0]:
            #            assert not clusters_overlap(cluster1[0], cluster2[0]),\
            #            """Clusters overlap: %s, and %s""" %\
            #            (cluster1[0][0].hit_id + str(get_cluster_range(cluster1[0])) +\
            #                    ' on the ' + get_strand_word(cluster1[0][0].hit_frame) + ' strand'\
            #             , cluster2[0][0].hit_id + str(get_cluster_range(cluster2[0])) +\
            #             ' on the ' + get_strand_word(cluster2[0][0].hit_frame) + ' strand')

            # Add clusters to complete list. 
            hsp_clusters = hsp_clusters + strand_hsp_clusters

    # Order clusters by ascending E-value.
    #print('\n\nPrinting E-values for tblastn hsp clusters:')
    #for h in hsp_clusters:
    #    print(min([y.evalue for y in h[0]]))
    hsp_clusters = sorted(hsp_clusters, key=lambda x: min([y.evalue for y in x[0]]))
    #print('\n\nPrinting E-values for tblastn hsp clusters after ranking:')
    #for h in hsp_clusters:
    #    print(min([y.evalue for y in h[0]]))
    #print('\n\n')

    ## Concatenate the hsp sequences in order from 5' to 3' as they appear
    ## in the nucleotide subject sequence.
    #hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))

    ## Make a description of where the hit sequence came from.
    ##seq_origin = ','.join(['[' + str(x.hit_range[0]) + '..' + str(x.hit_range[1])\
    ##    + ']' for x in proximate_hsp_objs])

    ## Define sequence coordinates as a list.
    #seq_coord = [[x.hit_range[0],x.hit_range[1]]\
    #    for x in proximate_hsp_objs]

    ## Generate a fasta sequence record.
    #hit_seq_record = SeqRecord(hit_sequence)
    #hit_seq_record.id = hit.id
    #descr_without_quotation_marks = hit.description.replace('\"', '')
    #hit_seq_record.description = '\"' + descr_without_quotation_marks + ' ' +\
    #str(seq_coord) + '\"'

    # Return list of HSPs.
    return hsp_clusters


def get_strand_word(hsp_hit_frame):
    """Take a frame (for a translation of a DNA sequence) and return 'plus' or
    'minus' depending on which strand the frame corresponds to.
    """
    if hsp_hit_frame > 0:
        return 'plus'
    else:
        return 'minus'


def get_cluster_range(hsp_list):
    """Take a list of SearchIO HSP objects an return a range from the furthest
    5' hsp end to the furthest 3' hsp end.
    """
    start = min([x.hit_range[0] for x in hsp_list])
    end = max([x.hit_range[1] for x in hsp_list])
    assert start < end
    return [start, end]


def clusters_overlap(cluster1, cluster2):
    """Take two lists of HSPs, and determine if they define overlapping ranges
    in the subject sequence.
    """
    # Get cluster ranges.
    cluster1_range = get_cluster_range(cluster1)
    cluster2_range = get_cluster_range(cluster2)

    overlap = check_if_two_hsp_ranges_overlap([cluster1_range, cluster2_range])
    ## Convert range lists to sets, and check for intersection of those sets.
    #x = set(range(cluster1_range[0], cluster1_range[1]))
    #y = set(range(cluster2_range[0], cluster2_range[1]))
    #intersect = x.intersection(y)

    #overlap = False
    #if intersect != set():
    #    overlap = True

    return overlap


def print_cluster(clusterplus, hit_num, cluster_num, num_dots, startend=None):
    """Print visualization of cluster.

    startend is a tuple containing the start and end position of the
    subsequence of the subject sequence to be visualized.
    """
    cluster = clusterplus[0]

    # Get cluster range.
    start = None
    end = None
    if startend is None: 
        cr = get_cluster_range(clusterplus[0] + clusterplus[1])
        start = cr[0]
        end = cr[1]
    else:
        start = startend[0]
        end = startend[1]

    subseqlen = end - start
    print('\tHit '+ str(hit_num) + ' HSP cluster ' + str(cluster_num) + ':')
    #print('\t\tcluster range: ' + str(cr))
    #print('\t\tsubseqlen: ' + str(subseqlen))
    #print('\t\t\t' + hit.id)
    print('\t' + 'HSP positions in subject sequence (1 dot = ' +\
    str(int(subseqlen / num_dots)) + ' bp):')
    print('\t' + str(start) + ' ' * (num_dots -1 - len(str(start))) + str(end))
    print('\t' + 'v' + ' ' * (num_dots -2) + 'v')
    print('\t' + '.' * num_dots)

    # Make a list of hsps.
    hsps = []
    for hsp in cluster:
        hsps.append(hsp)

    # Sort the HSPs.
    hsps2 = sorted(hsps, key=lambda x: x.hit_start)

    # Display the HSPs.
    for hsp in hsps2:
        prepend_dots = '.' * int(((hsp.hit_start - start)*num_dots)/(subseqlen))
        hsp_string = '#' * int(((hsp.hit_end - hsp.hit_start) * num_dots)/subseqlen)
        append_dots =  '.' * max([0, num_dots - len(prepend_dots) - len(hsp_string)])
        span_string = str(hsp.hit_start) + '..' + str(hsp.hit_end)
        strand = get_strand_word(hsp.hit_frame)
        info_string = ' ' + span_string + ', ' + strand + ', ' + str(hsp.evalue) 
        string = '\t' + prepend_dots + hsp_string + append_dots +\
        ' ' + info_string
        print(string)
    #for hsp in hsps2:
    #    string = '\t'
    #    prepend_dots = '.' * int(((hsp.hit_start - start)*num_dots)/(subseqlen))
    #    string = string + prepend_dots 
    #    span_string = str(hsp.hit_start) + ', ' + str(hsp.hit_end)
    #    string = string + span_string
    #    string = string + '.' * max([0, num_dots - len(prepend_dots) - len(span_string)])
    #    print(string)

    # Get list of hsps that are proximate to this cluster, but were
    # exluded due to overlap.
    unsorted_leftover_hsps = []
    for hsp in clusterplus[1]:
        unsorted_leftover_hsps.append(hsp)
    leftover_hsps = sorted(unsorted_leftover_hsps, key=lambda x: x.hit_start)
    
    # Display leftover HSPs that were excluded due to overlap.
    if len(leftover_hsps) > 0:
        t = 'HSPs excluded due to overlap with higher-scoring HSPs:'
        print('\t' + t + '=' * (num_dots - len(t)))
        for hsp in leftover_hsps:
            prepend_dots = '.' * int(((hsp.hit_start - start)*num_dots)/(subseqlen))
            hsp_string = '#' * int(((hsp.hit_end - hsp.hit_start) * num_dots)/subseqlen)
            append_dots =  '.' * max([0, num_dots - len(prepend_dots) - len(hsp_string)])
            span_string = str(hsp.hit_start) + '..' + str(hsp.hit_end)
            strand = get_strand_word(hsp.hit_frame)
            info_string = ' ' + span_string + ', ' + strand + ', ' + str(hsp.evalue) 
            string = '\t' + prepend_dots + hsp_string + append_dots +\
            ' ' + info_string
            print(string)

    print('\n')


def split_tblastn_hits_into_separate_genes(query_res_obj, max_gap):
    """Take a SearchIO QueryResult object and return a new object with hits
    split into groups of HSPs that represent distinct genes. This is important,
    because there may be multiple paralogous genes present in a single
    nucleotide subject sequence (such as a chromosome or scaffold).
    """
    # Print message.
    print('\n\tSearch program was tblastn.\n\tChecking number of distinct genes represented by HSPs.\n')

    # Copy the query result object.
    #query_res_obj2 = copy.deepcopy(query_res_obj)

    # Compile a list of all HSP clusters.
    # Display a simple visualization of HSP location.
    # List hits and HSPs in original object.
    num_dots = 150 
    all_hsp_clusters = []
    hit_num = 0
    for hit in query_res_obj:
        hit_num += 1
        print('\tQuery: ' + hit.query_id)
        print('\tHit '+ str(hit_num) + ': ' + hit.id + ' ' + hit.description)
        print('\t' + 'HSP positions in subject sequence (1 dot = ' +\
        str(int(hit.seq_len / num_dots)) + ' bp):')
        print('\t ' + '0' + ' ' * (num_dots -2) + str(hit.seq_len))
        print('\t ' + 'v' + ' ' * (num_dots -2) + 'v')
        print('\t ' + '.' * num_dots + ' ' + 'Query range:')

        # Make a list of hsps.
        hsps = []
        for hsp in hit:
            hsps.append(hsp)

        # Sort the HSPs.
        hsps2 = sorted(hsps, key=lambda x: x.hit_start)

        # Display the HSPs.
        for hsp in hsps2:
            string = '\t'
            sign = None
            if hsp.hit_frame > 0:
                sign = '+'
            elif hsp.hit_frame < 0:
                sign = '-'
            prepend_dots = '.' * int((hsp.hit_start*num_dots)/(hit.seq_len))
            string = string + sign + prepend_dots 
            span_string = str(hsp.hit_start) + ', ' + str(hsp.hit_end)
            string = string + span_string
            string = string + '.' * max([0, num_dots - len(prepend_dots) - len(span_string)])
            string = string + ' ' + str(hsp.query_range) #+ ' ' + str(hsp.evalue)
            print(string)
            #print(hsp.hit.seq)

        print('\n')

        # Generate an expanded list of hit objects.

        # Recursively find clusters of HSPs that likely represent different
        # genes, and return as a list of lists.
        hsp_clusters = get_hsp_clusters(hit, max_gap)
        all_hsp_clusters = all_hsp_clusters + hsp_clusters

        # Display HSPs in each cluster.
        cluster_num = 0
        for clusterplus in hsp_clusters:
            cluster = clusterplus[0]
            cluster_num += 1

            # Call function for printing visualization.
            print_cluster(clusterplus, hit_num, cluster_num, num_dots) #***

        ## ***Redundant?:
        ## Check that the clusters do not overlap with each other on the subject
        ## sequence.
        #for cluster1 in hsp_clusters:
        #    for cluster2 in hsp_clusters:
        #        if cluster1[0] != cluster2[0]:
        #            if clusters_overlap(cluster1[0], cluster2[0]):
        #                # Visualize overlapping clusters (for troubleshooting).
        #                startend = get_cluster_range(cluster1[0] + cluster2[0])
        #                print('Overlapping clusters:')
        #                print_cluster(cluster1,\
        #                        str(get_cluster_range(cluster1[0])),\
        #                        cluster_num, num_dots, startend)
        #                print_cluster(cluster2,\
        #                        str(get_cluster_range(cluster2[0])),\
        #                        cluster_num, num_dots, startend)

        #            ## Assert no overlap.
        #            #assert not clusters_overlap(cluster1[0], cluster2[0]),\
        #            #"""Clusters overlap: %s and %s""" %\
        #            #(cluster1[0][0].hit_id + str(get_cluster_range(cluster1[0])),\
        #            # cluster2[0][0].hit_id + str(get_cluster_range(cluster2[0])))

    ## Check that the clusters do not overlap with each other on the subject
    ## sequence.
    #for cluster1 in all_hsp_clusters:
    #    for cluster2 in all_hsp_clusters:
    #        if cluster1[0] != cluster2[0]:
    #            assert not clusters_overlap(cluster1[0], cluster2[0]),\
    #            """Clusters overlap: %s and %s""" %\
    #            (cluster1[0][0].hit_id + str(get_cluster_range(cluster1[0])),\
    #             cluster2[0][0].hit_id + str(get_cluster_range(cluster2[0])))

    # Sort HSPs according to E-value (the ranking may change because when
    # TBLASTN HSPs for the same scaffold sequence are split into those
    # representing potentially separate genes, then some may have higher
    # E-values).
    all_hsp_clusters.sort(key=lambda x: min([y.evalue for y in x[0]]))

    # Return the list of SearchIO HSP (not Hit) object clusters/lists.
    return all_hsp_clusters


def get_blastp_hit_seq_obj_and_coord(hit, max_gap):
    """Takes a SearchIO hit object (containing one or more HSPs) and returns a
    fasta sequence object likely to correspond to a single gene or the portion
    of the subject sequence that is actually homologous to the query sequence.
    """
    # Make a list of HSP objects within the given hit object.
    hsp_objs = []
    for hsp in hit:
        hsp_objs.append(hsp)

    ## Determine whether the top HSP is on the forward or reverse strand.
    top_hsp = hsp_objs[0]
    #fwd_strand = whether_fwd_strand(top_hsp) 

    ## Check that all the hsps are on the same strand.
    #hsps_on_different_strands = check_strands(hsp_objs, fwd_strand)

    ## If HSPs are on different strands, then make a reduced list of HSP objects
    ## composed of only those that are on the same strand as the top HSP.
    hsp_objs_red = hsp_objs 
    #if hsps_on_different_strands:
    #    hsp_objs_red = []
    #    for hsp in hsp_objs:
    #        x = whether_fwd_strand(hsp)
    #        # Only include in reduced list if on same strand as top HSP.
    #        if x == fwd_strand:
    #            hsp_objs_red.append(hsp)

    # (If more than one hsp)

    # Only include HSPs that are sufficiently proximate to the top HSP, and not
    # overlapping.
    proximate_hsp_objs = get_proximate_hsp_objs2(hsp_objs_red, top_hsp, max_gap)[0]

    # Concatenate the hsp sequences in order from 5' to 3' as they appear
    # in the nucleotide subject sequence.
    hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))

    # Make a description of where the hit sequence came from.
    #seq_origin = ','.join(['[' + str(x.hit_range[0]) + '..' + str(x.hit_range[1])\
    #    + ']' for x in proximate_hsp_objs])

    # Define sequence coordinates as a list.
    seq_coord = [[x.hit_range[0],x.hit_range[1]]\
        for x in proximate_hsp_objs]

    # Generate a fasta sequence record.
    hit_seq_record = SeqRecord(hit_sequence)
    hit_seq_record.id = hit.id
    descr_without_quotation_marks = hit.description.replace('\"', '')
    hit_seq_record.description = '\"' + descr_without_quotation_marks + ' ' +\
    str(seq_coord) + '\"'

    # Return the sequence record and sequence coordinates.
    return [hit_seq_record, seq_coord]




#def get_tblastn_hit_subseq_coord(hit, max_gap):
#    """Takes a SearchIO hit object (containing one or more HSPs) and returns a
#    coordinates for the subsequence(s) likely to correspond to the coding
#    sequence of a single gene (the most similar gene to the query in the
#    nucleotide subject sequence).
#    """
#    # Make a list of HSP objects within the given hit object.
#    hsp_objs = []
#    for hsp in hit:
#        hsp_objs.append(hsp)
#
#    # Determine whether the top HSP is on the forward or reverse strand.
#    top_hsp = hsp_objs[0]
#    fwd_strand = whether_fwd_strand(top_hsp) 
#
#    # Check that all the hsps are on the same strand.
#    hsps_on_different_strands = check_strands(hsp_objs, fwd_strand)
#
#    # If HSPs are on different strands, then make a reduced list of HSP objects
#    # composed of only those that are on the same strand as the top HSP.
#    hsp_objs_red = hsp_objs 
#    if hsps_on_different_strands:
#        hsp_objs_red = []
#        for hsp in hsp_objs:
#            x = whether_fwd_strand(hsp)
#            # Only include in reduced list if on same strand as top HSP.
#            if x == fwd_strand:
#                hsp_objs_red.append(hsp)
#
#    # (If more than one hsp)
#
#    # Only include HSPs that are sufficiently proximate to the top HSP, and not
#    # overlapping.
#    proximate_hsp_objs = get_proximate_hsp_objs(hsp_objs_red, top_hsp, max_gap)
#
#    # Concatenate the hsp sequences in order from 5' to 3' as they appear
#    # in the nucleotide subject sequence.
#    #hit_sequence = Seq(''.join([str(x.hit.seq).replace('-', '') for x in proximate_hsp_objs]))
#
#    # Make a description of where the hit sequence came from.
#    seq_origin = ','.join(['[' + str(x.hit_range[0]) + '..' + str(x.hit_range[1])\
#        + ']' for x in proximate_hsp_objs])
#
#    # Generate a fasta sequence record.
#    #hit_seq_record = SeqRecord(hit_sequence)
#    #hit_seq_record.id = hit.id
#    #descr_without_quotation_marks = hit.description.replace('\"', '')
#    #hit_seq_record.description = '\"' + descr_without_quotation_marks + ' ' + seq_origin + '\"'
#
#    #return hit_seq_record
#    print('seq_origin')
#    print(seq_origin)
#    return seq_origin


def parse_tblastn(tblastn_out_path, seq_out_path, max_gap):
    """Parse a given tblastn output file (XML) and write appropriate sequences
    corresponding to hits therein to a given sequence file path. Also, return
    relevant info for recording results. 
    """
    # Define input and output paths.
    outfile = tblastn_out_path
    outfile3 = seq_out_path

    # Parse tblastn output (using biopython "experimental code").
    fwd = open(outfile, 'U')
    fwd_recs = SearchIO.parse(fwd, 'blast-xml')
    rec_num = 0

    with open(outfile3, 'w') as o: 
        for record in fwd_recs:
            rec_num += 1
            hit_num = 0
            for hit in record:
                hit_num += 1
                hit_seq_obj = get_hit_seq_obj(hit, max_gap)

                # Write fasta string to output file.
                SeqIO.write(hit_seq_obj, o, 'fasta')


def search_scaffolds(query_file, subject_db, max_gap):
    """Uses functions to search a given nucleotide db with a given protein
    query.
    """
    # Define output for tblastn search.
    outfile = query_file.replace('.fa', '_tblastn_out.txt')
    outfile2 = query_file.replace('.fa', '_tblastn_out_text.txt')
    outfile3 = query_file.replace('.fa', '_tblastn_out_seqs.fa')

    # Run tblastn using given files.
    subprocess.call(["tblastn", "-query", query_file, "-db", subject_db, "-out",\
        outfile, "-outfmt", "5", "-evalue", str(0.05)])

    ## Do it again just to get text output (TEMPORARY).
    #subprocess.call(["tblastn", "-query", query_file, "-db", subject_db, "-out",\
    #    outfile2, "-evalue", str(0.05)])

    # Parse tblastn output file and write hit sequence to given file.
    parse_tblastn(outfile, outfile3, max_gap)



