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
"""For testing functions in module_nhmmer_scaffolds.py using unittest.

Usage:
    test_module_nhmmer_scaffolds.py

"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import unittest
import module_nhmmer_scaffolds

# Google "mocking" for testing

class TestNhmmerScaffolds(unittest.TestCase):

    def test_frameify_seq(self):
        inseq1 = 'ATGATGATG'
        hmmfrom1 = '1'
        result1 = 'ATGATGATG'

        self.assertEqual(module_nhmmer_scaffolds.frameify_seq(inseq1, hmmfrom1),
                result1)

        inseq2 = 'TGATGATG'
        hmmfrom2 = '2'
        result2 = 'ATGATG'

        self.assertEqual(module_nhmmer_scaffolds.frameify_seq(inseq2, hmmfrom2),
                result2)

        inseq3 = 'GATGATG'
        hmmfrom3 = '3'
        result3 = 'ATGATG'

        self.assertEqual(module_nhmmer_scaffolds.frameify_seq(inseq3, hmmfrom3),
                result3)

    #def test_check_for_duplicate_ranges(self):
    #    tuple_list1 = [(10,11),(10,11)]
    #    tuple_list2 = [(0,11),(10,11)]
    #    self.assertTrue(module_search_scaffolds.check_for_duplicate_ranges(tuple_list1))
    #    self.assertFalse(module_search_scaffolds.check_for_duplicate_ranges(tuple_list2))


    #def test_sort_left_ranges(self):
    #    range_list1 = [(1,10),(11,20)]
    #    result1 = [(11,20),(1,10)]
    #    self.assertEqual(module_search_scaffolds.sort_left_ranges(range_list1), result1)


    #def test_sort_right_ranges(self):
    #    range_list1 = [(11,20),(1,10)]
    #    result1 = [(1,10),(11,20)]
    #    self.assertEqual(module_search_scaffolds.sort_right_ranges(range_list1), result1)


    #def test_reduce_left_ranges(self):
    #    # Test for inclusion of ranges within maximum gap.
    #    sorted_left_ranges = [(51,70),(31,40),(11,20)]
    #    top_hsp_range = (81,100)
    #    max_gap = 15
    #    result1 = [(51,70),(31,40),(11,20)]        
    #    self.assertEqual(module_search_scaffolds.reduce_left_ranges(sorted_left_ranges,
    #        top_hsp_range, max_gap), result1)

    #    # Test for exclusion of ranges outside maximum gap.
    #    sorted_left_ranges = [(51,70),(31,40),(11,20)]
    #    top_hsp_range = (81,100)
    #    max_gap = 5
    #    result1 = []
    #    self.assertEqual(module_search_scaffolds.reduce_left_ranges(sorted_left_ranges,
    #        top_hsp_range, max_gap), result1)

    #    # Test for exclusion of overlapping ranges.
    #    sorted_left_ranges = [(51,82),(31,52),(11,32)]
    #    top_hsp_range = (81,100)
    #    max_gap = 15
    #    result1 = []
    #    self.assertEqual(module_search_scaffolds.reduce_left_ranges(sorted_left_ranges,
    #        top_hsp_range, max_gap), result1)

    #    # Test for correct measurement of maximum gap length.
    #    sorted_left_ranges = [(1,2)]
    #    top_hsp_range = (3,4)
    #    max_gap = 0
    #    result1 = [(1,2)]
    #    self.assertEqual(module_search_scaffolds.reduce_left_ranges(sorted_left_ranges,
    #        top_hsp_range, max_gap), result1)


    #def test_reduce_right_ranges(self):
    #    # Test for inclusion of ranges within maximum gap.
    #    sorted_right_ranges = [(21,30),(41,60),(71,90)]        
    #    top_hsp_range = (1,10)
    #    max_gap = 15
    #    result1 = [(21,30),(41,60),(71,90)]        
    #    self.assertEqual(module_search_scaffolds.reduce_right_ranges(sorted_right_ranges,
    #        top_hsp_range, max_gap), result1)

    #    # Test for exclusion of ranges outside maximum gap.
    #    sorted_right_ranges = [(21,30),(41,60),(71,90)]        
    #    top_hsp_range = (1,10)
    #    max_gap = 5
    #    result1 = []
    #    self.assertEqual(module_search_scaffolds.reduce_right_ranges(sorted_right_ranges,
    #        top_hsp_range, max_gap), result1)

    #    # Test for exclusion of overlapping ranges.
    #    sorted_right_ranges = [(21,42),(41,72),(71,90)]        
    #    top_hsp_range = (1,22)
    #    max_gap = 15
    #    result1 = []
    #    self.assertEqual(module_search_scaffolds.reduce_right_ranges(sorted_right_ranges,
    #        top_hsp_range, max_gap), result1)

    #    # Test for correct measurement of maximum gap length. 
    #    sorted_right_ranges = [(3,4)]        
    #    top_hsp_range = (1,2)
    #    max_gap = 0
    #    result1 = [(3,4)]
    #    self.assertEqual(module_search_scaffolds.reduce_right_ranges(sorted_right_ranges,
    #        top_hsp_range, max_gap), result1)


    #def test_get_proximate_hsp_ranges(self):
    #    #range_list1 = [(1,10),(11,20),(21,30),(31,40),(41,60),(51,70),(61,80),(71,90),(81,100),(91,110)]

    #    # Test whether it can handle ranges that do not overlap, and are spaced
    #    # closer than the maximum gap.
    #    range_list2 = [(1,10),(21,30),(41,60),(81,100),(121,140)]
    #    top_hsp_range1 = (41,60)
    #    result1 = [(1,10),(21,30),(41,60),(81,100),(121,140)]
    #    max_gap1 = 10

    #    self.assertEqual(module_search_scaffolds.get_proximate_hsp_ranges(range_list2,
    #        top_hsp_range1, max_gap1), result1)

    #    # Test whether it will properly handle gaps larger than the maximum.
    #    range_list2 = [(1,10),(21,30),(41,60),(81,100),(121,140)]
    #    top_hsp_range1 = (41,60)
    #    max_gap2 = 5
    #    result2 = [(41,60)]
    #    self.assertEqual(module_search_scaffolds.get_proximate_hsp_ranges(range_list2,
    #        top_hsp_range1, max_gap2), result2)


if __name__ == '__main__':
    unittest.main()