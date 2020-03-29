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
#!/usr/bin/env python3
"""Tests functions in module_search_scaffolds.py module using unittest.
"""

# Import necessary modules.
import unittest
import module_search_scaffolds

# Define a class with functions for performing tests.
class Testmodule_search_scaffolds(unittest.TestCase):

    def test_check_if_two_hsp_ranges_overlap(self):

        input_value1 = [[1, 5], [3, 7]]
        result_value1 = True

        self.assertEqual(module_search_scaffolds.check_if_two_hsp_ranges_overlap(input_value1), result_value1)

        input_value2 = [[1, 5], [7, 9]]
        result_value2 = False

        self.assertEqual(module_search_scaffolds.check_if_two_hsp_ranges_overlap(input_value2), result_value2)

        input_value3 = [[1, 5], [1, 5]]
        result_value3 = True

        self.assertEqual(module_search_scaffolds.check_if_two_hsp_ranges_overlap(input_value3), result_value3)


    def test_find_missed_hsp_ranges(self):

        a =  [(1,10), (20,30), (40,50), (60,70), (80,90), (19,21), (29,41), (49,61), (69,71)]
        b =  [        (20,30), (40,50), (60,70)         ]
        c =  [   (19,21), (29,41), (49,61), (69,71)     ]
        result_value1 = False

        self.assertEqual(module_search_scaffolds.find_missed_hsp_ranges(a,b,c), result_value1)

        a =  [(1,10), (20,30), (40,50), (60,70), (80,90), (19,21), (29,41), (49,61), (69,71)]
        b =  [        (20,30), (40,50), (60,70)         ]
        c =  [   (19,21),          (49,61), (69,71)     ]
        result_value2 = True

        self.assertEqual(module_search_scaffolds.find_missed_hsp_ranges(a,b,c), result_value2)

        a =  [(1,10), (20,30), (40,50), (60,70), (80,90), (19,21), (29,41), (49,61), (69,71), (32,35)]
        b =  [        (20,30), (40,50), (60,70)         ]
        c =  [   (19,21), (29,41), (49,61), (69,71)     ]
        result_value3 = True

        self.assertEqual(module_search_scaffolds.find_missed_hsp_ranges(a,b,c), result_value3)

        a =  [(8466018, 8466147), (8466353, 8466488), (8466585, 8466801),
                (8468034, 8468208), (8468456, 8468645), (8468717, 8468861),
                (8469013, 8469343)]
        b =  [(8466018, 8466147), (8466353, 8466488), (8466585, 8466801),
                                    (8468456, 8468645), (8468717, 8468861),
                (8469013, 8469343)]
        c =  [ ]
        result_value4 = True

        self.assertEqual(module_search_scaffolds.find_missed_hsp_ranges(a,b,c), result_value4)

        #a =  [(1,10), (20,30), (40,50), (60,70), (80,90)]
        #b =  [(1,10), (20,30), (40,50), (60,70), (80,90)]
        #c =  [(1,10), (20,30), (40,50), (60,70), (80,90)]
        #result_value4 = True

        #self.assertEqual(module_search_scaffolds.find_missed_hsp_ranges(a,b,c), result_value4)



if __name__ == '__main__':
    unittest.main()
