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
"""For testing functions in module_amoebae_get_datatype.py using unittest.
Usage:
    test_module_amoebae_get_datatype.py
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import unittest
from module_amoebae_get_datatype import get_datatype_for_sequence_string
from test_module_amoebae_get_datatype_testseqs import testseq1


class TestSearchScaffolds(unittest.TestCase):

    def test_get_datatype_for_sequence_string1(self):
        
        x = 'ATGCY'

        dbtype = get_datatype_for_sequence_string(x)

        self.assertEqual(dbtype, 'nucl')


    def test_get_datatype_for_sequence_string2(self):
        
        x = 'atgcy'

        dbtype = get_datatype_for_sequence_string(x)

        self.assertEqual(dbtype, 'nucl')


    def test_get_datatype_for_sequence_string3(self):
        
        x = 'ATGCV'

        dbtype = get_datatype_for_sequence_string(x)

        self.assertEqual(dbtype, 'prot')


    def test_get_datatype_for_sequence_string4(self):
        
        x = 'atgcv'

        dbtype = get_datatype_for_sequence_string(x)

        self.assertEqual(dbtype, 'prot')


    def test_get_datatype_for_sequence_string5(self):
        """Test whether 'Y' can be handled."""

        x = 'tccaaaaaatcgaaTTATYttattccccaccttcttttctcattttttga'

        dbtype = get_datatype_for_sequence_string(testseq1)

        self.assertEqual(dbtype, 'nucl')


    def test_get_datatype_for_sequence_stringX(self):

        x = testseq1

        dbtype = get_datatype_for_sequence_string(x)

        self.assertEqual(dbtype, 'nucl')


if __name__ == '__main__':
    unittest.main()