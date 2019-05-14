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
"""Script for calculating the number of possible sequence sets to use for
phylogenetic analysis of a collection of homologous sequences, given a total

number of possible homologous sequences to include, the number of
proteins/clades of interest, and the number of taxonomic groups that ideally
must be represented in each clade.

This is important, because not only are there many possible trees given a
specific sequence set, but also many possible sequence sets, and for analyses
of protein family evolution as for SNAREs, selection of an informative, if not
the most informative, sequence set is critical. So, if there is a huge number
of such sets, then an efficient algorithm is required to search that set space.

The formula for calculating the number of combinations of a given size from
elements of a set of a given size [ n!/(k!(n-k)!) ] was taken from wikipedia
(https://en.wikipedia.org/wiki/Combination).

Note: The numbers are generally large enough that one could not feasibly run
trees of all the sequence sets. The next question then would be whether it is
feasible to start with an initial set and simply remove sequences until you get
resolution. The problem is, you can't run a tree of all homologues at once
because there won't be any resolution, so you don't know which are really long
branches or not. So, you have to start with a tree of representatives, and
remove OR add sequences in some justifiable way. How many such trees are
possible? It would be a subset of all the possible ones, perhaps you could
estimate this number by reducing the total number of sequences entered into
this script by half or more and increasing the number of clades to be
represented. Even then, the number of possibilities would perhaps be near
1e+18. Still too many datasets to run individual phylogenetic analyses for.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import math

cmdln = sys.argv
total = int(cmdln[1])
clades = int(cmdln[2]) 
taxon_groups = int(cmdln[3]) 

# Calculate minimum sequence set size.
min_seq_set_size = clades * taxon_groups

# Calculate how many possible sequence sets.
summation = 0
for k in range(min_seq_set_size, total):
    summation += (math.factorial(total)/(math.factorial(k) * math.factorial(total - k)))

print('Total number of possible sequence sets to use for phylogenetic analysis:')
print(summation)

