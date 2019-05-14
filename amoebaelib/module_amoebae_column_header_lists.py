
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
# Define column header lists to be shared between functions in the
# module_amoebae_search module.

fwd_column_label_list = ['Query title',
                         'Query file',
                         'Query species (if applicable)',
                         'Query database name',
                         'Query accession (if applicable)',
                         'Query description',
                         'Query length',
                         'Subject database species (if applicable)',
                         'Subject database file',
                         'Forward search method',
                         'Forward hit rank',
                         'Forward hit score',
                         'Forward hit score difference from top hit score',
                         'Forward hit E-value (top HSP)',
                         'Forward hit E-value (top HSP) order of magnitude difference compared to top hit',
                         'Forward hit length',
                         'Forward hit length as a percentage of query length',
                         'Forward hit accession',
                         'Forward hit description',
                         'Forward hit sequence',
                         'Forward hit coordinates of subsequence(s) that align(s) to query',
                         'Forward hit description of subsequence(s) that align(s) to query',
                         'Forward hit subsequence(s) that align(s) to query',
                         'Proximity (bp) to end of subject sequence (if searching in nucleotide sequences)',
                         'Positive/redundant (+) or negative (-) hit based on E-value criterion'
                         ]


rev_column_label_list = ['Reverse search species (if applicable)',
                         'Reverse search database name',
                         'Redundant hit list applied',
                         'Top reverse hit description',
                         'Top reverse search hit E-value',
                         'Order of magnitude E-value difference between top hit and first non-redundant hit',
                         'Top reverse search hit score',
                         'Difference in score between top hit and first non-redundant hit', 'First non-redundant hit',
                         'Note',
                         'Positive (+) or negative (-) hit based on reverse search criteria'
                         ]

interp_column_label_list = ['Collective interpretation of reverse search results'
                            ]

phylo_class_column_label_list = [#'Model/backbone tree name',
                                 'Classification',
                                 'Second most likely classification',
                                 'AU topology test p-value for comparison with next most likely classification/topology',
                                 'All but top classification rejected?'
                                 ]

phylo_class_place_column_label_list = ['Model/backbone tree name',
                                     'Classification'
                                      ]
