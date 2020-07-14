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
"""This file contains definitions of variables used by scripts in the AMOEBAE
toolkit. 
"""
import os

class DataPaths:
    """Provides an object with directory and file paths as attributes based on
    a given main AMOEBAE data directory path.
    """
    def __init__(self, root_amoebae_data_dir):
        # Check that the input path exists as a directory.
        assert os.path.isdir(root_amoebae_data_dir), """Input AMOEBAE data
        directory path does not exist."""

        # Set path to directory containing genome data (predicted peptide sequences,
        # nucleotide scaffolds, etc.).
        self.dbdirpath =\
        os.path.join(root_amoebae_data_dir, 'Genomes')
        assert os.path.isfile(self.dbdirpath)
        
        # Set path to csv file with information about the relevant genome data.
        self.db_info_csv =\
        os.path.join(dbdirpath, '0_genome_info.csv')
        assert os.path.isfile(self.db_info_csv)
        
        # Set path to directory containing query files (.faa, .afaa, etc.).
        self.querydirpath =\
        os.path.join(root_amoebae_data_dir, 'Queries')
        assert os.path.isfile(self.querydirpath)
        
        # Set path to csv file with information about the relevant query files.
        self.query_info_csv =\
        os.path.join(querydirpath, '0_query_info.csv')
        assert os.path.isfile(self.query_info_csv)
        
        # Set path to directory containing reference tree (and alignment) files.
        self.model_dir_path =\
        os.path.join(root_amoebae_data_dir, 'Models')
        assert os.path.isfile(self.model_dir_path)
        
        # Set path to csv file with information about reference tree (and alignment)
        # files.
        self.model_info_csv =\
        os.path.join(model_dir_path, '0_models_info.csv')
        assert os.path.isfile(self.model_info_csv)


