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
"""This module is for making HMMs from alignment files in nexus format by first
converting them to aligned fasta files and then running hmmbuild.
"""

import os
import sys
import subprocess
from afa_to_nex import nex_to_afa, delete_extra_mesquite_lines


def nex_to_hmm(infilepath, outfilepath=None):
    """Takes a nexus file and outputs an hmm."""
    # Set output dir to input dir by default.
    if outfilepath is None:
        outfilepath = infilepath.rsplit('.', 1)[0] + '.hmm'
    
    # ?
    outdirpath = os.path.dirname(infilepath)

    # Delete extra blocks added by Mesquite, if present.
    delete_extra_mesquite_lines(infilepath)

    # Convert nex to afa.
    outafapath1 = os.path.join(outdirpath,
            os.path.basename(infilepath).replace('.nex', '_tmp.afa'))
    nex_to_afa(infilepath, outafapath1)

    # Replace '?' characters with 'X'.
    outafapath2 = os.path.join(outdirpath,
            os.path.basename(infilepath).replace('.nex', '.afa'))
    with open(outafapath1) as infh1, open(outafapath2, 'w') as o:
        for i in infh1:
            if not i.startswith('>'):
                # Replace problematic characters with 'X'
                o.write(i.replace('?', 'X'))

            else:
                o.write(i)

    # Call hmmbuild to make HMM.
    subprocess.call(["hmmbuild", outfilepath,
        outafapath2])

    # Delete intermediate afa file.
    os.remove(outafapath1)
    os.remove(outafapath2)


def afa_to_hmm(infilepath, outfilepath=None):
    """Takes a aligned fasta file and outputs an hmm."""
    # Set output dir to input dir by default.
    if outfilepath is None:
        outfilepath = infilepath.rsplit('.', 1)[0] + '.hmm'
    
    # ?
    outdirpath = os.path.dirname(infilepath)

    # Replace '?' characters with 'X'.
    temp_afa_1 = infilepath.rsplit('.', 1)[0] + '_temp1.afa' 
    with open(infilepath) as infh1, open(temp_afa_1, 'w') as o:
        for i in infh1:
            if not i.startswith('>'):
                # Replace problematic characters with 'X'
                o.write(i.replace('?', 'X'))

            else:
                o.write(i)

    # Call hmmbuild to make HMM.
    subprocess.call(["hmmbuild", outfilepath,
        temp_afa_1])

    # Delete intermediate afa file.
    os.remove(temp_afa_1)


