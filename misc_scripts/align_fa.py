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
"""This module is for calling alignment software (MUSCLE) to align fasta files
and save them as aligned fasta files.  This assumes that you want to use
default settings for MUSCLE.  By default, the input and output filepaths are
set to the current working directory, but these can be specified if the
align_all_fa function is imported by a script.  By default, it does all the .fa
files in a dir, but you can specify only one .fa file with the -i parameter.

Functions could be added for running different software, such as MAFFT.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))
import glob
import subprocess
import settings
from module_afa_to_nex import afa_to_nex
#import argparse

# # Set up argument parsing.
# parser = argparse.ArgumentParser(
#     description = """Outputs positive hits and sequences.""",
#     epilog = r"""epilog""")
# 
# # Positional arguments.
# parser.add_argument('indp', help='infilepath')
# 
# # Optional arguments.
# parser.add_argument('-p', '--program', help='Program (either muscle or mafft)',
#         default='muscle')
# parser.add_argument('-m', '--aamatrix', help='Amino acid matrix (NCBI format)',
#         default=None)
# 
# args = parser.parse_args()



def align_fa(infilepath, outfilepath, aamatrix_path, program='muscle'):
    """Calls MUSCLE or MAFFT."""
    # Run MUSCLE by default.
    if program is None:
        program='muscle'

    # Check that an alignment program was specified.
    #assert program == 'muscle' or program == 'mafft', 'Error: Must specify\
    #alignment program as either muscle or mafft.'

    if program == 'muscle':
        # Call MUSCLE with default options.
        subprocess.call(["muscle", "-in", infilepath, "-out", outfilepath])

    # Run MAFFT if specified using somewhat default parameters.
    elif program == 'mafft':
        # Get the matrix file path with user's home directory. 
        aamatrix_path = os.path.expanduser(aamatrix_path)
        with open(outfilepath, 'w') as o:
            # The --localpair option may not be optimal. There are two other
            # options suggested for accurate alignments in the mafft
            # documentation. The --genafpair option is what Surkont et al.
            # used.
            stdoutdata = subprocess.call(['mafft', '--maxiterate', '1000',
                '--ep', '0', '--genafpair', '--aamatrix', aamatrix_path, infilepath], stdout=o)

    # Run MAFFT with variable parameters (modify this code).
    elif program == 'altmafft':
        # Get the matrix file path with user's home directory. 
        aamatrix_path = os.path.expanduser(aamatrix_path)
        with open(outfilepath, 'w') as o:
            # The --localpair option may not be optimal. There are two other
            # options suggested for accurate alignments in the mafft
            # documentation.
            stdoutdata = subprocess.call(['mafft', '--maxiterate', '1000',
                '--localpair', '--aamatrix', aamatrix_path, infilepath], stdout=o)


def align_all_fa(indirpath=None, outdirpath=None):
    """Calls the align_fa_muscle function on all .fa files in a given
    directory, which is the cwd by default.  And, outputs to a given output
    directory, which is also cwd by default.
    """
    cwd = os.getcwd()
    if indirpath == None:
        indirpath = cwd
    if outdirpath == None:
        outdirpath = cwd
    for f in [x for x in glob.glob(os.path.join(indirpath, '*.fa'))]:
        outfilename = os.path.basename(f).replace('.fa', '.afa')
        align_fa(f, os.path.join(outdirpath, outfilename))

def align_one_fa(infilepath, outdirpath=None, program=None, aamatrix_path=None,
        conv_to_nex=None):
    """Calls the align_fa_muscle function on a .fa file.  And, outputs to a
    given output directory, which is the cwd by default.
    """
    if not infilepath.rsplit('.', 1)[1] == 'fa':
        print('\n*** Warning: The file specified does not have the extension\
fa.')
    if outdirpath == None:
        outdirpath = os.path.dirname(infilepath)
    outfilename = os.path.basename(infilepath).replace('.fa', '.afa')
    # Align with muscle with default settings and aamatrix.
    align_fa(infilepath, os.path.join(outdirpath, outfilename), aamatrix_path,\
            program)

    # Optionally convert output file to nex and delete afa.
    if conv_to_nex:
        o = os.path.join(outdirpath, outfilename)
        afa_to_nex(o, o.replace('.afa', '.nex'))
        os.remove(o)

if __name__ == '__main__':
    command_line_list = sys.argv
    indp = str(command_line_list[1])

    # Set amino acid scoring matrix path to use if mafft program called
    # (default is settings.blosum62_path).
    print('Matrix: ' + settings.align_fa_matrix)
    aamatrix_path = None
    if settings.align_fa_matrix == 'blosum62':
        aamatrix_path = settings.blosum62_path
    elif settings.align_fa_matrix == 'cc140':
        aamatrix_path = settings.cc140_path
    else:
        print('Error! matrix not specified.')
    print('Matrix file: ' + aamatrix_path)

    # Set alignment program to use (default is 'muscle', use 'altmafft' if
    # experimenting with different parameters).
    program = settings.align_fa_program 
    print('Program: ' + program + '\n')

    # Align the fasta file.
    align_one_fa(indp, None, program, aamatrix_path, True)

