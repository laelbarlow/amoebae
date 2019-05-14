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
"""This module contains functions for parsing and modifiying modification
numbers in file names.
"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
import glob
import re


def get_mod_num(infilename):
    """Takes a filename (string) and returns the modification number it
    contains as a string.
    """
    # Compile regular expression, and use it to search in the string.
    p = re.compile(r'_mod(\w+[.])+')
    x = p.search(infilename)

    # If the mod number was identified, assign it to a variable.
    mod_num = 'Modification number not identified in input file name.'
    if x is None:
        print('Modification number not identified in input file name.')   
    else:
        mod_num = infilename[x.start() + 4:x.end()]

    return  mod_num

def get_all_other_mod_nums(infilepath):
    """Takes a file path and returns a list of modification numbers for other
    versions of the file in the same directory.
    """
    # Get a list of file paths that are the same as the input except for the
    # modification number.
    filepath_without_mod_num = infilepath.rsplit('_mod', 1)[0]
    filepaths = glob.glob(filepath_without_mod_num + '*')

    # Check whether there are actually other versions.
    if len(filepaths) < 2:
        print('There are no other versions')
        pass

    # Get a list of file names from the list of file paths.
    filenames = []
    for f in filepaths:
        x = os.path.basename(f)
        filenames.append(x)

    # Get a list of modification numbers from list of file names.
    filemodnums = []
    for f in filenames:
        x = get_mod_num(f)
        filemodnums.append(x)

    # Remove the mod num of the input file from the list.
    filemodnums.remove(get_mod_num(os.path.basename(infilepath)))

    return filemodnums

def get_last_mod_num(mod_num_list):
    """Takes a list of mod numbers, sorts it numerically/alphabetically, and
    returns the last mod num.
    """
    return sorted(mod_num_list, key=last_value)[-1]

def last_value(last_mod_num):
    """Takes a mod num and returns the last numerical or alphabetical value it
    contains.
    """
    # Compile regular expression, and use it to search in the mod num.
    p = re.compile(r'\w+[.]$')
    x = p.search(last_mod_num)
    return last_mod_num[x.start():x.end() -1]

def get_next_mod_value(last_value):
    """Takes a numerical or alphabetical value and returns the next value in
    sequence.
    """
    if last_value.isnumeric():
        return str(int(last_value) + 1)
    elif last_value.isalpha():
        return str(chr(ord(last_value) + 1))

def replace_last_value(mod_num, next_mod_value):
    """Replaces the last value in the mod num with a new character (numeric
    or alphabetical).
    """
    # Compile regular expression, and use it to search in the mod num.
    p = re.compile(r'\w+[.]$')
    new_mod_num = p.sub(next_mod_value + '.', mod_num)
    return new_mod_num

def update_mod_num(infilepath, char_type='number'):
    """Appends the appropriate number or letter to the modification number.
    """
    mod_num = get_mod_num(infilepath)
    all_other_mod_nums = get_all_other_mod_nums(infilepath)

    # Identify any child modifications.
    child_mod_nums0 = []
    for num in all_other_mod_nums:
        if num.startswith(mod_num) and len(num) >= (len(mod_num) + 1):
            child_mod_nums0.append(num)

    # Remove extra characters from the child mod numbers, in a new list.
    child_mod_nums = []
    for childmod in child_mod_nums0:
        child_mod_nums.append(childmod[:len(mod_num) + 2])
    
    # If no child modifications, then add '1.' or 'A.' to mod num.
    new_mod_num = None
    if len(child_mod_nums) == 0:
        if char_type == 'number': 
            new_mod_num = mod_num + '1.'
        elif char_type == 'letter':
            new_mod_num = mod_num + 'A.'

    # If child mods, then determine what next child mod num should be.
    else:
        # Split the list of child mods into those ending in digits and those
        # ending in letters.
        number_ends = [x for x in child_mod_nums if\
                last_value(x).isnumeric()]
        letter_ends = [x for x in child_mod_nums if\
                last_value(x).isalpha()]

        # Determine what the new mod num should be.
        new_mod_num = None

        if char_type == 'number' and len(number_ends) is not 0:
            # Increase number value.
            last_mod_num = get_last_mod_num(number_ends)
            last_mod_value = last_value(last_mod_num)
            next_mod_value = get_next_mod_value(last_mod_value)
            new_mod_num = replace_last_value(last_mod_num, next_mod_value)

        elif char_type == 'letter' and len(letter_ends) is not 0:
            # Increase letter value.
            last_mod_num = get_last_mod_num(letter_ends)
            last_mod_value = last_value(last_mod_num)
            next_mod_value = get_next_mod_value(last_mod_value)
            new_mod_num = replace_last_value(last_mod_num, next_mod_value)

        else:
            # Make first child mod of numeric or alphabetic type.
            if char_type == 'number': 
                new_mod_num = mod_num + '1.'
            elif char_type == 'letter':
                new_mod_num = mod_num + 'A.'

    # Get the new file path by replacing the old mod num with the new one.    
    new_file_path = infilepath.replace(mod_num, new_mod_num, -1)
    
    return new_file_path

def update_mod_num_numeric(infilepath):
    """Takes a file path, and returns one with an updated modification number,
    either adding a '.1' if there are no child mods, or increasing the numeric
    value of the last child mod if there are child mods.
    """
    return update_mod_num(infilepath, 'number')

def update_mod_num_alphabetic(infilepath):
    """Takes a file path, and returns one with an updated modification number,
    either adding a '.A' if there are no child mods, or increasing the
    alphabetic value of the last child mod if there are any child mods.
    """
    return update_mod_num(infilepath, 'letter')
    
if __name__ == '__main__':
    command_line_list = sys.argv
    infp = str(command_line_list[1])
    print(update_mod_num_alphabetic(infp))
