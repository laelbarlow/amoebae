#/usr/bin/env python3
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
"""Changes module filenames and updates all cases where the names of the
modules are used in all relevant files.

***Note: Make sure to backup files before running this script!
"""
import os
import glob
import pprint
import subprocess

# Identify necessary name changes.
filepathdict = {}
namedict = {}
newnames = []
for module in glob.glob(os.path.join('../amoebaelib', '*.py')):
    basename = os.path.basename(module)
    dirname = os.path.dirname(module)
    if basename.startswith('module'):
        basename2 = basename.replace('module_amoebae_', '').replace('module_',
                '')

        # Add new name to list of new names.
        newnames.append(basename2)

        # Add strings to dict for replacements of module imports and calls.
        namedict[basename.rsplit('.', 1)[0]] = basename2.rsplit('.', 1)[0]
        
        # Define updated path for module file.
        newfilepath = os.path.join(dirname, basename2)

        # Add entry to dict of old and new file paths.
        filepathdict[module] = newfilepath


assert len(newnames) == len(list(set(newnames)))

print('\nDict of file path changes to make:')
pprint.pprint(filepathdict)

print('\nDict of string replacements to make within files:')
pprint.pprint(namedict)


# Change module filenames.
print('\n')
for oldpath in filepathdict.keys():
    command =  ['git', 'mv', oldpath, filepathdict[oldpath]]
    print(' '.join(command))
    #subprocess.call(command)


# Find instances within files that need to change (this is independent of what
# it is that needs to change).

# Define relevant files.
relevant_files = ['../amoebae'] \
                 + glob.glob(os.path.join('../amoebaelib', '*.py')) \
                 + glob.glob(os.path.join('../misc_scripts', '*.py'))


# Loop over relevant files and update contents to use new module filenames.
num_lines_changed = 0
for f in relevant_files:
    # Define temporary filename.
    f2 = f + '.tmp'
    # Open the file and a new temporary file.
    with open(f) as infh, open(f2, 'w') as o:
        for i in infh:
            newi = i
            for key in namedict.keys():
                if key in i:
                    newi = i.replace(key, namedict[key])
                    #print(i.strip())
                    #print(newi.strip())
                    #print('\n')
                    num_lines_changed += 1
                #if namedict[key] in i and not key in i:
                #    print(i.strip())
            # Write (updated) line to temporary file.
            #o.write(newi)
    # Remove original file.
    #os.remove(f)
    # Rename temporary file with original filename.
    #os.rename(f2, f)
print('\nNumber of lines changed:')
print(num_lines_changed)


# **** Potential problem: some functions have the same names as module filenames. 
