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
"""Module defining functions for testing import statements in the AMOEBAE
repository.
"""
import glob
import os
import sys
import subprocess
import random
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'amoebaelib'))
sys.path.append(os.path.dirname(sys.path[0]))

def find_and_test_all_import_statements_in_amoebae():
    """Find all import statements used in the AMOEBAE repository, write to a
    file, and run the file to see whether the imports work (i.e., whether all
    the modules are installed and accurately named on the current system).
    """

    # Compile a list of paths for files of interest with import statements.
    all_file_paths = \
          glob.glob(os.path.join(sys.path[0],'misc_scripts/*.py'))\
        + glob.glob(os.path.join(sys.path[0],'amoebaelib/*.py'))\
        + glob.glob(os.path.join(sys.path[0],'*.py'))\
        + glob.glob(os.path.join(sys.path[0],'amoebae'))
    
    # Define output file path.
    random_id_string = str(random.randint(1, 1000))
    output_path = os.path.join(sys.path[0],
            'get_nonredun_import_statments_for_amoebae_output_' +\
            random_id_string + '.py')
    
    # Compile a complete list of import statements.
    nonredun_list = []
    for f in all_file_paths:
        if f != output_path:
            with open(f) as infh:
                last_line_was_import = False
                last_line_ended_with_backslash = False
                for i in infh:
                    # Determine whether the line contains an import statement.
    
                    # If the last line was an import statement that ended with a
                    # backslash, then the current line is part of an import
                    # statement too.
                    if last_line_was_import and last_line_ended_with_backslash:
                        # Add line to list.
                        x = (i.strip(), os.path.basename(f))
                        if not x[0] in [y[0] for y in nonredun_list]: 
                            nonredun_list.append(x)
                        
                        # Current line is an import line, so swithch the variable
                        # to indicate for the next iteration of the loop.
                        last_line_was_import = True
    
                        if i.endswith('\\\n'):
                            last_line_ended_with_backslash = True
                        else:
                            last_line_ended_with_backslash = False
    
                    # Lines that are formatted like an import statement are import
                    # statements.
                    elif i.startswith('import ') or i.startswith('from ') and '. ' not in i:
                        # Add line to list.
                        x = (i.strip(), os.path.basename(f))
                        if not x[0] in [y[0] for y in nonredun_list]: 
                            nonredun_list.append(x)
    
                        # Current line is an import line, so swithch the variable
                        # to indicate for the next iteration of the loop.
                        last_line_was_import = True
    
                        if i.endswith('\\\n'):
                            last_line_ended_with_backslash = True
                        else:
                            last_line_ended_with_backslash = False
    
                    else:
                        last_line_was_import = False
    
    # Record the nonredundant import statements in an output file.
    with open(output_path, 'w') as o:
        o.write("""\
#!/usr/bin/env python3
import sys, os
sys.path.append(os.path.join(sys.path[0],'amoebaelib'))
sys.path.append(os.path.join(sys.path[0],'misc_scripts'))
""")
    
        print('\nNon-redundant list of import statements:\n')
        num = 0
        for i in nonredun_list:
            ## Exclude import statements that import modules from the amoebaelib
            ## library (may not want to do this). ****
            #if ' module_' not in i and '. ' not in i: 
            num += 1
            print(str(num) + '. ' + i[0] + '  # ' + i[1])
            #o.write('print(\"' + i.rstrip('\\') + '\")' + '\n')
            o.write(i[0] + '\n')

        # Write line to print success message.
        o.write('\nprint("All import statements ran successfully.")')
    
    
    # Run the output file to see whether all the modules could be imported.
    print('\nRunning output script to test import statements...')
    subprocess.call(['chmod', 'a+x', output_path])
    subprocess.call(['./' + os.path.basename(output_path)], cwd=sys.path[0])

    # Remove output file.
    os.remove(output_path)

    
