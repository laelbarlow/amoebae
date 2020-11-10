#!/usr/bin/env python3

import glob
import os
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt


# Count how many function definitions there are.
funct_def_count = 0
for d in ['amoebaelib', 'misc_scripts']:
    for f in glob.glob(os.path.join(d, '*.py')):
        with open(f) as infh:
            for i in infh:
                if i.startswith("def "):
                    #print(i.strip())
                    funct_def_count += 1
print("Total number of function definitions: " + str(funct_def_count))

# Count how many functions have tests in files in the tests directory.
tested_def_count = 0
for d in ['tests']:
    for f in glob.glob(os.path.join(d, '*.py')):
        with open(f) as infh:
            contents = infh.read()
            split_on_defs = contents.split()
            for i in infh:
                if i.startswith("def "):
                    #print(i.strip())
                    funct_def_count += 1

objects = ('Total Defs', 'Tested Defs')
y_pos = np.arange(len(objects))
performance = [funct_def_count, 0]

plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Def count')
plt.title('Test coverage')

plt.show()
