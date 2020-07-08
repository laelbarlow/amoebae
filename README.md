
<p align="center">
<img src="images/AMOEBAE_logo10.png" width="600">
</p>

## Quickstart

On **MacOS** or **Linux**, open a terminal window and enter the commands
described below. 

Clone this git repository:

    git clone https://github.com/laelbarlow/amoebae.git

Run the setup script (this will take about 5 minutes to complete). On
**MacOS**, and on **Linux** if Singularity is not installed, this will install
Vagrant and VirtualBox for managing and running Linux virtual machines, if not
already installed:

    cd amoebae
    /bin/bash setup.sh

Now open one of the tutorial jupyter notebooks using this script (this will
only work on a personal computer, see main documentation for alternatives).

    /bin/bash singularity_jupyter.sh


<!---
***Embedded video...

--->

## Introduction

**A**nalysis of **MO**lecular **E**volution with **BA**tch **E**ntry (AMOEBAE)
is a bioinformatics software toolkit primarily composed of scripts and modules
written in the Python3 language. AMOEBAE performs certain tasks involved in
identifying and classifying amino acid sequences. Specifically, AMOEBAE uses
Python packages including Biopython, ETE3, Pandas, and Matplotlib for setting
up, running, and summarizing analyses of molecular evolution using
bioinformatics software packages including MUSCLE, BLAST+, HMMER3, and IQ-TREE.
Applications include identifying and classifying predicted genes in novel
biological sequence data according to their evolutionary relationships with
homologues. AMOEBAE is open-source, and all dependencies are freely available.
However, please note that AMOEBAE is currently under development. [Lael D.
Barlow](https://scholar.google.com/citations?user=wohF-LAAAAAJ&hl=en&oi=ao) is the author. 


## Documentation

Instructions for automatically setting up a virtual environment
with all the dependencies are included in the main documentation file
here: [AMOEBAE_documentation.pdf](
https://github.com/laelbarlow/amoebae/blob/master/documentation/AMOEBAE_documentation.pdf). The whole process should take approximately 5 minutes (although more time may be necessary if you need to install Git or clear storage space). 

See the main documentation file for an overview and setup instructions:
[AMOEBAE_documentation.pdf](
https://github.com/laelbarlow/amoebae/blob/master/documentation/AMOEBAE_documentation.pdf).
Once you have set up an environment for running AMOEBAE, see the tutorials
which are in the form of Jupyter notebooks with code that runs analyses using
example data: [AMOEBAE
notebooks](https://github.com/laelbarlow/amoebae/tree/master/notebooks). For
documentation on the AMOEBAE code library, please see the html file(s), which
can be opened with your web browser: [Code
documentation](https://github.com/laelbarlow/amoebae/blob/master/documentation/code_documentation/html/index.html).


## Reporting problems

Please use the [issue tracker](https://github.com/laelbarlow/amoebae/issues) on
the GitHub webpage to report any problems you find in using AMOEBAE.


## How to cite AMOEBAE

Please cite the [AMOEBAE GitHub
repository](https://github.com/laelbarlow/amoebae) (or alternative permanent
repositories if relevant). Also, please remember to cite the software packages which
are dependencies of AMOEBAE (see main documentation).

## Acknowledgments

AMOEBAE was initially developed at the Dacks Laboratory at the University of
Alberta, and was supported by National Sciences and Engineering Council of
Canada (NSERC) Discovery grants RES0021028, RES0043758, and RES0046091 awarded
to Joel B. Dacks, as well as an NSERC Postgraduate Scholarship-Doctoral awarded
to Lael D. Barlow.

We acknowledge the support of the Natural Sciences and Engineering Research Council of Canada (NSERC).

Cette recherche a été financée par le Conseil de recherches en sciences naturelles et en génie du Canada (CRSNG).

<img src="images/NSERC_FIP_RGB.png">


Also, help with testing AMOEBAE has been kindly provided by Raegan T. Larson,
Shweta V. Pipalya, Kira More, and Kristína Záhonová.


## License

Copyright 2018 Lael D. Barlow

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


