## Quickstart

Open a terminal window and enter the following commands. This works on macOS.

Clone this git repository:

    git clone https://github.com/laelbarlow/amoebae.git

Run the setup script (this will take about 1 hour to complete):

    cd amoebae
    bash setup.sh

Now open one of the tutorial jupyter notebooks:

    bash run_env.sh


<!---
***Embedded video...

--->

## Introduction

**A**nalysis of **MO**lecular **E**volution with **BA**tch **E**ntry (AMOEBAE)
is a bioinformatics software toolkit composed of scripts and modules written in
the Python3 language. AMOEBAE performs certain tasks involved in identifying
and classifying amino acid sequences. Specifically, AMOEBAE uses Python
packages including Biopython, ETE3, pandas, and Matplotlib for setting up,
running, and summarizing analyses of molecular evolution using bioinformatics
software packages including MUSCLE, BLAST+, HMMER3, and IQ-Tree.  Applications
include identifying and classifying predicted genes in novel biological
sequence data according to their evolutionary relationships with homologues.
AMOEBAE is open-source, and all dependencies are freely available. However,
please note that AMOEBAE is currently under development. 


## Documentation

Setting up AMOEBAE is easy. Instructions for automatically setting up a virtual environment
with all the dependencies are included in the main documentation file
here: [AMOEBAE_documentation.pdf](
https://github.com/laelbarlow/amoebae/blob/master/documentation/AMOEBAE_documentation.pdf). The whole process takes approximately 1 hour. 

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
the GitHub webpage to report any problems you find in using this toolkit.


## How to cite AMOEBAE

Please cite the [GitHub webpage](https://github.com/laelbarlow/amoebae) (or alternative
permanent repositories if relevant). Also, remember to cite other software
packages which you may have also installed and used for your work via AMOEBAE
scripts.

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


