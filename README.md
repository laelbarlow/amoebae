

## Introduction

### Synopsis 
**A**nalysis of **MO**lecular **E**volution with **BA**tch **E**ntry (AMOEBAE)
is a bioinformatics software toolkit composed of scripts written in the Python3
language. AMOEBAE performs certain tasks involved in identifying and
classifying amino acid sequences. Specifically, AMOEBAE uses Python packages
including Biopython, ETE3, pandas, and Matplotlib for setting up, running, and
summarizing analyses of molecular evolution using bioinformatics software
packages including MUSCLE, BLAST+, HMMER3, and IQ-Tree.  Applications include
identifying and classifying predicted genes in novel biological sequence data
according to their evolutionary relationships with homologues.  AMOEBAE is
open-source, and all dependencies are freely available. However, please note
that AMOEBAE is currently under development and is not yet fully functional. 


## Setup

Setting up AMOEBAE is easy. Instructions for automatically installing a virtual environment
with all the dependencies are included in the main documentation file
here:[AMOEBAE_documentation.pdf](
https://github.com/laelbarlow/amoebae/blob/master/documentation/AMOEBAE_documentation.pdf). The whole process takes approximately 1 hour. 

<!---

1) Ensure that the following dependencies are installed and working on your
system:

- Python3 (Anaconda distribution works well on personal computers)

- pandas (packaged with the Anaconda Python3 distribution)

- Biopython (easily installed via Anaconda) ([Cock *et  al*., 2009](https://academic.oup.com/bioinformatics/article/25/11/1422/330687))

- exonerate (easily installed via Anaconda)

- ETE3 ([Huerta-Cepas *et  al*., 2016](https://academic.oup.com/mbe/article/33/6/1635/2579822))

- [gffutils](https://pythonhosted.org/gffutils/) 

- NCBI BLAST+ ([Camacho *et  al*., 2009](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-421))

- HMMer3 ([Eddy, 1998](https://doi.org/10.1093/bioinformatics/14.9.755))

- IQtree ([Nguyen *et  al*., 2015](https://academic.oup.com/mbe/article/32/1/268/2925592))

- MUSCLE ([Edgar, 2004](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-5-113))

-->


<!---
1) Install ([Docker](https://www.docker.com/products/docker-desktop) for
building virtual environments in which to install dependencies.

2) Download the amoebae source code, or clone the amoebae repository using git,
add the directory to your $PATH variable for access via the command line (as
well as the misc\_scripts directory, if you wish to access those files).

3) Save a copy of the settings.py.example file as settings.py, and customize
the relevant variable definitions therein.

4) Enter the following command in a terminal window, which should display a
list of commands that you can use with amoebae:
```
    amoebae -h
```

You can access further information about each listed command with the -h option
as follows:
```
    amoebae <commandnamehere> -h
```

-->


## Documentation

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


