

<p align="center">
<img src="images/AMOEBAE_logo10.png" width="550">
</p>


## Introduction

**A**nalysis of **MO**lecular **E**volution with **BA**tch **E**ntry (AMOEBAE)
is a customizable bioinformatics workflow for identifying homologues (and
potential orthologues) of genes of interest among a mid-size sampling of
genomes. This workflow is designed to be run on high-performance computing
(HPC) clusters and is executed via the
[SnakeMake](https://snakemake.readthedocs.io/en/stable/) workflow management
system. Code for steps in this workflow is written primarily in
[Python3](https://www.python.org/),
relying heavily on the [Biopython](https://biopython.org/) library, and apply
bioinformatics packages including
[BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/) and
[HMMER3](http://hmmer.org/) to input data files. AMOEBAE is open-source, and
all dependencies are freely available. [Lael D.
Barlow](https://scholar.google.com/citations?user=wohF-LAAAAAJ&hl=en&oi=ao) is
the author.

AMOEBAE is useful for certain mid-scale comparative genomics studies that
might otherwise require time-intensive and repetitive
manual/visual manipulation of data. Webservices such as those provided by
[NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
and [EMBL-EBI](https://www.ebi.ac.uk/Tools/hmmer/) provide a means to
readily investigate the evolution of one or a few genes via similarity searching,
and large-scale analysis workflows such as
[OrthoMCL](https://orthomcl.org/orthomcl/) and
[OrthoFinder](https://github.com/davidemms/OrthoFinder) attempt to rapidly perform orthology
prediction for all genes among several genomes. AMOEBAE addresses 
analyses which are too cumbersome to be performed via webservices or simple
scripts and yet require a level of detail and flexibility not offered by
large-scale analysis workflows. AMOEBAE is useful for analyzing the
distribution of homologues of up to approximately 30 genes/proteins among a
sampling of no more than approximately 100 eukaryotic genomes, especially when
follow-up with custom phylogenetic analysis is planned.

AMOEBAE serves this purpose by providing several unique features. The core
functionality of AMOEBAE is to run sequence similarity searches with multiple
algorithms, multiple queries, and multiple databases simultaneously and to
allow highly customizable implementation of reciprocal-best-hit search
strategies. The output includes detailed summaries of results in the form of a
spreadsheet and presence/absence plots. A particular advantage of AMOEBAE
compared to other workflows is its functionality for parsing results of TBLASTN
(which searches nucleotide sequences with peptide sequence queries) search
results. This allows rapid identification of High-scoring Segment Pair (HSP)
clusters at separate gene loci, automatic checking of those loci against
information in genome annotation files, and systematic use of the
[Exonerate](http://europepmc.org/article/MED/15713233) package where possible
for obtaining exon predictions. In addition, AMOEBAE provides many
    options which can be tailored to the specific genes/proteins being
    analyzed. Despite the complexity of this workflow, analyses performed using
    AMOEBAE can be reproduced via
    [SnakeMake](https://snakemake.readthedocs.io/en/stable/).

The output files include a plot of the number of identified
homologues (potential orthologues) of several genes across several genomes, as
well as a spreadsheet in CSV format providing a detailed summary of search
results. 

Here's a diagram of the steps in the overall workflow:

<p align="center">
<img src="images/example_workflow_diagram.png" width="500">
</p>


Here's an example coulson plot output by the workflow:

<p align="center">
<img src="images/example_coulson_plot.png" width="500">
</p>

## Protocol

See the [workflow
protocol](documentation/workflow_protocol.md) for instructions and
guidelines for running the AMOEBAE workflow.


## Questions and bug reporting

Please use the [issue tracker](https://github.com/laelbarlow/amoebae/issues) on
the GitHub webpage to report any problems you encounter while using AMOEBAE.


## How to cite AMOEBAE

Please cite the [AMOEBAE GitHub
repository](https://github.com/laelbarlow/amoebae) (or alternative permanent
repositories if relevant). 


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
Shweta V. Pipalya, Kira More, Kristína Záhonová, and Emily K. Herman.


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


