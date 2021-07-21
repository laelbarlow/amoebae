

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

Unique feature: Allows/requires manual selection of valid reverse search hits
for interpretation of results.

Here's a diagram of the steps in the overall workflow:

<p align="center">
<img src="images/example_workflow_diagram.png" width="500">
</p>

Here's an example coulson plot output by the workflow.

<p align="center">
<img src="images/example_coulson_plot.png" width="500">
</p>


## Requirements 

The following setup procedure should work on most Linux HPC clusters. This can
also be run on Linux or MacOS personal computers, but this is generally not
recommended due to requirements of storage (~30GB or more) and computation time. 

This workflow has minimal essential dependencies for installation, which are
all widely used in the life sciences (see below for installation instructions):
- Python version 3.
- [Conda](https://docs.conda.io/en/latest/miniconda.html) **or**
  [Singularity](https://sylabs.io/docs/)
  version 3.6+.


## Installation

These instructions are for setting up and running AMOEBAE via the
[SnakeMake](https://snakemake.readthedocs.io/en/stable/) command-line
interface, which is well-documented and provides the flexibility to run AMOEBAE
in a wide variety of systems (further automations via custom scripts and
Makefiles is possible, but will likely be very user-specific in design). 

1. If you do not already have the conda package and environment manager
   installed, this may be installed using the latest version of Miniconda3 from
   the [conda website](https://docs.conda.io/en/latest/miniconda.html). You may
   need help from a system administrator to set this up properly on a HPC cluster.
   Some HPC cluster administrators do not allow use of conda at all, but in
   such cases conda can still be used to install software in
   [Singularity](https://sylabs.io/docs/) containers (this process is automated
   by snakemake). If singularity is not already installed on your system, you
   will need to request that you system administrator install it for you. Also,
   if you need to use singularity, then you will need to add
   `--use-singularity` to most of the snakemake commands described below.

2. Install snakemake and a few other packages. If you do not already have the
   snakemake workflow manager installed then you will need to install it. If
   you have conda installed, then you can follow the [snakemake installation
   instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
   on the snakemake website. If conda cannot be installed, then snakemake can
   be installed in a Python virtual environment. You will need a few additional
   dependencies: [cookiecutter](https://github.com/cookiecutter/cookiecutter)
   and [graphviz](https://graphviz.org/).

   These can be installed at the same time as snakemake. For example, after
   installing [mamba](https://github.com/mamba-org/mamba) using conda, you can
   run the following command to create an environment for running snakemake: 
   ```
   mamba create -c conda-forge -c bioconda \
       -n snakemake \
           snakemake \
           cookiecutter \
           graphviz
    ```

   An alternative method, if you cannot use conda on your system, is to
   install snakemake in a [Python virtual
   environment](https://docs.python.org/3/library/venv.html).
   ```
   python3 -m venv ~/amoebae_python_env
   source ~/amoebae_python_env/bin/activate
   pip install --upgrade pip
   pip install \
       snakemake \
       cookiecutter \
       graphviz
   ```


3. If running AMOEBAE on an HPC cluster (recommended) then you will need to generate
   cluster configuration files so that snakemake knows how to submit jobs
   appropriately on your system. This is described in the [snakemake
   documentation on
   profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).
   To set up a profile, select the appropriate template profile for your system
   from the [snakemake profiles](https://github.com/snakemake-profiles/doc),
   and follow the relevant setup instructions provided. This will write new
   files to a new hidden directory in your filesystem. For example, if your HPC
   cluster uses PBS-TORQUE, and you choose to name the new snakemake
   profile "pbs-torque", then cluster configuration files will be written to a
   directory with the path `~/.config/snakemake/pbs-torque`. For example:
   ```
   conda activate snakemake
   mkdir -p ~/.config/snakemake
   cd ~/.config/snakemake
   cookiecutter https://github.com/Snakemake-Profiles/pbs-torque.git
   chmod a+x pbs-torque/*.py
   cd -
   ```
   To edit the cluster configuration in the snakemake profile (if necessary),
   edit the `cluster_config.yaml` or `cluster.yaml` file (depends on what type
   of profile, and what you chose to name the files) using your favourite text
   editor. Details of the configuration will depend on the job scheduler and
   resources available on your system, and can be modified at any time. Your
   system administrator will have provided job submission guidelines and
   example job submission commands/scripts, and these may help guide you in
   adapting the cluster configuration files to your system. To use a snakemake
   profile, you will need to append an option to all snakemake commands to
   indicate the name of the profile (for example, `--profile pbs-torque`). 

4. If you are **not** running AMOEBAE on an HPC cluster, for example if you are
   running it on a personal computer, then snakemake profiles will not be
   necessary. In this case, simply activate your previously defined conda
   environment for snakemake, and you will ready to run amoebae or other
   snakemake workflows:
   ```
   conda activate snakemake
   ```

5. Finally, clone the AMOEBAE code repository into an appropriate directory.
    ```
    git clone --branch dev2 https://github.com/laelbarlow/amoebae.git
    ```

## Running the workflow

With example (default) input files, this workflow should take between 30 and 60
minutes to run, depending on resource availability on your system. If running
on an HPC cluster, it may be useful to use
[tmux](https://github.com/tmux/tmux/wiki) or
[nohup](https://www.gnu.org/software/coreutils/manual/html_node/nohup-invocation.html)
to prevent snakemake processes from being interrupted. 

1. Collect genome/proteome/transcriptome FASTA files to be searched:
    - Copy the example genomes.csv file in the config subdirectory.
        ```
        cd amoebae
        cp config/example_genomes.csv config/genomes.csv
        ```
    - If you want to use example genome files, leave this copy as it
      is, and proceed to the next step.
    - Otherwise, modify the `config/genomes.csv` file by adding information
      about each predicted peptide FASTA files (.faa), nucleotide FASTA files
      (.fna), and/or GFF3 annotation files (.gff3) of interest.
        - In the "FASTA header delimiter" column, enter the text character that
          separates sequence IDs from other elements of the FASTA header. In the
          case of FASTA files from NCBI for example, this is usually a space
          character.
        - In the "Sequence ID position" column, enter the position of the sequence
          ID in a list resulting from splitting the whole FASTA header on the
          character defined in the "FASTA header delimiter" column. Importantly,
          counting starts from zero, so if the sequence header starts with the
          sequence ID (as in the case of most NCBI FASTA files), then the value in
          this column should be "0".
        - If the FASTA files are to be downloaded from a website, enter the URL in
          the "Location" column. Otherwise, it will be assumed that the files have
          been copied to the `resources/local_db_files` directory.
        - Note: If you use a spreadsheet editor such as Excel, then make sure to save the modified .csv file with UTF-8 encoding (plain text).
        - If you wish to search in any local FASTA files (instead of downloading
          directly from [NCBI](https://www.ncbi.nlm.nih.gov/)), copy those to the
          `resources/local_db_files` directory (in addition to listing them in the
          genomes.csv file).

2. Collect query sequence FASTA files:
    - Copy the example queries.csv file in the config subdirectory.
        ```
        cp config/example_queries.csv config/queries.csv
        ```
    - If you want to use example query files, leave this copy as it is, and
      proceed to the next step.
    - Otherwise, modify the `config/queries.csv` file by adding information
      about each query sequence.
        - There are only two columns: "Filename" and "Sequence ID". Any given
          filename may appear on multiple rows, but each row must contain a
          unique sequence ID, *which must be a valid NCBI accession number*, or a blank
          space in the "Sequence ID" column.
        - If a single sequence ID is associated with a given filename, then a
          single-sequence FASTA file will be generated by downloading the
          corresponding sequence from NCBI.
        - If multiple sequence IDs are listed for a given filename, then a
          multiple-sequence FASTA file will be generated by downloading each of
          the corresponding sequences from NCBI and writing them to a FASTA
          file with the specified name.
        - If no sequence IDs are provided for a given filename, then the filename will
          be assumed to correspond to an existing local file in the
          `resources/local_query_files` directory. 
        - So, if you wish to search using local query (FASTA) files, copy them
          to the `resources/local_query_files` directory (in addition to
          listing their filenames in the queries.csv file). 
        - Importantly, query filenames must start with a query title which is
          followed by the first underscore. This query title is just the name
          of the proteins you are searching for. For example, "AP1beta" for the
          adaptor protein complex 1 beta subunit. Multiple query files can have
          the same query title, allowing you to search for the same proteins
          with multiple queries.

3. Specify an appropriate reference sequence file, to query in
   reverse searches (this is the second set of searches in reciprocal-best-hit
   sequence similarity searching). 
    - Copy the example file in the config subdirectory.
        ```
        cp config/example_reference_db_list.txt \
           config/reference_db_list.txt
        ```
    - Again, if you simply want to use the example file (*Arabidopsis thaliana*
      amino acid sequences), then proceed to the next step.
    - To use different files, modify the `config/reference_db_list.txt`
      file, adding/removing filenames (one filename per line). Well-annotated
      reference genomes, such as those of *Arabidopsis thaliana*,
      *Saccharomyces cerevisiae*, or *Homo sapiens* may be most useful for this
      purpose.

4. Execute initial workflow steps to download (if necessary) and format
   sequence data and generate lists of potential reference orthologues (for
   interpreting reverse searches). 
    ```
    snakemake get_ref_seqs -j 100 --use-conda --profile pbs-torque
    ```

2. Select relevant reference sequences for interpreting reverse search results.
    - If you ran searches just with the example files, then you can use the
      example reference sequence selection file by copying it in the `config`
      subdirectory.
        ```
        cp config/example_Ref_seqs_1_manual_selections.csv \
           config/Ref_seqs_1_manual_selections.csv
        ```
    - Otherwise, copy the relevant file from the results subdirectory.
        ```
        cp results/Ref_seqs_1_auto_selections.csv \
           config/Ref_seqs_1_manual_selections.csv
        ```
    - The purpose of this `config/Ref_seqs_1_manual_selections.csv` file is to
      identify all sequences in the reference genome which are expected to be
      retrieved as top hits by sequences of interest from other genomes. Each
      row in this file corresponds to a sequence in the reference genome
      retrieved with one of the queries in a BLASTP sequence similarity search.
      Each reference sequence is identified as either accepted or 
      unaccepted as a top hit in reverse searches.
    - If you are using example data, proceed to the next step. Otherwise,
      modify values in the 5th column of the
      `config/Ref_seqs_1_manual_selections.csv` file. In this column, '+'
      indicates inclusion of a sequence as a representative of the query used,
      and '-' indicates that the sequence is too distantly related to the query
      to be relevant.

3. Configure the organization of output plots.
    - First, copy the relavent example configuration files.
        ```
        cp config/example_coulson_plot_organization.csv \
           config/coulson_plot_organization.csv

        cp config/example_output_plot_row_order.txt \
           config/output_plot_row_order.txt
        ```
    - Then re-arrange or re-write as necessary. The
      `coulson_plot_organization.csv` file contains two columns: The first is
      the name assigned to a column in the output coulson plot, and the second
      lists query titles for queries assigned to each column (multiple query
      titles can be associated with the same column name). The
      `output_plot_row_order.txt` file simply lists the file names of the FASTA
      files searched (including any nucleotide FASTA files), in the order that
      you want them to appear in the plots.

4. Execute remainder of the workflow.
    ```
    snakemake -j 100 --use-conda --profile pbs-torque
    ```

5. View results summary files and positive hit sequence alignments at these
   paths. For convenience, you may wish to download these files using an SFTP
   client such as [Cyberduck](https://cyberduck.io/download/) or
   [FileZilla](https://filezilla-project.org/):
    ```
    results/plot_coulson_both.pdf
    results/fwd_srchs_1_rev_srch_1_interp_with_ali_col_nonredun.csv
    results/fwd_srchs_1_rev_srch_1_interp_with_ali_col_nonredun_fasta_ali_files
    ```


Note that these results require careful interpretation, and in most cases
re-analysis with modified parameters will be necessary as well as follow up
with additional methods such as phylogenetic analysis.

6. Generate a report of results in HTML format which can be opened in a web browser.
    ```
    snakemake --cores 1 --report results/amoebae_report.html
    ```


## Customizing parameters

For more advanced customization, refer to the Snakemake documentation and
modify the Snakefile (snakemake workflow definition file): 
    ```
    vim workflow/Snakefile    
    ```
For customizing shell commands in workflow rules (steps) that run the `amoebae`
script, refer to the [AMOEBAE command-line interface
documentation](https://github.com/laelbarlow/amoebae/blob/master/documentation/AMOEBAE_documentation.pdf)
for available options. 


## Questions and bug reporting

Please use the [issue tracker](https://github.com/laelbarlow/amoebae/issues) on
the GitHub webpage to report any problems you encounter while using AMOEBAE.


## How to cite AMOEBAE

Please cite the [AMOEBAE GitHub
repository](https://github.com/laelbarlow/amoebae) (or alternative permanent
repositories if relevant). Also, please remember to cite the software packages which
are dependencies of AMOEBAE.


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


