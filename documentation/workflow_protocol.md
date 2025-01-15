
# AMOEBAE workflow protocol

## Requirements 

The following setup procedure should work on most Linux or MacOS
systems. Please ensure that you have sufficient storage, as this protocol
will generate files totalling ~30GB or more in size.


## Installation

### Local servers and personal computers

These instructions are for setting up and running AMOEBAE via the
[SnakeMake](https://snakemake.readthedocs.io/en/stable/) command-line
interface, which is well-documented and provides the flexibility to run AMOEBAE
in a wide variety of systems. 

1. Install Pixi according to the instructions on the [Pixi
website](https://pixi.sh/latest/). Pixi will install Conda packages without the
need for you to install Conda itself (or Mamba, etc.). This approach is new as
of November 2024, and addresses the problems including those related to
installing dependencies on Apple Silicon MacOS systems.

2. Clone the AMOEBAE code repository into an appropriate directory.
    ```
    git clone https://github.com/laelbarlow/amoebae.git my_amoebae_directory
    ```
    - Note: Replace `my_amoebae_directory` with the name of the directory where
      you want to clone AMOEBAE.

3. Navigate into the new directory, and install dependencies listed in the new `pixi.toml` file:
    ```
    cd my_amoebae_directory
    pixi install
    ```

4. Install `blast`, `muscle`, and `exonerate` globally using Pixi, because
    this is operating-system specific.

    First, make sure you don't already have blast, muscle, or exonerate globally installed:
    ```
    which blastp
    which muscle
    which exonerate
    ```
    If you do, then you should may need to uninstall these before proceeding.
    
    Then, if you are not installing on an Apple Silicon MacOS system, then
    simply do this:
    ```
    pixi global install -c conda-forge -c bioconda blast muscle=3.8 exonerate
    ```

    If you are using Apple Silicon MacOS, then install Rosetta 2 (if you haven't already):
    ```
    /usr/sbin/softwareupdate --install-rosetta --agree-to-license
    ```

    Then use the `--platform` option to install the osx64 binaries like this:
    ```
    pixi global install --platform osx-64 -c conda-forge -c bioconda blast muscle=3.8 exonerate
    ```

5. **Note**: In contrast to previous versions of this protocol, now with Pixi we
    prefix all snakemake commands with `pixi run` so that the workflow is run in
    the Conda environment created using Pixi, and no longer use the
    `--use-conda` option. For example, `pixi run snakemake`.

### High Performance Computing Clusters

See the [instructions for using AMOEBAE on Digital Research Alliance of Canada clusters](./drac.md).

## Running the workflow

With example (default) input files, this workflow should take between 30 and 60
minutes to run, depending on resource availability on your system. If running
on an HPC cluster, it may be useful to use
[tmux](https://github.com/tmux/tmux/wiki) or
[nohup](https://www.gnu.org/software/coreutils/manual/html_node/nohup-invocation.html)
to prevent snakemake processes from being interrupted. 

1. Collect genome/proteome/transcriptome files to be searched:
    - Note: This includes reference genomes that will be used for reverse
      searches.
    - Copy the [example genomes.csv
      file](https://github.com/laelbarlow/amoebae/blob/126fd9dfaff6b21165cda0631bd77d994b0a69aa/config/example_genomes.csv)
      in the config subdirectory.
        ```
        cd amoebae
        cp config/example_genomes.csv config/genomes.csv
        ```
    - If you want to use example genome files, leave this `config/genomes.csv`
      file as it is, and proceed to the next step.
    - Otherwise, define the files you wish to search in.  
        - Input FASTA files for searching may be predicted peptide FASTA files
          (.faa) and/or nucleotide FASTA files (.fna). For any genomic nucleotide FASTA
          files used, you may also include an associated [General Feature
          Format Version 3
          (GFF3)](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
          annotation file (.gff3), which defines where genes are located in the
          genomic nucleotide sequences. These GFF3 files are usually provided
          with genomic data. 
        - Files that are publicly accessible online at specific URLs (not
          password protected) can be listed in the `config/genomes.csv` file
          with the relevant URLs, and they will be downloaded automatically
          when you run the workflow (see below).
        - Files that are not accessible online, such as those for novel
          sequence data, must be copied to the `resources/local_db_files`
          directory. However, these will not be searched unless they are also
          listed in the `config/genomes.csv` file. 
        - To define a file for searching in the `config/genomes.csv` file, open
          the genomes file in a spreadsheet editor or plain text editor, and enter
          the following information (each row contains information about a
          single file):
            - In the "Filename" column, enter the filename.
                - Formatting of filenames is important, as information from the
                  filenames will be used to identify results of searches in the
                  specific files and to identify which FASTA and/or GFF3 files
                  correspond to the same genome.
                - Filenames should be simply the name of the organism that was
                  sequenced with underscores instead of spaces. For example,
                  for the yeast **Saccharomyces cerevisiae**, you may wish to
                      simply provide the binomial species name, but
                      strain names can be included as well.
                  ```
                  Saccharomyces_cerevisiae_S288C.faa
                  Saccharomyces_cerevisiae_S288C.fna
                  Saccharomyces_cerevisiae_S288C.gff3
                  ```
                - All filenames listed in this column must end with one of the
                  following filename extensions as appropriate:`.faa`, `.fna`,
                  or `.gff3`. 
                - Files that have differing filenames (minus the filename
                  extension) will be assumed to be for different genomes. 
                - Filenames listed here are not the file paths. Files that are
                  downloaded by the workflow will be saved in an appropriate
                  results subdirectory with the specified filename, and local
                  files will be identified in the `resources/local_db_files`
                  directory based on the filename specified here.
            - In the "FASTA header delimiter" column, enter the text character that
              separates sequence IDs from other elements of the FASTA header. In the
              case of FASTA files from NCBI for example, this is usually a space
              character.
            - In the "Sequence ID position" column, enter the position of the sequence
              ID in a list resulting from splitting the whole FASTA header on the
              character defined in the "FASTA header delimiter" column. Importantly,
              counting starts from zero, so if the sequence header starts with the
              sequence ID (as in the case of most NCBI FASTA files), then the value in
              this column should be "0". The purpose of this is for extracting
              appropriate sequence IDs which can be listed in your result
              tables.
            - If the file is to be downloaded from a website: Enter the file
              compression type in the "Compression type" column.  Currently,
              acceptable compression types are gzip or none (if the file for
              download is not compressed then just leave this cell empty).
              Also, enter the URL for the file in the "Location" column (see
              example URLs for files on NCBI).
            - If the file is a local file, then leave the cells in the
              "Compression type" and "Location" columns empty for this row. See
              the `config/example_local_genomes.csv` file for an example (each
              row must contain the same number of commas for this to be
              interpreted properly as CSV format).
            - Note: If you use a spreadsheet editor such as Microsoft Excel, then
              make sure to save the modified file with UTF-8 encoding (plain
              text) and in CSV format with the filename extension `.csv`.
        

2. Collect query sequence FASTA files:
    - Copy the [example queries.csv
      file](https://github.com/laelbarlow/amoebae/blob/126fd9dfaff6b21165cda0631bd77d994b0a69aa/config/example_queries.csv)
      in the config subdirectory.
        ```
        cp config/example_queries.csv config/queries.csv
        ```
    - If you want to use example query files, leave this copy as it is, and
      proceed to the next step.
    - Otherwise, define the queries you want to search with. 
        - This is similar to how genome files are specified (see above).
        - For simplicity, it is best to select query sequences from the genomes that you plan
          to use as reference genomes (see below).
        - Input query files must contain peptide (amino acid) sequences, and
          may be in either single-FASTA (for BLAST searches) or aligned
          multi-FASTA format (for profile searches with HMMer). These must be
          in valid FASTA format, which includes a sequence header prefixed by
          `>`.
        - Sequences that are accessible online via the [NCBI Protein
          database](https://www.ncbi.nlm.nih.gov/protein/) can be listed in the
          `config/queries.csv` file with the relevant [NCBI protein
          accessions](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#AccessionB),
          and they will be downloaded automatically when you run the workflow
          (see below).
        - Files containing sequences that are not easily downloadable must be
          copied to the `resources/local_query_files` directory. However, these
          will not be searched unless they are also listed in the
          `config/queries.csv` file. 
        - To define a query sequence/file for searching in the
          `config/queries.csv` file, open the queries file in a spreadsheet
          editor or plain text editor, and enter the following information
          (each row contains information about a single query sequence or an
          alignment):
            - In the "Filename" column, enter the filename. 
                - Formatting of filenames is important, as information from the filenames
                  will be used to identify results of searches with these specific query
                  files.  

                - Importantly, query filenames must start with a query title
                  which is followed by the first underscore. This query title
                  is just the name of the proteins you are searching for. For
                  example, "AP1beta" for the adaptor protein complex 1 beta
                  subunit. 
                - Multiple query files can have the same query title, allowing
                  you to search for the same proteins with multiple queries
                  (without duplication in the final results).
                - For example, a single-FASTA file containing a peptide sequence of an
                  Arabidopsis orthologue of the protein Rab2 could be used to search with
                  BLASTP, while a multi-FASTA alignment file with Rab2 orthologue sequences
                  from Arabidopsis and other embryophytes could be used for profile
                  searches for Rab2 orthologues using HMMer3. These files could be named as
                  follows:
                  ```
                  Rab2_Athaliana_NP_193449.1.faa,
                  Rab2_Embryophyta.afaa,
                  ```
                - All filenames listed in this column must end with one of the
                  following filename extensions as appropriate:`.faa` for
                  single-FASTA peptide sequences, and `.afaa` for multi-FASTA
                  peptide alignments. 
                - Only the filename is necessary, not a file path.
            - In the "Sequence ID" column, enter the NCBI protein sequence
              accession when relevant.
                - If you wish to have the query sequence downloaded
                  automatically from NCBI, then include a sequence ID here.
                  These *must be a valid NCBI protein accessions*. 
                - If instead the sequence or sequences are in a local file in
                  the `resources/local_query_files` directory, then simply
                  leave this cell empty.
        - Note: It is also possible to list a filename (with the extension
          `.faa`) in multiple rows, each with different NCBI accessions (as in
          the example `config/queries.csv` file). In this case all the relevant
          sequences will be downloaded and written to a single file and aligned
          to make an HMM for profile searching with HMMer. However, this is
          unlikely to be the best way to construct alignments in most cases.



3. Specify appropriate reference genomes to query in
   reverse searches.
   	- Protein FASTA file(s) (with `.faa` extensions) selected here will be used in the second set of searches 
   	  performed during reciprocal-best-hit sequence similarity searching. 
    - Copy the example file in the config subdirectory.
        ```
        cp config/example_reference_db_list.txt \
           config/reference_db_list.txt
        ```
    - Again, if you simply want to use the example file (*Rhizophagus
      irregularis* amino acid sequences), then proceed to the next step.
    - To use different files, modify the `config/reference_db_list.txt`
      file accordingly.
        - List filenames with `.faa` extensions, one filename per line.
        - Any filenames listed here must chosen from among those listed in your 
          `config/genomes.csv` file (see above).  
        - Note: Well-annotated reference genomes, such as those of *Arabidopsis
          thaliana*, *Saccharomyces cerevisiae*, or *Homo sapiens* may be most
          useful for this purpose. These are available from NCBI (see the
          example `config/genomes.csv` file).


4. Select relevant reference sequences for interpreting reverse search results.
    - If you have already run amoebae from this directory, and
      have some previous result files in the `amoebae/results` subdirectory,
      then you will need to remove your `amoebae/results` directory (or move it
      elsewhere) before running the workflow again with snakemake commands
      described below in this step and subsequent steps. 
    - Perform sequence similarity searches in your reference genomes using your
      queries. Here, several initial workflow steps will be run to download (if
      necessary) and format sequence data. 
      ```
      pixi run snakemake get_ref_seqs 
      ```
    - If you ran searches just with the example files, then you can use the
      example reference sequence selection file by copying it in the `config`
      subdirectory (this defines appropriate *Rhizophagus* sequences as
      orthologs of the example queries).
        ```
        cp config/example_Ref_seqs_1_manual_selections.csv \
           config/Ref_seqs_1_manual_selections.csv
        ```
    - Otherwise, copy the relevant file from the results subdirectory.
        ```
        cp results/Ref_seqs_1_auto_predictions.csv \
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
    - Notes:
        - Occasionally, temporary disruptions of NCBI servers will prevent
          download of files (queries or genomes). If the workflow fails to
          download FASTA files, check the FTP download URL or NCBI accession,
          as applicable, and try running the workflow again. 
        - If your `results/Ref_seqs_1_auto_predictions.csv` file is missing
          search results that you expected it to contain, then you may be using
          the wrong protein FASTA file for one or more reference genomes. Check
          that you have specified the relevant URLs in the `config/genomes.csv`
          file and/or correctly named the relevant FASTA files in your
          `resources/local_db_files` directory.
        - AMOEBAE re-ranks all hmmsearch (HMMer) hits according to the best 1
          domain E-value instead of the full sequence E-value, which can be
          useful when searching among sequences with multiple homologous
          domains. See the [HMMer user
          guide](http://eddylab.org/software/hmmer/Userguide.pdf) for further
          information.
        - If you do another analysis with the same reference genome and query set,
          then you can re-use the same `config/Ref_seqs_1_manual_selections.csv`
          file.


5. Configure the organization of output plots.
    - First, copy the relevent configuration files from the results directory.
        ```
        cp results/coulson_plot_organization.csv \
           config/coulson_plot_organization.csv

        cp results/db_list.txt \
           config/output_plot_row_order.txt
        ```
    - Then re-arrange or re-write these configuration files as necessary. The
      `coulson_plot_organization.csv` file contains two columns: The first is
      the name assigned to a column in the output coulson plot, and the second
      lists query titles for queries assigned to each column (multiple query
      titles can be associated with the same column name). The
      `output_plot_row_order.txt` file simply lists the file names of the FASTA
      files searched (including any nucleotide FASTA files), in the order that
      you want them to appear in the plots. For examples, see
      `config/example_coulson_plot_organization.csv` and
      `config/example_output_plot_row_order.txt`.


6. Execute remainder of the workflow to perform searches in your genomes of
   interest.
    - Note: This is where searches in all the FASTA files listed in your
      `config/genomes.csv` file are performed.
    - Run the following command:
      ```
      pixi run snakemake 
      ```
    - This will take minutes to hours to run, depending on the input data,
      computing resources, etc.
    - View results summary files and positive hit sequence alignments at these
      paths. 
      ```
      results/plot_coulson_both.pdf
      results/fwd_srchs_1_rev_srch_1_interp_with_ali_col_nonredun.csv
      results/fwd_srchs_1_rev_srch_1_interp_with_ali_col_nonredun_fasta_ali_files
      ```
    - You may also wish to generate a report of results in HTML format which
      can be opened in a web browser. This can be done with the following
      command:
      ```
      pixi run snakemake --report results/amoebae_report.html
      ```
    - If you are running the workflow on an HPC cluster, for convenience you may
      wish to download these files using an SFTP client such as
      [Cyberduck](https://cyberduck.io/download/) or
      [FileZilla](https://filezilla-project.org/).
    - Note that these results require careful interpretation, and in most cases
      re-analysis with modified parameters will be necessary as well as follow
      up with additional methods such as phylogenetic analysis.


## Customizing parameters

For more advanced customization, refer to the Snakemake documentation and
modify the Snakefile (snakemake workflow definition file): 
    ```
    vim workflow/Snakefile    
    ```
For customizing shell commands in workflow rules (steps) that run the `amoebae`
script, refer to the [AMOEBAE command-line interface
documentation](amoebae_commands.pdf) for available options. 


