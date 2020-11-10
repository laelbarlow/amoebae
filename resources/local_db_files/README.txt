
Purpose of this directory:

    To store local genome/proteome/transcriptome (FASTA) and annotation (GFF3)
    files for input to the workflow, which are not easily downloaded from
    online databases. 

Instructions:

    Copy relevant FASTA files (and any relevant GFF3 files) to this directory. 

    Formatting of filenames is important, as information from the filenames
    will be used to identify results of searches in the specific files and to
    identify which GFF3 files correspond to which nucleotide FASTA files.
    Filenames should be simply the name of the organism that was sequenced with
    underscores instead of spaces. For example, for the yeast Saccharomyces
    cerevisiae, you may wish to simply provide the binomial species name as
    follows, but strain names could be included as well.
        
        Saccharomyces_cerevisiae.faa
        Saccharomyces_cerevisiae.fna
        Saccharomyces_cerevisiae.gff3

    Upon running the workflow, these will be detected and added to the data in
    the results directory in a usable format.

    Remember to clearly identify the sources of data files in any publications.
