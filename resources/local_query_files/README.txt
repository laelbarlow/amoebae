
Purpose of this directory:

    To store local query files in single- or multi-FASTA format. This is useful
    in the case where sequences are not easily retrieved from online databases
    (and aligned) programmatically.

Instructions:

    Copy relevant FASTA files to this directory. 

    Formatting of filenames is important, as information from the filenames
    will be used to identify results of searches with these specific query
    files.  Filenames can have one or more fields (separated by underscores),
    with the first and most important field being the query title. Multiple
    query files may have the same query title, thereby helping to organize
    searching for orthologues or homologues of a single gene/protein with
    multiple queries (and multiple search methods) without duplication of final
    results. For example, a single-FASTA file containing a peptide sequence of
    an Arabidopsis orthologue of the protein Rab2 could be used to search with
    BLASTP, while a multi-FASTA alignment file with Rab2 orthologue sequences
    from Arabidopsis and other embryophytes could be used for profile searches
    for Rab2 orthologues using HMMer3. These files could be named as follows:

        Rab2_Athaliana_NP_193449.1.faa

        Rab2_Embryophyta.afaa

    Upon running the workflow, these will be detected and added to the data in
    the results directory in a usable format.

    Remember to clearly identify the sources of sequences/query files in
    publications.
