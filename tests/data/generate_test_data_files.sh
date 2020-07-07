
# Script for generating input files for tests.

# Generate BLASTP output files.

# Make FASTA file blastable.
makeblastdb -in blastp_database_1.faa -dbtype prot

# Run blastp search with text output.
blastp -query blastp_query_1.faa -db blastp_database_1.faa -out blastp_query_1_database_1_output.txt
# Run blastp search with XML output.
blastp -query blastp_query_1.faa -db blastp_database_1.faa -out blastp_query_1_database_1_output.xml -outfmt 5
