#/usr/bin/env bash

echo "\\subsection{amoebae}" > amoebae_help_output.tex
echo "\\begin{lstlisting}" >> amoebae_help_output.tex
amoebae -h >> amoebae_help_output.tex
echo "\\end{lstlisting}" >> amoebae_help_output.tex

CMDLIST="mkdatadir
         setup_hmmdb
         add_to_dbs
         list_dbs
         add_to_queries    
         list_queries      
         get_redun_hits    
         setup_fwd_srch    
         run_fwd_srch      
         sum_fwd_srch      
         setup_rev_srch    
         run_rev_srch      
         sum_rev_srch      
         interp_srchs      
         find_redun_seqs   
         plot              
         csv_to_fasta      
         check_depend      
         check_imports
         regen_genome_info"     


for CMD in $CMDLIST;
do 
    echo "" >> amoebae_help_output.tex
    echo "\\subsection{amoebae ${CMD//_/$'\_'}}" >> amoebae_help_output.tex
    echo "\\begin{lstlisting}" >> amoebae_help_output.tex
    amoebae $CMD -h >> amoebae_help_output.tex
    echo "\\end{lstlisting}" >> amoebae_help_output.tex
done
