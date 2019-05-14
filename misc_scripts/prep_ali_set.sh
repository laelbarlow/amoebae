#!/bin/bash

cedarcount=`ls -d *_cedar/ 2>/dev/null | wc -l`
reducednexcount=`ls -l *reduced*.nex 2>/dev/null | wc -l`

if [ $cedarcount != 0 ]
then 
    printf "Visualizing trees"
    printf ""
    for X in *cedar; do amoebae visualize_tree $X iqtree --root_taxon\
        Klebsormidium --highlight_paralogues; done &&\
    for X in *cedar; do amoebae auto_prune $X; done &&\
    for X in *cedar/*pdf; do open $X; done
elif [ $reducednexcount != 0 ]
then
    printf "Realigning reduced alignments"
    printf ""
    for X in *reduced*.nex; do realign_nex.py $X; done &&\
    for X in *realigned.nex; do mask_nex.py $X; done &&\
    for X in *mask.nex; do trim_nex.py $X; done &&\
    for X in *trim.nex; do codenames_nex.py $X; done
else
    printf "Aligning initial FASTA files"
    printf ""
    for X in *.fa; do align_fa.py $X; done &&\
    for X in *.nex; do mask_nex.py $X; done &&\
    for X in *mask.nex; do trim_nex.py $X; done &&\
    for X in *trim.nex; do codenames_nex.py $X; done
fi

