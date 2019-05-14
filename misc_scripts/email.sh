#!/bin/bash
##This is a script for sending an email when a job finishes. 
##These commands set up the Grid Environment for your job:
## Job finished.
#PBS -N Job finished.
##You can specify the node but is best to leave this to the Queue manager to decide
##using -l nodes=1,and the Walltime(maximum executing time in HH:MM:SS)
#PBS -l nodes=1,walltime=1:00:00
#PBS -q default
#PBS -M lael.barlow@gmail.com
#PBS -m abe
##print Node Hostname the time and date handy for trouble shooting
echo $HOSTNAME
date
