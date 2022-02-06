#!/bin/bash

# script to take a whole lot of assemblies through for sscaffold stats

perl ~/Scripts/scaffold_stats.pl -f All_De_Novo/Guppy/sks048/sks048.contigs.fasta All_De_Novo/Guppy/sks048/sks048_clean.contigs.fasta All_De_Novo/Qcat/sks048/sks048.contigs.fasta All_De_Novo/Qcat/sks048/sks048_clean.contigs.fasta All_De_Novo/MaxGuppy/sks048/sks048.contigs.fasta All_De_Novo/MaxGuppy/sks048/sks048_clean.contigs.fasta All_De_Novo/MaxQcat/sks048/sks048.contigs.fasta All_De_Novo/MaxQcat/sks048/sks048_clean.contigs.fasta -t 1000 1000000 -c 10000 > Assembly_QC/ScaffStats/sks048_all.txt
