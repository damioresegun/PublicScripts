#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N PkCanu_stats
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk
#Script to take through files through assembly stats. It can be infinitely expanded to take in different fasta files of different types. It is currently made for aseembly outputs of Canu and Flye. 
## needs the input folder to be the folder holding the folders for the different isolates. 
###CHANGE FILE NAME AND IF IT FLYE OR CANU!!!!!!
#input files
INP="$HOME/All_De_Novo/VsPk/Canu"
OUTP="$HOME/All_Stats/postVsPkAssembly/Canu/"
mkdir -p $OUTP
#load needed modules
module load assembly-stats/gitv0_003a372
## canu assembly stats -- comment/uncomment block where necessary
for i in INP/*
do
echo "The folder you are working on is ${i}"
    fname=$(basename "$i")    
    echo $fname
	assembly-stats -t $INP/${fname}/${fname}.contigs.fasta >> $OUTP/VsPkassembly.tsv
done
#convert canu assem to csv
#cat $OUTP/VsPkassembly.txt | tr -s '[:blank:]' ',' > $OUTP/VsPkassembly.csv
## fly assembly stats -- comment/uncomment block where necessary
# for i in $INP/*
# do
# ls $i
# echo "The folder you are working on is $i"
    # fname=$(basename "$i")    
    # echo $fname
	# assembly-stats $INP/${fname}/assembly.fasta >> $OUTP/flyeVsPk_assembly.txt
# done
# #convert flye assem to csv
# cat $OUTP/flyeVsPk_assembly.txt | tr -s '[:blank:]' ',' > $OUTP/flyeVsPk_assembly.csv

echo "All done"
