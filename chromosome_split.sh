#!/bin/bash
# echo "Your current directory is: $(pwd)"
# cwd=$(pwd)
# echo "Your current directory is set as your working directory"
# mkdir -p chromosomes
# cd chromosomes
# read -p "What is the input file to split. Input the filename.extension... `echo $'\n> '`" file
# input=${cwd}/${file}
# echo "Your input file is ${input}"
# #split the chromosome based on the > pattern seen on the first line
# csplit -f chr -s -z ${input} '/>/' '{*}'
# #rename the new files with the chromosome names and save as fasta files
# for i in chr* ; do \
  # head -n 1 ${i}
  # n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  # mv ${i} "${i}_{${n}}.fa" ; \
 # done

input=$1 #$HOME/All_CompanionGenomes/VsPk/Canu
output=$2 #$HOME/quickCanu
mkdir -p $output
for f in $input/*
do
fna=$(basename "$f")
fname=${fna%.*}
mkdir -p $output/${fname}
cd $output/${fname}
csplit -f chr -s -z ${f} '/>/' '{*}'

#rename the new files with the chromosome names and save as fasta files
for i in chr* ; do
  head -n 1 ${i}
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ;
  #mv ${i} "${i}_{${n}}.fa" ;
  mv ${i} "${n}.fa" ;
 done
 done

