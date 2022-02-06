#!/bin/bash

input=$HOME/FinalGenomes
output=$HOME/FinalGenomes
#mkdir -p $output
for f in $input/*
do
fna=$(basename "$f")
#fname=${fna%.*}
#mkdir -p $output/${fna}
cd $output/${fna}
for f in $output/${fna}/*
do
nme=$(basename "$f")
numu=${nme%.*}
if [[ $numu != *_00 ]];
then
cat $input/${fna}/${nme} >> $output/${fna}/Full_${fna}.fasta
else
echo "Something has gone wrong"
fi
 done
done
