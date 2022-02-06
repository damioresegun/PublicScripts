#!/bin/bash
#folder where you want to delete something
folder="$HOME/All_Aligned"
# just change the file extension and pattern to be looking for
cd $folder
for i in ${folder}/*
do
fname=$(basename $i)
echo $fname
cd $fname
rm ${fname}VsErnest.bam
rm ${fname}VsErnest_sorted.fastq.gz
cd $folder
done
