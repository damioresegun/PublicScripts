#!/bin/bash
# Script to take fastq files through assembly stats and then take the text files for be changed to csv for opening in excel. The yext files will be deleted after to save space.

input="$HOME/All_AdapterRemoved"
output="$HOME/All_Stats/postAdapterRemoved_isolates"
#load necessary modules
module load assembly-stats/gitv0_003a372
#state the current working directory
hre=$(pwd)
echo $hre
mkdir -p $output
for f in ${input}/*
do
ls $f
#get the filename
file=$(basename $f)
echo $file
file=${file%.*}
echo $file
#cd $f
#assembly-stats it
assembly-stats *.fastq >> $output/$file.txt
echo "assembly stats done"
# convert to csv
cat $output/$file.txt | tr -s '[:blank:]' ',' > $output/$file.csv
#rm assembly-stats output
rm $output/$file.txt
echo "${file} done"
cd $input
done
exit

