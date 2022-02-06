#!/bin/bash 
### script to carry out batch fastqc of different fastq files in recursive directories of the working directory
#set working directory
wdr="$HOME/FastQCheck"
cd $wdr
echo "First working directory"
pwd
# for every folder in the working directory, change into that folder
for fldr in $wdr/*;
do
cd $fldr
echo "you have changed to the barcode folder"
pwd
# for every fastq file in each child folder, do the fastqc and save in the fastqc folder in the child folder
for i in $fldr/*.fastq; do
echo "now stating the fastq files"
mkdir ${fldr}/fastqc/
ls $i
$HOME/Tools/FastQC/fastqc -t 32 -o ${fldr}/fastqc/ $i
cd $wdr
done
echo "All done"
done
