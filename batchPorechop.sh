#!/bin/bash
#Script to take fastq files through porechop based on the names they already have. Hence these may need renaming manually after this. 

input="$HOME/GPU_Demultiplexed"
output="$HOME/All_GPU_AdapterRemoved"
THREADS="16"

set -e
#load necessary modules
module load python/3.6.4
for folder in $input/*
do
	#ls $folder
	echo "It begins"
	echo $folder
    echo "The experiment you are working on is $folder"
    fname=$(basename "$folder")    
    #echo $fname
    echo
	mkdir -p $output/$fname
	for file in $folder/*
	do
	ls $file
    echo "Beginning adapter removal for $file..."
	fasa=$(basename "$file")
    echo $fasa
	fasna=${fasa%.*}
	echo $fasna
	porechop -i $file -o $output/$fname/$fasa --verbosity 2 --threads $THREADS
done
echo
done
echo "All adapters removed and the files are saved in $output"
echo "You may need to manually change the name of these files to your isolates..."
exit