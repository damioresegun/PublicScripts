#!/bin/bash
cPATH=/home/dami/Videos/Mapping
for g in $(find $cPATH -name *_genome); do
  echo $g
  cd $g
  fileNamefull=$(basename "$g")
	echo $fileNamefull
   isoName=$(echo ${fileNamefull} |sed s/\_genome//g)
    echo $isoName
    run="rm test_file"
    fg="rm Total_Sample_CNV_out.txt"
    ft="rm putput_file.txt"
    echo $run
    eval $run
    eval $fg
    eval $ft
done