#!/bin/bash
cPATH=/home/dami/Mapping/
for g in $(find $cPATH -name *_genome); do
  echo $g
  cd $g
  fileNamefull=$(basename "$g")
	echo $fileNamefull
   isoName=$(echo ${fileNamefull} |sed s/\_genome//g)
   echo $isoName
   run="Rscript ~/CopyCat/subscripts/coverage_grapher.R $g $isoName Pf3D7_12_v3 46788 56805 /home/dami/CopyCat/"
   echo $run
   eval $run
   mv ~/ERR*.pdf ~/FinalGraphs
done
