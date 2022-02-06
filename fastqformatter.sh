#!/bin/bash
# Script to carry out formatting of fastq files that are being directed to. 
# Rationale: often fastq files might have some odd encoding, line ends etc. This script seeeks those out and attempts to fix it by calling the tool script convert_fq_to_fa.py.
## has to be run in the Flye or Canu conda environment as it needs bioawk and seqkit! 
inp="$HOME/All_Isolate_PkMapped"
outr="$HOME/All_Isolate_PkMapped_Formatted/"
mkdir -p $outr
module load python/3.6.4
for i in $inp/*; 
do
ls $i
echo
echo "The file you are working on is $i"
    fname=$(basename "$i")    
    echo $fname
    #nmr=${fname%.*}
	nmr=${fname%VsPk*_*.*}
    echo "$nmr"
    # format the file
    python3 $HOME/Tools/convert_fq_to_fa.py -i $i -o $outr/${nmr}_reformat.fastq
	# remove empty reads
	bioawk -cfastx 'length($seq) > 1 {print "@"$name"\n"$seq"\n+\n"$qual}' $outr/${nmr}_reformat.fastq > $outr/${nmr}_empty0.fastq
	#remove duplicates
	rm $outr/${nmr}_reformat.fastq
    seqkit rename $outr/${nmr}_empty0.fastq > $outr/${nmr}_reformat.fastq
	rm $outr/${nmr}_empty0.fastq
    echo "done"
    
 done 
