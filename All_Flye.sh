#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N Pk_Flye

#$ -M dro@st-andrews.ac.uk

### Entire adapter removed data
Aligned="/storage/home/users/dro/All_Isolate_PkMapped_Formatted"
THREADS="32"
#number of polishing iterations
#ITER="4"
#cd $Aligned

for file in $Aligned/*
do
  ls $file
    echo "The file you are working on is $file"
    fname=$(basename "$file")    
    echo $fname
    #nmr=${fname%VsHuman*_*.*}
	nmr=${fname%_reformat.*}
    echo "You are ${nmr}"
    mkdir -p $HOME/All_De_Novo/VsPk/Flye/${nmr}
	#mkdir -p $HOME/All_De_Novo/VsPk/Flye/Meta/${nmr}
    dest="$HOME/All_De_Novo/VsPk/Flye/${nmr}"
	#dest2="$HOME/All_De_Novo/VsPk/Flye/Meta/${nmr}"
    #cd ${Aligned}/${fname}
    # flye --nano-raw ${file} --genome-size 25m --out-dir ${dest} --threads $THREADS -i $ITER
	# flye --nano-raw ${file} --genome-size 25m --out-dir ${dest2} --threads $THREADS -i $ITER --meta
	flye --nano-raw ${file} --genome-size 25m --out-dir ${dest} --threads $THREADS
	#flye --nano-raw ${file} --genome-size 25m --out-dir ${dest2} --threads $THREADS --meta
    cd $Aligned
done
#canu gnuplotTested=true -p sks339 \ "${nmr}VsHuman.fastq"
#-d /storage/home/users/dro/All_De_Novo/newAssem/sks339 genomeSize=24m -nanopore-raw 'sks339VsHuman_sorted_unmapped.fastq.gz' \
#-useGrid=False -maxMemory=120 \
#-maxThreads=12
