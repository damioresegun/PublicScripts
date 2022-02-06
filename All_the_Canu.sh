#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N PkCanu
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk
# Uncomment for subsect of adapter removed data
#Aligned="/storage/home/users/dro/AssembledIsolates_Reads/Guppy/quick"
### Entire adapter removed data
Aligned="$HOME/All_Isolate_PkMapped_Formatted"
THREADS="32"
#cd $Aligned
#files=$(cat $HOME/Isolates_To_Be_Assembled.txt)

for file in $Aligned/*
do
  ls $file
    echo "The file you are working on is $file"
    fname=$(basename "$file")    
    echo $fname
    #nmr=${fname%VsHuman*_*.*}
	nmr=${fname%_reformat.*}
    echo "You are ${nmr}"
    mkdir -p $HOME/All_De_Novo/VsPk/Canu/${nmr}
    #mkdir -p $HOME/DenovoGuppy/${nmr}
    dest="$HOME/All_De_Novo/VsPk/Canu/${nmr}"
    #dest="$HOME/DenovoGuppy/${nmr}"
    #cd ${Aligned}/${fname}
    canu -p ${nmr} -d ${dest} genomeSize=25m -nanopore-raw ${file} -useGrid=False -maxMemory=120 -maxThreads=$THREADS -corOutCoverage=999
    #canu -p ${nmr} -d ${dest} genomeSize=24m -nanopore-raw ${file} -useGrid=False -maxMemory=120 -maxThreads=16
    cd $Aligned
done
#canu gnuplotTested=true -p sks339 \ "${nmr}VsHuman.fastq"
#-d /storage/home/users/dro/All_De_Novo/newAssem/sks339 genomeSize=24m -nanopore-raw 'sks339VsHuman_sorted_unmapped.fastq.gz' \
#-useGrid=False -maxMemory=120 \
#-maxThreads=12


