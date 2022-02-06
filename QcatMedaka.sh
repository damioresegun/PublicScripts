#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N QcatMedaka
#$ -M dro@st-andrews.ac.uk
#Take shit through medaka. Use on phylo AND in MEDAKA CONDA

rawsPath="$HOME/AssembledIsolates_Reads/Qcat"
raconOut="$HOME/All_Racon/Qcat"
outputPath="$HOME/All_Medaka/Qcat"
QcAssem="$HOME/Assembly_QC/AssemStats"
mkdir -p $outputPath
mkdir -p $QcAssem
module load assembly-stats/gitv0_003a372

for i in $raconOut/* ; 
do 
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/${fname}
#cd $outputPath/${fname}
#no more than 8 threads
medaka_consensus -i $rawsPath/${fname}*.fastq -d $raconOut/${fname}/${fname}_iter4.fasta -o $outputPath/${fname} -t 8
#rename the output
mv $outputPath/${fname}/consensus.fasta $outputPath/${fname}/${fname}_consensus.fasta
# take some assembly stats of this
assembly-stats $outputPath/${fname}/${fname}_consensus.fasta >> ${QcAssem}/QcatAssemStats.txt
echo "Isolate ${fname} done"
done
echo "All done"
