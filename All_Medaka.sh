#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N VsPk_CanuMedaka
#$ -r y
#$ -j y
## add $ after first # to activate these
#$ -q centos7.q
#$ -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk

#Take shit through medaka. Use on phylo AND in MEDAKA CONDA

#add the path to the raw fastq files used for de-novo
rawsPath="$HOME/All_Isolate_PkMapped_Formatted/"
# add path to the racon outputs
raconOut="$HOME/All_Racon/VsPk/Canu"
#add path to output for medaka
outputPath="$HOME/All_Medaka/VsPk/Canu"
# add path to place assembly stats
QcAssem="$HOME/All_AssemblyQC/Medaka/VsPk/Canu"
#use no more than 8 threads
THREADS="8"
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
medaka_consensus -i $rawsPath/${fname}*.fastq -d $raconOut/${fname}/${fname}_iter4.fasta -o $outputPath/${fname} -t $THREADS
#rename the output
mv $outputPath/${fname}/consensus.fasta $outputPath/${fname}/${fname}_consensus.fasta
# take some assembly stats of this
assembly-stats $outputPath/${fname}/${fname}_consensus.fasta >> ${QcAssem}/MaxGuppyAssemStats.txt
echo "Isolate ${fname} done"
done
echo "All done"
