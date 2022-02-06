#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N AssemQC
## add $ after first # to activate these
#$ -M dro@st-andrews.ac.uk
#Script to take assemblies through for QC-ing. This can be adjusted for raw, clean and polished assemblies. All that need to be changed is the input and output directories. The script will align the reads vs the assembly, sort the file, then take the stats of the alignment. It will also take in the assembly-stats of the assembly files, put it in txt file and then a CSV.
# intermediate files will be deleted i.e. SAM files and unsorted BAM files. Sorted BAM files will be renamed to just .bam rather than _sorted.bam

# For future implementations:
#QUAST, BUSCO, Scaffoldingstats.pl?

set -e
#set input and output folders
rawReads="$HOME/AssembledIsolates_Reads/Guppy"
assembly="$HOME/All_De_Novo/Guppy"
outSave="$HOME/Assembly_QC/Guppy"

#inputFol="$HOME/All_AdapterRemoved/Guppy"
#fastqSave="$HOME/All_Isolate_HumanUnMapped/Guppy"

mkdir -p $outSave
#mkdir -p $fastqSave
#gupsta=$(dirname $inputFol)
#cd $gupsta
#state the current working directory
hre=$(pwd)
echo $hre
#load necessary modules
module load assembly-stats/gitv0_003a372
module load samtools/1.6
for f in ${rawReads}/*
do
#get the filename
file=$(basename $f)
echo $file
#remove the extension
fname=${file%VsHuman*_*.*}
echo $fname
echo
#set the alignment name to the path and the file
alname=${outSave}/${fname}
echo $alname
#align to human reference
minimap2 -ax map-ont $assembly/$fname/${fname}_clean.contigs.fasta $rawReads/${fname}*.fastq > ${alname}_readsVsCleanAssem.sam -t 16
# now convert to bam, sort and delete intermediate files
echo "alignment done"
samtools view -bS ${alname}_readsVsCleanAssem.sam > ${alname}_readsVsCleanAssem.bam
echo
echo "bam conversion done"
samtools sort -o ${alname}_readsVsCleanAssem_sorted.bam ${alname}_readsVsCleanAssem.bam
samtools flagstat ${alname}_readsVsCleanAssem_sorted.bam >> ${outSave}/Guppy_readsVsCleanAssem_stats.txt
#delete the alignment files
rm ${alname}_readsVsCleanAssem.sam
rm ${alname}_readsVsCleanAssem.bam
rm ${alname}_readsVsCleanAssem_sorted.bam

echo "all done"
#cd $hre
done
#convert to csv
#cat ${outSave}/Guppy_readsVsCleanAssem_stats.txt | tr -s '[:blank:]' ',' > ${outSave}/Guppy_readsVsCleanAssem_stats.csv
