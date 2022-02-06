#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -q all.q
#$ -l hostname="node6"
#$ -N All_DeN_CanuRacon
#$ -M dro@st-andrews.ac.uk

# Script to take through data for racon polishing. Will carry out racon polishing through 4 rounds in order to take this through to Medaka. 

#output to a log file -- only if not using queue system
#exec 1>All_DeN_FlyeRacon.out 2>&1

contigsPath="$HOME/All_De_Novo/All_DeNoVo/Canu"
rawsPath="$HOME/All_AdapterRemoved"
outputPath="$HOME/All_Racon/All_DeNoVo/Canu"
alignPath="$HOME/All_Reads_Vs_Assembly/All_DeNoVo/Clean_Canu"
iterVread="$HOME/All_Reads_Vs_Assembly/Racon/All_DeNoVo/Canu"
raconPath="$HOME/Tools/racon/build/bin/racon"
THREADS="16"
#set -e
#spacePath="$HOME/Tools/SSPACE-LongRead_v1-1"
#perlPath="/usr/bin/perl"
mkdir -p $iterVread
mkdir -p $outputPath
for i in $contigsPath/* ;
do
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/${fname}
cd $outputPath/${fname}
# Do iteration one for each isolate
#canu
$raconPath $rawsPath/${fname}*.fastq ${alignPath}/${fname}_readsVsCleanAssem.sam $contigsPath/${fname}/${fname}_clean.contigs.fasta -t $THREADS > ${fname}_iter1.fasta

#flye
# $raconPath $rawsPath/${fname}*.fastq ${alignPath}/${fname}_readsVsCleanAssem.sam $contigsPath/${fname}/clean_assembly.fasta -t $THREADS > ${fname}_iter1.fasta

echo "Iteration one done for ${fname}"
pwd
echo "#######################################################################################"
echo
echo "Moving to iteration two"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter1.fasta ${rawsPath}/${fname}*.fastq > ${iterVread}/${fname}_readsVsRacon1.sam -t $THREADS
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}*.fastq ${iterVread}/${fname}_readsVsRacon1.sam ${fname}_iter1.fasta -t $THREADS > ${fname}_iter2.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${iterVread}/${fname}_readsVsRacon1.sam > ${iterVread}/${fname}_readsVsRacon1.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${iterVread}/${fname}_readsVsRacon1_sorted.bam ${iterVread}/${fname}_readsVsRacon1.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    rm ${iterVread}/${fname}_readsVsRacon1.sam
    rm ${iterVread}/${fname}_readsVsRacon1.bam
    mv ${iterVread}/${fname}_readsVsRacon1_sorted.bam ${iterVread}/${fname}_readsVsRacon1.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Iteration two done for ${fname}" 
echo "#######################################################################################"
echo
echo "Moving to iteration three"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter2.fasta ${rawsPath}/${fname}*.fastq > ${iterVread}/${fname}_readsVsRacon2.sam -t $THREADS
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}*.fastq ${iterVread}/${fname}_readsVsRacon2.sam ${fname}_iter2.fasta -t $THREADS > ${fname}_iter3.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${iterVread}/${fname}_readsVsRacon2.sam > ${iterVread}/${fname}_readsVsRacon2.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${iterVread}/${fname}_readsVsRacon2_sorted.bam ${iterVread}/${fname}_readsVsRacon2.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    rm ${iterVread}/${fname}_readsVsRacon2.sam
    rm ${iterVread}/${fname}_readsVsRacon2.bam
    mv ${iterVread}/${fname}_readsVsRacon2_sorted.bam ${iterVread}/${fname}_readsVsRacon2.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Iteration three done for ${fname}" 
echo "#######################################################################################"
echo
echo "Moving to iteration four"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter3.fasta ${rawsPath}/${fname}*.fastq > ${iterVread}/${fname}_readsVsRacon3.sam -t $THREADS
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}*.fastq ${iterVread}/${fname}_readsVsRacon3.sam ${fname}_iter3.fasta -t $THREADS > ${fname}_iter4.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${iterVread}/${fname}_readsVsRacon3.sam > ${iterVread}/${fname}_readsVsRacon3.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${iterVread}/${fname}_readsVsRacon3_sorted.bam ${iterVread}/${fname}_readsVsRacon3.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    rm ${iterVread}/${fname}_readsVsRacon3.sam
    rm ${iterVread}/${fname}_readsVsRacon3.bam
    mv ${iterVread}/${fname}_readsVsRacon3_sorted.bam ${iterVread}/${fname}_readsVsRacon3.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Iteration four done for ${fname}" 
done
