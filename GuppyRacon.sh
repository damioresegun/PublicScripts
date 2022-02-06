#!/bin/bash
#$ -cwd
#exec 1 2>&1 | tee GupRacon_log.out

#output to a log file -- only if not using queue system
#exec 1>Racon_log.out 2>&1

#echo "$(date) : part 1 - start" >&3
# Script to take through data for racon polishing. Will carry out racon polishing through 4 rounds in order to take this through to Medaka. 

contigsPath="$HOME/All_De_Novo/Guppy"
rawsPath="$HOME/All_Isolate_HumanUnMapped/Guppy"
outputPath="$HOME/All_Racon/Guppy"
alignPath="$HOME/All_Reads_Vs_Assembly/iter1/Guppy"
iterVread="$HOME/All_Reads_Vs_Assembly/Racon/Guppy"
raconPath="$HOME/Tools/racon/build/bin/racon"
#spacePath="$HOME/Tools/SSPACE-LongRead_v1-1"
#perlPath="/usr/bin/perl"
mkdir -p $iterVread
for i in $contigsPath/* ;
do
fname=$(basename "$i")
echo "You are working on the ${fname} folder"
mkdir -p $outputPath/${fname}_iter1
cd $outputPath/${fname}_iter1
# Do iteration one for each isolate
$raconPath $rawsPath/${fname}.fastq ${alignPath}/${fname}_readsVsCleanAssem.sam $contigsPath/${fname}/${fname}_clean.contigs.fasta -t 16 > ${fname}_iter1.fasta
echo "Iteration one done for ${fname}"
pwd
echo "#######################################################################################"
echo
echo "Moving to iteration two"
echo "#######################################################################################"
# aligning output of iteration 1 to the raw reads
minimap2 -ax map-ont ${fname}_iter1.fasta ${rawsPath}/${fname}.fastq > ${iterVread}/${fname}_readsVsRacon1.sam
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}.fastq ${iterVread}/${fname}_readsVsRacon1.sam ${fname}_iter1.fasta -t 16 > ${fname}_iter2.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view -bS ${iterVread}/${fname}_readsVsRacon1.sam > ${iterVread}/${fname}_readsVsRacon1.bam
    # sort the bam file
    samtools sort -o ${iterVread}/${fname}_readsVsRacon1_sorted.bam ${iterVread}/${fname}_readsVsRacon1.bam
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
minimap2 -ax map-ont ${fname}_iter2.fasta ${rawsPath}/${fname}.fastq > ${iterVread}/${fname}_readsVsRacon2.sam
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}.fastq ${iterVread}/${fname}_readsVsRacon2.sam ${fname}_iter2.fasta -t 16 > ${fname}_iter3.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view -bS ${iterVread}/${fname}_readsVsRacon2.sam > ${iterVread}/${fname}_readsVsRacon2.bam
    # sort the bam file
    samtools sort -o ${iterVread}/${fname}_readsVsRacon2_sorted.bam ${iterVread}/${fname}_readsVsRacon2.bam
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
minimap2 -ax map-ont ${fname}_iter3.fasta ${rawsPath}/${fname}.fastq > ${iterVread}/${fname}_readsVsRacon3.sam
    echo
    echo "Alignment done"
$raconPath $rawsPath/${fname}.fastq ${iterVread}/${fname}_readsVsRacon3.sam ${fname}_iter3.fasta -t 16 > ${fname}_iter4.fasta
# cleaning up. converting and sorting sam file and zipping it
    #convert sam file to bam file
    samtools view -bS ${iterVread}/${fname}_readsVsRacon3.sam > ${iterVread}/${fname}_readsVsRacon3.bam
    # sort the bam file
    samtools sort -o ${iterVread}/${fname}_readsVsRacon3_sorted.bam ${iterVread}/${fname}_readsVsRacon3.bam
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