#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -N VsPkFlye_GenomeBlobing
#$ -r y
#$ -j y
## add $ after first # to activate these
#$ -q all.q
#$ -l hostname="node8"
#$ -M dro@st-andrews.ac.uk
# script to take your assemblies through for contamination detection using minimap2, samtools, blast and blobtools.
# REQUIRES A TEXT FILE WITH THE FILE PATHS FOR THE RAQ FASTQ FILES OF ISOLATES THAT ASSEMBLED!!!
# Will need to change the different commands depending on if Canu or Flye!!!!
set -e
## NOTE: RUN THIS IN THE IN THE BLOBTOOL CONDA ENVIRONMENT!!!!!!!!!!
# set folder for inputs and outputs
export PATH=/shelf/apps/ncbi-blast-2.7.1+/bin/:$PATH
export BLASTDB=/shelf/public/blastntnr/blastDatabases
## inputs
#change these below as appropriate
data="$HOME/InputFiles/Pk_Flye_AssembledIsolates.txt"
assembly="$HOME/All_CompanionGenomes/VsPk/Flye"
##outputs
alignOutput="$HOME/All_Reads_Vs_Genome/VsPk/Flye"
blastOutput="$HOME/All_GenomeBLAST/VsPk/Flye"
blobOutput="$HOME/All_GenomeBlobTools/VsPk/Flye"

THREADS="16"

## make the output folders if not existing
mkdir -p $alignOutput
mkdir -p $blastOutput
mkdir -p $blobOutput

while read -r line;
#for folder in $assembly/*
do
    echo "The isolate you are working on is $line"
    fname=$(basename "$line")    
    echo $fname
    echo
	#get the isolate name
	#iso=${fname%VsHuman.*}
	iso=${fname%_*.*}
	echo $iso
    echo "Beginning alignment of assembly to reads..."
    echo
    # minimap2 for alignment against the raw reads
    minimap2 -ax map-ont ${assembly}/${iso}.fasta ${line} > ${alignOutput}/${iso}_readsVsGenome.sam -t $THREADS
    echo
    echo "Alignment done. Moving on to sam file conversion and sorting..."
    echo
    #convert sam file to bam file
    samtools view --thread $THREADS -bS ${alignOutput}/${iso}_readsVsGenome.sam > ${alignOutput}/${iso}_readsVsGenome.bam
    # sort the bam file
    samtools sort --thread $THREADS -o ${alignOutput}/${iso}_readsVsGenome_sorted.bam ${alignOutput}/${iso}_readsVsGenome.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    #rm ${alignOutput}/${iso}_readsVsGenome.sam
    rm ${alignOutput}/${iso}_readsVsGenome.bam
    mv ${alignOutput}/${iso}_readsVsGenome_sorted.bam ${alignOutput}/${iso}_readsVsGenome.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Moving on to Blasting..."
    echo
	#blast
    blastn -task megablast -query ${assembly}/${iso}.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${iso}_vs_nt.out -num_threads $THREADS
    echo
    echo "Blast is finished... Continuing to blobtools"
    echo
    echo "Creating a JSON file of the data"
    mkdir -p ${blobOutput}/${iso}
	
	#canu blobtools
    blobtools create -i ${assembly}/${iso}.fasta -b ${alignOutput}/${iso}_readsVsGenome.bam -t ${blastOutput}/${iso}_vs_nt.out -o ${blobOutput}/${iso}/
    echo
    echo "Creating blobtools view file..."
    blobtools view -i ${blobOutput}/${iso}/blobDB.json -o ${blobOutput}/${iso}/
    echo
    echo "Blobtools table created. Plotting blobplots..."
    echo
    blobtools plot -i ${blobOutput}/${iso}/blobDB.json -o ${blobOutput}/${iso}/
    echo
    echo "Blobplot complete. Drawing covplot..."
    echo
    blobtools covplot -i ${blobOutput}/${iso}/blobDB.json -c ${blobOutput}/${iso}/${iso}_readsVsGenome.bam.cov -o ${blobOutput}/${iso}/ --max 1e03
    echo
    echo "Covplots done"
    done < ${data}
echo "This iteration is complete. Please look at your results to manually remove contamination. Run this script once contamination has been removed to confirm no more contamination is present"

