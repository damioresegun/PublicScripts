#!/bin/bash
#$ -cwd
#$ -m e
#$ -V
#$ -r y
#$ -j y
#$ -N All_DeN_Flye_iter1Blob
## add $ after first # to activate these
#$ -q all.q
#$ -l hostname="node8"
#$ -M dro@st-andrews.ac.uk

#Script to take data through part of decontamination after blobtools. Should go through the blobDB text file to keep only entries that have Apicomplexa, no-hit, undef and # thus theoretically keeping only the phylum for Apicomplexa and other undefined sequences. After this, the script will get the names of the remaining sequences and put them in a list to then use to filter the contigs file. The resulting contig file will then be taken through blobtools again for confirmation of decontamination before being taken forward for polishing
# script to take your assemblies through for contamination detection using minimap2, samtools, blast and blobtools after decontamination. It is currently made for your first decontamintation iteration so just change numbers as necessary
## NOTE: RUN THIS IN THE IN THE BLOBTOOL CONDA ENVIRONMENT!!!!!!!!!!
export PATH=/shelf/apps/ncbi-blast-2.7.1+/bin/:$PATH
export BLASTDB=/shelf/public/blastntnr/blastDatabases
# set folder for inputs and outputs
## inputs
rawReads="$HOME/All_AdapterRemoved"
assembly="$HOME/All_De_Novo/All_DeNoVo/Flye"
THREADS="16"
##outputs
alignOutput="$HOME/All_Reads_Vs_Assembly/All_DeNoVo/Clean_Flye"
blastOutput="$HOME/All_BLAST/All_DeNoVo/Clean_Flye"
blobOutput="$HOME/All_BlobTools/All_DeNoVo/Clean_Flye"

## make the output folders if not existing
mkdir -p $alignOutput
mkdir -p $blastOutput
mkdir -p $blobOutput

for folder in $assembly/*
do
    echo "The isolate you are working on is $folder"
    fname=$(basename "$folder")    
    #echo $fname
    echo
    echo "Beginning alignment of assembly to reads..."
    echo
    # minimap2 for alignment against the raw reads
	#canu
    #minimap2 -ax map-ont ${assembly}/${fname}/${fname}_clean.contigs.fasta ${rawReads}/${fname}*.fastq > ${alignOutput}/${fname}_readsVsCleanAssem.sam -t $THREADS
	#flye
	minimap2 -ax map-ont ${assembly}/${fname}/assembly.fasta ${rawReads}/${fname}*.fastq > ${alignOutput}/${fname}_readsVsCleanAssem.sam -t $THREADS
    echo
    echo "Alignment done. Moving on to sam file conversion and sorting..."
    echo
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${alignOutput}/${fname}_readsVsCleanAssem.sam > ${alignOutput}/${fname}_readsVsCleanAssem.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${alignOutput}/${fname}_readsVsCleanAssem_sorted.bam ${alignOutput}/${fname}_readsVsCleanAssem.bam
    # clean up
    echo "Deleting intermediate alignment files..."
	# sam files are needed for racon. so they are not deleted here
    #rm ${alignOutput}/${fname}_readsVsCleanAssem.sam
    rm ${alignOutput}/${fname}_readsVsCleanAssem.bam
    mv ${alignOutput}/${fname}_readsVsCleanAssem_sorted.bam ${alignOutput}/${fname}_readsVsCleanAssem.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Moving on to Blasting..."
    echo
	#canu
    #blastn -task megablast -query ${assembly}/${fname}/${fname}_clean.contigs.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${fname}_vs_nt.out -num_threads $THREADS
	#flye
	blastn -task megablast -query ${assembly}/${fname}/clean_assembly.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${fname}_vs_nt.out -num_threads $THREADS
    echo
    echo "Blast is finished... Continuing to blobtools"
    echo
    echo "Creating a JSON file of the data"
    mkdir -p ${blobOutput}/${fname}
	#canu
    #blobtools create -i ${assembly}/${fname}/${fname}_clean.contigs.fasta -b ${alignOutput}/${fname}_readsVsCleanAssem.bam -t ${blastOutput}/${fname}_vs_nt.out -o ${blobOutput}/${fname}/
	#flye
	blobtools create -i ${assembly}/${fname}/clean_assembly.fasta -b ${alignOutput}/${fname}_readsVsCleanAssem.bam -t ${blastOutput}/${fname}_vs_nt.out -o ${blobOutput}/${fname}/
    echo
    echo "Creating blobtools view file..."
    blobtools view -i ${blobOutput}/${fname}/blobDB.json -o ${blobOutput}/${fname}/
    echo
    echo "Blobtools table created. Plotting blobplots..."
    echo
    blobtools plot -i ${blobOutput}/${fname}/blobDB.json -o ${blobOutput}/${fname}/
    echo
    echo "Blobplot complete. Drawing covplot..."
    echo
    blobtools covplot -i ${blobOutput}/${fname}/blobDB.json -c ${blobOutput}/${fname}/${fname}_readsVsCleanAssem.bam.cov -o ${blobOutput}/${fname}/ --max 1e03
    echo
    echo "Covplots done"
    done
echo "This iteration is complete. Please look at your results to manually remove contamination. Run this script once contamination has been removed to confirm no more contamination is present"
