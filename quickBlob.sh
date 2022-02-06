#$ -cwd
#$ -m e
#$ -V
#$ -N All_DeN_CanuBlobing2
## add $ after first # to activate these
# -q centos7.q
# -l hostname="phylo"
#$ -M dro@st-andrews.ac.uk
# script to take your assemblies through for contamination detection using minimap2, samtools, blast and blobtools.
#set -e
## NOTE: RUN THIS IN THE IN THE BLOBTOOL CONDA ENVIRONMENT!!!!!!!!!!
# set folder for inputs and outputs
export PATH=/shelf/apps/ncbi-blast-2.7.1+/bin/:$PATH
export BLASTDB=/shelf/public/blastntnr/blastDatabases
## inputs
rawReads="$HOME/All_AdapterRemoved"
assembly="$HOME/All_De_Novo/All_DeNoVo/Canu/quicky"
##outputs
alignOutput="$HOME/All_Reads_Vs_Assembly/All_DeNovo/Canu"
blastOutput="$HOME/All_BLAST/All_DeNovo/Canu"
blobOutput="$HOME/All_BlobTools/All_DeNovo/Canu"

THREADS="32"

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
	#canu
    # minimap2 for alignment against the raw reads
    minimap2 -ax map-ont ${assembly}/${fname}/${fname}.contigs.fasta ${rawReads}/${fname}*.fastq > ${alignOutput}/${fname}_readsVsAssem.sam -t $THREADS
	#flye
	#minimap2 -ax map-ont ${assembly}/${fname}/assembly.fasta ${rawReads}/${fname}*.fastq > ${alignOutput}/${fname}_readsVsAssem.sam -t $THREADS
    echo
    echo "Alignment done. Moving on to sam file conversion and sorting..."
    echo
    #convert sam file to bam file
    samtools view --threads $THREADS -bS ${alignOutput}/${fname}_readsVsAssem.sam > ${alignOutput}/${fname}_readsVsAssem.bam
    # sort the bam file
    samtools sort --threads $THREADS -o ${alignOutput}/${fname}_readsVsAssem_sorted.bam ${alignOutput}/${fname}_readsVsAssem.bam
    # clean up
    echo "Deleting intermediate alignment files..."
    #rm ${alignOutput}/${fname}_readsVsAssem.sam
    rm ${alignOutput}/${fname}_readsVsAssem.bam
    mv ${alignOutput}/${fname}_readsVsAssem_sorted.bam ${alignOutput}/${fname}_readsVsAssem.bam
    echo
    echo "SAM to BAM conversion done. Sorting of BAM file done."
    echo 
    echo "Moving on to Blasting..."
    echo
	#canu
    blastn -task megablast -query ${assembly}/${fname}/${fname}.contigs.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${fname}_vs_nt.out -num_threads $THREADS
	#flye
	#blastn -task megablast -query ${assembly}/${fname}/assembly.fasta -db nt -outfmt '6 qseqid staxids bitscore std scomnames sscinames sblastnames sskingdoms stitle' -evalue 1e-20 -out ${blastOutput}/${fname}_vs_nt.out -num_threads $THREADS
    echo
    echo "Blast is finished... Continuing to blobtools"
    echo
    echo "Creating a JSON file of the data"
    mkdir -p ${blobOutput}/${fname}
	#canu
    blobtools create -i ${assembly}/${fname}/${fname}.contigs.fasta -b ${alignOutput}/${fname}_readsVsAssem.bam -t ${blastOutput}/${fname}_vs_nt.out -o ${blobOutput}/${fname}/
	#flye
	#blobtools create -i ${assembly}/${fname}/assembly.fasta -b ${alignOutput}/${fname}_readsVsAssem.bam -t ${blastOutput}/${fname}_vs_nt.out -o ${blobOutput}/${fname}/
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
    blobtools covplot -i ${blobOutput}/${fname}/blobDB.json -c ${blobOutput}/${fname}/${fname}_readsVsAssem.bam.cov -o ${blobOutput}/${fname}/ --max 1e03
    echo
    echo "Covplots done"
    done
echo "This iteration is complete. Please look at your results to manually remove contamination. Run this script once contamination has been removed to confirm no more contamination is present"
