#!/bin/bash
#$ -cwd
#$ -V
#$ -m e
#$ -N AligningStep
#$ -M dro@st-andrews.ac.uk
#for aligning and sorting SAM files to BAM files for all sequences within the input folder and to the output folder
# intermediate files will be deleted i.e. SAM files and unsorted BAM files. Sorted BAM files will be renamed to just .bam rather than _sorted.bam
##NOTE: MIGHT HAVE TO RUN THIS IN A SAMTOOLS CONDA ENVIRONMENT
set -e
#set input and output folders
inputFol="$HOME/All_AdapterRemoved"
outSave="$HOME/All_Aligned"
fastqSave="$HOME/All_Isolate_HumanUnMapped"
THREADS="16"
statsFol="$HOME/All_Stats/"
# load necessary modules
module load samtools/1.6
mkdir -p $fastqSave
mkdir -p $statsFol/postAlign_Extraction
#gupsta=$(dirname $inputFol)
#cd $gupsta
#state the current working directory
hre=$(pwd)
echo $hre
#load necessary modules
module load assembly-stats/gitv0_003a372
#make an unmapped folder
mkdir -p $outSave/Unmapped
unmpd=$outSave/Unmapped
#mkdir -p $unmpd/VsLapp
#lapr=$unmpd/VsLapp
for f in ${inputFol}/*
do
#get the filename
file=$(basename $f)
echo $file
#remove the extension
fname="${file%.*}"
echo $fname
echo
#set the alignment name to the path and the file
alname=${outSave}/${fname}
echo $alname
#align to human reference
minimap2 -ax map-ont $HOME/Index/HumanGRCh38_patch12_genomic.mmi $inputFol/${fname}.fastq > ${alname}VsHuman.sam -t $THREADS
# now convert to bam, sort and delete intermediate files
echo "alignment done"
samtools view --threads $THREADS -bS ${alname}VsHuman.sam > ${alname}VsHuman.bam
echo
echo "bam conversion done"
samtools sort --threads $THREADS -o ${alname}VsHuman_sorted.bam ${alname}VsHuman.bam 
rm ${alname}VsHuman.sam
rm ${alname}VsHuman.bam
mv ${alname}VsHuman_sorted.bam ${alname}VsHuman.bam
echo ${fname} >> All_FlagstatMappedVsHuman_stats.txt
samtools flagstat --threads $THREADS ${alname}VsHuman.bam >> ${statsFol}/postAlign_Extraction/All_FlagstatMappedVsHuman_stats.txt
# extract unmapped reads
cd $unmpd
pwd
samtools view --threads $THREADS -f 4 -b ${alname}VsHuman.bam > ${unmpd}/${fname}VsHuman_unmapped.bam
echo "unmapped extraction done"
#convert unmapped reads to fastq
bedtools bamtofastq -i ${unmpd}/${fname}VsHuman_unmapped.bam -fq ${unmpd}/${fname}VsHuman_unmapped.fastq
echo
#check the stats
assembly-stats ${unmpd}/${fname}VsHuman_unmapped.fastq > ${statsFol}/postAlign_Extraction/${fname}_VsHumanUnmappedStats.txt
#move fastqs to their own folder
mv ${unmpd}/${fname}VsHuman_unmapped.fastq $fastqSave

#align against P knowlesi reference
#minimap2 -ax map-ont $HOME/Index/PkLappRefSeq_index.mmi ${unmpd}/${fname}VsHuman_unmapped.fastq > $lapr/${fname}VsLapp.sam
#samtools view -bS $lapr/${fname}VsLapp.sam > $lapr/${fname}VsLapp.bam
#echo
#echo "bam conversion done"
#samtools sort -o $lapr/${fname}VsLapp_sorted.bam $lapr/${fname}VsLapp.bam
#rm $lapr/${fname}VsLapp.sam
#rm $lapr/${fname}VsLapp.bam
#mv $lapr/${fname}VsLapp_sorted.bam $lapr/${fname}VsLapp.bam
#samtools flagstat $lapr/${fname}VsLapp.bam >> $lapr/${fname}VsLapp_stats.txt
echo "all done"
#cd $hre
done
cat ${statsFol}/postAlign_Extraction/*_VsHumanUnmappedStats.txt >> ${statsFol}/postAlign_Extraction/All_VsHumanUnmappedStats.txt
cat ${statsFol}/postAlign_Extraction/All_VsHumanUnmappedStats.txt | tr -s '[:blank:]' ',' > ${statsFol}/postAlign_Extraction/All_VsHumanUnmappedStats.csv
rm ${statsFol}/postAlign_Extraction/All_VsHumanUnmappedStats.txt
exit
